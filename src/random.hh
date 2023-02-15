/*

 random.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions
 for cosmological simulations

 Copyright (C) 2010-23 by  Oliver Hahn

*/

//... for testing purposes.............
// #define DEGRADE_RAND1
// #define DEGRADE_RAND2
//.....................................

#pragma once

#define DEF_RAN_CUBE_SIZE 32

#include <fstream>
#include <algorithm>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "general.hh"
#include "mesh.hh"
#include "mg_operators.hh"
#include "constraints.hh"
// #include "convolution_kernel.hh"
#include "density_grid.hh"

class RNG_plugin
{
protected:
	config_file *pcf_;									//!< pointer to config_file from which to read parameters
	const refinement_hierarchy *prefh_; //!< pointer to refinement hierarchy structure containing the grid sizes
public:
	explicit RNG_plugin(config_file &cf) //, const refinement_hierarchy& refh )
			: pcf_(&cf)											 //, prefh_( & refh )
	{
	}
	virtual ~RNG_plugin() {}
	virtual bool is_multiscale() const = 0;
	virtual void fill_grid(int level, DensityGrid<real_t> &R) = 0;
	virtual void initialize_for_grid_structure(const refinement_hierarchy &refh) = 0;
};

struct RNG_plugin_creator
{
	virtual RNG_plugin *create(config_file &cf) const = 0;
	virtual ~RNG_plugin_creator() {}
};

std::map<std::string, RNG_plugin_creator *> &
get_RNG_plugin_map();

void print_RNG_plugins(void);

template <class Derived>
struct RNG_plugin_creator_concrete : public RNG_plugin_creator
{
	//! register the plugin by its name
	RNG_plugin_creator_concrete(const std::string &plugin_name)
	{
		get_RNG_plugin_map()[plugin_name] = this;
	}

	//! create an instance of the plugin
	RNG_plugin *create(config_file &cf) const //, const refinement_hierarchy& refh ) const
	{
		return new Derived(cf); //, refh );
	}
};

typedef RNG_plugin RNG_instance;
RNG_plugin *select_RNG_plugin(config_file &cf); //, const refinement_hierarchy& refh );

/*!
 * @brief encapsulates all things for multi-scale white noise generation
 */
template <typename T>
class random_number_generator
{
protected:
	config_file *pcf_;
	// const refinement_hierarchy * prefh_;
	RNG_plugin *generator_;
	int levelmin_, levelmax_;

public:
	//! constructor
	explicit random_number_generator(config_file &cf)
			: pcf_(&cf) //, prefh_( &refh )
	{
		levelmin_ = pcf_->get_value<int>("setup", "levelmin");
		levelmax_ = pcf_->get_value<int>("setup", "levelmax");
		generator_ = select_RNG_plugin(cf);
	}

	//! destructor
	~random_number_generator()
	{
	}

	//! initialize_for_grid_structure
	void initialize_for_grid_structure(const refinement_hierarchy &refh)
	{
		generator_->initialize_for_grid_structure(refh);
	}

	//! load random numbers to a new array
	template <typename array>
	void load(array &A, int ilevel)
	{
		generator_->fill_grid(ilevel, A);
	}
};

using noise_generator = random_number_generator<real_t>;
