// This file is part of MUSIC
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2010-2024 by Oliver Hahn
// 
// monofonIC is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// monofonIC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

//... for testing purposes.............
// #define DEGRADE_RAND1
// #define DEGRADE_RAND2
//.....................................

#pragma once

#define DEF_RAN_CUBE_SIZE 32

#include <fstream>
#include <algorithm>
#include <map>
#include <array>
#include <cstdio>
#include <string>
#include <vector>
#include <unistd.h>

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

	//! cache of the white noise per level, so it is generated only once per run
	struct cached_noise
	{
		std::array<size_t, 3> n;
		std::array<int, 3> o;
		size_t size;
		std::vector<T> data; // used when the cache is kept in memory
		std::string fname;   // used when the cache is kept on disk
	};
	std::string cache_mode_; // "disk" (default), "memory" or "no"
	std::string cache_dir_;
	std::map<int, cached_noise> noise_cache_;

public:
	//! constructor
	explicit random_number_generator(config_file &cf)
			: pcf_(&cf) //, prefh_( &refh )
	{
		levelmin_ = pcf_->get_value<int>("setup", "levelmin");
		levelmax_ = pcf_->get_value<int>("setup", "levelmax");
		generator_ = select_RNG_plugin(cf);

		cache_mode_ = pcf_->get_value_safe<std::string>("random", "cache_noise", "disk");
		cache_dir_ = pcf_->get_value_safe<std::string>("random", "cache_dir", ".");
		if (pcf_->get_value_safe<bool>("setup", "WhiteNoiseOnly", false))
			cache_mode_ = "no"; // every level is generated only once anyway
		if (cache_mode_ != "disk" && cache_mode_ != "memory" && cache_mode_ != "no")
		{
			music::elog << "Unknown [random] cache_noise = " << cache_mode_ << " (use disk, memory or no)" << std::endl;
			throw std::runtime_error("Unknown value for [random] cache_noise");
		}
		if (cache_mode_ == "disk")
			music::ilog << "White noise is generated once per level and cached on disk in '" << cache_dir_ << "'" << std::endl;
		else if (cache_mode_ == "memory")
			music::ilog << "White noise is generated once per level and cached in memory" << std::endl;
		else
			music::ilog << "White noise is not cached and is regenerated for every field" << std::endl;
	}

	//! destructor, removes cache files
	~random_number_generator()
	{
		for (auto &c : noise_cache_)
			if (!c.second.fname.empty())
				std::remove(c.second.fname.c_str());
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
		if (cache_mode_ != "no")
		{
			auto it = noise_cache_.find(ilevel);
			if (it != noise_cache_.end() && it->second.n == std::array<size_t, 3>{A.nx_, A.ny_, A.nz_} &&
					it->second.o == std::array<int, 3>{A.ox_, A.oy_, A.oz_} && it->second.size == A.data_.size())
			{
				double wtime = get_wtime();
				if (cache_mode_ == "memory")
					A.data_ = it->second.data;
				else
				{
					std::ifstream ifs(it->second.fname, std::ios::binary);
					ifs.read(reinterpret_cast<char *>(A.data_.data()), A.data_.size() * sizeof(T));
					if (!ifs)
						throw std::runtime_error("Could not read white noise cache file " + it->second.fname);
				}
				music::ilog << "Reusing cached white noise for level " << ilevel << " (took " << get_wtime() - wtime << " s)" << std::endl;
				return;
			}
		}

		generator_->fill_grid(ilevel, A);

		if (cache_mode_ == "no")
			return;

		cached_noise c{{A.nx_, A.ny_, A.nz_}, {A.ox_, A.oy_, A.oz_}, A.data_.size(), {}, {}};
		if (cache_mode_ == "memory")
			c.data = A.data_;
		else
		{
			c.fname = cache_dir_ + "/music_noise_cache_" + std::to_string(getpid()) + "_level_" + std::to_string(ilevel) + ".bin";
			std::ofstream ofs(c.fname, std::ios::binary | std::ios::trunc);
			ofs.write(reinterpret_cast<const char *>(A.data_.data()), A.data_.size() * sizeof(T));
			if (!ofs)
				throw std::runtime_error("Could not write white noise cache file " + c.fname + " (set [random] cache_dir or cache_noise)");
		}
		if (noise_cache_.count(ilevel) && !noise_cache_[ilevel].fname.empty() && noise_cache_[ilevel].fname != c.fname)
			std::remove(noise_cache_[ilevel].fname.c_str());
		noise_cache_[ilevel] = std::move(c);
	}
};

using noise_generator = random_number_generator<real_t>;
