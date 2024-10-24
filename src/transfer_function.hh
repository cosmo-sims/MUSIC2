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

#pragma once

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <complex>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>

#include <Numerics.hh>
#include <config_file.hh>
#include <cosmology_parameters.hh>

enum tf_type
{
	delta_matter,
	delta_cdm,
	delta_baryon,
	theta_matter,
	theta_cdm,
	theta_baryon,
	delta_bc,
	theta_bc,
	delta_matter0,
	delta_cdm0,
	delta_baryon0,
	theta_matter0,
	theta_cdm0,
	theta_baryon0,
};

#define GSL_INTEGRATION_ERR 1e-5

//! Abstract base class for transfer functions
/*!
 This class implements a purely virtual interface that can be
 used to derive instances implementing various transfer functions.
 */
class transfer_function_plugin
{
public:
	config_file *pcf_;													//!< pointer to config_file from which to read parameters
	const cosmology::parameters &cosmo_params_; //!< cosmological parameters are stored here
	bool tf_distinct_;													//!< bool if density transfer function is distinct for baryons and DM
	bool tf_withvel_;														//!< bool if also have velocity transfer functions
	bool tf_withtotal0_;												//!< have the z=0 spectrum for normalisation purposes
	bool tf_velunits_;													//!< velocities are in velocity units (km/s)
	bool tf_isnormalised_;											//!< assume that transfer functions come already correctly normalised and need be re-normalised to a specified value

	//! constructor
	transfer_function_plugin(config_file &cf,const cosmology::parameters &cosmo_params)
			: pcf_(&cf), cosmo_params_(cosmo_params), tf_distinct_(false), tf_withvel_(false), tf_withtotal0_(false), tf_velunits_(false), tf_isnormalised_(false)
	{
	}

	//! destructor
	virtual ~transfer_function_plugin(){};

	//! initialise, i.e. prepare data for later usage
	virtual void intialise(void) {}

	//! compute value of transfer function at waven umber
	virtual double compute(double k, tf_type type) const = 0;

	//! return maximum wave number allowed
	virtual double get_kmax(void) const = 0;

	//! return minimum wave number allowed
	virtual double get_kmin(void) const = 0;

	//! return if density transfer function is distinct for baryons and DM
	bool tf_is_distinct(void) const
	{
		return tf_distinct_;
	}

	//! return if we also have velocity transfer functions
	bool tf_has_velocities(void) const
	{
		return tf_withvel_;
	}

	//! return if we also have a z=0 transfer function for normalisation
	bool tf_has_total0(void) const
	{
		return tf_withtotal0_;
	}

	//! return if velocity returned is in velocity or in displacement units
	bool tf_velocity_units(void) const
	{
		return tf_velunits_;
	}
};

//! Implements abstract factory design pattern for transfer function plug-ins
struct transfer_function_plugin_creator
{
	//! create an instance of a transfer function plug-in
	virtual std::unique_ptr<transfer_function_plugin> create(config_file &cf, const cosmology::parameters& cp) const = 0;

	//! destroy an instance of a plug-in
	virtual ~transfer_function_plugin_creator() {}
};

//! Write names of registered transfer function plug-ins to stdout
std::map<std::string, transfer_function_plugin_creator *> &get_transfer_function_plugin_map();
void print_transfer_function_plugins(void);

//! Concrete factory pattern for transfer function plug-ins
template <class Derived>
struct transfer_function_plugin_creator_concrete : public transfer_function_plugin_creator
{
	//! register the plug-in by its name
	transfer_function_plugin_creator_concrete(const std::string &plugin_name)
	{
		get_transfer_function_plugin_map()[plugin_name] = this;
	}

	//! create an instance of the plug-in
	std::unique_ptr<transfer_function_plugin> create(config_file &cf, const cosmology::parameters& cp) const
	{
		return std::make_unique<Derived>(cf,cp);
	}
};

typedef transfer_function_plugin transfer_function;

std::unique_ptr<transfer_function_plugin> select_transfer_function_plugin(config_file &cf, const cosmology::parameters &cp);

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

//! k-space transfer function
class TransferFunction_k
{
public:
	static transfer_function *ptf_;
	static real_t nspec_;
	double pnorm_, sqrtpnorm_;
	static tf_type type_;

	TransferFunction_k(tf_type type, transfer_function *tf, real_t nspec, real_t pnorm)
			: pnorm_(pnorm)
	{
		ptf_ = tf;
		nspec_ = nspec;
		sqrtpnorm_ = sqrt(pnorm_);
		type_ = type;

		std::string fname("input_powerspec.txt");
		if (type == delta_cdm || type == delta_matter)
		{
			std::ofstream ofs(fname.c_str());
			double kmin = log10(tf->get_kmin()), kmax = log10(tf->get_kmax());
			double dk = (kmax - kmin) / 300.;

			ofs << "# The power spectrum definition is smaller than CAMB by a factor 8 pi^3."
					<< std::endl;

			if (tf->tf_is_distinct())
			{
				ofs << "#"
						<< std::setw(15) << "k [h/Mpc]"
						<< std::setw(16) << "P_cdm"
						<< std::setw(16) << "P_theta_cdm"
						<< std::setw(16) << "P_bar"
						<< std::setw(16) << "P_vbar"
						<< std::setw(16) << "P_total"
						<< std::setw(16) << "P_vtotal"
						<< std::endl;

				for (int i = 0; i < 300; ++i)
				{
					double k = pow(10.0, kmin + i * dk);
					ofs << std::setw(16) << k
							<< std::setw(16) << pow(sqrtpnorm_ * pow(k, 0.5 * nspec_) * ptf_->compute(k, delta_cdm), 2)
							<< std::setw(16) << pow(sqrtpnorm_ * pow(k, 0.5 * nspec_) * ptf_->compute(k, theta_cdm), 2)
							<< std::setw(16) << pow(sqrtpnorm_ * pow(k, 0.5 * nspec_) * ptf_->compute(k, delta_baryon), 2)
							<< std::setw(16) << pow(sqrtpnorm_ * pow(k, 0.5 * nspec_) * ptf_->compute(k, theta_baryon), 2)
							<< std::setw(16) << pow(sqrtpnorm_ * pow(k, 0.5 * nspec_) * ptf_->compute(k, delta_matter), 2)
							<< std::setw(16) << pow(sqrtpnorm_ * pow(k, 0.5 * nspec_) * ptf_->compute(k, theta_matter), 2)
							<< std::endl;
				}
			}
			else
			{
				ofs << "#"
						<< std::setw(16) << "k [h/Mpc]"
						<< std::setw(16) << "P_cdm"
						<< std::setw(16) << "P_theta_cdm"
						<< std::setw(16) << "P_total"
						<< std::endl;

				for (int i = 0; i < 300; ++i)
				{
					double k = pow(10.0, kmin + i * dk);
					ofs << std::setw(16) << k
							<< std::setw(16) << pow(sqrtpnorm_ * pow(k, 0.5 * nspec_) * ptf_->compute(k, delta_cdm), 2)
							<< std::setw(16) << pow(sqrtpnorm_ * pow(k, 0.5 * nspec_) * ptf_->compute(k, theta_cdm), 2)
							<< std::setw(16) << pow(sqrtpnorm_ * pow(k, 0.5 * nspec_) * ptf_->compute(k, delta_matter), 2)
							<< std::endl;
				}
			}
			if (ofs.fail()) {
				std::string message = "Could not write to file: " + fname;
				music::elog << message << std::endl;
				throw std::runtime_error(message);
			}
			ofs.close();
		}
	}

	inline real_t compute(real_t k) const
	{
		return sqrtpnorm_ * pow(k, 0.5 * nspec_) * ptf_->compute(k, type_);
	}
};

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
