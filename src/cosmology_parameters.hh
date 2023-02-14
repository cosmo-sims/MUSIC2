// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2020 by Oliver Hahn
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

#include <map>
#include <string>
#include <cmath>
#include <physical_constants.hh>
#include <config_file.hh>
#include <general.hh>

namespace cosmology_new
{
    //! structure for cosmological parameters
    class parameters
    {
    public:
        using defaultmmap_t = std::map<std::string,std::map<std::string,real_t>>;

    private:
        std::map<std::string, double> pmap_;  //!< All parameters are stored here as key-value pairs

        static defaultmmap_t default_pmaps_;  //!< Holds pre-defined parameter sets, see src/cosmology_parameters.cc 

    public:
        //! get routine for cosmological parameter key-value pairs
        double get(const std::string &key) const
        {
            auto it = pmap_.find(key);
            if (it == pmap_.end())
            {
                auto errmsg = std::string("Cosmological parameter \'") + key + std::string("\' does not exist in internal list.");
                music::elog << errmsg << std::endl;
                throw std::runtime_error(errmsg.c_str());
            }
            return it->second;
        }

        //! set routine for cosmological parameter key-value pairs
        void set(const std::string &key, const double value)
        {
            auto it = pmap_.find(key);
            if (it != pmap_.end())
            {
                pmap_[key] = value;
            }
            else
            {
                auto errmsg = std::string("Cosmological parameter \'") + key + std::string("\' does not exist in internal list. Needs to be defaulted before it can be set!");
                music::elog << errmsg << std::endl;
                throw std::runtime_error(errmsg.c_str());
            }
        }

        //! shortcut get routine for cosmological parameter key-value pairs through bracket operator
        inline double operator[](const std::string &key) const { return this->get(key); }

        //! default constructor does nothing
        parameters() {}

        //! default copy constructor
        parameters(const parameters &) = default;

        //! main constructor for explicit construction from input config file
        explicit parameters( config_file &cf )
        {
            music::ilog << "-------------------------------------------------------------------------------" << std::endl;
            
            if( cf.get_value_safe<std::string>("cosmology","ParameterSet","none") == std::string("none"))
            {
                //-------------------------------------------------------------------------------
                // read parameters from config file
                //-------------------------------------------------------------------------------

                auto defaultp = default_pmaps_.find("Planck2018EE+BAO+SN")->second;
            
                // CMB
                pmap_["Tcmb"] = cf.get_value_safe<double>("cosmology", "Tcmb", defaultp["Tcmb"]);
                pmap_["YHe"] = cf.get_value_safe<double>("cosmology", "YHe", defaultp["YHe"]);

                // H0
                pmap_["h"] = cf.get_value_safe<double>("cosmology", "H0", defaultp["h"]*100.0) / 100.0;

                // primordial and normalisation
                if(!cf.contains_key("cosmology/n_s"))
                    pmap_["n_s"] = cf.get_value_safe<double>("cosmology", "nspec", defaultp["n_s"]);
                else
                    pmap_["n_s"] = cf.get_value_safe<double>("cosmology", "n_s", defaultp["n_s"]);

                pmap_["A_s"] = cf.get_value_safe<double>("cosmology", "A_s", cf.contains_key("cosmology/sigma_8")? -1.0 : defaultp["A_s"]);
                pmap_["k_p"] = cf.get_value_safe<double>("cosmology", "k_p", defaultp["k_p"]);

                pmap_["sigma_8"] = cf.get_value_safe<double>("cosmology", "sigma_8", -1.0);
                
                // baryon and non-relativistic matter content
                pmap_["Omega_b"] = cf.get_value_safe<double>("cosmology", "Omega_b", defaultp["Omega_b"]);
                pmap_["Omega_m"] = cf.get_value_safe<double>("cosmology", "Omega_m", defaultp["Omega_m"]);

                // massive neutrino species
                pmap_["m_nu1"] = cf.get_value_safe<double>("cosmology", "m_nu1", defaultp["m_nu1"]);
                pmap_["m_nu2"] = cf.get_value_safe<double>("cosmology", "m_nu2", defaultp["m_nu2"]);
                pmap_["m_nu3"] = cf.get_value_safe<double>("cosmology", "m_nu3", defaultp["m_nu3"]);
                int N_nu_massive = int(this->get("m_nu1") > 1e-9) + int(this->get("m_nu2") > 1e-9) + int(this->get("m_nu3") > 1e-9);;
                
                // number ultrarelativistic neutrinos
                pmap_["N_ur"] = cf.get_value_safe<double>("cosmology", "N_ur", 3.046 - N_nu_massive);
            
                // dark energy
                pmap_["Omega_DE"] = cf.get_value_safe<double>("cosmology", "Omega_L", defaultp["Omega_DE"]);
                pmap_["w_0"] = cf.get_value_safe<double>("cosmology", "w_0", defaultp["w_0"]);
                pmap_["w_a"] = cf.get_value_safe<double>("cosmology", "w_a", defaultp["w_a"]);

            }else{

                //-------------------------------------------------------------------------------
                // load predefined parameter set
                //-------------------------------------------------------------------------------
                auto pset_name = cf.get_value<std::string>("cosmology","ParameterSet");
                auto it = default_pmaps_.find( pset_name );
                if( it == default_pmaps_.end() ){
                    music::elog << "Unknown cosmology parameter set \'" << pset_name << "\'!" << std::endl;
                    music::ilog << "Valid pre-defined sets are: " << std::endl;
                    for( auto ii : default_pmaps_ ){
                        music::ilog << "  " << ii.first << std::endl;
                    }
                    throw std::runtime_error("Invalid value for cosmology/ParameterSet");
                }else{
                    music::ilog << "Loading cosmological parameter set \'" << it->first << "\'..." << std::endl;
                }

                // copy pre-defined set in
                for( auto entry : it->second ){
                    pmap_[entry.first] = entry.second;
                }
            }

            //-------------------------------------------------------------------------------
            // derived parameters
            //-------------------------------------------------------------------------------
            double h = this->get("h");
            pmap_["H0"] = 100.0 * h;

            // massive neutrinos
            pmap_["N_nu_massive"] = int(this->get("m_nu1") > 1e-9) + int(this->get("m_nu2") > 1e-9) + int(this->get("m_nu3") > 1e-9);
            const double sum_m_nu = this->get("m_nu1") + this->get("m_nu2") + this->get("m_nu3");
            pmap_["Omega_nu_massive"] = sum_m_nu / (93.14 * h * h); // Omega_nu_m = \sum_i m_i / (93.14 eV h^2)

            // calculate energy density in ultrarelativistic species from Tcmb and Neff
            // photons
            pmap_["Omega_gamma"] = 4 * phys_const::sigma_SI / std::pow(phys_const::c_SI, 3) * std::pow(this->get("Tcmb"), 4.0) 
                    / phys_const::rhocrit_h2_SI / (this->get("h") * this->get("h"));
            // massless neutrinos
            pmap_["Omega_nu_massless"] = this->get("N_ur") * this->get("Omega_gamma") * 7. / 8. * std::pow(4. / 11., 4. / 3.);
            // total relativistic
            pmap_["Omega_r"] = this->get("Omega_gamma") + this->get("Omega_nu_massless");

            // compute amount of cold dark matter as the rest
            pmap_["Omega_c"] = this->get("Omega_m") - this->get("Omega_b") - this->get("Omega_nu_massive");

            if (cf.get_value_safe<bool>("cosmology", "ZeroRadiation", false))
            {
                pmap_["Omega_r"] = 0.0;
            }

            pmap_["f_b"] = this->get("Omega_b") / this->get("Omega_m");
            pmap_["f_c"] = 1.0 - this->get("f_b"); // this means we add massive neutrinos to CDM here

#if 1
            // assume zero curvature, take difference from dark energy
            pmap_["Omega_DE"] += 1.0 - this->get("Omega_m") - this->get("Omega_DE") - this->get("Omega_r");
            // Omega_DE += 1.0 - Omega_m - Omega_DE - Omega_r;
            pmap_["Omega_k"] = 0.0;
#else
            // allow for curvature
            Omega_k = 1.0 - Omega_m - Omega_DE - Omega_r;
#endif

            pmap_["dplus"] = 0.0;
            pmap_["pnorm"] = 0.0;
            pmap_["sqrtpnorm"] = 0.0;
            pmap_["vfact"] = 0.0;

            music::ilog << "Cosmological parameters are: " << std::endl;
            music::ilog << " h        = " << std::setw(16) << this->get("h");
            if( this->get("A_s") > 0.0 )
              music::ilog << "A_s      = " << std::setw(16) << this->get("A_s");
            else
              music::ilog << "sigma_8  = " << std::setw(16) << this->get("sigma_8");
            music::ilog << "n_s     = " << std::setw(16) << this->get("n_s") << std::endl;
            music::ilog << " Omega_c  = " << std::setw(16) << this->get("Omega_c")  << "Omega_b  = " << std::setw(16) << this->get("Omega_b") << "Omega_m = " << std::setw(16) << this->get("Omega_m") << std::endl;
            music::ilog << " Omega_r  = " << std::setw(16) << this->get("Omega_r")  << "Omega_nu = " << std::setw(16) << this->get("Omega_nu_massive") << "∑m_nu   = " << sum_m_nu << "eV" << std::endl;
            music::ilog << " Omega_DE = " << std::setw(16) << this->get("Omega_DE") << "w_0      = " << std::setw(16) << this->get("w_0")      << "w_a     = " << std::setw(16) << this->get("w_a") << std::endl;
            //music::ilog << " Omega_k  = " << 1.0 - this->get("Omega_m") - this->get("Omega_r") - this->get("Omega_DE") << std::endl;
            if (this->get("Omega_r") > 0.0)
            {
                music::wlog << " Radiation enabled, using Omega_r=" << this->get("Omega_r") << " internally for backscaling." << std::endl;
                music::wlog << " Make sure your sim code supports this, otherwise set [cosmology] / ZeroRadiation=true." << std::endl;
            }
        }

        //! print the names of all pre-defined parameter sets
        void print_available_sets( void ) const
        {
            for( auto ii : default_pmaps_ ){
                music::ilog << "\t\'" << ii.first << "\'" << std::endl;
            }
        }
    };

    void print_ParameterSets( void );
} // namespace cosmology