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

#include <array>
#include "math/vec.hh"

#include <cosmology_parameters.hh>
#include <physical_constants.hh>
#include <transfer_function.hh>
#include <math/ode_integrate.hh>
#include <logger.hh>

#include <math/interpolate.hh>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

namespace cosmology
{

/*!
 * @class cosmology::calculator
 * @brief provides functions to compute cosmological quantities
 *
 * This class provides member functions to compute cosmological quantities
 * related to the Friedmann equations and linear perturbation theory, it also
 * provides the functionality to work with back-scaled cosmological fields
 */
class calculator
{
public:
    //! data structure to store cosmological parameters
    cosmology::parameters cosmo_param_;

    //! pointer to an instance of a transfer function plugin
    std::unique_ptr<transfer_function_plugin> transfer_function_;

private:
    static constexpr double REL_PRECISION = 1e-10;
    interpolated_function_1d<true,true,false> D_of_a_, f_of_a_, a_of_D_;
    double Dnow_, Dplus_start_, Dplus_target_, astart_, atarget_;

    double m_n_s_, m_sqrtpnorm_;

    //! wrapper for GSL adaptive integration routine, do not use if many integrations need to be done as it allocates and deallocates memory
    //! set to 61-point Gauss-Kronrod and large workspace, used for sigma_8 normalisation
    real_t integrate(double (*func)(double x, void *params), double a, double b, void *params) const
    {
        constexpr size_t wspace_size{100000};

        double result{0.0};
        double error{0.0};

        gsl_function F;
        F.function = func;
        F.params = params;
        
        auto errh = gsl_set_error_handler_off();
        gsl_integration_workspace *wspace = gsl_integration_workspace_alloc(wspace_size);
        gsl_integration_qag(&F, a, b, 0, REL_PRECISION, wspace_size, GSL_INTEG_GAUSS61, wspace, &result, &error);
        gsl_integration_workspace_free(wspace);
        gsl_set_error_handler(errh);

        if (error / result > REL_PRECISION)
            music::wlog << "no convergence in function 'integrate', rel. error=" << error / result << std::endl;

        return static_cast<real_t>(result);
    }

    //! compute the linear theory growth factor D+ by solving the single fluid ODE, returns tables D(a), f(a)
    /*!
	 * @param tab_a reference to STL vector for values of a at which table entries exist
     * @param tab_D reference to STL vector for values D(a) with a from tab_a
     * @param tab_f reference to STL vector for values f(a)=dlog D / dlog a with a from tab_a
     */
    void compute_growth( std::vector<double>& tab_a, std::vector<double>& tab_D, std::vector<double>& tab_f )
    {
        using v_t = vec_t<3,double>;

        // set ICs, very deep in radiation domination
        const double a0 = 1e-10;
        const double D0 = a0;
        const double Dprime0 = 2.0 * D0 * H_of_a(a0) / std::pow(phys_const::c_SI, 2);
        const double t0 = 1.0 / (a0 * H_of_a(a0));

        v_t y0({a0, D0, Dprime0});

        // set up integration
        double dt = 1e-9;
        double dtdid, dtnext;
        const double amax = 2.0;

        v_t yy(y0);
        double t = t0;
        const double eps = 1e-10;

        const double Omega_m( cosmo_param_["Omega_m"] ), H0( cosmo_param_["H0"] );

        while (yy[0] < amax)
        {
            // RHS of ODEs
            auto rhs = [&](double t, v_t y) -> v_t {
                auto a = y[0];
                auto D = y[1];
                auto Dprime = y[2];
                v_t dy;
                // da/dtau = a^2 H(a)
                dy[0] = a * a * H_of_a(a);
                // d D/dtau
                dy[1] = Dprime;
                // d^2 D / dtau^2
                dy[2] = -a * H_of_a(a) * Dprime + 3.0 / 2.0 * Omega_m * std::pow(H0, 2) * D / a;
                return dy;
            };

            // scale by predicted value to get approx. constant fractional errors
            v_t yyscale = yy.abs() + dt * rhs(t, yy).abs();
            
            // call integrator
            ode_integrate::rk_step_qs(dt, t, yy, yyscale, rhs, eps, dtdid, dtnext);

            tab_a.push_back(yy[0]);
            tab_D.push_back(yy[1]);
            tab_f.push_back(yy[2]); // temporarily store D' in table

            dt = dtnext;
        }

        // compute f, before we stored here D'
        for (size_t i = 0; i < tab_a.size(); ++i)
        {
            tab_f[i] = tab_f[i] / (tab_a[i] * H_of_a(tab_a[i]) * tab_D[i]);
            tab_D[i] = tab_D[i];
            tab_a[i] = tab_a[i];
        }
    }

public:
    
    calculator() = delete;
    
    calculator(const calculator& c) = delete;

    //! constructor for a cosmology calculator object
    /*!
	 * @param acosmo a cosmological parameters structure
	 * @param pTransferFunction pointer to an instance of a transfer function object
	 */
    explicit calculator(config_file &cf)
        : cosmo_param_(cf), astart_( 1.0/(1.0+cf.get_value<double>("setup","zstart")) ),
            atarget_( 1.0/(1.0+cf.get_value_safe<double>("cosmology","ztarget",0.0)) )
    {
        // pre-compute growth factors and store for interpolation
        std::vector<double> tab_a, tab_D, tab_f;
        this->compute_growth(tab_a, tab_D, tab_f);
        D_of_a_.set_data(tab_a,tab_D);
        f_of_a_.set_data(tab_a,tab_f);
        a_of_D_.set_data(tab_D,tab_a);
        Dnow_ = D_of_a_(1.0);

        Dplus_start_  = D_of_a_( astart_ ) / Dnow_;
        Dplus_target_ = D_of_a_( atarget_ ) / Dnow_;

        music::ilog << "Linear growth factors: D+_target = " << Dplus_target_ << ", D+_start = " << Dplus_start_ << std::endl;
        music::ilog << "-------------------------------------------------------------------------------" << std::endl;

        // set up transfer functions and compute normalisation
        transfer_function_ = select_transfer_function_plugin(cf, cosmo_param_);
        transfer_function_->intialise();
        if( !transfer_function_->tf_isnormalised_ ){
            cosmo_param_.set("pnorm", this->compute_pnorm_from_sigma8()*Dplus_start_*Dplus_start_ );
        }else{
            // WARNING: we do explicit back-scaling here, which is notably different from monofonIC
            cosmo_param_.set("pnorm", Dplus_start_*Dplus_start_/Dplus_target_/Dplus_target_);
            auto sigma8 = this->compute_sigma8()*Dplus_start_/Dplus_target_;
            music::ilog << "Measured sigma_8 for given PS normalisation is " <<  sigma8 << std::endl;
        }
        cosmo_param_.set("sqrtpnorm", std::sqrt(cosmo_param_["pnorm"]));

        // if (!transfer_function_->tf_is_distinct())
        //     music::wlog << " - WARNING: The selected transfer function does not support" << std::endl
        //                 << "            distinct amplitudes for baryon and DM fields!" << std::endl
        //                 << "            Perturbation amplitudes will be identical!" << std::endl;

        music::ilog << std::setw(32) << std::left << " . TF supports distinct CDM+baryons"
                    << " : " << (transfer_function_->tf_is_distinct() ? "yes" : "no") << std::endl;
        music::ilog << std::setw(32) << std::left << " . TF maximum wave number"
                    << " : " << transfer_function_->get_kmax() << " h/Mpc" << std::endl;

        m_n_s_ = cosmo_param_["n_s"];
        m_sqrtpnorm_ = cosmo_param_["sqrtpnorm"];
    }

    ~calculator() { }

    //! Write out a correctly scaled power spectrum at time a
    void write_powerspectrum(real_t a, std::string fname) const
    {
        // const real_t Dplus0 = this->get_growth_factor(a);

        if (CONFIG::MPI_task_rank == 0)
        {
            double kmin = std::max(1e-4, transfer_function_->get_kmin());

            // write power spectrum to a file
            std::ofstream ofs(fname.c_str());
            std::stringstream ss;
            ss << " ,ap=" << a << "";
            ofs << "# " << std::setw(18) << "k [h/Mpc]"
                << std::setw(20) << ("P_dtot(k,a=ap)")
                << std::setw(20) << ("P_dcdm(k,a=ap)")
                << std::setw(20) << ("P_dbar(k,a=ap)")
                << std::setw(20) << ("P_tcdm(k,a=ap)")
                << std::setw(20) << ("P_tbar(k,a=ap)")
                << std::setw(20) << ("P_dtot(k,a=1)")
                << std::setw(20) << ("P_dcdm(k,a=1)")
                << std::setw(20) << ("P_dbar(k,a=1)")
                << std::setw(20) << ("P_tcdm(k,a=1)")
                << std::setw(20) << ("P_tbar(k,a=1)")
                << std::endl;
            for (double k = kmin; k < transfer_function_->get_kmax(); k *= 1.01)
            {
                ofs << std::setw(20) << std::setprecision(10) << k
                    << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_matter)*Dplus_start_, 2.0)
                    << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_cdm)*Dplus_start_, 2.0)
                    << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_baryon)*Dplus_start_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_matter)*Dplus_start_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_cdm)*Dplus_start_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_baryon)*Dplus_start_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, theta_cdm)*Dplus_start_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, theta_baryon)*Dplus_start_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_matter0)* Dplus_start_ / Dplus_target_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_cdm0)* Dplus_start_ / Dplus_target_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, delta_baryon0)* Dplus_start_ / Dplus_target_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, theta_cdm0)* Dplus_start_ / Dplus_target_, 2.0)
                    // << std::setw(20) << std::setprecision(10) << std::pow(this->get_amplitude(k, theta_baryon0)* Dplus_start_ / Dplus_target_, 2.0)
                    << std::endl;
            }
        }
        music::ilog << "Wrote power spectrum at a=" << a << " to file \'" << fname << "\'" << std::endl;
    }

    //! Write out a correctly scaled power spectrum at starting time
    void write_transfer( std::string fname ) const
    {
        // const real_t Dplus0 = this->get_growth_factor(a);

        if (CONFIG::MPI_task_rank == 0)
        {
            double kmin = std::max(1e-4, transfer_function_->get_kmin());

            // write power spectrum to a file
            std::ofstream ofs(fname.c_str());
            std::stringstream ss;
            ss << " ,ap=" << astart_ << "";
            ofs << "# " << std::setw(18) << "k [h/Mpc]"
                << std::setw(20) << ("delta_c(k,a=ap)")
                << std::setw(20) << ("delta_b(k,a=ap)")
                << std::setw(20) << ("delta_m(k,a=ap)")
                << std::setw(20) << ("delta_bc(k,a=ap)")
                << std::endl;
            double fb = cosmo_param_["f_b"], fc = cosmo_param_["f_c"];
            for (double k = kmin; k < transfer_function_->get_kmax(); k *= 1.01)
            {
                const double dm  = this->get_amplitude(k, delta_matter) * Dplus_start_ / Dplus_target_;
                const double dbc = this->get_amplitude(k, delta_bc);
                const double db  = dm + fc * dbc;
                const double dc  = dm - fb * dbc;
                const double tm  = this->get_amplitude(k, delta_matter) * Dplus_start_ / Dplus_target_;
                const double tbc = this->get_amplitude(k, theta_bc);
                const double tb  = dm + fc * dbc;
                const double tc  = dm - fb * dbc;
                
                ofs << std::setw(20) << std::setprecision(10) << k
                    << std::setw(20) << std::setprecision(10) << dc
                    << std::setw(20) << std::setprecision(10) << db
                    << std::setw(20) << std::setprecision(10) << dm
                    << std::setw(20) << std::setprecision(10) << dbc + 2 * tbc * (std::sqrt( Dplus_target_ / Dplus_start_ ) - 1.0)
                    << std::setw(20) << std::setprecision(10) << tc / std::pow( Dplus_start_ / Dplus_target_, 0.5 )
                    << std::setw(20) << std::setprecision(10) << tb / std::pow( Dplus_start_ / Dplus_target_, 0.5 )
                    << std::setw(20) << std::setprecision(10) << tm / std::pow( Dplus_start_ / Dplus_target_, 0.5 )
                    << std::setw(20) << std::setprecision(10) << tbc / std::pow( Dplus_start_ / Dplus_target_, 0.5 )
                    << std::endl;
            }
        }
        music::ilog << "Wrote input transfer functions at a=" << astart_ << " to file \'" << fname << "\'" << std::endl;
    }

    const cosmology::parameters &get_parameters(void) const noexcept
    {
        return cosmo_param_;
    }

    //! return the value of the Hubble function H(a) = dloga/dt 
    inline double H_of_a(double a) const noexcept
    {
        double HH2 = 0.0;
        HH2 += cosmo_param_["Omega_r"] / (a * a * a * a);
        HH2 += cosmo_param_["Omega_m"] / (a * a * a);
        HH2 += cosmo_param_["Omega_k"] / (a * a);
        HH2 += cosmo_param_["Omega_DE"] * std::pow(a, -3. * (1. + cosmo_param_["w_0"] + cosmo_param_["w_a"])) * exp(-3. * (1.0 - a) * cosmo_param_["w_a"]);
        return cosmo_param_["H0"] * std::sqrt(HH2);
    }

    //! Computes the linear theory growth factor D+, normalised to D+(a=1)=1
    real_t get_growth_factor(real_t a) const noexcept
    {
        return D_of_a_(a) / Dnow_;
    }

    //! Computes the inverse of get_growth_factor
    real_t get_a( real_t Dplus ) const noexcept
    {
        return a_of_D_( Dplus * Dnow_ );
    }

    //! Computes the linear theory growth rate f
    /*! Function computes (by interpolating on precalculated table)
     *   f = dlog D+ / dlog a
     */
    real_t get_f(real_t a) const noexcept
    {
        return f_of_a_(a);
    }

    //! Compute the factor relating particle displacement and velocity
    /*! Function computes
     *  vfac = a * (H(a)/h) * dlogD+ / dlog a 
     */
    real_t get_vfact(real_t a) const noexcept
    {
        return f_of_a_(a) * a * H_of_a(a) / cosmo_param_["h"];
    }

    //! Integrand for the sigma_8 normalization of the power spectrum
    /*! Returns the value of the primordial power spectrum multiplied with 
     the transfer function and the window function of 8 Mpc/h at wave number k */
    static double dSigma8(double k, void *pParams)
    {
        cosmology::calculator *pcc = reinterpret_cast<cosmology::calculator *>(pParams);

        const double x = k * 8.0;
        const double w = (x < 0.001)? 1.0-0.1*x*x : 3.0 * (std::sin(x) - x * std::cos(x)) / (x * x * x);
            
        static double nspect = (double)pcc->cosmo_param_["n_s"];
        double tf = pcc->transfer_function_->compute(k, delta_matter);

        //... no growth factor since we compute at z=0 and normalize so that D+(z=0)=1
        return k * k * w * w * pow((double)k, (double)nspect) * tf * tf;
    }

    //! Integrand for the sigma_8 normalization of the power spectrum
    /*! Returns the value of the primordial power spectrum multiplied with 
	 the transfer function and the window function of 8 Mpc/h at wave number k */
    static double dSigma8_0(double k, void *pParams)
    {
        cosmology::calculator *pcc = reinterpret_cast<cosmology::calculator *>(pParams);

        const double x = k * 8.0;
        const double w = (x < 0.001)? 1.0-0.1*x*x : 3.0 * (std::sin(x) - x * std::cos(x)) / (x * x * x);

        static double nspect = static_cast<double>(pcc->cosmo_param_["n_s"]);
        double tf = pcc->transfer_function_->compute(k, delta_matter0);

        //... no growth factor since we compute at z=0 and normalize so that D+(z=0)=1
        return k * k * w * w * std::pow(k, nspect) * tf * tf;
    }

    //! Computes the amplitude of a mode from the power spectrum
    /*! Function evaluates the supplied transfer function transfer_function_
	 * and returns the amplitude of fluctuations at wave number k (in h/Mpc) back-scaled to z=z_start
	 * @param k wave number at which to evaluate
     * @param type one of the species: {delta,theta}_{matter,cdm,baryon,neutrino}
	 */
    inline real_t get_amplitude( const real_t k, const tf_type type) const
    {
        return std::pow(k, 0.5 * m_n_s_) * transfer_function_->compute(k, type) * m_sqrtpnorm_;
    }

    //! Compute amplitude of the back-scaled delta_bc mode, with decaying velocity v_bc included or not (in which case delta_bc=const)
    inline real_t get_amplitude_delta_bc( const real_t k, bool withvbc ) const
    {
        const real_t Dratio = Dplus_target_ / Dplus_start_;
        const real_t dbc = transfer_function_->compute(k, delta_bc) + (withvbc? 2 * transfer_function_->compute(k, theta_bc) * (std::sqrt(Dratio) - 1.0) : 0.0);
        // need to multiply with Dplus_target since sqrtpnorm rescales like that
        return std::pow(k, 0.5 * m_n_s_) * dbc * (m_sqrtpnorm_ * Dplus_target_);
    }

    //! Compute amplitude of the back-scaled relative velocity theta_bc mode if withvbc==true, otherwise return zero
    inline real_t get_amplitude_theta_bc( const real_t k, bool withvbc ) const
    {
        const real_t Dratio = Dplus_target_ / Dplus_start_;
        const real_t tbc = transfer_function_->compute(k, theta_bc) * std::sqrt(Dratio);
        // need to multiply with Dplus_target since sqrtpnorm rescales like that
        return withvbc ? std::pow(k, 0.5 * m_n_s_) * tbc * (m_sqrtpnorm_ * Dplus_target_) : 0.0;
    }


    //! Computes the normalization for the power spectrum
    /*!
	 * integrates the power spectrum to fix the normalization to that given
	 * by the sigma_8 parameter
	 */
    real_t compute_sigma8(void)
    {
        real_t sigma0, kmin, kmax;
        kmax = transfer_function_->get_kmax();
        kmin = transfer_function_->get_kmin();

        if (!transfer_function_->tf_has_total0())
            sigma0 = 4.0 * M_PI * integrate(&dSigma8, static_cast<double>(kmin), static_cast<double>(kmax), this);
        else{
            sigma0 = 4.0 * M_PI * integrate(&dSigma8_0, static_cast<double>(kmin), static_cast<double>(kmax), this);
        }

        return std::sqrt(sigma0);
    }

    //! Computes the normalization for the power spectrum
    /*!
	 * integrates the power spectrum to fix the normalization to that given
	 * by the sigma_8 parameter
	 */
    real_t compute_pnorm_from_sigma8(void)
    {
        auto measured_sigma8 = this->compute_sigma8();
        return cosmo_param_["sigma_8"] * cosmo_param_["sigma_8"] / (measured_sigma8  * measured_sigma8);
    }
};

//! compute the jeans sound speed
/*! given a density in g/cm^-3 and a mass in g it gives back the sound
 *  speed in cm/s for which the input mass is equal to the jeans mass
 *  @param rho density 
 *  @param mass mass scale
 *  @returns jeans sound speed
 */
// inline double jeans_sound_speed(double rho, double mass)
// {
//     const double G = 6.67e-8;
//     return pow(6.0 * mass / M_PI * std::sqrt(rho) * std::pow(G, 1.5), 1.0 / 3.0);
// }

} // namespace cosmology
