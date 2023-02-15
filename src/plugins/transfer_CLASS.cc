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

#ifdef USE_CLASS

#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <sstream>

#include <ClassEngine.hh>

#include <general.hh>
#include <config_file.hh>
#include <transfer_function.hh>
// #include <ic_generator.hh>

#include <math/interpolate.hh>

class transfer_CLASS_plugin : public transfer_function_plugin
{
private:
  
  using transfer_function_plugin::cosmo_params_;

  interpolated_function_1d<true, true, false> delta_c_, delta_b_, delta_n_, delta_m_, theta_c_, theta_b_, theta_n_, theta_m_;
  interpolated_function_1d<true, true, false> delta_c0_, delta_b0_, delta_n0_, delta_m0_, theta_c0_, theta_b0_, theta_n0_, theta_m0_;

  double zstart_, ztarget_, astart_, atarget_, kmax_, kmin_, h_, tnorm_;

  ClassParams pars_;
  std::unique_ptr<ClassEngine> the_ClassEngine_;
  std::ofstream ofs_class_input_;

  template <typename T>
  void add_class_parameter(std::string parameter_name, const T parameter_value)
  {
    pars_.add(parameter_name, parameter_value);
    ofs_class_input_ << parameter_name << " = " << parameter_value << std::endl;
  }

  //! Set up class parameters from MUSIC cosmological parameters
  void init_ClassEngine(void)
  {
    //--- general parameters ------------------------------------------
    add_class_parameter("z_max_pk", std::max(std::max(zstart_, ztarget_),199.0)); // use 1.2 as safety
    add_class_parameter("P_k_max_h/Mpc", std::max(2.0,kmax_));
    add_class_parameter("output", "dTk,vTk");
    add_class_parameter("extra metric transfer functions","yes");
    // add_class_parameter("lensing", "no");

    //--- choose gauge ------------------------------------------------
    // add_class_parameter("extra metric transfer functions", "yes");
    add_class_parameter("gauge", "synchronous");

    //--- cosmological parameters, densities --------------------------
    add_class_parameter("h", cosmo_params_.get("h"));

    add_class_parameter("Omega_b", cosmo_params_.get("Omega_b"));
    add_class_parameter("Omega_cdm", cosmo_params_.get("Omega_c"));
    add_class_parameter("Omega_k", cosmo_params_.get("Omega_k"));
    add_class_parameter("Omega_fld", 0.0);
    add_class_parameter("Omega_scf", 0.0);


    // add_class_parameter("fluid_equation_of_state","CLP");
    // add_class_parameter("w0_fld", -1 );
    // add_class_parameter("wa_fld", 0. );
    // add_class_parameter("cs2_fld", 1);

    //--- massive neutrinos -------------------------------------------
#if 0
    //default off
    // add_class_parameter("Omega_ur",0.0);
    add_class_parameter("N_ur", cosmo_params_.get("N_ur"));
    add_class_parameter("N_ncdm", 0);

#else
    
    add_class_parameter("N_ur", cosmo_params_.get("N_ur"));
    add_class_parameter("N_ncdm", cosmo_params_.get("N_nu_massive"));
    if( cosmo_params_.get("N_nu_massive") > 0 ){
      std::stringstream sstr;
      if( cosmo_params_.get("m_nu1") > 1e-9 ) sstr << cosmo_params_.get("m_nu1");
      if( cosmo_params_.get("m_nu2") > 1e-9 ) sstr << ", " << cosmo_params_.get("m_nu2");
      if( cosmo_params_.get("m_nu3") > 1e-9 ) sstr << ", " << cosmo_params_.get("m_nu3");
      add_class_parameter("m_ncdm", sstr.str().c_str());
    }
    
    // change above to enable
    //add_class_parameter("omega_ncdm", 0.0006451439);
    //add_class_parameter("m_ncdm", "0.4");
    //add_class_parameter("T_ncdm", 0.71611);
#endif

    //--- cosmological parameters, primordial -------------------------
    add_class_parameter("P_k_ini type", "analytic_Pk");

    if( cosmo_params_.get("A_s") > 0.0 ){
      add_class_parameter("A_s", cosmo_params_.get("A_s"));
    }else{
      add_class_parameter("sigma8", cosmo_params_.get("sigma_8"));
    }
    add_class_parameter("n_s", cosmo_params_.get("n_s"));
    add_class_parameter("alpha_s", 0.0);
    add_class_parameter("T_cmb", cosmo_params_.get("Tcmb"));
    add_class_parameter("YHe", cosmo_params_.get("YHe"));

    // additional parameters
    add_class_parameter("reio_parametrization", "reio_none");

    // precision parameters
    add_class_parameter("k_per_decade_for_pk", 100);
    add_class_parameter("k_per_decade_for_bao", 100);
    add_class_parameter("compute damping scale", "yes");
    add_class_parameter("tol_perturb_integration", 1.e-8);
    add_class_parameter("tol_background_integration", 1e-9);

    // high precision options from cl_permille.pre:
    // precision file to be passed as input in order to achieve at least percent precision on scalar Cls
    add_class_parameter("hyper_flat_approximation_nu", 7000.);
    add_class_parameter("transfer_neglect_delta_k_S_t0", 0.17);
    add_class_parameter("transfer_neglect_delta_k_S_t1", 0.05);
    add_class_parameter("transfer_neglect_delta_k_S_t2", 0.17);
    add_class_parameter("transfer_neglect_delta_k_S_e", 0.13);
    add_class_parameter("delta_l_max", 1000);

    int class_verbosity = 0;

    add_class_parameter("background_verbose", class_verbosity);
    add_class_parameter("thermodynamics_verbose", class_verbosity);
    add_class_parameter("perturbations_verbose", class_verbosity);
    add_class_parameter("transfer_verbose", class_verbosity);
    add_class_parameter("primordial_verbose", class_verbosity);
    add_class_parameter("spectra_verbose", class_verbosity);
    add_class_parameter("nonlinear_verbose", class_verbosity);
    add_class_parameter("lensing_verbose", class_verbosity);
    add_class_parameter("output_verbose", class_verbosity);

    // output parameters, only needed for the control CLASS .ini file that we output
    std::stringstream zlist;
    if (ztarget_ == zstart_)
      zlist << ztarget_ << ((ztarget_!=0.0)? ", 0.0" : "");
    else
      zlist << std::max(ztarget_, zstart_) << ", " << std::min(ztarget_, zstart_) << ", 0.0";
    add_class_parameter("z_pk", zlist.str());

    music::ilog << "Computing transfer function via ClassEngine..." << std::endl;
    double wtime = get_wtime();

    the_ClassEngine_ = std::make_unique<ClassEngine>(pars_, false);

    wtime = get_wtime() - wtime;
    music::ilog << "CLASS took " << wtime << " s." << std::endl;
  }

  //! run ClassEngine with parameters set up
  void run_ClassEngine(double z, std::vector<double> &k, std::vector<double> &dc, std::vector<double> &tc, std::vector<double> &db, std::vector<double> &tb,
                       std::vector<double> &dn, std::vector<double> &tn, std::vector<double> &dm, std::vector<double> &tm)
  {
    k.clear(); 
    dc.clear(); db.clear(); dn.clear(); dm.clear();
    tc.clear(); tb.clear(); tn.clear(); tm.clear();
    
    the_ClassEngine_->getTk(z, k, dc, db, dn, dm, tc, tb, tn, tm);

    const double h  = cosmo_params_.get("h");

    for (size_t i = 0; i < k.size(); ++i)
    {
      // convert to 'CAMB' format, since we interpolate loglog and
      // don't want negative numbers...
      auto ik2 = 1.0 / (k[i] * k[i]) * h * h;
      dc[i] = -dc[i] * ik2;
      db[i] = -db[i] * ik2;
      dn[i] = -dn[i] * ik2;
      dm[i] = -dm[i] * ik2;
      tc[i] = -tc[i] * ik2;
      tb[i] = -tb[i] * ik2;
      tn[i] = -tn[i] * ik2;
      tm[i] = -tm[i] * ik2;
    }
  }

public:
  explicit transfer_CLASS_plugin(config_file &cf, const cosmology::parameters& cosmo_params)
      : transfer_function_plugin(cf,cosmo_params)
  {
    this->tf_isnormalised_ = true;

    ofs_class_input_.open("input_class_parameters.ini", std::ios::trunc);

    // all cosmological parameters need to be passed through the_cosmo_calc
    
    ztarget_ = pcf_->get_value_safe<double>("cosmology", "ztarget", 0.0);
    atarget_ = 1.0 / (1.0 + ztarget_);
    zstart_ = pcf_->get_value<double>("setup", "zstart");
    astart_ = 1.0 / (1.0 + zstart_);

    h_ = cosmo_params_["h"];
    
    if (cosmo_params_["A_s"] > 0.0) {
      music::ilog << "CLASS: Using A_s=" << cosmo_params_["A_s"] << " to normalise the transfer function." << std::endl;
    }else{
      double sigma8 = cosmo_params_["sigma_8"];
      if( sigma8 < 0 ){
        throw std::runtime_error("Need to specify either A_s or sigma_8 for CLASS plugin...");
      }
      music::ilog << "CLASS: Using sigma8_ =" << sigma8<< " to normalise the transfer function." << std::endl;
    }

    // determine highest k we will need for the resolution selected
    double lbox = pcf_->get_value<double>("setup", "boxlength");
    int levelmax = pcf_->get_value<int>("setup", "levelmax");
    double dx = lbox / (1<<levelmax);
    kmax_ = std::max(20.0, M_PI /dx  * sqrt(3) * 2.0); // 200% of spatial diagonal, or k=20h Mpc-1

    // initialise CLASS and get the normalisation
    this->init_ClassEngine();
    double A_s_ = the_ClassEngine_->get_A_s(); // this either the input one, or the one computed from sigma8
    
    // compute the normalisation to interface with MUSIC
    double k_p = cosmo_params["k_p"] / cosmo_params["h"];
    tnorm_ = std::sqrt(2.0 * M_PI * M_PI * A_s_ * std::pow(1.0 / k_p, cosmo_params["n_s"] - 1) / std::pow(2.0 * M_PI, 3.0));

    // compute the transfer function at z=0 using CLASS engine
    std::vector<double> k, dc, tc, db, tb, dn, tn, dm, tm;
    this->run_ClassEngine(0.0, k, dc, tc, db, tb, dn, tn, dm, tm);

    delta_c0_.set_data(k, dc);
    theta_c0_.set_data(k, tc);
    delta_b0_.set_data(k, db);
    theta_b0_.set_data(k, tb);
    delta_n0_.set_data(k, dn);
    theta_n0_.set_data(k, tn);
    delta_m0_.set_data(k, dm);
    theta_m0_.set_data(k, tm);

     // compute the transfer function at z=z_target using CLASS engine
    this->run_ClassEngine(ztarget_, k, dc, tc, db, tb, dn, tn, dm, tm);
    delta_c_.set_data(k, dc);
    theta_c_.set_data(k, tc);
    delta_b_.set_data(k, db);
    theta_b_.set_data(k, tb);
    delta_n_.set_data(k, dn);
    theta_n_.set_data(k, tn);
    delta_m_.set_data(k, dm);
    theta_m_.set_data(k, tm);

    kmin_ = k[0];
    kmax_ = k.back();

    music::ilog << "CLASS table contains k = " << this->get_kmin() << " to " << this->get_kmax() << " h Mpc-1." << std::endl;

    tf_distinct_ = true;
    tf_withvel_ = true;
    tf_withtotal0_ = true;
  }

  ~transfer_CLASS_plugin()
  {
  }

  inline double compute(double k, tf_type type) const
  {
    k *= h_;

    if (k < kmin_ || k > kmax_)
    {
      return 0.0;
    }

    real_t val(0.0);
    switch (type)
    {
      // values at ztarget:
    case delta_matter:
      val = delta_m_(k); break;
    case delta_cdm:
      val = delta_c_(k); break;
    case delta_baryon:
      val = delta_b_(k); break;
    case theta_matter:
      val = theta_m_(k); break;
    case theta_cdm:
      val = theta_c_(k); break;
    case theta_baryon:
      val = theta_b_(k); break;
    case delta_bc:
      val = delta_b_(k)-delta_c_(k); break;
    case theta_bc:
      val = theta_b_(k)-theta_c_(k); break;

      // values at zstart:
    case delta_matter0:
      val = delta_m0_(k); break;
    case delta_cdm0:
      val = delta_c0_(k); break;
    case delta_baryon0:
      val = delta_b0_(k); break;
    case theta_matter0:
      val = theta_m0_(k); break;
    case theta_cdm0:
      val = theta_c0_(k); break;
    case theta_baryon0:
      val = theta_b0_(k); break;
    default:
      throw std::runtime_error("Invalid type requested in transfer function evaluation");
    }
    return val * tnorm_;
  }

  inline double get_kmin(void) const { return kmin_ / h_; }
  inline double get_kmax(void) const { return kmax_ / h_; }
};

namespace
{
transfer_function_plugin_creator_concrete<transfer_CLASS_plugin> creator("CLASS");
}

#endif // USE_CLASS
