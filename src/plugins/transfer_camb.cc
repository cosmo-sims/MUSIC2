/*

 transfer_camb.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions for cosmological simulations

 Copyright (C) 2010  Oliver Hahn

*/

#include "transfer_function.hh"

const double tiny = 1e-30;

class transfer_CAMB_plugin : public transfer_function_plugin {

private:
  std::string m_filename_Pk, m_filename_Tk;
  std::vector<double> m_tab_k, m_tab_Tk_tot, m_tab_Tk_cdm, m_tab_Tk_baryon;
  std::vector<double> m_tab_Tvk_tot, m_tab_Tvk_cdm, m_tab_Tvk_baryon;
  gsl_interp_accel *acc_tot, *acc_cdm, *acc_baryon;
  gsl_interp_accel *acc_vtot, *acc_vcdm, *acc_vbaryon;
  gsl_spline *spline_tot, *spline_cdm, *spline_baryon;
  gsl_spline *spline_vtot, *spline_vcdm, *spline_vbaryon;

  double m_kmin, m_kmax, m_Omega_b, m_Omega_m, m_zstart;
  unsigned m_nlines;

  bool m_linbaryoninterp;

  void read_table(void) {

    m_nlines = 0;
    m_linbaryoninterp = false;

#ifdef WITH_MPI
    if (MPI::COMM_WORLD.Get_rank() == 0) {
#endif
      LOGINFO("Reading tabulated transfer function data from file \n    \'%s\'",
              m_filename_Tk.c_str());

      LOGINFO("CAUTION: make sure that this transfer function \n\t has been output for z=%f!",m_zstart);

      std::string line;
      std::ifstream ifs(m_filename_Tk.c_str());

      if (!ifs.good())
        throw std::runtime_error("Could not find transfer function file \'" +
                                 m_filename_Tk + "\'");

      m_tab_k.clear();
      m_tab_Tk_tot.clear();
      m_tab_Tk_cdm.clear();
      m_tab_Tk_baryon.clear();
      m_tab_Tvk_tot.clear();
      m_tab_Tvk_cdm.clear();    //>[150609SH: add]
      m_tab_Tvk_baryon.clear(); //>[150609SH: add]

      m_kmin = 1e30;
      m_kmax = -1e30;
      std::ofstream ofs("dump_transfer.txt");

      while (!ifs.eof()) {
        getline(ifs, line);
        if (ifs.eof()) break;

        // OH: ignore line if it has a comment:
        if (line.find("#") != std::string::npos) continue;

        std::stringstream ss(line);

        double k, Tkc, Tkb, Tktot, Tkvtot, Tkvc, Tkvb, dummy;

        ss >> k;
        ss >> Tkc;   // cdm
        ss >> Tkb;   // baryon
        ss >> dummy; // photon
        ss >> dummy; // nu
        ss >> dummy; // mass_nu
        ss >> Tktot; // total
        ss >> dummy; // no_nu
        ss >> dummy; // total_de
        ss >> dummy; // Weyl
        ss >> Tkvc;  //>[150609SH: add] // v_cdm
        ss >> Tkvb;  //>[150609SH: add] // v_b
        ss >> dummy; // v_b-v_cdm

        if( ss.bad() || ss.fail() ){
          LOGERR("Error reading the transfer function file (corrupt or not in expected format)!");
          throw std::runtime_error("Error reading transfer function file \'" +
                                   m_filename_Tk + "\'");
        }

	if( m_Omega_b < 1e-6 ) Tkvtot = Tktot;
	else Tkvtot=((m_Omega_m-m_Omega_b)*Tkvc+m_Omega_b*Tkvb)/m_Omega_m; //MvD

        m_linbaryoninterp |= Tkb < 0.0 || Tkvb < 0.0;

        m_tab_k.push_back(log10(k));

        m_tab_Tk_tot.push_back(Tktot);
        m_tab_Tk_baryon.push_back(Tkb);
        m_tab_Tk_cdm.push_back(Tkc);
        m_tab_Tvk_tot.push_back(Tkvtot);
        m_tab_Tvk_baryon.push_back(Tkvb);
        m_tab_Tvk_cdm.push_back(Tkvc);

        ++m_nlines;

        if (k < m_kmin)
          m_kmin = k;
        if (k > m_kmax)
          m_kmax = k;
      }

      for (size_t i = 0; i < m_tab_k.size(); ++i) {
        m_tab_Tk_tot[i] = log10(m_tab_Tk_tot[i]);
        m_tab_Tk_cdm[i] = log10(m_tab_Tk_cdm[i]);
        m_tab_Tvk_cdm[i] = log10(m_tab_Tvk_cdm[i]);
        m_tab_Tvk_tot[i] = log10(m_tab_Tvk_tot[i]);

        if (!m_linbaryoninterp) {
          m_tab_Tk_baryon[i] = log10(m_tab_Tk_baryon[i]);
          m_tab_Tvk_baryon[i] = log10(m_tab_Tvk_baryon[i]);
        }
      }

      ifs.close();

      LOGINFO("Read CAMB transfer function table with %d rows", m_nlines);

      if (m_linbaryoninterp)
        LOGINFO("Using log-lin interpolation for baryons\n    (TF is not "
                "positive definite)");

#ifdef WITH_MPI
    }

    unsigned n = m_tab_k.size();
    MPI::COMM_WORLD.Bcast(&n, 1, MPI_UNSIGNED, 0);

    if (MPI::COMM_WORLD.Get_rank() > 0) {
      m_tab_k.assign(n, 0);
      m_tab_Tk_tot.assign(n, 0);
      m_tab_Tk_cdm.assign(n, 0);
      m_tab_Tk_baryon.assign(n, 0);
      m_tab_Tvk_tot.assign(n, 0);
      m_tab_Tvk_cdm.assign(n, 0);    //>[150609SH: add]
      m_tab_Tvk_baryon.assign(n, 0); //>[150609SH: add]
    }

    MPI::COMM_WORLD.Bcast(&m_tab_k[0], n, MPI_DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&m_tab_Tk_tot[0], n, MPI_DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&m_tab_Tk_cdm[0], n, MPI_DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&m_tab_Tk_baryon[0], n, MPI_DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&m_tab_Tvk_tot[0], n, MPI_DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&m_tab_Tvk_cdm[0], n, MPI_DOUBLE,
                          0); //>[150609SH: add]
    MPI::COMM_WORLD.Bcast(&m_tab_Tvk_baryon[0], n, MPI_DOUBLE,
                          0); //>[150609SH: add]

#endif
  }

public:
  transfer_CAMB_plugin(config_file &cf)
  : transfer_function_plugin(cf)
  {
    m_filename_Tk = pcf_->getValue<std::string>("cosmology", "transfer_file");
    m_Omega_m=cf.getValue<double>("cosmology","Omega_m"); //MvD
    m_Omega_b=cf.getValue<double>("cosmology","Omega_b"); //MvD
    m_zstart =cf.getValue<double>("setup","zstart"); //MvD

    read_table();

    acc_tot = gsl_interp_accel_alloc();
    acc_cdm = gsl_interp_accel_alloc();
    acc_baryon = gsl_interp_accel_alloc();
    acc_vtot = gsl_interp_accel_alloc();    //>[150609SH: add]
    acc_vcdm = gsl_interp_accel_alloc();    //>[150609SH: add]
    acc_vbaryon = gsl_interp_accel_alloc(); //>[150609SH: add]

    spline_tot = gsl_spline_alloc(gsl_interp_cspline, m_tab_k.size());
    spline_cdm = gsl_spline_alloc(gsl_interp_cspline, m_tab_k.size());
    spline_baryon = gsl_spline_alloc(gsl_interp_cspline, m_tab_k.size());
    spline_vtot = gsl_spline_alloc(gsl_interp_cspline, m_tab_k.size());
    spline_vcdm =
        gsl_spline_alloc(gsl_interp_cspline, m_tab_k.size()); //>[150609SH: add]
    spline_vbaryon =
        gsl_spline_alloc(gsl_interp_cspline, m_tab_k.size()); //>[150609SH: add]

    gsl_spline_init(spline_tot, &m_tab_k[0], &m_tab_Tk_tot[0], m_tab_k.size());
    gsl_spline_init(spline_cdm, &m_tab_k[0], &m_tab_Tk_cdm[0], m_tab_k.size());
    gsl_spline_init(spline_baryon, &m_tab_k[0], &m_tab_Tk_baryon[0],
                    m_tab_k.size());
    gsl_spline_init(spline_vtot, &m_tab_k[0], &m_tab_Tvk_tot[0],
                    m_tab_k.size()); //>[150609SH: add]
    gsl_spline_init(spline_vcdm, &m_tab_k[0], &m_tab_Tvk_cdm[0],
                    m_tab_k.size()); //>[150609SH: add]
    gsl_spline_init(spline_vbaryon, &m_tab_k[0], &m_tab_Tvk_baryon[0],
                    m_tab_k.size()); //>[150609SH: add]

    tf_distinct_ = true; // [150612SH: different density between CDM v.s. Baryon]
    tf_withvel_  = true; // [150612SH: using velocity transfer function]
  }

  ~transfer_CAMB_plugin() {
    gsl_spline_free(spline_tot);
    gsl_spline_free(spline_cdm);
    gsl_spline_free(spline_baryon);
    gsl_spline_free(spline_vtot);
    gsl_spline_free(spline_vcdm);    //>[150609SH: add]
    gsl_spline_free(spline_vbaryon); //>[150609SH: add]

    gsl_interp_accel_free(acc_tot);
    gsl_interp_accel_free(acc_cdm);
    gsl_interp_accel_free(acc_baryon);
    gsl_interp_accel_free(acc_vtot);
    gsl_interp_accel_free(acc_vcdm);    //>[150609SH: add]
    gsl_interp_accel_free(acc_vbaryon); //>[150609SH: add]
  }

  // linear interpolation in log-log
  inline double extrap_right(double k, const tf_type &type) {
    int n = m_tab_k.size() - 1, n1 = n - 1;

    double v1(1.0), v2(1.0);

    double lk = log10(k);
    double dk = m_tab_k[n] - m_tab_k[n1];
    double delk = lk - m_tab_k[n];

    switch (type) {
    case cdm:
      v1 = m_tab_Tk_cdm[n1];
      v2 = m_tab_Tk_cdm[n];
      return pow(10.0, (v2 - v1) / dk * (delk) + v2);
    case baryon:
      v1 = m_tab_Tk_baryon[n1];
      v2 = m_tab_Tk_baryon[n];
      if (m_linbaryoninterp)
        return std::max((v2 - v1) / dk * (delk) + v2, tiny);
      return pow(10.0, (v2 - v1) / dk * (delk) + v2);
    case vtotal: //>[150609SH: add]
      v1 = m_tab_Tvk_tot[n1];
      v2 = m_tab_Tvk_tot[n];
      return pow(10.0, (v2 - v1) / dk * (delk) + v2);
    case vcdm: //>[150609SH: add]
      v1 = m_tab_Tvk_cdm[n1];
      v2 = m_tab_Tvk_cdm[n];
      return pow(10.0, (v2 - v1) / dk * (delk) + v2);
    case vbaryon: //>[150609SH: add]
      v1 = m_tab_Tvk_baryon[n1];
      v2 = m_tab_Tvk_baryon[n];
      if (m_linbaryoninterp)
        return std::max((v2 - v1) / dk * (delk) + v2, tiny);
      return pow(10.0, (v2 - v1) / dk * (delk) + v2);
    case total:
      v1 = m_tab_Tk_tot[n1];
      v2 = m_tab_Tk_tot[n];
      return pow(10.0, (v2 - v1) / dk * (delk) + v2);
    default:
      throw std::runtime_error(
          "Invalid type requested in transfer function evaluation");
    }

    return 0.0;
  }

  inline double compute(double k, tf_type type) {
    // use constant interpolation on the left side of the tabulated values
    if (k < m_kmin) {
      switch (type) {
      case cdm:
        return pow(10.0, m_tab_Tk_cdm[0]);
      case baryon:
        if (m_linbaryoninterp)
          return m_tab_Tk_baryon[0];
        return pow(10.0, m_tab_Tk_baryon[0]);
      case vtotal:
        return pow(10.0, m_tab_Tvk_tot[0]);
      case vcdm:
        return pow(10.0, m_tab_Tvk_cdm[0]);
      case vbaryon:
        if (m_linbaryoninterp)
          return m_tab_Tvk_baryon[0];
        return pow(10.0, m_tab_Tvk_baryon[0]);
      case total:
        return pow(10.0, m_tab_Tk_tot[0]);
      default:
        throw std::runtime_error(
            "Invalid type requested in transfer function evaluation");
      }
    }
    // use linear interpolation on the right side of the tabulated values
    else if (k > m_kmax)
      return extrap_right(k, type);

    double lk = log10(k);
    switch (type) {
    case cdm:
      return pow(10.0, gsl_spline_eval(spline_cdm, lk, acc_cdm));
    case baryon:
      if (m_linbaryoninterp)
        return gsl_spline_eval(spline_baryon, lk, acc_baryon);
      return pow(10.0, gsl_spline_eval(spline_baryon, lk, acc_baryon));
    case vtotal:
      return pow(10.0, gsl_spline_eval(spline_vtot, lk, acc_vtot)); //MvD
    case vcdm:
      return pow(10.0, gsl_spline_eval(spline_vcdm, lk, acc_vcdm));
    case vbaryon:
      if (m_linbaryoninterp)
        return gsl_spline_eval(spline_vbaryon, lk, acc_vbaryon);
      return pow(10.0, gsl_spline_eval(spline_vbaryon, lk, acc_vbaryon));
    case total:
      return pow(10.0, gsl_spline_eval(spline_tot, lk, acc_tot));
    default:
      throw std::runtime_error(
          "Invalid type requested in transfer function evaluation");
    }
  }

  inline double get_kmin(void) { return pow(10.0, m_tab_k[1]); }

  inline double get_kmax(void) {
    return pow(10.0, m_tab_k[m_tab_k.size() - 2]);
  }
};

namespace {
transfer_function_plugin_creator_concrete<transfer_CAMB_plugin>
    creator("camb_file");
}
