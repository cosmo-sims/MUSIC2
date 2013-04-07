/*
 
 transfer_camb.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
*/

#include "transfer_function.hh"

class transfer_CAMB_plugin : public transfer_function_plugin
{
	
private:
	std::string m_filename_Pk, m_filename_Tk;
	std::vector<double> m_tab_k, m_tab_Tk_tot, m_tab_Tk_cdm, m_tab_Tk_baryon;
	gsl_interp_accel *acc_tot, *acc_cdm, *acc_baryon;
	gsl_spline *spline_tot, *spline_cdm, *spline_baryon;
	
    double m_kmin, m_kmax;
    unsigned m_nlines;
	
	void read_table( void ){
#ifdef WITH_MPI
		if( MPI::COMM_WORLD.Get_rank() == 0 ){
#endif
			std::cerr 
			<< " - reading tabulated transfer function data from file \n"
			<< "    \'" << m_filename_Tk << "\'\n";
			
			std::string line;
			std::ifstream ifs( m_filename_Tk.c_str() );
			
			if(! ifs.good() )
				throw std::runtime_error("Could not find transfer function file \'"+m_filename_Tk+"\'");
			
			m_tab_k.clear();
			m_tab_Tk_tot.clear();
			m_tab_Tk_cdm.clear();
			m_tab_Tk_baryon.clear();
          
            m_kmin = 1e30;
            m_kmax = -1e30;
            m_nlines = 0;
			
			while( !ifs.eof() ){
				getline(ifs,line);
				
				if(ifs.eof()) break;
				
				std::stringstream ss(line);
				
				double k, Tkc, Tkb, Tkg, Tkr, Tknu, Tktot;
				ss >> k;
				ss >> Tkc;
				ss >> Tkb;
				ss >> Tkg;
				ss >> Tkr;
				ss >> Tknu;
				ss >> Tktot;
              
                if( k < m_kmin ) m_kmin = k;
                if( k > m_kmax ) m_kmax = k;
				
				m_tab_k.push_back( log10(k) );
				
				m_tab_Tk_tot.push_back( log10(Tktot) );
				m_tab_Tk_baryon.push_back( log10(Tkb) );
				m_tab_Tk_cdm.push_back( log10(Tkc) );
                ++m_nlines;
			}
			
			ifs.close();
			
			
			
			
#ifdef WITH_MPI
		}
		
		unsigned n=m_tab_k.size();
		MPI::COMM_WORLD.Bcast( &n, 1, MPI_UNSIGNED, 0 );
		
		if( MPI::COMM_WORLD.Get_rank() > 0 ){
			m_tab_k.assign(n,0);
			m_tab_Tk_tot.assign(n,0);
			m_tab_Tk_cdm.assign(n,0);
			m_tab_Tk_baryon.assign(n,0);

		}
		
		MPI::COMM_WORLD.Bcast( &m_tab_k[0],  n, MPI_DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( &m_tab_Tk_tot[0], n, MPI_DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( &m_tab_Tk_cdm[0], n, MPI_DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( &m_tab_Tk_baryon[0], n, MPI_DOUBLE, 0 );

#endif
		
	}
	
public:
	transfer_CAMB_plugin( config_file& cf )
	: transfer_function_plugin( cf )
	{
		m_filename_Tk = pcf_->getValue<std::string>("cosmology","transfer_file");
		
		read_table( );
		
		acc_tot = gsl_interp_accel_alloc();
		acc_cdm = gsl_interp_accel_alloc();
		acc_baryon = gsl_interp_accel_alloc();
		
		
		spline_tot = gsl_spline_alloc( gsl_interp_cspline, m_tab_k.size() );
		spline_cdm = gsl_spline_alloc( gsl_interp_cspline, m_tab_k.size() );
		spline_baryon = gsl_spline_alloc( gsl_interp_cspline, m_tab_k.size() );
		
		gsl_spline_init (spline_tot, &m_tab_k[0], &m_tab_Tk_tot[0], m_tab_k.size() );
		gsl_spline_init (spline_cdm, &m_tab_k[0], &m_tab_Tk_cdm[0], m_tab_k.size() );
		gsl_spline_init (spline_baryon, &m_tab_k[0], &m_tab_Tk_baryon[0], m_tab_k.size() );
		
		tf_distinct_ = true;
		tf_withvel_  = false;
	}
	
	~transfer_CAMB_plugin()
	{
		gsl_spline_free (spline_tot);
		gsl_spline_free (spline_cdm);
		gsl_spline_free (spline_baryon);
		
		gsl_interp_accel_free (acc_tot);
		gsl_interp_accel_free (acc_cdm);
		gsl_interp_accel_free (acc_baryon);
	}

    // linear interpolation in log-log
    inline double extrap_right( double k, const tf_type& type )
    {
      double v1(1.0), v2(1.0);
      
      int n=m_tab_k.size()-1, n1=n-1;
      switch( type )
      {
        case cdm:
          v1 = m_tab_Tk_cdm[n1];
          v2 = m_tab_Tk_cdm[n];
          break;
        case baryon:
          v1 = m_tab_Tk_baryon[n1];
          v2 = m_tab_Tk_baryon[n];
          break;
        case vcdm:
        case vbaryon:
        case total:
          v1 = m_tab_Tk_tot[n1];
          v2 = m_tab_Tk_tot[n];
          break;
          
        default:
          throw std::runtime_error("Invalid type requested in transfer function evaluation");
      }
      
      double lk = log10(k);
      double dk = m_tab_k[n]-m_tab_k[n1];
      double delk = lk-m_tab_k[n];
      
      return pow(10.0,(v2-v1)/dk*(delk)+v2);
    }
    
	inline double compute( double k, tf_type type )
    {
	    // use constant interpolation on the left side of the tabulated values
        if( k < m_kmin )
        {
          if( type == cdm )
            return pow(10.0,m_tab_Tk_cdm[0]);
          
          else if( type == baryon )
            return pow(10.0,m_tab_Tk_baryon[0]);
          
          return pow(10.0,m_tab_Tk_tot[0]);
          
        }
        // use linear interpolation on the right side of the tabulated values
        else if( k>m_kmax )
          return extrap_right( k, type );
          
      
        double lk = log10(k);
		if( type == cdm )
			return pow(10.0, gsl_spline_eval (spline_cdm, lk, acc_cdm) );

		if( type == baryon )
			return pow(10.0, gsl_spline_eval (spline_baryon, lk, acc_baryon) );
		
		return pow(10.0, gsl_spline_eval (spline_tot, lk, acc_tot) );
	}
	
	inline double get_kmin( void ){
		return pow(10.0,m_tab_k[1]);
	}
	
	inline double get_kmax( void ){
		return pow(10.0,m_tab_k[m_tab_k.size()-2]);
	}
	
};

namespace{
	transfer_function_plugin_creator_concrete< transfer_CAMB_plugin > creator("camb_file");
}


