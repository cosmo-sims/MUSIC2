/*
 
 transfer_camb.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#include "transfer_function.hh"

class transfer_MUSIC_plugin : public transfer_function_plugin
{
	
private:
	std::string m_filename_Pk, m_filename_Tk;
	std::vector<double> m_tab_k, m_tab_Tk_tot, m_tab_Tk_cdm, m_tab_Tk_baryon, m_tab_Tvk_cdm, m_tab_Tvk_baryon;
	gsl_interp_accel *acc_dtot, *acc_dcdm, *acc_dbaryon, *acc_vcdm, *acc_vbaryon;
	gsl_spline *spline_dtot, *spline_dcdm, *spline_dbaryon, *spline_vcdm, *spline_vbaryon;
	
	
	
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
			m_tab_Tvk_cdm.clear();
			m_tab_Tvk_baryon.clear();
			
			const double zero = 1e-10;
			
			while( !ifs.eof() ){
				getline(ifs,line);
				
				if(ifs.eof()) break;
				
				std::stringstream ss(line);
				
				double k, Tkc, Tkb, Tktot, Tkvc, Tkvb;
				ss >> k;
				ss >> Tktot;
				ss >> Tkc;
				ss >> Tkb;
				ss >> Tkvc;
				ss >> Tkvb;
				
				Tktot = std::max(zero,Tktot);
				Tkc   = std::max(zero,Tkc);
				Tkb   = std::max(zero,Tkb);
				Tkvc  = std::max(zero,Tkvc);
				Tkvb  = std::max(zero,Tkvb);
				
				m_tab_k.push_back( log10(k) );
				
				m_tab_Tk_tot.push_back( log10(Tktot) );
				m_tab_Tk_baryon.push_back( log10(Tkb) );
				m_tab_Tk_cdm.push_back( log10(Tkc) );
				m_tab_Tvk_cdm.push_back( log10(Tkvc) );
				m_tab_Tvk_baryon.push_back( log10(Tkvb) );
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
			m_tab_Tvk_cdm.assign(n,0);
			m_tab_Tvk_baryon.assign(n,0);
		}
		
		MPI::COMM_WORLD.Bcast( &m_tab_k[0],  n, MPI_DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( &m_tab_Tk_tot[0], n, MPI_DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( &m_tab_Tk_cdm[0], n, MPI_DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( &m_tab_Tk_baryon[0], n, MPI_DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( &m_tab_Tvk_cdm[0], n, MPI_DOUBLE, 0 );
		MPI::COMM_WORLD.Bcast( &m_tab_Tvk_baryon[0], n, MPI_DOUBLE, 0 );
#endif
		
	}
	
public:
	transfer_MUSIC_plugin( config_file& cf )
	: transfer_function_plugin( cf )
	{
		m_filename_Tk = pcf_->getValue<std::string>("cosmology","transfer_file");
		
		read_table( );
		
		acc_dtot = gsl_interp_accel_alloc();
		acc_dcdm = gsl_interp_accel_alloc();
		acc_dbaryon = gsl_interp_accel_alloc();
		acc_vcdm = gsl_interp_accel_alloc();
		acc_vbaryon = gsl_interp_accel_alloc();
		
		spline_dtot = gsl_spline_alloc( gsl_interp_cspline, m_tab_k.size() );
		spline_dcdm = gsl_spline_alloc( gsl_interp_cspline, m_tab_k.size() );
		spline_dbaryon = gsl_spline_alloc( gsl_interp_cspline, m_tab_k.size() );
		spline_vcdm = gsl_spline_alloc( gsl_interp_cspline, m_tab_k.size() );
		spline_vbaryon = gsl_spline_alloc( gsl_interp_cspline, m_tab_k.size() );
		
		gsl_spline_init (spline_dtot, &m_tab_k[0], &m_tab_Tk_tot[0], m_tab_k.size() );
		gsl_spline_init (spline_dcdm, &m_tab_k[0], &m_tab_Tk_cdm[0], m_tab_k.size() );
		gsl_spline_init (spline_dbaryon, &m_tab_k[0], &m_tab_Tk_baryon[0], m_tab_k.size() );
		gsl_spline_init (spline_vcdm, &m_tab_k[0], &m_tab_Tvk_cdm[0], m_tab_k.size() );
		gsl_spline_init (spline_vbaryon, &m_tab_k[0], &m_tab_Tvk_baryon[0], m_tab_k.size() );
		
		tf_distinct_ = true;
		tf_withvel_  = true;
	}
	
	~transfer_MUSIC_plugin()
	{
		gsl_spline_free (spline_dtot);
		gsl_spline_free (spline_dcdm);
		gsl_spline_free (spline_dbaryon);
		gsl_spline_free (spline_vcdm);
		gsl_spline_free (spline_vbaryon);
		
		gsl_interp_accel_free (acc_dtot);
		gsl_interp_accel_free (acc_dcdm);
		gsl_interp_accel_free (acc_dbaryon);
		gsl_interp_accel_free (acc_vcdm);
		gsl_interp_accel_free (acc_vbaryon);
	}
	
	inline double extrap_left( double k, const tf_type& type ) 
	{
		if( k<1e-8 )
			return 1.0;
		
		double v1(1.0), v2(1.0);
		switch( type )
		{
			case cdm:
				v1 = m_tab_Tk_cdm[0];
				v2 = m_tab_Tk_cdm[1];
				break;
			case baryon:
				v1 = m_tab_Tk_baryon[0];
				v2 = m_tab_Tk_baryon[1];
				break;
			case vcdm:
				v1 = m_tab_Tvk_cdm[0];
				v2 = m_tab_Tvk_cdm[1];
				break;
			case vbaryon:
				v1 = m_tab_Tvk_baryon[0];
				v2 = m_tab_Tvk_baryon[1];
				break;
			case total: 
				v1 = m_tab_Tk_tot[0];
				v2 = m_tab_Tk_tot[1];
				break;
				
			default:
				throw std::runtime_error("Invalid type requested in transfer function evaluation");
		}
		
		double lk = log10(k);
		double dk = m_tab_k[1]-m_tab_k[0];
		double delk = lk-m_tab_k[0];
		
		//double xi = (v2-v1)/dk;
		return pow(10.0,(v2-v1)/dk*(delk)+v1);
	}
	
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
				v1 = m_tab_Tvk_cdm[n1];
				v2 = m_tab_Tvk_cdm[n];
				break;
			case vbaryon:
				v1 = m_tab_Tvk_baryon[n1];
				v2 = m_tab_Tvk_baryon[n];
				break;
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
		
		//double xi = (v2-v1)/dk;
		return pow(10.0,(v2-v1)/dk*(delk)+v2);
	}
	
	inline double compute( double k, tf_type type ){
		
		double lk = log10(k);
		
		//if( lk<m_tab_k[1])
		//	return 1.0;
		
		//if( lk>m_tab_k[m_tab_k.size()-2] );
		//	return m_tab_Tk_cdm[m_tab_k.size()-2]/k/k;
		
		if( k<get_kmin() )
			return extrap_left(k, type );
		
		if( k>get_kmax() )
			return extrap_right(k,type );
		
		
		switch( type )
		{
			case cdm:
				return pow(10.0, gsl_spline_eval (spline_dcdm, lk, acc_dcdm) );
			case baryon:
				return pow(10.0, gsl_spline_eval (spline_dbaryon, lk, acc_dbaryon) );
			case vcdm:
				return pow(10.0, gsl_spline_eval (spline_vcdm, lk, acc_vcdm) );
			case vbaryon:
				return pow(10.0, gsl_spline_eval (spline_vbaryon, lk, acc_vbaryon) );
			case total: 
				return pow(10.0, gsl_spline_eval (spline_dtot, lk, acc_dtot) );
				
			default:
				throw std::runtime_error("Invalid type requested in transfer function evaluation");
		}
		
		return 1.0;
	}
	
	inline double get_kmin( void ){
		return pow(10.0,m_tab_k[0]);
	}
	
	inline double get_kmax( void ){
		return pow(10.0,m_tab_k[m_tab_k.size()-1]);
	}
	
};

namespace{
	transfer_function_plugin_creator_concrete< transfer_MUSIC_plugin > creator("music");
}


