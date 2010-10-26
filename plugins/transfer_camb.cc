/*
 
 transfer_camb.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 */

#include "transfer_function.hh"

class transfer_CAMB_plugin : public transfer_function_plugin
{
	
private:
	//Cosmology m_Cosmology;
	
	std::string m_filename_Pk, m_filename_Tk;
	std::vector<double> m_tab_k, m_tab_Tk_tot, m_tab_Tk_cdm, m_tab_Tk_baryon;
	
	Spline_interp *m_psinterp;
	gsl_interp_accel *acc_tot, *acc_cdm, *acc_baryon;
	gsl_spline *spline_tot, *spline_cdm, *spline_baryon;
	
	
	
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
				
				m_tab_k.push_back( log10(k) );
				
				m_tab_Tk_tot.push_back( log10(Tktot) );
				m_tab_Tk_baryon.push_back( log10(Tkb) );
				m_tab_Tk_cdm.push_back( log10(Tkc) );
				
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
	transfer_CAMB_plugin( config_file& cf )//Cosmology aCosm, std::string filename_Tk, TFtype iwhich )
	: transfer_function_plugin( cf )//, m_filename_Tk( filename_Tk ), m_psinterp( NULL )
	{
		m_filename_Tk = pcf_->getValue<std::string>("cosmology","transfer_file");
		
		read_table( );
		
		acc_tot = gsl_interp_accel_alloc();
		acc_cdm = gsl_interp_accel_alloc();
		acc_baryon = gsl_interp_accel_alloc();
		
		
		spline_tot = gsl_spline_alloc( gsl_interp_akima, m_tab_k.size() );
		spline_cdm = gsl_spline_alloc( gsl_interp_akima, m_tab_k.size() );
		spline_baryon = gsl_spline_alloc( gsl_interp_akima, m_tab_k.size() );
		
		/*spline_tot = gsl_spline_alloc( gsl_interp_linear, m_tab_k.size() );
		spline_cdm = gsl_spline_alloc( gsl_interp_linear, m_tab_k.size() );
		spline_baryon = gsl_spline_alloc( gsl_interp_linear, m_tab_k.size() );*/
		
		gsl_spline_init (spline_tot, &m_tab_k[0], &m_tab_Tk_tot[0], m_tab_k.size() );
		gsl_spline_init (spline_cdm, &m_tab_k[0], &m_tab_Tk_cdm[0], m_tab_k.size() );
		gsl_spline_init (spline_baryon, &m_tab_k[0], &m_tab_Tk_baryon[0], m_tab_k.size() );
		
		tf_distinct_ = true;
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
	
	inline double compute( double k, tf_type type ){
		
		double lk = log10(k);
		
		//if( lk<m_tab_k[1])
		//	return 1.0;
		
		//if( lk>m_tab_k[m_tab_k.size()-2] );
		//	return m_tab_Tk_cdm[m_tab_k.size()-2]/k/k;
		
		if( type == total )
			return pow(10.0, gsl_spline_eval (spline_tot, lk, acc_tot) );
		if( type == cdm )
			return pow(10.0, gsl_spline_eval (spline_cdm, lk, acc_cdm) );
		
		return pow(10.0, gsl_spline_eval (spline_baryon, lk, acc_baryon) );
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


