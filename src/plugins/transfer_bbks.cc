/*
 
 transfer_bbks.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#include "transfer_function.hh"

//! Implementation of class TransferFunction_BBKS for the BBKS transfer function 
/*!
 This class implements the analytical fit to the matter transfer
 function by Bardeen, Bond, Kaiser & Szalay (BBKS).
 ( see Bardeen et al. (1986) )
 */
class transfer_bbks_plugin : public transfer_function_plugin{
private:
	double      m_Gamma;
	
public:
	//! Constructor
	/*!
	 \param aCosm Structure of type Cosmology carrying the cosmological parameters
	 \param bSugiyama flag whether the Sugiyama (1995) correction shall be applied (default=true)
	 */
	transfer_bbks_plugin( config_file& cf )
    : transfer_function_plugin( cf )
	{  
		double Omega0 = cosmo_.Omega_m;
		double FreeGamma = -1.0;
		
		bool bSugiyama(true);
		
		try{
			bSugiyama= pcf_->getValue<bool>( "cosmology", "sugiyama_corr" );
		}catch(...){
			throw std::runtime_error("Error in \'tranfer_bbks_plugin\': need to specify \'[cosmology]/sugiyama_corr = [true/false]");
		}
		
		FreeGamma = pcf_->getValueSafe<double>( "cosmology", "gamma", FreeGamma );
		
		if( FreeGamma <= 0.0 ){
			m_Gamma = Omega0*0.01*cosmo_.H0;
			if( bSugiyama )
				m_Gamma *= exp(-cosmo_.Omega_b*(1.0+sqrt(2.0*0.01*cosmo_.H0)/Omega0));
		}else
			m_Gamma = FreeGamma;
		
		tf_distinct_ = false;
		tf_withvel_  = false;
		
	}
	
	//! computes the value of the BBKS transfer function for mode k (in h/Mpc)
	inline double compute( double k, tf_type type ){
		double q, f1, f2;
		
		if(k < 1e-7 )
			return 1.0;
		
		q = k/(m_Gamma);
		f1 = log(1.0 + 2.34*q)/(2.34*q);
		f2 = 1.0 + q*(3.89 + q*(259.21 + q*(162.771336 + q*2027.16958081)));
		
		return f1/sqrt(sqrt(f2));
		
	}
	
	inline double get_kmin( void ){
		return 1e-4;
	}
	
	inline double get_kmax( void ){
		return 1.e4;
	}
};


namespace{
	transfer_function_plugin_creator_concrete< transfer_bbks_plugin > creator("bbks");
}

