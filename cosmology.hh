/*
 
 cosmology.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#ifndef _COSMOLOGY_HH
#define _COSMOLOGY_HH


#include "transfer_function.hh"
#include "mesh.hh"
#include "general.hh"

/*!
 * @class CosmoCalc
 * @brief provides functions to compute cosmological quantities
 *
 * This class provides member functions to compute cosmological quantities
 * related to the Friedmann equations and linear perturbation theory
 */
class CosmoCalc
{
public:
	//! data structure to store cosmological parameters
	Cosmology m_Cosmology;
	
	//! pointer to an instance of a transfer function plugin
	transfer_function_plugin *m_pTransferFunction;
	
	
	//! constructor for a cosmology calculator object
	/*!
	 * @param acosmo a cosmological parameters structure
	 * @param pTransferFunction pointer to an instance of a transfer function object
	 */
	 
	CosmoCalc( const Cosmology acosmo, transfer_function_plugin *pTransferFunction )
	{
		m_Cosmology = acosmo;
		m_pTransferFunction = pTransferFunction;
	}
	
	//! returns the amplitude of amplitude of the power spectrum
	/*!
	 * @param k the wave number in h/Mpc
	 * @param a the expansion factor of the universe
	 * @returns power spectrum amplitude for wave number k at time a
	 */
	inline real_t Power( real_t k, real_t a ){
		real_t m_Dplus    = CalcGrowthFactor( a );
		real_t m_DplusOne = CalcGrowthFactor( 1.0 );
		real_t m_pNorm = ComputePNorm( 1e4 );
		m_Dplus    /= m_DplusOne;
		m_DplusOne = 1.0;
		real_t scale = m_Dplus/m_DplusOne;
		return m_pNorm*scale*scale*TransferSq(k)*pow((double)k,(double)m_Cosmology.nspect);
	}
	
	//! integrand function for Calc_fPeebles
	/*!
	 * @sa Calc_fPeebles
	 */
	inline static double fy( double a, void *Params )
	{
		Cosmology *cosm = (Cosmology*)Params;
		double y = cosm->Omega_m*(1.0/a-1.0) + cosm->Omega_L*(a*a-1.0) + 1.0;
		return 1.0/pow(y,1.5);
	}
	
	//! calculates d log D+/d log a
	/*! this version follows the Peebles (TBD: add citation)
	 *  formula to compute Bertschinger's vfact
	 */
	inline real_t Calc_fPeebles( real_t a )
	{
		real_t y = m_Cosmology.Omega_m*(1.0/a-1.0) + m_Cosmology.Omega_L*(a*a-1.0) + 1.0;
		real_t fact = integrate( &fy, 1e-6, a, (void*)&m_Cosmology );
		return (m_Cosmology.Omega_L*a*a-0.5*m_Cosmology.Omega_m/a)/y - 1.0 + a*fy(a,(void*)&m_Cosmology)/fact;
	}
	
	//! Computes the linear theory growth factor D+
	/*! Function integrates over member function GrowthIntegrand */
	inline real_t CalcGrowthFactor( real_t a )
	{ 
		real_t eta =  sqrt((double)(m_Cosmology.Omega_r/a/a+m_Cosmology.Omega_m/a+m_Cosmology.Omega_L*a*a
								  +1.0-m_Cosmology.Omega_m-m_Cosmology.Omega_L));
		
		real_t integral = integrate( &GrowthIntegrand, 0.0, a, (void*)&m_Cosmology );
		return eta/a*integral;
	}
    
    inline static double InvHubble( double a, void *Params )
    {
        Cosmology *cosm = (Cosmology*)Params;
		double eta = 1.0/(100.0*sqrt((double)(cosm->Omega_m/(a*a*a)+cosm->Omega_L
								   +(1.0-cosm->Omega_m-cosm->Omega_L)/(a*a))));
        return eta;
    }
    
    inline real_t CalcTFac( real_t a0, real_t a1 )
    {
        real_t fact = integrate( &InvHubble, a0, a1, (void*)&m_Cosmology );
        return 1.0/fact;
    }
    
	//! Integrand used by function CalcGrowthFactor to determine the linear growth factor D+
	inline static double GrowthIntegrand( double a, void *Params )
	{
		Cosmology *cosm = (Cosmology*)Params;
		double eta = sqrt((double)(cosm->Omega_r/a/a+cosm->Omega_m/a+cosm->Omega_L*a*a
								   +1.0-cosm->Omega_m-cosm->Omega_L));
		return 2.5/(eta*eta*eta);
	}
    
    //! Computes the linear theory growth factor D+
	/*! Function integrates over member function GrowthIntegrand */
	inline real_t CalcGrowthFactor_Matter( real_t a )
	{ 
		real_t eta =  sqrt((double)(m_Cosmology.Omega_m/a+m_Cosmology.Omega_L*a*a
                                    +1.0-m_Cosmology.Omega_m-m_Cosmology.Omega_L));
		
		real_t integral = integrate( &GrowthIntegrand_Matter, 0.0, a, (void*)&m_Cosmology );
		return eta/a*integral;
	}
    
	//! Integrand used by function CalcGrowthFactor to determine the linear growth factor D+
	inline static double GrowthIntegrand_Matter( double a, void *Params )
	{
		Cosmology *cosm = (Cosmology*)Params;
		double eta = sqrt((double)(cosm->Omega_m/a+cosm->Omega_L*a*a
								   +1.0-cosm->Omega_m-cosm->Omega_L));
		return 2.5/(eta*eta*eta);
	}
	
    real_t ComputeVFactnew( real_t a ){
     return Calc_fPeebles(a)*sqrt(m_Cosmology.Omega_m/a+m_Cosmology.Omega_L*a*a+1.0-m_Cosmology.Omega_m-m_Cosmology.Omega_L)*100.0;
    }
     
     real_t ComputedDdt( real_t a )
     {
     return Calc_fPeebles(a);
     }
    
	//! Compute the factor relating particle displacement and velocity
	real_t ComputeVFact( real_t a ){
		real_t fomega, dlogadt, eta;
		real_t Omega_k = 1.0 - m_Cosmology.Omega_m - m_Cosmology.Omega_L;
		
		real_t Dplus = CalcGrowthFactor( a );
		
		eta     = sqrt( (double)(m_Cosmology.Omega_r/a/a+m_Cosmology.Omega_m/a+ m_Cosmology.Omega_L*a*a + Omega_k ));
		fomega  = (2.5/Dplus-1.5*m_Cosmology.Omega_m/a-Omega_k)/eta/eta;
		dlogadt = a*eta;
		
		//... /100.0 since we would have to multiply by H0 to convert
		//... the displacement to velocity units. But displacement is
		//... in Mpc/h, and H0 in units of h is 100.
        
        //std::cerr << "vfact1 = " << fomega * dlogadt/a *100.0 << "\n";
        //std::cerr << "vfact2 = " << ComputeVFactnew(a) << "\n";
        
		return fomega * dlogadt/a *100.0;
	}
	
	
	
	
	//! Integrand for the sigma_8 normalization of the power spectrum
	/*! Returns the value of the primordial power spectrum multiplied with 
	 the transfer function and the window function of 8 Mpc/h at wave number k */
	static double dSigma8( double k, void *Params )
	{
		if( k<=0.0 )
			return 0.0f;
		
		transfer_function *ptf = (transfer_function *)Params;
		
		double x = k*8.0;
		double w = 3.0*(sin(x)-x*cos(x))/(x*x*x);
		static double nspect = (double)ptf->cosmo_.nspect;
		
		double tf = ptf->compute(k, total);
		
		//... no growth factor since we compute at z=0 and normalize so that D+(z=0)=1
		return k*k * w*w * pow((double)k,(double)nspect) * tf*tf;
		
	}
    
    //! Integrand for the sigma_8 normalization of the power spectrum
	/*! Returns the value of the primordial power spectrum multiplied with 
	 the transfer function and the window function of 8 Mpc/h at wave number k */
	static double dSigma8_0( double k, void *Params )
	{
		if( k<=0.0 )
			return 0.0f;
		
		transfer_function *ptf = (transfer_function *)Params;
		
		double x = k*8.0;
		double w = 3.0*(sin(x)-x*cos(x))/(x*x*x);
		static double nspect = (double)ptf->cosmo_.nspect;
		
		double tf = ptf->compute(k, total0);
		
		//... no growth factor since we compute at z=0 and normalize so that D+(z=0)=1
		return k*k * w*w * pow((double)k,(double)nspect) * tf*tf;
		
	}
	
	
	//! Computes the square of the transfer function
	/*! Function evaluates the supplied transfer function m_pTransferFunction
	 * and returns the square of its value at wave number k
	 * @param k wave number at which to evaluate the transfer function
	 */
	inline real_t TransferSq( real_t k ){
		//.. parameter supplied transfer function
		real_t tf1 = m_pTransferFunction->compute(k, total);
		return tf1*tf1;
	}
	
	
	//! Computes the normalization for the power spectrum
	/*!
	 * integrates the power spectrum to fix the normalization to that given
	 * by the sigma_8 parameter
	 */
	real_t ComputePNorm( real_t kmax )
	{
		real_t sigma0, kmin;
		kmax = m_pTransferFunction->get_kmax();//m_Cosmology.H0/8.0;
		kmin = m_pTransferFunction->get_kmin();//0.0;
        
        if( !m_pTransferFunction->tf_has_total0() )
            sigma0 = 4.0 * M_PI * integrate( &dSigma8, (double)kmin, (double)kmax, (void*)m_pTransferFunction );
		else
            sigma0 = 4.0 * M_PI * integrate( &dSigma8_0, (double)kmin, (double)kmax, (void*)m_pTransferFunction );
		
        return m_Cosmology.sigma8*m_Cosmology.sigma8/sigma0;
	}
    
    
	
    
	//! integrand function for comoving line element
	inline static double dline( double a, void *Params )
	{
		Cosmology *cosm = (Cosmology*)Params;
		double y = (cosm->Omega_m + cosm->Omega_L*a*a*a)*a;
		
		return 1./sqrt(y);
	}
	
	//! compute necessary initial velocity kick to prevent motion of zoom region 
	real_t ComputeVelocityCompensation( double astart, double afinal )
	{
		double s = integrate( &dline, astart, afinal, (void*)&m_Cosmology );
		//std::cerr << "s      = " << s << std::endl;
		//std::cerr << "D(z_f) = " << CalcGrowthFactor(afinal) << std::endl;
		//std::cerr << "D(z_i) = " << CalcGrowthFactor(astart) << std::endl;
		return m_Cosmology.H0 * CalcGrowthFactor(afinal)/CalcGrowthFactor(astart)/s;
	}
	
	real_t ComputeVelocityCompensation_2LPT( double astart, double afinal )
	{
		double s = integrate( &dline, astart, afinal, (void*)&m_Cosmology );
		//std::cerr << "s      = " << s << std::endl;
		//std::cerr << "D(z_f) = " << CalcGrowthFactor(afinal) << std::endl;
		//std::cerr << "D(z_i) = " << CalcGrowthFactor(astart) << std::endl;
		return m_Cosmology.H0 * pow(CalcGrowthFactor(afinal)/CalcGrowthFactor(astart),2.0)/s;
	}
	

	
};


//! compute the jeans sound speed
/*! given a density in g/cm^-3 and a mass in g it gives back the sound
 *  speed in cm/s for which the input mass is equal to the jeans mass
 *  @param rho density 
 *  @param mass mass scale
 *  @returns jeans sound speed
 */
inline double jeans_sound_speed( double rho, double mass )
{
	const double G = 6.67e-8;
	return pow( 6.0*mass/M_PI*sqrt(rho)*pow(G,1.5), 1.0/3.0 );
}

//! computes the density from the potential using the Laplacian
void compute_Lu_density( const grid_hierarchy& u, grid_hierarchy& fnew, unsigned order=4 );

//! computes the 2nd order density perturbations using also off-diagonal terms in the potential Hessian 
void compute_LLA_density( const grid_hierarchy& u, grid_hierarchy& fnew, unsigned order=4 );

//! computes the source term for the 2nd order perturbations in the displacements
void compute_2LPT_source( const grid_hierarchy& u, grid_hierarchy& fnew, unsigned order=4 );

void compute_2LPT_source_FFT( config_file& cf_, const grid_hierarchy& u, grid_hierarchy& fnew );


#endif // _COSMOLOGY_HH

