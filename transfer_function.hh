/*
 
 transfer_function.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
*/

#ifndef __TRANSFERFUNCTION_HH
#define __TRANSFERFUNCTION_HH

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

#include "Numerics.hh"
#include "config_file.hh"

enum tf_type{
	total, cdm, baryon, vcdm, vbaryon
};


//! Abstract base class for transfer functions
/*!
 This class implements a purely virtual interface that can be
 used to derive instances implementing various transfer functions.
 */ 
class transfer_function_plugin{
public:
	Cosmology cosmo_;		//!< cosmological parameter, read from config_file
	config_file *pcf_;		//!< pointer to config_file from which to read parameters
	bool tf_distinct_;		//!< bool if density transfer function is distinct for baryons and DM
	bool tf_withvel_;		//!< bool if also have velocity transfer functions
	
public:
	
	//! constructor
	transfer_function_plugin( config_file& cf ) 
	: pcf_( &cf )
	{
		real_t zstart;
		zstart				= pcf_->getValue<real_t>( "setup", "zstart" );
		cosmo_.astart		= 1.0/(1.0+zstart);
		cosmo_.Omega_b		= pcf_->getValue<real_t>( "cosmology", "Omega_b" );
		cosmo_.Omega_m		= pcf_->getValue<real_t>( "cosmology", "Omega_m" );
		cosmo_.Omega_L		= pcf_->getValue<real_t>( "cosmology", "Omega_L" );
		cosmo_.H0			= pcf_->getValue<real_t>( "cosmology", "H0" );
		cosmo_.sigma8		= pcf_->getValue<real_t>( "cosmology", "sigma_8" );
		cosmo_.nspect		= pcf_->getValue<real_t>( "cosmology", "nspec" );
	}
	
	//! destructor
	virtual ~transfer_function_plugin(){ };
	
	//! compute value of transfer function at waven umber
	virtual double compute( double k, tf_type type) = 0;

	//! return maximum wave number allowed
	virtual double get_kmax( void ) = 0;
	
	//! return minimum wave number allowed
	virtual double get_kmin( void ) = 0;
	
	//! return if density transfer function is distinct for baryons and DM
	bool tf_is_distinct( void )
	{	return tf_distinct_;	}
	
	//! return if we also have velocity transfer functions
	bool tf_has_velocities( void )
	{	return tf_withvel_; }
};


//! Implements abstract factory design pattern for transfer function plug-ins
struct transfer_function_plugin_creator
{
	//! create an instance of a transfer function plug-in
	virtual transfer_function_plugin * create( config_file& cf ) const = 0;
	
	//! destroy an instance of a plug-in 
	virtual ~transfer_function_plugin_creator() { }
};

//! Write names of registered transfer function plug-ins to stdout
std::map< std::string, transfer_function_plugin_creator *>& get_transfer_function_plugin_map();
void print_transfer_function_plugins( void );

//! Concrete factory pattern for transfer function plug-ins
template< class Derived >
struct transfer_function_plugin_creator_concrete : public transfer_function_plugin_creator
{
	//! register the plug-in by its name 
	transfer_function_plugin_creator_concrete( const std::string& plugin_name )
	{
		get_transfer_function_plugin_map()[ plugin_name ] = this;
	}
	
	//! create an instance of the plug-in 
	transfer_function_plugin * create( config_file& cf ) const
	{
		return new Derived( cf );
	}
};

typedef transfer_function_plugin transfer_function;

transfer_function_plugin *select_transfer_function_plugin( config_file& cf );


/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

//! k-space transfer function
class TransferFunction_k
{
public:
	static transfer_function *ptf_;
	static real_t nspec_;
	double dplus_, pnorm_, sqrtpnorm_;
	static tf_type type_;
	
	TransferFunction_k( tf_type type, transfer_function *tf, real_t nspec, real_t pnorm, real_t dplus )
	: dplus_(dplus), pnorm_(pnorm)
	{
		ptf_ = tf;
		nspec_ = nspec;
		sqrtpnorm_ = sqrt( pnorm_ );
		type_ = type;

		std::string fname("input_powerspec_cdm.txt");
		
		if( type == baryon )
			fname = "input_powerspec_baryon.txt";
		if( type == total )
			fname = "input_powerspec_total.txt";
		
		
		if( type == cdm || type == baryon || type == total )
		{
			std::ofstream ofs(fname.c_str());
			double kmin=-3, kmax=3, dk=(kmax-kmin)/100.;
			
			for( int i=0; i<100; ++i )
			{ 
				double k = pow(10.0,kmin+i*dk);
				ofs << std::setw(16) << k 
			    << std::setw(16) << pow(dplus_*sqrtpnorm_*pow(k,0.5*nspec_)*ptf_->compute(k,type_),2)
			    << std::endl;
			}
			
		}


	}
	
	inline real_t compute( real_t k ) const
	{
		return dplus_*sqrtpnorm_*pow(k,0.5*nspec_)*ptf_->compute(k,type_);
	}
};


/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

#define NZERO_Q
typedef std::complex<double> complex;
class TransferFunction_real
{
	
public:
	
	double Tr0_;
	real_t Tmin_, Tmax_, Tscale_;
	real_t rneg_, rneg2_, kny_;
	static transfer_function *ptf_;
	static real_t nspec_;
	
protected:
	
	real_t krgood( real_t mu, real_t q, real_t dlnr, real_t kr )
	{
		double krnew = kr;
		complex cdgamma, zm, zp;
		double arg, iarg, xneg, xpos, y;
		gsl_sf_result g_a, g_p;
		
		xpos = 0.5*(mu+1.0+q);
		xneg = 0.5*(mu+1.0-q);
		y = M_PI/(2.0*dlnr);
		zp=complex(xpos,y);
		zm=complex(xneg,y);
		
		gsl_sf_lngamma_complex_e (zp.real(), zp.imag(), &g_a, &g_p);
		zp=std::polar(exp(g_a.val),g_p.val);
		real_t zpa = g_p.val;

		gsl_sf_lngamma_complex_e (zm.real(), zm.imag(), &g_a, &g_p);
		zm=std::polar(exp(g_a.val),g_p.val);
		real_t zma = g_p.val;
		
		arg=log(2.0/kr)/dlnr+(zpa+zma)/M_PI;
		iarg=(real_t)((int)(arg + 0.5));
		
		if( arg!=iarg )
			krnew=kr*exp((arg-iarg)*dlnr);
		
		return krnew;
	}
	
	void transform( real_t pnorm, real_t dplus, unsigned N, real_t q, std::vector<double>& rr, std::vector<double>& TT )
	{
		const double mu = 0.5;
		double qmin = 1.0e-7, qmax = 1.0e+7;
		
		q = 0.0;
		
		//N = 16384;
		N = 1<<12;
		
#ifdef NZERO_Q
		//q=0.4;
		q=0.4;
#endif
		
		double kmin = qmin, kmax=qmax;
		double rmin = qmin, rmax = qmax;
		double k0 = exp(0.5*(log(kmax)+log(kmin)));
		double r0 = exp(0.5*(log(rmax)+log(rmin)));
		double L = log(rmax)-log(rmin);
		double k0r0 = k0*r0;
		double dlnk = L/N, dlnr = L/N;
		
		double sqrtpnorm = sqrt(pnorm);
		
		double dir = 1.0;
		
		double fftnorm = 1.0/N;
		
		fftw_complex *in, *out;
		
		in = new fftw_complex[N];
		out = new fftw_complex[N];
		
		//... perform anti-ringing correction from Hamilton (2000)
		k0r0 = krgood( mu, q, dlnr, k0r0 );

		std::ofstream ofsk("transfer_k.txt");
		double sum_in = 0.0;
		

		
		for( unsigned i=0; i<N; ++i )
		{
			double k = k0*exp(((int)i - (int)N/2+1) * dlnk);
			double T = ptf_->compute( k, type_ );
#ifdef FFTW3
			
			in[i][0] = dplus*sqrtpnorm*T*pow(k,0.5*nspec_)*pow(k,1.5-q);
			in[i][1] = 0.0;
			
			sum_in += in[i][0];	
			ofsk << std::setw(16) << k <<std::setw(16) << in[i][0] << std::setw(16) << T << std::endl;
#else
			
			in[i].re = dplus*sqrtpnorm*T*pow(k,0.5*nspec_)*pow(k,1.5-q);
			in[i].im = 0.0;
			
			sum_in += in[i].re;
			ofsk << std::setw(16) << k <<std::setw(16) << in[i].re << std::setw(16) << T << std::endl;
#endif
			
		}
		ofsk.close();
		
#ifdef FFTW3
		fftw_plan p,ip;
		p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		ip = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(p);
#else
		fftw_plan p,ip;
		p = fftw_create_plan(N, FFTW_FORWARD, FFTW_ESTIMATE);
		ip = fftw_create_plan(N, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_one(p, in, out);
#endif
		
		//... compute the Hankel transform by convolution with the Bessel function
		for( unsigned i=0; i<N; ++i )
		{
			int ii=i;
			if( ii > (int)N/2 )
				ii -= N;
			
#ifndef NZERO_Q
			double y=ii*M_PI/L;
			complex zp((mu+1.0)*0.5,y);
			gsl_sf_result g_a, g_p;
			gsl_sf_lngamma_complex_e(zp.real(), zp.imag(), &g_a, &g_p);
			
			double arg = 2.0*(log(2.0/k0r0)*y+g_p.val);
			complex cu = complex(out[i].re,out[i].im)*std::polar(1.0,arg);
			out[i].re = cu.real()*fftnorm;
			out[i].im = cu.imag()*fftnorm;
			
#else		
			complex x(dir*q, (double)ii*2.0*M_PI/L);
			gsl_sf_result g_a, g_p;
			
			complex g1, g2, garg, U, phase;						
			complex twotox = pow(complex(2.0,0.0),x);
			
			/////////////////////////////////////////////////////////
			//.. evaluate complex Gamma functions
			
			garg = 0.5*(mu+1.0+x);
			gsl_sf_lngamma_complex_e (garg.real(), garg.imag(), &g_a, &g_p);
			g1 = std::polar(exp(g_a.val),g_p.val);

			
			garg = 0.5*(mu+1.0-x);
			gsl_sf_lngamma_complex_e (garg.real(), garg.imag(), &g_a, &g_p);
			g2 = std::polar(exp(g_a.val),g_p.val);

			/////////////////////////////////////////////////////////
			//.. compute U
			
			if( (fabs(g2.real()) < 1e-19 && fabs(g2.imag()) < 1e-19) )
			{
				//std::cerr << "Warning : encountered possible singularity in TransferFunction_real::transform!\n";
				g1 = 1.0; g2 = 1.0;
			}
			
			
			U = twotox * g1 / g2;
			phase = pow(complex(k0r0,0.0),complex(0.0,2.0*M_PI*(double)ii/L));
			
#ifdef FFTW3
			complex cu = complex(out[i][0],out[i][1])*U*phase*fftnorm;
			
			out[i][0] = cu.real();
			out[i][1] = cu.imag();
#else
			complex cu = complex(out[i].re,out[i].im)*U*phase*fftnorm;
			
			out[i].re = cu.real();
			out[i].im = cu.imag();
#endif
			/*if( (out[i].re != out[i].re)||(out[i].im != out[i].im) )
			{	std::cerr << "NaN @ i=" << i << ", U= " << U << ", phase = " << phase << ", g1 = " << g1 << ", g2 = " << g2 << std::endl;
				std::cerr << "mu+1+q = " << mu+1.0+q << std::endl;
				//break;
			}*/
			
#endif

		}
		
#ifdef FFTW3
		fftw_execute(ip);
#else
		fftw_one(ip, out, in);
#endif
		
		rr.assign(N,0.0);
		TT.assign(N,0.0);
		
		r0 = k0r0/k0;
		
		for( unsigned i=0; i<N; ++i )
		{
			int ii = i;
			ii -= N/2-1;
			double r = r0*exp(-ii*dlnr);
			rr[N-i-1] = r;
#ifdef FFTW3
			TT[N-i-1] = 4.0*M_PI* sqrt(M_PI/2.0) *  in[i][0]*pow(r,-(1.5+q));
#else
			TT[N-i-1] = 4.0*M_PI* sqrt(M_PI/2.0) *  in[i].re*pow(r,-(1.5+q));
#endif
		}
		
		
		
		{
			std::string fname;
			if(type_==total) fname = "transfer_real_total.txt";
			if(type_==cdm) fname = "transfer_real_cdm.txt";
			if(type_==baryon) fname = "transfer_real_baryon.txt";
			if(type_==vcdm) fname = "transfer_real_vcdm.txt";
			if(type_==vbaryon) fname = "transfer_real_vbaryon.txt";
			
			std::ofstream ofs(fname.c_str());//"transfer_real.txt");
							  
			for( unsigned i=0; i<N; ++i )
			{
				int ii = i;
				ii -= N/2-1;
				
				double r = r0*exp(-ii*dlnr);//r0*exp(ii*dlnr);
#ifdef FFTW3
				double T = 4.0*M_PI* sqrt(M_PI/2.0) *  in[i][0]*pow(r,-(1.5+q));
				ofs << r << "\t\t" << T << "\t\t" << in[i][1] << std::endl;				
#else
				double T = 4.0*M_PI* sqrt(M_PI/2.0) *  in[i].re*pow(r,-(1.5+q));
				ofs << r << "\t\t" << T << "\t\t" << in[i].im << std::endl;
#endif
			}
		}
		
		delete[] in;
		delete[] out;
		fftw_destroy_plan(p);
		fftw_destroy_plan(ip);
	}
	std::vector<real_t> m_xtable,m_ytable,m_dytable;
	double m_xmin, m_xmax, m_dx, m_rdx;
	static tf_type type_;
	
	
public:
	
	TransferFunction_real( double boxlength, int nfull, tf_type type, transfer_function *tf, 
						   real_t nspec, real_t pnorm, real_t dplus, real_t rmin, real_t rmax, real_t knymax, unsigned nr )
	{
		real_t q = 0.8;
		
		ptf_	= tf;
		nspec_	= nspec;
		type_	= type;
		kny_	= knymax;
		
		
		/*****************************************************************/
		//... compute the FFTlog transform of the k^n T(k) kernel

		std::vector<double> r,T;
		transform( pnorm, dplus, nr, q, r, T );
		
		gsl_set_error_handler_off ();
		
		
		/*****************************************************************/
		//... compute T(r=0) by 3D k-space integration
		{
			const double REL_PRECISION=1.e-5;
			
			gsl_integration_workspace * ws = gsl_integration_workspace_alloc (1000);
			gsl_function ff;
			
			double kmin = 2.0*M_PI/boxlength;
			double kmax = nfull*M_PI/boxlength;

			//... integrate 0..kmax
			double a[6];
			a[3] = 0.0;
			a[4] = kmax;
			  
			ff.function = &call_x;
			ff.params = reinterpret_cast<void*> (a);
			double res, err, res2, err2;
			
			gsl_integration_qags( &ff, a[3], a[4], 0.0, REL_PRECISION, 1000, ws, &res, &err );
			
			if( err/res > REL_PRECISION )
				std::cerr << " - Warning: no convergence in \'TransferFunction_real\', rel. error=" << err/res << std::endl;
			
			//... integrate 0..kmin
			a[3] = 0.0;
			a[4] = kmin;
			gsl_integration_qags( &ff, a[3], a[4], 0.0, REL_PRECISION, 1000, ws, &res2, &err2 );
			
			if( err2/res2 > 10*REL_PRECISION )
				std::cerr << " - Warning: no convergence in \'TransferFunction_real\', rel. error=" << err2/res2 << std::endl;

			gsl_integration_workspace_free ( ws );

			//.. get kmin..kmax
			res -= res2;
			//.. *8 because we only integrated one octant
			res *= 8.0*dplus*sqrt(pnorm);
			Tr0_ = res;
		}
		
		/*****************************************************************/
		//... store as table for spline interpolation
		
		gsl_interp_accel *accp;
		gsl_spline *splinep;
		
		std::vector<double> xsp, ysp;
		
		for( unsigned i=0; i<r.size(); ++i )
		{
			if( r[i] > rmin && r[i] < rmax )
			{
				xsp.push_back( log10(r[i]) );
				ysp.push_back( T[i]*r[i]*r[i] );
			}
			
		}
		
		accp = gsl_interp_accel_alloc ();

		//... spline interpolation is only marginally slower here
		splinep = gsl_spline_alloc (gsl_interp_akima, xsp.size() );

		//... set up everything for spline interpolation
		gsl_spline_init (splinep, &xsp[0], &ysp[0], xsp.size() );
		
		//.. build lookup table using spline interpolation
		m_xmin = log10(rmin);
		m_xmax = log10(rmax);
		m_dx   = (m_xmax-m_xmin)/nr;
		m_rdx  = 1.0/m_dx;
		
		for(unsigned i=0; i<nr; ++i )
		{
			m_xtable.push_back( m_xmin+i*m_dx );
			m_ytable.push_back( gsl_spline_eval(splinep, (m_xtable.back()), accp) );
		}
		
		for(unsigned i=0; i<nr-1; ++i )
		{
			real_t dy,dr;
			dy = m_ytable[i+1]/pow(m_xtable[i+1],2)-m_ytable[i]/pow(m_xtable[i],2);
			dr = pow(10.0,m_xtable[i+1])-pow(10.0,m_xtable[i]);
			m_dytable.push_back(dy/dr);
		}
		
		gsl_spline_free (splinep);
		gsl_interp_accel_free (accp);
	}
	
	
	static double call_wrapper( double k, void *arg )
	{
		double *a = (double*)arg;

		double T = ptf_->compute( k, type_ );
		
		return 4.0*M_PI*a[0]*T*pow(k,0.5*nspec_)*k*k;
	}
	
	
	
	static double call_x( double kx, void *arg )
	{
		gsl_integration_workspace * wx = gsl_integration_workspace_alloc (1000);
		
		double *a = (double*)arg;
		double kmin = a[3], kmax = a[4];
		
		a[0] = kx;
		
		gsl_function FX;
		FX.function = &call_y;
		FX.params = reinterpret_cast<void*> (a);
		
		double resx, errx;
		gsl_integration_qags( &FX, kmin, kmax, 1e-5, 1e-5, 1000, wx, &resx, &errx );
							 
		gsl_integration_workspace_free (wx);
		
		return resx;
	}
	
	static double call_y( double ky, void *arg )
	{
		gsl_integration_workspace * wy = gsl_integration_workspace_alloc (1000);
		
		double *a = (double*)arg;
		double kmin = a[3], kmax = a[4];
		
		a[1] = ky;
		
		gsl_function FY;
		FY.function = &call_z;
		FY.params = reinterpret_cast<void*> (a);
		
		double resy, erry;
		gsl_integration_qags( &FY, kmin, kmax, 1e-5, 1e-5, 1000, wy, &resy, &erry );
		
		gsl_integration_workspace_free (wy);
		
		return resy;
	}
	
	static double call_z( double kz, void *arg )
	{
		double *a = (double*)arg;
		double kx = a[0], ky = a[1];
		
		double kk = sqrt(kx*kx+ky*ky+kz*kz);
		double T = ptf_->compute( kk, type_ );
		
		return pow(kk,0.5*nspec_)*T;
		
	}
	
	~TransferFunction_real()
	{ }

	inline real_t get_grad( real_t r2 ) const
	{
		double r = 0.5*fast_log10(r2);
		double ii = (r-m_xmin)*m_rdx;
		int i = (int)ii;
		
		i=std::max(0,i);
		i=std::min(i, (int)m_xtable.size()-2);
		
		return (real_t)m_dytable[i];
	}
	
	//... fast version
	inline real_t compute_real( real_t r2 ) const
	{
		const double EPS = 1e-8;
		const double Reps2 = EPS*EPS;
		
		if( r2 <Reps2 )
			return Tr0_;
		
		//double r = 0.5*log10(r2);
		double r = 0.5*fast_log10(r2);
		
		double ii = (r-m_xmin)*m_rdx;
		int i = (int)ii;
		
		i=std::max(0,i);
		i=std::min(i, (int)m_xtable.size()-2);
		
		double y1,y2;
		y1 = m_ytable[i];
		y2 = m_ytable[i+1];
		
		//divide by r**2 because r^2 T is tabulated
		return (real_t)((y1 + (y2-y1)*(ii-(double)i))/r2);
	}
};


#endif
