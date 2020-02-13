/*
 This file is part of MUSIC -
 a tool to generate initial conditions for cosmological simulations
 
 Copyright (C) 2008-12  Oliver Hahn, ojha@gmx.de
 
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

#ifndef __TRANSFERFUNCTION_HH
#define __TRANSFERFUNCTION_HH

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdexcept>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>

#include "Numerics.hh"
#include "general.hh"

#include <complex>

#define NZERO_Q

typedef std::complex<double> complex;

//! Abstract base class for transfer functions
/*!
    This class implements a purely virtual interface that can be
    used to derive instances implementing various transfer functions.
*/ 
class TransferFunction{
public:
	Cosmology m_Cosmology;
	
public:
  
	TransferFunction( Cosmology acosm ) : m_Cosmology( acosm ) { };
	virtual double compute( double k ) = 0;
	virtual ~TransferFunction(){ };
	virtual double get_kmax( void ) = 0;
	virtual double get_kmin( void ) = 0;
};

class TransferFunction_real
{
	
public:
	gsl_interp_accel *accp, *accn;
	gsl_spline *splinep, *splinen;
	double Tr0_, Tmin_, Tmax_, Tscale_;
	double rneg_, rneg2_;
	static TransferFunction *ptf_;
	static double nspec_;
	
protected:
	
	double krgood( double mu, double q, double dlnr, double kr )
	{
		double krnew = kr;
		complex cdgamma, zm, zp;
		double arg, iarg, xm, xp, y;
		gsl_sf_result g_a, g_p;
		
		xp = 0.5*(mu+1.0+q);
		xm = 0.5*(mu+1.0-q);
		y = M_PI/(2.0*dlnr);
		zp=complex(xp,y);
		zm=complex(xm,y);
		
		gsl_sf_lngamma_complex_e (zp.real(), zp.imag(), &g_a, &g_p);
		zp=std::polar(exp(g_a.val),g_p.val);
		double zpa = g_p.val;

		gsl_sf_lngamma_complex_e (zm.real(), zm.imag(), &g_a, &g_p);
		zm=std::polar(exp(g_a.val),g_p.val);
		double zma = g_p.val;
		
		arg=log(2.0/kr)/dlnr+(zpa+zma)/M_PI;
		iarg=(double)((int)(arg + 0.5));
		
		if( arg!=iarg )
			krnew=kr*exp((arg-iarg)*dlnr);
		
		return krnew;
	}
	
	void transform( double pnorm, double dplus, unsigned N, double q, std::vector<double>& rr, std::vector<double>& TT )
	{
		const double mu = 0.5;
		double qmin = 1.0e-6, qmax = 1.0e+6;
		
		q = 0.0;
		
		N = 16384;
		
#ifdef NZERO_Q
		//q = 0.4;
		q = 0.2;
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
		
		fftw_complex in[N], out[N];
		fftw_plan p,ip;
		
		//... perform anti-ringing correction from Hamilton (2000)
		k0r0 = krgood( mu, q, dlnr, k0r0 );

		std::ofstream ofsk("transfer_k.txt");
		double sum_in = 0.0;
		for( unsigned i=0; i<N; ++i )
		{
			
			double k = k0*exp(((int)i - (int)N/2+1) * dlnk);
			//double k = k0*exp(((int)i - (int)N/2) * dlnk);
			//double k = k0*exp(ii * dlnk);
			
			//... some constants missing ...//
			in[i].re = dplus*sqrtpnorm*ptf_->compute( k )*pow(k,0.5*nspec_)*pow(k,1.5-q);
			in[i].im = 0.0;
			
			sum_in += in[i].re;
			ofsk << std::setw(16) << k <<std::setw(16) << in[i].re << std::endl;
		}
		ofsk.close();
		
		
		p = fftw_create_plan(N, FFTW_FORWARD, FFTW_ESTIMATE);
		ip = fftw_create_plan(N, FFTW_BACKWARD, FFTW_ESTIMATE);
		
		//fftw_one(p, in, out);
		fftw_one(p, in, out);
		
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
			//complex x(dir*q, (double)ii*2.0*M_PI/L);
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
			
			complex cu = complex(out[i].re,out[i].im)*U*phase*fftnorm;
			
			out[i].re = cu.real();
			out[i].im = cu.imag();
			
			if( (out[i].re != out[i].re)||(out[i].im != out[i].im) )
			{	std::cerr << "NaN @ i=" << i << ", U= " << U << ", phase = " << phase << ", g1 = " << g1 << ", g2 = " << g2 << std::endl;
				std::cerr << "mu+1+q = " << mu+1.0+q << std::endl;
				//break;
			}
			
#endif

		}
			
		/*out[N/2].im = 0.0;
		out[N/2+1].im = 0.0;
		out[N/2+1].re = out[N/2].re;
		out[N/2].im = 0.0;*/
		
		fftw_one(ip, out, in);
		
		rr.assign(N,0.0);
		TT.assign(N,0.0);
		
		r0 = k0r0/k0;
		
		for( unsigned i=0; i<N; ++i )
		{
			int ii = i;
			ii -= N/2-1;
			//ii -= N/2;
			//if( ii>N/2)
			//	ii-=N;
			
			
			
			double r = r0*exp(-ii*dlnr);
			rr[N-i-1] = r;
			TT[N-i-1] = 4.0*M_PI* sqrt(M_PI/2.0) *  in[i].re*pow(r,-(1.5+q));
			
			//TT[N-i-1] = 4.0*M_PI* sqrt(M_PI/2.0) *  in[i].re*exp( -dir*(q+1.5)*ii*dlnr +q*log(k0r0))/r0;
			
			//rr[i] = r;
			//TT[i] = 4.0*M_PI* sqrt(M_PI/2.0) *  in[i].re*pow(r,-(1.5+q));
			
		}
		
		
		{
			std::ofstream ofs("transfer_real_new.txt");
			for( unsigned i=0; i<N; ++i )
			{
				int ii = i;
				ii -= N/2-1;
				
				double r = r0*exp(-ii*dlnr);//r0*exp(ii*dlnr);
				double T = 4.0*M_PI* sqrt(M_PI/2.0) *  in[i].re*pow(r,-(1.5+q));
				ofs << r << "\t\t" << T << "\t\t" << in[i].im << std::endl;
			}
		}
		

		fftw_destroy_plan(p);
		fftw_destroy_plan(ip);
	}
	
public:
	TransferFunction_real( TransferFunction *tf, double nspec, double pnorm, double dplus, double rmin, double rmax, double knymax, unsigned nr )
	{
				
		ptf_ = tf;
		nspec_ = nspec;
	
		double q = 0.8;
		
		std::vector<double> r,T,xp,yp,xn,yn;
		
		transform( pnorm, dplus, nr, q, r, T );
		
		//... determine r=0 zero component by integrating up to the Nyquist frequency
		gsl_integration_workspace * wp; 
		gsl_function F;
		wp = gsl_integration_workspace_alloc(20000);
		F.function = &call_wrapper;
		double par[2]; par[0] = dplus*sqrt(pnorm); //par[1] = M_PI/kny;
		F.params = (void*)par;
		double error;
		
		//#warning factor of sqrt(1.5) needs to be adjusted for non-equilateral boxes
		//.. need a factor sqrt( 2*kny^2_x + 2*kny^2_y + 2*kny^2_z )/2 = sqrt(3/2)kny (in equilateral case)
		gsl_integration_qag (&F, 0.0, sqrt(1.5)*knymax, 0, 1e-8, 20000, GSL_INTEG_GAUSS21, wp, &Tr0_, &error); 
		//Tr0_ = 0.0;
		gsl_integration_workspace_free(wp);
				
		
		for( unsigned i=0; i<r.size(); ++i )
		{
			// spline positive and negative part separately
			/*if( T[i] > 0.0 )
			{
				xp.push_back( 2.0*log10(r[i]) );
				yp.push_back( log10(T[i]) );
				rneg_ = r[i];
				rneg2_ = rneg_*rneg_;
			}else {
				xn.push_back( 2.0*log10(r[i]) );
				yn.push_back( log10(-T[i]) );
			}*/
			
			
			if( r[i] > rmin && r[i] < rmax )
			{
				xp.push_back( 2.0*log10(r[i]) );
				yp.push_back( log10(fabs(T[i])) );
				xn.push_back( 2.0*log10(r[i]) );
				if( T[i] >= 0.0 ) 
					yn.push_back( 1.0 );
				else
					yn.push_back( -1.0 );
				
				
				//ofs << std::setw(16) << xp.back() << std::setw(16) << yp.back() << std::endl;
			}
			
		}
		

		
		
		
		accp = gsl_interp_accel_alloc ();
		accn = gsl_interp_accel_alloc ();
		
		//... spline interpolation is only marginally slower here
		splinep = gsl_spline_alloc (gsl_interp_cspline, xp.size() );
		splinen = gsl_spline_alloc (gsl_interp_cspline, xn.size() );

		//... set up everything for spline interpolation
		gsl_spline_init (splinep, &xp[0], &yp[0], xp.size() );
		gsl_spline_init (splinen, &xn[0], &yn[0], xn.size() );		
		

		
		
		{
			double dlogr = (log10(rmax)-log10(rmin))/100;
			std::ofstream ofs("transfer_splinep.txt");			
			
			for( int i=0; i< 100; ++i ) 
			{
				double r = rmin*pow(10.0,i*dlogr);
				ofs << std::setw(16) << r << std::setw(16) << compute_real(r*r) << std::endl;
			}
		}
		
	}
	
	static double call_wrapper( double k, void *arg )
	{
		double *a = (double*)arg;
		return 4.0*M_PI*a[0]*ptf_->compute( k )*pow(k,0.5*nspec_)*k*k;
	}
	
	~TransferFunction_real()
	{
		gsl_spline_free (splinen);
		gsl_interp_accel_free (accn);
		gsl_spline_free (splinep);
		gsl_interp_accel_free (accp);

	}
	
	inline double compute_real( double r2 ) const
	{
		const double EPS = 1e-8;
		const double Reps2 = EPS*EPS;
		
		if( r2 <Reps2 )
			return Tr0_;
		double q;
		/*if( r2 < rneg2_ )
			q = pow(10.0,gsl_spline_eval (splinep, log10(r2), accp));
		else
			q = -pow(10.0,gsl_spline_eval(splinen, log10(r2), accn));*/
		
		double logr2 = log10(r2);
		q = pow(10.0,gsl_spline_eval(splinep, logr2, accp));
		double sign = 1.0;
		if( gsl_spline_eval(splinen, logr2, accn) < 0.0 )
			sign = -1.0;
		return q*sign;
	}
};

#endif
