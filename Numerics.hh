/*
 
 numerics.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#ifndef __NUMERICS_HH
#define __NUMERICS_HH

#ifdef WITH_MPI
  #ifdef MANNO
    #include <mpi.h>
  #else
    #include <mpi++.h>
  #endif
#endif

#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include <vector>
#include <algorithm>
#include "general.hh"

real_t integrate( double (* func) (double x, void * params), double a, double b, void *params=NULL);

typedef __attribute__((__may_alias__)) int aint;

inline float fast_log2 (float val)
{
	//if( sizeof(int) != sizeof(float) )
	//	throw std::runtime_error("fast_log2 will fail on this system!!");
	aint * const    exp_ptr = reinterpret_cast <aint *> (&val);
	aint            x = *exp_ptr;
	const int      log_2 = ((x >> 23) & 255) - 128;
	x &= ~(255 << 23);
	x += 127 << 23;
	*exp_ptr = x;
	
	val = ((-1.0f/3) * val + 2) * val - 2.0f/3;   // (1)
	
	return (val + log_2);
} 

inline float fast_log (const float &val)
{
	return (fast_log2 (val) * 0.69314718f);
} 

inline float fast_log10 (const float &val)
{
	return (fast_log2 (val) * 0.3010299956639812f);
} 

inline unsigned locate( const double x, const std::vector<double> vx )
{
	long unsigned ju,jm,jl;
	bool ascnd=(vx[vx.size()-1]>=vx[0]);
	jl = 0;
	ju = vx.size()-1;
	while( ju-jl > 1 ) {
		jm = (ju+jl)>>1;
		if( (x >= vx[jm]) == ascnd )
			jl = jm;
		else
			ju = jm;
	}
	return std::max((long unsigned)0,std::min((long unsigned)(vx.size()-2),(long unsigned)jl));
}


inline real_t linint( const double x, const std::vector<double>& xx, const std::vector<double>& yy )
{
	unsigned i = locate(x,xx);

	if( x<xx[0] )
		return yy[0];
	if( x>=xx[xx.size()-1] )
		return yy[yy.size()-1]; 
	double a  = 1.0/(xx[i+1]-xx[i]);
	double dy = (yy[i+1]-yy[i])*a;
	double y0 = (yy[i]*xx[i+1]-xx[i]*yy[i+1])*a;
  return dy*x+y0;
}


#endif


