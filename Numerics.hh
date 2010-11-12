/*
 
 numerics.hh - This file is part of MUSIC -
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



struct Base_interp
{
	int n, mm, jsav, cor, dj; 
	const double *xx, *yy; 
	
	Base_interp(std::vector<double> &x, const double *y, int m)
	: n(x.size()), mm(m), jsav(0), cor(0), dj(0), xx(&x[0]), yy(y) 
	{ 
		dj = std::min(1,(int)pow((double)n,0.25));
	}

	double interp(double x) {
		int jlo = cor ? hunt(x) : locate(x); 
		return rawinterp(jlo,x);
	}
	
	virtual ~Base_interp()
	{ }
	
	int locate(const double x);
	int hunt(const double x);

	double virtual rawinterp(int jlo, double x) = 0;

};




class Spline_interp : public Base_interp
{
protected:
	std::vector<double> y2;
	
	void sety2( const double *xv, const double *yv, double yp1, double ypn)
	{
		int i,k; 
		double p,qn,sig,un; 
		//int n=y2.size(); 
		std::vector<double> u(n-1,0.0); 
		
		if (yp1 > 0.99e99)
			y2[0]=u[0]=0.0; 
		else {
			y2[0] = -0.5; 
			u[0]=(3.0/(xv[1]-xv[0]))*((yv[1]-yv[0])/(xv[1]-xv[0])-yp1);
		} 
		
		for (i=1;i<n-1;i++) {
			sig=(xv[i]-xv[i-1])/(xv[i+1]-xv[i-1]); 
			p=sig*y2[i-1]+2.0; 
			y2[i]=(sig-1.0)/p; 
			u[i]=(yv[i+1]-yv[i])/(xv[i+1]-xv[i]) - (yv[i]-yv[i-1])/(xv[i]-xv[i-1]); 
			u[i]=(6.0*u[i]/(xv[i+1]-xv[i-1])-sig*u[i-1])/p;
		} 
		
		if (ypn > 0.99e99)
			qn=un=0.0; 
		else {
			qn=0.5; 
			un=(3.0/(xv[n-1]-xv[n-2]))*(ypn-(yv[n-1]-yv[n-2])/(xv[n-1]-xv[n-2]));
		} 
		
		y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0); 
		
		for (k=n-2;k>=0;k--)
			y2[k]=y2[k]*y2[k+1]+u[k];
	}
	
	
public:
	
	Spline_interp( std::vector<double> &xv, std::vector<double> &yv, double yp1=1e99, double ypn=1.e99 )
	: Base_interp( xv, &yv[0], xv.size() ), y2( xv.size(), 0.0 )
	{
		sety2( &xv[0], &yv[0], yp1, ypn );
	}
	
	virtual ~Spline_interp()
	{ }

	double rawinterp( int jl, double x )
	{
		int klo=jl,khi=jl+1; 
		double y,h,b,a; 
		
		h=xx[khi]-xx[klo]; 
		
		if (fabs(h) < 1e-10) throw("Bad input to routine splint"); 
		
		a=(xx[khi]-x)/h;
		b=(x-xx[klo])/h;
		y=a*yy[klo]+b*yy[khi]+((a*a*a-a)*y2[klo] +(b*b*b-b)*y2[khi])*(h*h)/6.0; 
		
		return y;
	}

};



#if 1
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
#else
inline unsigned locate( const real_t x, const std::vector<double> vx )
{
	unsigned i=0;
	while( vx[i]<x && (unsigned)i< vx.size() ) ++i;
  if(i>0)i=i-1;
	return i;
}
#endif


inline real_t linint( const double x, const std::vector<double>& xx, const std::vector<double>& yy )
{
	unsigned i = locate(x,xx);
/*	if(i==xx.size()-1) //--i;
		return 0.0;
	if( x<xx[0] )
		return 0.0;*/
	if( x<xx[0] )
		return yy[0];
	if( x>=xx[xx.size()-1] )
		return yy[yy.size()-1]; 
	double a  = 1.0/(xx[i+1]-xx[i]);
	double dy = (yy[i+1]-yy[i])*a;
	double y0 = (yy[i]*xx[i+1]-xx[i]*yy[i+1])*a;
  return dy*x+y0;
}
/*
inline float linint( float x, const std::vector<float>& xx, const std::vector<float>& yy )
{

        int i=0;
        while( xx[i]<x && (unsigned)i< xx.size() ) ++i;
        if(i>0)i=i-1;
        float dy = (yy[i+1]-yy[i])/(xx[i+1]-xx[i]);
        float y0 = (yy[i]*xx[i+1]-xx[i]*yy[i+1])/(xx[i+1]-xx[i]);
        return dy*x+y0;
}
*/


#endif


