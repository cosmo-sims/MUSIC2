/*
 
 numerics.cc - This file is part of MUSIC -
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

#ifdef WITH_MPI
  #ifdef MANNO
    #include <mpi.h>
  #else
    #include <mpi++.h>
  #endif
#endif
#include <iostream>
#include "Numerics.hh"


#ifndef REL_PRECISION
#define REL_PRECISION 1.e-5
#endif

int Base_interp::locate(const double x)
{
	int ju,jm,jl; 
	if (n < 2 || mm < 2 || mm > n) throw("locate size error");
	bool ascnd=(xx[n-1] >= xx[0]); 
	jl=0; 
	ju=n-1; 
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1; 
		if ((x >= xx[jm]) == ascnd)
			jl=jm; 
		else
			ju=jm;
	}
	cor = abs(jl-jsav) > dj ? 0 : 1; 
	jsav = jl; 
	return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
}

int Base_interp::hunt(const double x)
{
	int jl=jsav, jm, ju, inc=1; 
	if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
	
	bool ascnd=(xx[n-1] >= xx[0]); 
	
	if (jl < 0 || jl > n-1) {
		jl=0;
		ju=n-1; 
	} else {
		if ((x >= xx[jl]) == ascnd) {
			for (;;) {
				ju = jl + inc; 
				if (ju >= n-1) { 
					ju = n-1; 
					break;
				} else if ((x < xx[ju]) == ascnd) break;
				else {
					jl = ju;
					inc += inc;
				}
			}
		} else {
			ju = jl; 
			for (;;) {
				jl = jl - inc;
				if (jl <= 0) { jl = 0; break;}
				else if ((x >= xx[jl]) == ascnd) break;
				else {
					ju = jl; 
					inc += inc;
				}
			}
		}
	}
	
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1; 
		if ((x >= xx[jm]) == ascnd)
			jl=jm; 
		else
			ju=jm;
		
	}
	
	cor = abs(jl-jsav) > dj ? 0 : 1; 
	jsav = jl; 
	
	return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
}

#if 1
real_t integrate( double (* func) (double x, void * params), double a, double b, void *params )
{
  gsl_function F;
  F.function = func;
  F.params = params;//NULL;

  double result;
  double error;
	//size_t neval;

	
	gsl_set_error_handler_off ();
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(100000);
	gsl_integration_qag( &F, a, b, 0, REL_PRECISION, 100000, 6, w, &result, &error );
	
	//gsl_integration_qags( &F, a, b, 0, REL_PRECISION, 1000, w, &result, &error );
	//gsl_integration_qng( &F, a, b, 0, REL_PRECISION, &result, &error, &neval );

	//gsl_integration_qags( &F, a, b, 0, REL_PRECISION, 10000, w, &result, &error );
	
  gsl_integration_workspace_free(w);

	gsl_set_error_handler(NULL);

  if( error/result > REL_PRECISION )
    std::cerr << " - Warning: no convergence in function 'integrate', rel. error=" << error/result << std::endl;

  return (real_t)result;
}
#else

real_t integrate( double (* func) (double x, void * params), double a, double b, void *params )
{
	unsigned nn = 1000;
	double la = log10(a), lb = log10(b), dlk = (lb-la)/(nn-1);
	
	double sum = 0.0;
	for( unsigned i=1; i<nn; ++i )
	{
		
		double xr = pow(10.0, la+i*dlk );
		double xl = pow(10.0, la+(i-1)*dlk );
		
		
		sum += (xr-xl)*func(0.5*(xl+xr),params);
		
	}
	return sum;
}

#endif
