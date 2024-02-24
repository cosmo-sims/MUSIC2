// This file is part of MUSIC
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2010-2024 by Oliver Hahn
// 
// monofonIC is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// monofonIC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

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

real_t integrate( double (* func) (double x, void * params), double a, double b, void *params )
{
	gsl_function F;
	F.function = func;
	F.params = params;

	double result;
	double error;

	
	gsl_set_error_handler_off ();
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(100000);
	gsl_integration_qag( &F, a, b, 0, REL_PRECISION, 100000, 6, w, &result, &error );
	
	
	gsl_integration_workspace_free(w);

	gsl_set_error_handler(NULL);

	if( error/result > REL_PRECISION )
		std::cerr << " - Warning: no convergence in function 'integrate', rel. error=" << error/result << std::endl;

	return (real_t)result;
}
