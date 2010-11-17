/*
 
 numerics.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
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
