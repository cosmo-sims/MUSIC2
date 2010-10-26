/*
 
 general.hh - This file is part of MUSIC -
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

#ifndef __GENERAL_HH
#define __GENERAL_HH

#include "log.hh"

#include <omp.h>

#ifdef WITH_MPI
  #ifdef MANNO
    #include <mpi.h>
  #else
    #include <mpi++.h>
  #endif
#else
#include <time.h>
#endif

/*#ifdef SINGLE_PRECISION
 #ifdef WITH_MPI
  #include "srfftw_mpi.h"
  #define MPI_TREAL MPI_FLOAT
 #else
  #include "srfftw.h"
 #endif
  typedef float real_t;
#else
 #ifdef WITH_MPI
  #include "drfftw_mpi.h"
  #define MPI_TREAL MPI_DOUBLE
 #else
  #include "drfftw.h"
 #endif
  typedef double real_t;
#endif*/
#ifdef FFTW3
	#include <fftw3.h>
	#if defined(SINGLE_PRECISION)
	typedef float fftw_real;
	#else
	typedef double fftw_real;
	#endif

#else
	#if defined(SINGLE_PRECISION) and not defined(SINGLETHREAD_FFTW)
	#include <srfftw.h>
	#include <srfftw_threads.h>
	#elif defined(SINGLE_PRECISION) and defined(SINGLETHREAD_FFTW)
	#include <srfftw.h>
	#elif not defined(SINGLE_PRECISION) and not defined(SINGLETHREAD_FFTW)
	#include <drfftw.h>
	#include <drfftw_threads.h>
	#elif not defined(SINGLE_PRECISION) and defined(SINGLETHREAD_FFTW)
	#include <drfftw.h>
	#endif
#endif

#ifdef SINGLE_PRECISION
	typedef float real_t;
#else
	typedef double real_t;
#endif


#include <vector>


#include "mesh.hh"
typedef GridHierarchy<real_t> grid_hierarchy;
typedef MeshvarBnd<real_t> meshvar_bnd;
typedef Meshvar<real_t> meshvar;

#include "random.hh"
typedef random_numbers<real_t> rand_nums;
typedef random_number_generator< rand_nums,real_t> rand_gen;


//! compute square of argument
template< typename T >
inline T SQR( T a ){
  return a*a;
}

//! compute cube of argument
template< typename T >
inline T CUBE( T a ){
  return a*a*a;
}

//! compute 4th power of argument
template< typename T >
inline T POW4( T a ){
	return SQR(SQR(a));
  //return a*a*a*a;
}


//! structure for cosmological parameters
typedef struct cosmology{
  double 
    Omega_m,		//!< baryon+dark matter density
    Omega_b,		//!< baryon matter density
    Omega_L,		//!< dark energy density
    H0,				//!< Hubble constant
    nspect,			//!< long-wave spectral index (scale free is nspect=1) 
    sigma8,			//!< power spectrum normalization
	//Gamma,		//!< shape parameter (of historical interest, as a free parameter)
    //fnl,			//!< non-gaussian contribution parameter
	//w0,			//!< dark energy equation of state parameter (not implemented, i.e. =1 at the moment)
	//wa,			//!< dark energy equation of state parameter (not implemented, i.e. =1 at the moment)
	dplus,			//!< linear perturbation growth factor
	pnorm,			//!< actual power spectrum normalisation factor
	vfact,			//!< velocity<->displacement conversion factor in Zel'dovich approx.
	WDMmass,		//!< Warm DM particle mass
	WDMg_x,			//!< Warm DM particle degrees of freedom
	astart;			//!< expansion factor a for which to generate initial conditions
}Cosmology;

//! basic box/grid/refinement structure parameters
typedef struct {
	unsigned levelmin, levelmax;
	double boxlength;
	std::vector<unsigned> offx,offy,offz,llx,lly,llz;
}Parameters;

//! measure elapsed wallclock time
inline double time_seconds( void )
{
  #ifdef WITH_MPI
    return MPI_Wtime();
  #else
    return ((double) clock()) / CLOCKS_PER_SEC;
  #endif
}

#endif
