/*
 
 general.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
*/

#ifndef __GENERAL_HH
#define __GENERAL_HH

#include "log.hh"

#include <cassert>
#include "omp.h"

#ifdef WITH_MPI
  #ifdef MANNO
    #include <mpi.h>
  #else
    #include <mpi++.h>
  #endif
#else
#include <time.h>
#endif

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


#ifdef FFTW3
	#define RE(x) ((x)[0])
	#define IM(x) ((x)[1])
#else
	#define RE(x) ((x).re)
	#define IM(x) ((x).im)
#endif

#if defined(FFTW3) && defined(SINGLE_PRECISION)
#define fftw_complex fftwf_complex
#endif



#include <vector>

#include "config_file.hh"
//#include "mesh.hh"



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
    Omega_r,        //!< photon + relativistic particle density
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
	
	cosmology( config_file cf )
	{
		double zstart = cf.getValue<double>( "setup", "zstart" );
		
		astart		= 1.0/(1.0+zstart);
		Omega_b		= cf.getValue<double>( "cosmology", "Omega_b" );
		Omega_m		= cf.getValue<double>( "cosmology", "Omega_m" );
		Omega_L		= cf.getValue<double>( "cosmology", "Omega_L" );
        Omega_r     = cf.getValueSafe<double>( "cosmology", "Omega_r", 8.3e-5 );
		H0			= cf.getValue<double>( "cosmology", "H0" );
		sigma8		= cf.getValue<double>( "cosmology", "sigma_8" );
		nspect		= cf.getValue<double>( "cosmology", "nspec" );
		WDMg_x		= cf.getValueSafe<double>( "cosmology", "WDMg_x", 1.5 );
		WDMmass		= cf.getValueSafe<double>( "cosmology", "WDMmass", 0.0 );
		
		dplus			= 0.0;
		pnorm			= 0.0;
		vfact			= 0.0;
	}
	
	cosmology( void )
	{
		
	}
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


inline bool is_number(const std::string& s)
{
	for (unsigned i = 0; i < s.length(); i++)
		if (!std::isdigit(s[i])&&s[i]!='-' )
			return false;
	
	return true;
}


#endif
