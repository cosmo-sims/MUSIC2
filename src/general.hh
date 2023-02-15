/*
 
 general.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
*/

#pragma once

#include <logger.hh>
#include <config_file.hh>
#include <memory>
#include <cassert>
#include <omp.h>
#include <complex>

#include <fftw3.h>

// include CMake controlled configuration settings
#include "cmake_config.hh"

#if defined(USE_PRECISION_FLOAT)
  using real_t = float;
  using complex_t = fftwf_complex;
  #define FFTW_PREFIX fftwf
#elif defined(USE_PRECISION_DOUBLE)
  using real_t = double;
  using complex_t = fftw_complex;
  #define FFTW_PREFIX fftw
#elif defined(USE_PRECISION_LONGDOUBLE)
  using real_t = long double;
  using complex_t = fftwl_complex;
  #define FFTW_PREFIX fftwl
#endif

using ccomplex_t = std::complex<real_t>;

#define FFTW_GEN_NAME_PRIM(a, b) a##_##b
#define FFTW_GEN_NAME(a, b) FFTW_GEN_NAME_PRIM(a, b)
#define FFTW_API(x) FFTW_GEN_NAME(FFTW_PREFIX, x)

using fftw_plan_t = FFTW_GEN_NAME(FFTW_PREFIX, plan);

#define RE(x) ((x)[0])
#define IM(x) ((x)[1])

#include <vector>
#include <array>

namespace CONFIG
{
// extern int MPI_thread_support;
extern int MPI_task_rank;
extern int MPI_task_size;
extern bool MPI_ok;
// extern bool MPI_threads_ok;
extern bool FFTW_threads_ok;
extern int num_threads;
} // namespace CONFIG

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
/*
typedef struct cosmology_old{
  double 
    Omega_m,		//!< baryon+dark matter density
    Omega_b,		//!< baryon matter density
    Omega_DE,		//!< dark energy density (cosmological constant or parameterised)
    Omega_r,            //!< photon + relativistic particle density
    Omega_k,            //!< curvature density
    H0,			//!< Hubble constant in km/s/Mpc
    nspect,		//!< long-wave spectral index (scale free is nspect=1) 
    sigma8,		//!< power spectrum normalization
    w_0,                //!< dark energy equation of state parameter 1: w = w0 + a * wa
    w_a,                //!< dark energy equation of state parameter 2: w = w0 + a * wa

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
    double zstart = cf.get_value<double>( "setup", "zstart" );
    
    astart     	= 1.0/(1.0+zstart);
    Omega_b    	= cf.get_value<double>( "cosmology", "Omega_b" );
    Omega_m    	= cf.get_value<double>( "cosmology", "Omega_m" );
    Omega_DE    = cf.get_value<double>( "cosmology", "Omega_L" );
    w_0         = cf.get_value_safe<double>( "cosmology", "w0", -1.0 );
    w_a         = cf.get_value_safe<double>( "cosmology", "wa", 0.0 );	
    
    Omega_r     = cf.get_value_safe<double>( "cosmology", "Omega_r", 0.0 ); // no longer default to nonzero (8.3e-5)
    Omega_k     = 1.0 - Omega_m - Omega_DE - Omega_r;

    H0	       	= cf.get_value<double>( "cosmology", "H0" );
    sigma8     	= cf.get_value<double>( "cosmology", "sigma_8" );
    nspect      = cf.get_value<double>( "cosmology", "nspec" );
    WDMg_x     	= cf.get_value_safe<double>( "cosmology", "WDMg_x", 1.5 );
    WDMmass    	= cf.get_value_safe<double>( "cosmology", "WDMmass", 0.0 );
    
    dplus      	= 0.0;
    pnorm      	= 0.0;
    vfact      	= 0.0;
  }
  
  cosmology( void )
  {
    
  }
}Cosmology;*/

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
