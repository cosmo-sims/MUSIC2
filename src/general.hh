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

#if defined(_OPENMP)
#include <omp.h>
inline double get_wtime()
{
  return omp_get_wtime();
}
#else
#include <ctime>
inline double get_wtime()
{
  return std::clock() / double(CLOCKS_PER_SEC);
}
#endif

// //! basic box/grid/refinement structure parameters
// typedef struct {
// 	unsigned levelmin, levelmax;
// 	double boxlength;
// 	std::vector<unsigned> offx,offy,offz,llx,lly,llz;
// }Parameters;

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
