/*
 
 constraints.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#ifndef __CONSTRAINTS_HH
#define __CONSTRAINTS_HH

#include <vector>
#include <complex>

#include <gsl/gsl_linalg.h>

#include "general.hh"
#include "config_file.hh"
#include "transfer_function.hh"
#include "cosmology.hh"


//! matrix class serving as a gsl wrapper
class matrix
{
protected:
	gsl_matrix * m_;
	//double *data_;
	size_t M_, N_;
	
public:
	matrix( size_t M, size_t N )
	: M_(M), N_(N)
	{
		m_ = gsl_matrix_alloc(M_,N_);
	}
	
	matrix( size_t N )
	: M_(N), N_(N)
	{
		m_ = gsl_matrix_alloc(M_,N_);
	}
	
	matrix( const matrix& o )
	{
		M_ = o.M_;
		N_ = o.N_;
		m_ = gsl_matrix_alloc(M_,N_);
		gsl_matrix_memcpy(m_, o.m_ );
	}
	
	~matrix()
	{
		gsl_matrix_free( m_ );
	}
	
	double& operator()( size_t i, size_t j )
	{	return *gsl_matrix_ptr( m_, i, j );	}
	
	const double& operator()( size_t i, size_t j ) const
	{	return *gsl_matrix_const_ptr( m_, i, j );	}
	
	matrix& operator=( const matrix& o )
	{
		gsl_matrix_free( m_ );
		
		M_ = o.M_;
		N_ = o.N_;
		m_ = gsl_matrix_alloc(M_,N_);
		gsl_matrix_memcpy(m_, o.m_ );
		return *this;
	}
	
	
	matrix& invert()
	{
		if( M_!=N_ )
			throw std::runtime_error("Attempt to invert a non-square matrix!");
		
		int s;
		gsl_matrix* im = gsl_matrix_alloc(M_,N_);
		
		gsl_permutation * p = gsl_permutation_alloc (M_);
		gsl_linalg_LU_decomp( m_, p, &s );
		gsl_linalg_LU_invert( m_, p, im );
		
		gsl_matrix_memcpy(m_, im);
		
		gsl_permutation_free(p);
		gsl_matrix_free(im);
		return *this;
	}
};


//! class to impose constraints on the white noise field (van de Weygaert & Bertschinger 1996)
class constraint_set
{
	
public:
	enum constr_type{ halo, peak };
	
protected:
	
	struct constraint{
		constr_type type;
		double x,y,z;
		double gx,gy,gz;
		double Rg, Rg2;
		double gRg, gRg2;
		double sigma;
	};
	
	config_file *pcf_;
	std::vector<constraint> cset_;
	transfer_function *ptf_;
	CosmoCalc *pccalc_;
	Cosmology *pcosmo_;
	double dplus0_;
	unsigned constr_level_;
	
	
	inline std::complex<double> eval_constr( size_t icon, double kx, double ky, double kz )
	{
		double re, im, kdotx, k2;
		
		kdotx = cset_[icon].gx*kx+cset_[icon].gy*ky+cset_[icon].gz*kz;
		k2    = kx*kx+ky*ky+kz*kz;
		
		re  = im = exp(-k2*cset_[icon].gRg2/2.0);
		re *= cos( kdotx );
		im *= sin( kdotx );
		
		return std::complex<double>(re,im);
	}
	
	
#if defined(FFTW3) && defined(SINGLE_PRECISION)
	
	//! apply constraints to the white noise
	void wnoise_constr_corr( double dx, size_t nx, size_t ny, size_t nz, std::vector<double>& g0, matrix& cinv, fftwf_complex* cw );
	
	//! measure sigma for each constraint in the unconstrained noise
	void wnoise_constr_corr( double dx, fftwf_complex* cw, size_t nx, size_t ny, size_t nz, std::vector<double>& g0 );
	
#else
	//! apply constraints to the white noise
	void wnoise_constr_corr( double dx, size_t nx, size_t ny, size_t nz, std::vector<double>& g0, matrix& cinv, fftw_complex* cw );
	
	//! measure sigma for each constraint in the unconstrained noise
	void wnoise_constr_corr( double dx, fftw_complex* cw, size_t nx, size_t ny, size_t nz, std::vector<double>& g0 );
	
#endif
	
	//! compute the covariance between the constraints
	void icov_constr( double dx, size_t nx, size_t ny, size_t nz, matrix& cij );
	
	
public:
	
	
	//! constructor 
	constraint_set( config_file& cf, transfer_function *ptf );
	
	//! destructor
	~constraint_set()
	{
		delete pccalc_;
		delete pcosmo_;
	}
	
	
	template< typename rng >
	void apply( unsigned ilevel, int x0[], int lx[], rng* wnoise )
	{
		if( cset_.size() == 0 || constr_level_ != ilevel )
			return;
		
		unsigned nlvl = 1<<ilevel;
		double boxlength = pcf_->getValue<double>("setup","boxlength");
		
		//... compute constraint coordinates for grid
		for( size_t i=0; i<cset_.size(); ++i )
		{
			cset_[i].gx = cset_[i].x * (double)nlvl;
			cset_[i].gy = cset_[i].y * (double)nlvl;
			cset_[i].gz = cset_[i].z * (double)nlvl;
			cset_[i].gRg = cset_[i].Rg/boxlength * (double)nlvl;
			cset_[i].gRg2 = cset_[i].gRg*cset_[i].gRg;
			
			if(cset_[i].gRg > 0.5*lx[0])
				LOGWARN("Constraint %d appears to be too large scale",i);
		}
		
		
		std::vector<double> g0;
		
//		unsigned levelmax = pcf_->getValue<unsigned>("setup","levelmax");
		unsigned levelmin = pcf_->getValue<unsigned>("setup","levelmin_TF");
		
		bool bperiodic = ilevel==levelmin;
		double dx = pcf_->getValue<double>("setup","boxlength")/(1<<ilevel);
		

		LOGINFO("Computing constrained realization...");		
		
		if( bperiodic )
		{
			//... we are operating on the periodic coarse grid
			size_t nx = lx[0], ny = lx[1], nz = lx[2], nzp = nz+2;
			fftw_real * w = new fftw_real[nx*ny*nzp];
			
			
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
			fftwf_complex * cw = reinterpret_cast<fftwf_complex*> (w);
			fftwf_plan	p  = fftwf_plan_dft_r2c_3d( nx, ny, nz, w, cw, FFTW_ESTIMATE),
						ip = fftwf_plan_dft_c2r_3d( nx, ny, nz, cw, w, FFTW_ESTIMATE);
	#else
			fftw_complex * cw = reinterpret_cast<fftw_complex*> (w);
			fftw_plan	p  = fftw_plan_dft_r2c_3d( nx, ny, nz, w, cw, FFTW_ESTIMATE),
						ip = fftw_plan_dft_c2r_3d( nx, ny, nz, cw, w, FFTW_ESTIMATE);
	#endif
#else
			fftw_complex * cw = reinterpret_cast<fftw_complex*> (w);
			rfftwnd_plan p	= rfftw3d_create_plan( nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE),
						 ip = rfftw3d_create_plan( nx, ny, nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
#endif
			
			double fftnorm = 1.0/sqrt(nx*ny*nz);
			
			#pragma omp parallel for
			for( int i=0; i<(int)nx; i++ )
				for( int j=0; j<(int)ny; j++ )
					for( int k=0; k<(int)nz; k++ )
					{
						size_t q = ((size_t)i*ny+(size_t)j)*nzp+(size_t)k;
						w[q] = (*wnoise)((x0[0]+i)%nx,(x0[1]+j)%ny,(x0[2]+k)%nz)*fftnorm;
					}
			
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
			fftwf_execute( p );
	#else
			fftw_execute( p );
	#endif
#else
#ifndef SINGLETHREAD_FFTW		
			rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), p, w, NULL );
#else
			rfftwnd_one_real_to_complex( p, w, NULL );
#endif
#endif
			wnoise_constr_corr( dx, cw, nx, ny, nz, g0 );
			
			matrix c(2,2);
			icov_constr( dx, nx, ny, nz, c );
			
			
			wnoise_constr_corr( dx, nx, ny, nz, g0, c, cw );
			
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
			fftwf_execute( ip );
	#else
			fftw_execute( ip );
	#endif
#else
#ifndef SINGLETHREAD_FFTW		
			rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), ip, cw, NULL );
#else
			rfftwnd_one_complex_to_real( ip, cw, NULL );
#endif
#endif
			
			#pragma omp parallel for
			for( int i=0; i<(int)nx; i++ )
				for( int j=0; j<(int)ny; j++ )
					for( int k=0; k<(int)nz; k++ )
					{
						size_t q = ((size_t)i*ny+(size_t)j)*nzp+(size_t)k;
						(*wnoise)((x0[0]+i),(x0[1]+j),(x0[2]+k)) = w[q]*fftnorm;
					}
			
			LOGINFO("Applied constraints to level %d.",ilevel);
			
						
			delete[] w;
			
			
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
			fftwf_destroy_plan(p);
	#else
			fftw_destroy_plan(p);
	#endif
#else
			fftwnd_destroy_plan(p);
#endif
		}else{
			
			//... we are operating on a refinement grid, not necessarily the finest
			
			size_t nx = lx[0], ny = lx[1], nz = lx[2], nzp = nz+2;
			fftw_real * w = new fftw_real[nx*ny*nzp];
			
			
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
			fftwf_complex * cw = reinterpret_cast<fftwf_complex*> (w);
			fftwf_plan	p  = fftwf_plan_dft_r2c_3d( nx, ny, nz, w, cw, FFTW_ESTIMATE),
						ip = fftwf_plan_dft_c2r_3d( nx, ny, nz, cw, w, FFTW_ESTIMATE);
	#else
			fftw_complex * cw = reinterpret_cast<fftw_complex*> (w);
			fftw_plan	p  = fftw_plan_dft_r2c_3d( nx, ny, nz, w, cw, FFTW_ESTIMATE),
						ip = fftw_plan_dft_c2r_3d( nx, ny, nz, cw, w, FFTW_ESTIMATE);
	#endif
#else
			fftw_complex * cw = reinterpret_cast<fftw_complex*> (w);
			rfftwnd_plan p	= rfftw3d_create_plan( nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE),
			ip = rfftw3d_create_plan( nx, ny, nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
#endif
			
			double fftnorm = 1.0/sqrt(nx*ny*nz);
			
			int il = nx/4, ir = 3*nx/4, jl=ny/4, jr = 3*ny/4, kl = nz/4, kr = 3*nz/4;
			
			#pragma omp parallel for
			for( int i=0; i<(int)nx; i++ )
				for( int j=0; j<(int)ny; j++ )
					for( int k=0; k<(int)nz; k++ )
					{
						size_t q = ((size_t)i*ny+(size_t)j)*nzp+(size_t)k;
						
						if( i>=il && i<ir && j>=jl && j<jr && k>=kl && k<kr )
							w[q] = (*wnoise)((x0[0]+i),(x0[1]+j),(x0[2]+k))*fftnorm;
						else
							w[q] = 0.0;

					}
			
			int nlvl05 = 1<<(ilevel-1);
			int xs = nlvl05-x0[0], ys = nlvl05-x0[1], zs = nlvl05-x0[2];
			
			for( size_t i=0; i<cset_.size(); ++i )
			{
				cset_[i].gx -= xs;
				cset_[i].gy -= ys;
				cset_[i].gz -= zs;
			}
			
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
			fftwf_execute( p );
	#else
			fftw_execute( p );
	#endif
#else
#ifndef SINGLETHREAD_FFTW		
			rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), p, w, NULL );
#else
			rfftwnd_one_real_to_complex( p, w, NULL );
#endif
#endif
			wnoise_constr_corr( dx, cw, nx, ny, nz, g0 );
			
			matrix c(2,2);
			icov_constr( dx, nx, ny, nz, c );
			
			
			wnoise_constr_corr( dx, nx, ny, nz, g0, c, cw );
			
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
			fftwf_execute( ip );
	#else
			fftw_execute( ip );
	#endif
#else
#ifndef SINGLETHREAD_FFTW		
			rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), ip, cw, NULL );
#else
			rfftwnd_one_complex_to_real( ip, cw, NULL );
#endif
#endif
			
			#pragma omp parallel for
			for( int i=0; i<(int)nx; i++ )
				for( int j=0; j<(int)ny; j++ )
					for( int k=0; k<(int)nz; k++ )
					{
						size_t q = ((size_t)i*ny+(size_t)j)*nzp+(size_t)k;
						if( i>=il && i<ir && j>=jl && j<jr && k>=kl && k<kr )
							(*wnoise)((x0[0]+i),(x0[1]+j),(x0[2]+k)) = w[q]*fftnorm;
					}
			

			LOGINFO("Applied constraints to level %d.",ilevel);	
			
			delete[] w;
			
			
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
			fftwf_destroy_plan(p);
	#else
			fftw_destroy_plan(p);
	#endif
#else
			fftwnd_destroy_plan(p);
#endif
			
		}
		
	}
	
};


#endif // __CONSTRAINTS_HH
