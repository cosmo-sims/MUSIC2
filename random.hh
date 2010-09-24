/*
 
 random.hh - This file is part of MUSIC -
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

//#define DEGRADE_RAND1
//#define DEGRADE_RAND2

#ifndef __RANDOM_HH
#define __RANDOM_HH

#include <fstream>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>

#ifdef SINGLE_PRECISION
#ifdef SINGLETHREAD_FFTW
#include <srfftw.h>
#else
#include <srfftw_threads.h>
#endif
#else
#ifdef SINGLETHREAD_FFTW
#include <drfftw.h>
#else
#include <drfftw_threads.h>
#endif
#endif

#include "constraints.hh"
#include "mesh.hh"
#include "mg_operators.hh"

/*!
 * @brief encapsulates all things random number generator related
 */
template< typename T >
class random_numbers
{
public:
	unsigned 
		res_,		//!< resolution of the full mesh
		cubesize_,	//!< size of one independent random number cube
		ncubes_;	//!< number of random number cubes to cover the full mesh
	long baseseed_;	//!< base seed from which cube seeds are computed 
	
	//! vector of 3D meshes (the random number cubes) with random numbers
	std::vector< Meshvar<T>* > rnums_;	
	
protected:
	
	//! fills a subcube with random numbers
	double fill_cube( int i, int j, int k)
	{
		
		gsl_rng	*RNG = gsl_rng_alloc( gsl_rng_mt19937 );
		
		i = (i+ncubes_)%ncubes_;
		j = (j+ncubes_)%ncubes_;
		k = (k+ncubes_)%ncubes_;
		
		long icube = (i*ncubes_+j)*ncubes_+k;
		long cubeseed = baseseed_+icube; //... each cube gets its unique seed
		
		gsl_rng_set( RNG, cubeseed );
		
		if( rnums_[icube] == NULL )
			rnums_[icube] =  new Meshvar<T>( cubesize_, 0, 0, 0 );
		
		double mean = 0.0;
		
		for( int ii=0; ii<(int)cubesize_; ++ii )
			for( int jj=0; jj<(int)cubesize_; ++jj )
				for( int kk=0; kk<(int)cubesize_; ++kk )
				{	
					(*rnums_[icube])(ii,jj,kk) = gsl_ran_ugaussian_ratio_method( RNG );
					mean += (*rnums_[icube])(ii,jj,kk);
				}
		
		gsl_rng_free( RNG );
		
		return mean/(cubesize_*cubesize_*cubesize_);
	}
	
	//! subtract a constant from an entire cube
	void subtract_from_cube( int i, int j, int k, double val )
	{
		i = (i+ncubes_)%ncubes_;
		j = (j+ncubes_)%ncubes_;
		k = (k+ncubes_)%ncubes_;
		
		long icube = (i*ncubes_+j)*ncubes_+k;
		
		for( int ii=0; ii<(int)cubesize_; ++ii )
			for( int jj=0; jj<(int)cubesize_; ++jj )
				for( int kk=0; kk<(int)cubesize_; ++kk )
					(*rnums_[icube])(ii,jj,kk) -= val;
		
	}
	
	//! copy random numbers from a cube to a full grid array
	template< class C >
	void copy_cube( int i, int j, int k, C& dat )
	{
		int offi, offj, offk;
		
		offi = i*cubesize_;
		offj = j*cubesize_;
		offk = k*cubesize_;
		
		i = (i+ncubes_)%ncubes_;
		j = (j+ncubes_)%ncubes_;
		k = (k+ncubes_)%ncubes_;
		
		long icube = (i*ncubes_+j)*ncubes_+k;
		
		for( int ii=0; ii<(int)cubesize_; ++ii )
			for( int jj=0; jj<(int)cubesize_; ++jj )
				for( int kk=0; kk<(int)cubesize_; ++kk )
					dat(offi+ii,offj+jj,offk+kk) = (*rnums_[icube])(ii,jj,kk);
	}
	
	//! free the memory associated with a subcube
	void free_cube( int i, int j, int k ) 
	{

		i = (i+ncubes_)%ncubes_;
		j = (j+ncubes_)%ncubes_;
		k = (k+ncubes_)%ncubes_;
		long icube = (i*ncubes_+j)*ncubes_+k;
		
		if( rnums_[icube] != NULL )
		{
			delete rnums_[icube];
			rnums_[icube] = NULL;
		}
	}
	
	//! initialize member variables and allocate memory
	void initialize( void )
	{
		
		ncubes_ = std::max((int)((double)res_/cubesize_),1);
		if( res_ < cubesize_ )
		{	
			ncubes_ = 1;
			cubesize_ = res_;
		}
		
		std::cout << " - Generating random numbers w/ sample cube size of " << cubesize_ << std::endl;
		
		rnums_.assign( ncubes_*ncubes_*ncubes_, NULL );
	}
	
	//! fill a cubic subvolume of the full grid with random numbers
	double fill_subvolume( int *i0, int *n )
	{
		int i0cube[3], ncube[3];
		
		i0cube[0] = (int)((double)i0[0]/cubesize_);
		i0cube[1] = (int)((double)i0[1]/cubesize_);
		i0cube[2] = (int)((double)i0[2]/cubesize_);
		
		ncube[0] = (int)((double)(i0[0]+n[0])/cubesize_+1.0)-i0cube[0];
		ncube[1] = (int)((double)(i0[1]+n[1])/cubesize_+1.0)-i0cube[1];
		ncube[2] = (int)((double)(i0[2]+n[2])/cubesize_+1.0)-i0cube[2];
		
		double mean = 0.0;
		
		std::cerr << i0[0] << ", " << i0[1] << ", " << i0[2] << "\n";
		std::cerr << n[0] << ", " << n[1] << ", " << n[2] << "\n";
		std::cerr << i0cube[0] << ", " << i0cube[1] << ", " << i0cube[2] << "\n";
		std::cerr << ncube[0] << ", " << ncube[1] << ", " << ncube[2] << "\n";
		
		#pragma omp parallel for reduction(+:mean)
		for( int i=i0cube[0]; i<i0cube[0]+ncube[0]; ++i )
			for( int j=i0cube[1]; j<i0cube[1]+ncube[1]; ++j )
				for( int k=i0cube[2]; k<i0cube[2]+ncube[2]; ++k )
				{
					int ii(i),jj(j),kk(k);
					
					ii = (ii+ncubes_)%ncubes_;
					jj = (jj+ncubes_)%ncubes_;
					kk = (kk+ncubes_)%ncubes_;
					
					mean += fill_cube(ii, jj, kk);
				}
		return mean/(ncube[0]*ncube[1]*ncube[2]);
	}
	
	//! fill an entire grid with random numbers
	double fill_all( void )
	{
		double sum = 0.0;
		
		#pragma omp parallel for reduction(+:sum)
		for( int i=0; i<(int)ncubes_; ++i )
			for( int j=0; j<(int)ncubes_; ++j )
				for( int k=0; k<(int)ncubes_; ++k )
				{
					int ii(i),jj(j),kk(k);
					
					ii = (ii+ncubes_)%ncubes_;
					jj = (jj+ncubes_)%ncubes_;
					kk = (kk+ncubes_)%ncubes_;
					
					sum+=fill_cube(ii, jj, kk);
				}
		
		//... subtract mean
		#pragma omp parallel for reduction(+:sum)
		for( int i=0; i<(int)ncubes_; ++i )
			for( int j=0; j<(int)ncubes_; ++j )
				for( int k=0; k<(int)ncubes_; ++k )
				{
					int ii(i),jj(j),kk(k);
					
					ii = (ii+ncubes_)%ncubes_;
					jj = (jj+ncubes_)%ncubes_;
					kk = (kk+ncubes_)%ncubes_;
					subtract_from_cube(ii,jj,kk,sum/(ncubes_*ncubes_*ncubes_));
				}
		
		/////////////////////////////////////////////////////

#if defined(DEGRADE_RAND1)
		
		{
			std::cerr << " - degrading field for 1 level...(" << res_ << ")\n";
			//unsigned ixoff=51,iyoff=51,izoff=51;
			//unsigned nx=52, ny=52, nz=52;
			unsigned ixoff=102, iyoff=102, izoff=102;
			unsigned nx=104, ny=104, nz=104;
			
#pragma omp parallel for
			for( unsigned ix=0; ix<res_; ix+=2 )
				for( unsigned iy=0; iy<res_; iy+=2 )
					for( unsigned iz=0; iz<res_; iz+=2 )
					{
						if( ix>=2*ixoff && ix < 2*ixoff+nx 
						   && iy>=2*iyoff && iy < 2*iyoff+ny
						   && iz>=2*izoff && iz < 2*izoff+nz )
						{
							continue;
						}
						
						double avg = 0.125*((*this)(ix,iy,iz)+(*this)(ix+1,iy,iz)+(*this)(ix,iy+1,iz)+(*this)(ix,iy,iz+1)+
											(*this)(ix+1,iy+1,iz)+(*this)(ix+1,iy,iz+1)+(*this)(ix,iy+1,iz+1)+(*this)(ix+1,iy+1,iz+1));
						
						
						
						(*this)(ix,iy,iz) = avg;
						(*this)(ix+1,iy,iz) = avg;
						(*this)(ix,iy+1,iz) = avg;
						(*this)(ix,iy,iz+1) = avg;
						(*this)(ix+1,iy+1,iz) = avg;
						(*this)(ix+1,iy,iz+1) = avg;
						(*this)(ix,iy+1,iz+1) = avg;
						(*this)(ix+1,iy+1,iz+1) = avg;
						
					}
			
			
		}
		
#elif defined(DEGRADE_RAND2)
		
		{
			std::cerr << " - degrading field for 2 level...(" << res_ << ")\n";
			
			unsigned ixoff2=102,iyoff2=102,izoff2=102;
			unsigned nx2=52, ny2=52, nz2=52;
			
			unsigned ixoff1=86,iyoff1=86,izoff1=86;
			unsigned nx1=168,ny1=168,nz1=168;

			
			//unsigned ixoff2=51,iyoff2=51,izoff2=51;
			//unsigned nx2=52, ny2=52, nz2=52;
			
			//unsigned ixoff1=34,iyoff1=34,izoff1=34;
			//unsigned nx1=120,ny1=120,nz1=120;
			
			//unsigned ixoff2=84, iyoff2=84, izoff2=84;
			//unsigned nx2=90, ny2=90, nz2=90;
			
			//unsigned ixoff1=100, iyoff1=100, izoff1=100;
			//unsigned nx1=112, ny1=112, nz1=112;
			
			
			/*unsigned ixoff2=100, iyoff2=100, izoff2=100;
			 unsigned nx2=112, ny2=112, nz2=112;
			 
			 unsigned ixoff1=84, iyoff1=84, izoff1=84;
			 unsigned nx1=180, ny1=180, nz1=180;*/
			
#pragma omp parallel for
			for( unsigned ix=0; ix<res_; ix+=2 )
				for( unsigned iy=0; iy<res_; iy+=2 )
					for( unsigned iz=0; iz<res_; iz+=2 )
					{
						if( ix>=2*ixoff2 && ix < 2*ixoff2+nx2 
						   && iy>=2*iyoff2 && iy < 2*iyoff2+ny2
						   && iz>=2*izoff2 && iz < 2*izoff2+nz2 )
						{
							continue;
						}
						
						double avg = 0.125*((*this)(ix,iy,iz)+(*this)(ix+1,iy,iz)+(*this)(ix,iy+1,iz)+(*this)(ix,iy,iz+1)+
											(*this)(ix+1,iy+1,iz)+(*this)(ix+1,iy,iz+1)+(*this)(ix,iy+1,iz+1)+(*this)(ix+1,iy+1,iz+1));
						
						
						
						(*this)(ix,iy,iz) = avg;
						(*this)(ix+1,iy,iz) = avg;
						(*this)(ix,iy+1,iz) = avg;
						(*this)(ix,iy,iz+1) = avg;
						(*this)(ix+1,iy+1,iz) = avg;
						(*this)(ix+1,iy,iz+1) = avg;
						(*this)(ix,iy+1,iz+1) = avg;
						(*this)(ix+1,iy+1,iz+1) = avg;
						
					}
			
#pragma omp parallel for
			for( unsigned ix=0; ix<res_; ix+=4 )
				for( unsigned iy=0; iy<res_; iy+=4 )
					for( unsigned iz=0; iz<res_; iz+=4 )
					{
						if( ix>=2*ixoff1 && ix < 2*ixoff1+nx1
						   && iy>=2*iyoff1 && iy < 2*iyoff1+ny1
						   && iz>=2*izoff1 && iz < 2*izoff1+nz1 )
						{
							continue;
						}
						double avg = 0.0;//0.125*((*this)(ix,iy,iz)+(*this)(ix+1,iy,iz)+(*this)(ix,iy+1,iz)+(*this)(ix,iy,iz+1)+
										 //(*this)(ix+1,iy+1,iz)+(*this)(ix+1,iy,iz+1)+(*this)(ix,iy+1,iz+1)+(*this)(ix+1,iy+1,iz+1));
						for( int i=0; i<4; ++i )
							for( int j=0; j<4; ++j )
								for( int k=0; k<4; ++k )
									avg += (*this)(ix+i,iy+j,iz+k);
						avg /=4.*4.*4.;
						
						for( int i=0; i<4; ++i )
							for( int j=0; j<4; ++j )
								for( int k=0; k<4; ++k )
									(*this)(ix+i,iy+j,iz+k) = avg;
					}
			
			
		}
#endif
		
		
		/////////////////////////////////////////////////////
		
		return sum/(ncubes_*ncubes_*ncubes_);
	}

public:
	template< class C >
	double fill_all( C& dat )
	{
		double sum = 0.0;
		
		std::cerr << "CHECK\n";
		
		
		#pragma omp parallel for reduction(+:sum)
		for( int i=0; i<(int)ncubes_; ++i )
			for( int j=0; j<(int)ncubes_; ++j )
				for( int k=0; k<(int)ncubes_; ++k )
				{
					int ii(i),jj(j),kk(k);
					
					ii = (ii+ncubes_)%ncubes_;
					jj = (jj+ncubes_)%ncubes_;
					kk = (kk+ncubes_)%ncubes_;
					
					sum+=fill_cube(ii, jj, kk);
					copy_cube(ii,jj,kk,dat);
					free_cube(ii, jj, kk);
				}
		
		return sum/(ncubes_*ncubes_*ncubes_);
	}
	
	
	void print_allocated( void )
	{
		unsigned ncount = 0, ntot = rnums_.size();
		for( int i=0; i<rnums_.size(); ++i )
			if( rnums_[i]!=NULL ) ncount++;
		
		std::cerr << " -> " << ncount << " of " << ntot << " random number cubes currently allocated\n";
	}
	
public:
	
	//! constructor
	random_numbers( unsigned res, unsigned cubesize, long baseseed, int *x0, int *lx )
	: res_( res ), cubesize_( cubesize ), ncubes_( 1 ), baseseed_( baseseed )
	{
		std::cout << " - Generating random numbers (1) with seed " << baseseed << std::endl;
		initialize();
		fill_subvolume( x0, lx );
	}
	
	
	//... create constrained new field
	random_numbers( random_numbers<T>& rc, unsigned cubesize, long baseseed, bool kspace=false, int *x0_=NULL, int *lx_=NULL, bool zeromean=true )
	: res_( 2*rc.res_ ), cubesize_( cubesize ), ncubes_( 1 ), baseseed_( baseseed )
	{
		

		//double mean = 0.0;
		
		initialize();
		
		int x0[3],lx[3];
		if( x0_==NULL || lx_==NULL )
		{	
			for(int i=0;i<3;++i ){
				x0[i]=0;
				lx[i]=res_;
			}
			fill_all();
			
		}
		else
		{
			for(int i=0;i<3;++i ){
				x0[i]=x0_[i];
				lx[i]=lx_[i];
			}
			fill_subvolume( x0, lx );
		}

			
			
		if( kspace )
		{
			
			std::cout << " - Generating a constrained random number set with seed " << baseseed << "\n"
					  << "    using coarse mode replacement...\n";
#if 1
			int nx=lx[0],ny=lx[1],nz=lx[2],nxc=lx[0]/2,nyc=lx[1]/2,nzc=lx[2]/2;
			
			//std::cerr << "nx = " << nx << ", nxc = " << nxc << std::endl;
			//std::cerr << "off = " << x0[0] << ", " << x0[1] << ", " << x0[2] << std::endl;
			
			fftw_real 
				*rcoarse = new fftw_real[nxc*nyc*(nzc+2)], 
				*rfine = new fftw_real[nx*ny*(nz+2)];
			
			fftw_complex
				*ccoarse = reinterpret_cast<fftw_complex*> (rcoarse),
				*cfine = reinterpret_cast<fftw_complex*> (rfine);
			
			rfftwnd_plan 
				pc	= rfftw3d_create_plan( nxc, nyc, nzc, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE),
				pf	= rfftw3d_create_plan( nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE),
				ipf	= rfftw3d_create_plan( nx, ny, nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
			
			#pragma omp parallel for
			for( int i=0; i<nx; i++ )
				for( int j=0; j<ny; j++ )
					for( int k=0; k<nz; k++ )
					{
						unsigned q = (i*ny+j)*(nz+2)+k;
						rfine[q] = (*this)(x0[0]+i,x0[1]+j,x0[2]+k);
					}
			
			#pragma omp parallel for
			for( int i=0; i<nxc; i++ )
				for( int j=0; j<nyc; j++ )
					for( int k=0; k<nzc; k++ )
					{
						unsigned q = (i*nyc+j)*(nzc+2)+k;
						rcoarse[q] = rc(x0[0]/2+i,x0[1]/2+j,x0[2]/2+k);
					}
			
		#ifndef SINGLETHREAD_FFTW		
			rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), pc, rcoarse, NULL );
			rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), pf, rfine, NULL );
		#else
			rfftwnd_one_real_to_complex( pc, rcoarse, NULL );
			rfftwnd_one_real_to_complex( pf, rfine, NULL );
		#endif
			
			
			double fftnorm = 1.0/(nx*ny*nz);
			
			#pragma omp parallel for
			for( int i=0; i<nxc; i++ )
				for( int j=0; j<nyc; j++ )
					for( int k=0; k<nzc/2+1; k++ )
					{
						int ii(i),jj(j),kk(k);
						
						if( i > nxc/2 ) ii += nx/2;
						if( j > nyc/2 ) jj += ny/2;
						
						unsigned qc,qf;
						qc = (i*nyc+j)*(nzc/2+1)+k;
						qf = (ii*ny+jj)*(nz/2+1)+kk;
					
						cfine[qf].re = sqrt(8.0)*ccoarse[qc].re;
						cfine[qf].im = sqrt(8.0)*ccoarse[qc].im;
					}
			
			delete[] rcoarse;
			
			#pragma omp parallel for
			for( int i=0; i<nx; i++ )
				for( int j=0; j<ny; j++ )
					for( int k=0; k<nz/2+1; k++ )
					{
						unsigned q = (i*ny+j)*(nz/2+1)+k;
						cfine[q].re *= fftnorm;
						cfine[q].im *= fftnorm;
					}
			
		#ifndef SINGLETHREAD_FFTW		
			rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), ipf, cfine, NULL );
		#else
			rfftwnd_one_complex_to_real( ipf, cfine, NULL );
		#endif
			
			for( int i=0; i<nx; i++ )
				for( int j=0; j<ny; j++ )
					for( int k=0; k<nz; k++ )
					{
						unsigned q = (i*ny+j)*(nz+2)+k;
						(*this)(x0[0]+i,x0[1]+j,x0[2]+k) = rfine[q];
					}
			
			delete[] rfine;
			
			/*for( int i=x0[0],ii=x0[0]/2; ii<lx[0]; i+=2,++ii )
				for( int j=x0[1],jj=x0[1]/2; jj<lx[1]; j+=2,++jj )
					for( int k=x0[2],kk=x0[2]/2; kk<lx[2]; k+=2,++kk )
					{
						double locmean = 0.125*((*this)(i,j,k)+(*this)(i+1,j,k)+(*this)(i,j+1,k)+(*this)(i,j,k+1)+
												(*this)(i+1,j+1,k)+(*this)(i+1,j,k+1)+(*this)(i,j+1,k+1)+(*this)(i+1,j+1,k+1));
						
						rc(ii,jj,kk) = sqrt(8)*locmean;
					}
			 */
#endif
		}
#warning Applying also Hoffman Ribak
		else
		{
			std::cout << " - Generating a constrained random number set with seed " << baseseed << "\n"
			          << "    using Hoffman-Ribak constraints...\n";
			
			double fac = 1./sqrt(8.0);
			
			for( int i=x0[0],ii=x0[0]/2; ii<lx[0]; i+=2,++ii )
				for( int j=x0[1],jj=x0[1]/2; jj<lx[1]; j+=2,++jj )
					for( int k=x0[2],kk=x0[2]/2; kk<lx[2]; k+=2,++kk )
					{
						double topval = rc(ii,jj,kk);
						double locmean = 0.125*((*this)(i,j,k)+(*this)(i+1,j,k)+(*this)(i,j+1,k)+(*this)(i,j,k+1)+
												(*this)(i+1,j+1,k)+(*this)(i+1,j,k+1)+(*this)(i,j+1,k+1)+(*this)(i+1,j+1,k+1));
						double dif = fac*topval-locmean;
						
						(*this)(i,j,k) += dif;
						(*this)(i+1,j,k) += dif;
						(*this)(i,j+1,k) += dif;
						(*this)(i,j,k+1) += dif;
						(*this)(i+1,j+1,k) += dif;
						(*this)(i+1,j,k+1) += dif;
						(*this)(i,j+1,k+1) += dif;
						(*this)(i+1,j+1,k+1) += dif;
						
						
						/*dif = topval;
						
						(*this)(i,j,k) = dif;
						(*this)(i+1,j,k) = dif;
						(*this)(i,j+1,k) = dif;
						(*this)(i,j,k+1) = dif;
						(*this)(i+1,j+1,k) = dif;
						(*this)(i+1,j,k+1) = dif;
						(*this)(i,j+1,k+1) = dif;
						(*this)(i+1,j+1,k+1) = dif;*/
					}
			
			/*for( int i=0; i<res_; ++i )
				for( int j=0; j<res_; ++j )				
					for( int k=0; k<res_; ++k )
					{
						(*this)(i,j,k) *=10;
					}*/
		}
	}
	
	//! constructor
	random_numbers( unsigned res, unsigned cubesize, long baseseed, bool zeromean=true )
	: res_( res ), cubesize_( cubesize ), ncubes_( 1 ), baseseed_( baseseed )
	{
		std::cout << " - Generating random numbers (2) with seed " << baseseed << std::endl;
		
		double mean = 0.0;
		initialize();
		mean = fill_all();
		
		if( true )//zeromean )
		{
			mean = 0.0;
			for(unsigned i=0; i<res_; ++i )
				for( unsigned j=0; j<res_; ++j )
					for( unsigned k=0; k<res_; ++k )
						mean += (*this)(i,j,k);
			mean *= 1.0/(res_*res_*res_);
			
			for(unsigned i=0; i<res_; ++i )
				for( unsigned j=0; j<res_; ++j )
					for( unsigned k=0; k<res_; ++k )
						(*this)(i,j,k) = (*this)(i,j,k) - mean;
		}
		
	}
	
	//! constructor
	random_numbers( unsigned res, std::string randfname )
	: res_( res ), cubesize_( res ), ncubes_(1)
	{
		rnums_.push_back( new Meshvar<T>( res, 0, 0, 0 ) );
		
		std::ifstream ifs(randfname.c_str(), std::ios::binary);
		if( !ifs )
			throw std::runtime_error(std::string("Could not open random number file \'")+randfname+std::string("\'!"));
		
		unsigned vartype;
		unsigned nx,ny,nz,blksz;
		int iseed;
		long seed;
		//ifs.read( (char*)&vartype, sizeof(unsigned) );
		
		//... read header .../
		ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
		
		if( blksz != 4*sizeof(int) )
		   throw std::runtime_error("corrupt random number file");
		   
		ifs.read( reinterpret_cast<char*> (&nx), sizeof(unsigned) );
		ifs.read( reinterpret_cast<char*> (&ny), sizeof(unsigned) );
		ifs.read( reinterpret_cast<char*> (&nz), sizeof(unsigned) );
		ifs.read( reinterpret_cast<char*> (&iseed), sizeof(int) );
		seed = (long)iseed;
		
		if( nx!=res_ || ny!=res_ || nz!=res_ )
		{	
			char errmsg[128];
			sprintf(errmsg,"White noise file dimensions do not match level dimensions: %ux%ux%u vs. %u**3",nx,ny,nz,res_);
			throw std::runtime_error(errmsg);
			
		}
		
		ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
		
		//... read data ...//
		
		ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
		if( blksz == nx*ny*sizeof(float) )
			vartype = 4;
		else if( blksz == nx*ny*sizeof(double) )
			vartype = 8;
		else
			throw std::runtime_error("corrupt random number file");
		
		ifs.seekg(-sizeof(int),std::ios::cur);
			
		std::vector<float> in_float;
		std::vector<double> in_double;
		
		std::cout << " - Random number file \'" << randfname << "\'\n"
				  << "   contains " << nx*ny*nz << " numbers. Reading..." << std::endl;
		
		double sum = 0.0, sum2 = 0.0;
		unsigned count = 0;
		
		if( vartype == 4 )
		{
			for( int ii=0; ii<(int)nz; ++ii )
			{
				ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
				if( blksz != nx*ny*sizeof(float) )
					throw std::runtime_error("corrupt random number file");
				
				in_float.assign(nx*ny,0.0f);
				ifs.read( (char*)&in_float[0], nx*ny*sizeof(float) );
				
				for( int jj=0,q=0; jj<(int)ny; ++jj )
					for( int kk=0; kk<(int)nx; ++kk ){
						sum += in_float[q];
						sum2 += in_float[q]*in_float[q];
						++count;

						(*rnums_[0])(kk,jj,ii) = -in_float[q++];
					}
				ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
				if( blksz != nx*ny*sizeof(float) )
					throw std::runtime_error("corrupt random number file");
				
			}
		}
		else if( vartype == 8 )
		{
			for( int ii=0; ii<(int)nz; ++ii )
			{
				ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
				if( blksz != nx*ny*sizeof(double) )
					throw std::runtime_error("corrupt random number file");
				
				in_double.assign(nx*ny,0.0f);				
				ifs.read( (char*)&in_double[0], nx*ny*sizeof(double) );
				
				for( int jj=0,q=0; jj<(int)ny; ++jj )
					for( int kk=0; kk<(int)nx; ++kk )
					{
						sum += in_double[q];
						sum2 += in_double[q]*in_double[q];
						++count;
						(*rnums_[0])(kk,jj,ii) = -in_double[q++];
					}
				ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
				if( blksz != nx*ny*sizeof(double) )
					throw std::runtime_error("corrupt random number file");
				
			}
		}
		
		double mean, var;
		mean = sum/count;
		var = sum2/count-mean*mean;
		
		std::cout << " - Random numbers in file have \n"
				  << "     mean = " << mean << " and var = " << var << std::endl;
		
	}
	
	//... create constrained new field
	/*random_numbers( random_numbers<T>& rc, unsigned cubesize, long baseseed, bool kspace=false, bool zeromean=true )
	: res_( 2*rc.res_ ), cubesize_( cubesize ), ncubes_( 1 ), baseseed_( baseseed )
	{
		double mean = 0.0;
		initialize();
		mean = fill_all();
		
		if( zeromean )
		{
			for(unsigned i=0; i<res_; ++i )
				for( unsigned j=0; j<res_; ++j )
					for( unsigned k=0; k<res_; ++k )
						(*this)(i,j,k) -= mean;
		}
		
		if( kspace )
		{
		}
		else
		{
			std::cout << " - Generating a constrained random number set with seed " << baseseed << "...\n";
			
			double fac = 1./sqrt(8.0);
			
			for( unsigned i=0,ii=0; i<res_; i+=2,++ii )
				for( unsigned j=0,jj=0; j<res_; j+=2,++jj )
					for( unsigned k=0,kk=0; k<res_; k+=2,++kk )
					{
						double topval = rc(ii,jj,kk);
						double locmean = 0.125*((*this)(i,j,k)+(*this)(i+1,j,k)+(*this)(i,j+1,k)+(*this)(i,j,k+1)+
												(*this)(i+1,j+1,k)+(*this)(i+1,j,k+1)+(*this)(i,j+1,k+1)+(*this)(i+1,j+1,k+1));
						double dif = fac*topval-locmean;
						
						(*this)(i,j,k) += dif;
						(*this)(i+1,j,k) += dif;
						(*this)(i,j+1,k) += dif;
						(*this)(i,j,k+1) += dif;
						(*this)(i+1,j+1,k) += dif;
						(*this)(i+1,j,k+1) += dif;
						(*this)(i,j+1,k+1) += dif;
						(*this)(i+1,j+1,k+1) += dif;
						
					}
			
		}
	}*/

	//... copy construct by averaging down
	explicit random_numbers( /*const*/ random_numbers <T>& rc, bool kdegrade = true )
	{
		//if( res > rc.m_res || (res/rc.m_res)%2 != 0 )
		//			throw std::runtime_error("Invalid restriction in random number container copy constructor.");
		
		double sum = 0.0, sum2 = 0.0;
		unsigned count = 0;
		
		if( kdegrade )
		{
			std::cout << " - Generating a coarse white noise field by k-space degrading\n";
			//... initialize properties of container		
			res_		= rc.res_/2;
			cubesize_	= res_;
			ncubes_		= 1;
			baseseed_	= -2;
			
			if( sizeof(fftw_real)!=sizeof(T) )
				throw std::runtime_error("type mismatch with fftw_real in k-space averaging");
			
			fftw_real 
				*rfine = new fftw_real[rc.res_*rc.res_*2*(rc.res_/2+1)],
				*rcoarse = new fftw_real[res_*res_*2*(res_/2+1)];
			
			fftw_complex
				*ccoarse = reinterpret_cast<fftw_complex*> (rcoarse),
				*cfine = reinterpret_cast<fftw_complex*> (rfine);
			
			int nx(rc.res_), ny(rc.res_), nz(rc.res_), nxc(res_), nyc(res_), nzc(res_);
			
			rfftwnd_plan 
				pf	= rfftw3d_create_plan( nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE),
				ipc	= rfftw3d_create_plan( nxc, nyc, nzc, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
			
			for( int i=0; i<nx; i++ )
				for( int j=0; j<ny; j++ )
					for( int k=0; k<nz; k++ )
					{
						unsigned q = (i*ny+j)*(nz+2)+k;
						rfine[q] = rc(i,j,k);
					}
			
			
#ifndef SINGLETHREAD_FFTW		
			rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), pf, rfine, NULL );
#else
			rfftwnd_one_real_to_complex( pf, rfine, NULL );
#endif
			
			double fftnorm = 1.0/(nxc*nyc*nzc);
			
			for( int i=0; i<nxc; i++ )
				for( int j=0; j<nyc; j++ )
					for( int k=0; k<nzc/2+1; k++ )
					{
						int ii(i),jj(j),kk(k);
						
						if( i > nxc/2 ) ii += nx/2;
						if( j > nyc/2 ) jj += ny/2;
						
						unsigned qc,qf;
						
						qc = (i*nyc+j)*(nzc/2+1)+k;
						qf = (ii*ny+jj)*(nz/2+1)+kk;
						
						ccoarse[qc].re = 1.0/sqrt(8.0)*cfine[qf].re*fftnorm;
						ccoarse[qc].im = 1.0/sqrt(8.0)*cfine[qf].im*fftnorm;
					}
			
			delete[] rfine;
			
#ifndef SINGLETHREAD_FFTW		
			rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), ipc, ccoarse, NULL );
#else
			rfftwnd_one_complex_to_real( ipf, cfine, NULL );
#endif
			rnums_.push_back( new Meshvar<T>( res_, 0, 0, 0 ) );
			
			for( int i=0; i<nxc; i++ )
				for( int j=0; j<nyc; j++ )
					for( int k=0; k<nzc; k++ )
					{
						unsigned q = (i*nyc+j)*(nzc+2)+k;
						(*rnums_[0])(i,j,k) = rcoarse[q];
						sum += (*rnums_[0])(i,j,k);
						sum2+= (*rnums_[0])(i,j,k) * (*rnums_[0])(i,j,k);
						++count;
					}

			delete[] rcoarse;
			
			rfftwnd_destroy_plan(pf);
			rfftwnd_destroy_plan(ipc);
			
		}
		else
		{
			std::cout << " - Generating a coarse white noise field by averaging\n";
			if( rc.rnums_.size() == 1 )
			{
				//... initialize properties of container
				res_		= rc.res_/2;
				cubesize_	= res_;
				ncubes_		= 1;
				baseseed_	= -2;
				
				//... use restriction to get consistent random numbers on coarser grid
				mg_straight gop;
				rnums_.push_back( new Meshvar<T>( res_, 0, 0, 0 ) );		
				gop.restrict( *rc.rnums_[0], *rnums_[0] );
				
				
				for( unsigned i=0; i< rnums_[0]->size(0); ++i )
					for( unsigned j=0; j< rnums_[0]->size(1); ++j )
						for( unsigned k=0; k< rnums_[0]->size(2); ++k )
						{
							(*rnums_[0])(i,j,k) *= sqrt(8); //.. maintain that var(delta)=1
							sum += (*rnums_[0])(i,j,k);
							sum2+= (*rnums_[0])(i,j,k) * (*rnums_[0])(i,j,k);
							++count;
						}
			}
			else 
			{
				//... initialize properties of container
				res_		= rc.res_/2;
				cubesize_	= res_;
				ncubes_		= 1;
				baseseed_	= -2;
				
				rnums_.push_back( new Meshvar<T>( res_, 0, 0, 0 ) );
				double fac = 1.0/sqrt(8);
				
				for( unsigned i=0,ii=0; i<rc.res_; i+=2,++ii )
					for( unsigned j=0,jj=0; j<rc.res_; j+=2,++jj )
						for( unsigned k=0,kk=0; k<rc.res_; k+=2,++kk )
						{
							(*rnums_[0])(ii,jj,kk) = fac * 
							( rc(i,j,k)+rc(i+1,j,k)+rc(i,j+1,k)+rc(i,j,k+1)+
							 rc(i+1,j+1,k)+rc(i+1,j,k+1)+rc(i,j+1,k+1)+rc(i+1,j+1,k+1));
							
							sum += (*rnums_[0])(ii,jj,kk);
							sum2+= (*rnums_[0])(ii,jj,kk) * (*rnums_[0])(ii,jj,kk);
							++count;
						}
				
			}
		}
		
		double rmean, rvar;
		rmean = sum/count;
		rvar = sum2/count-rmean*rmean;
		
		std::cout << " - Restricted random numbers have\n"
				  << "       mean = " << rmean << ", var = " << rvar << std::endl;
	}
	
	
	~random_numbers()
	{
		for( unsigned i=0; i<rnums_.size(); ++i )
			if( rnums_[i] != NULL )
				delete rnums_[i];
		rnums_.clear();
	}
	
	T& operator()( int i, int j, int k )
	{
		int ic, jc, kc, is, js, ks;
		
		if( ncubes_ == 0 )
			throw std::runtime_error("random_numbers: internal error, not properly initialized");
		
		//... determine cube
		ic = (int)((double)i/cubesize_ + ncubes_) % ncubes_;
		jc = (int)((double)j/cubesize_ + ncubes_) % ncubes_;
		kc = (int)((double)k/cubesize_ + ncubes_) % ncubes_;
		
		long icube = (ic*ncubes_+jc)*ncubes_+kc;
		
		if( rnums_[ icube ] == NULL )
		{	
			//... cube has not been precomputed. fill now with random numbers
			//std::cerr << "*";
			rnums_[ icube ] = new Meshvar<T>( cubesize_, 0, 0, 0 );
			fill_cube(ic, jc, kc);
		}
		
		//... determine cell in cube
		is = (i - ic * cubesize_ + cubesize_) % cubesize_;
		js = (j - jc * cubesize_ + cubesize_) % cubesize_;
		ks = (k - kc * cubesize_ + cubesize_) % cubesize_;
		
		return (*rnums_[ icube ])(is,js,ks);
	}
	
	
};

#define DEF_RAN_CUBE_SIZE	32

template< typename rng, typename T >
class random_number_generator
{
protected:
	config_file * pcf_;
	refinement_hierarchy * prefh_;
	constraint_set	constraints;
	
	int levelmin_, levelmax_;
	int levelmin_seed_;
	std::vector<long>			rngseeds_;
	std::vector<std::string>	rngfnames_;
	
	unsigned					ran_cube_size_;
	
	bool is_number(const std::string& s)
	{
		for (unsigned i = 0; i < s.length(); i++)
			if (!std::isdigit(s[i])&&s[i]!='-' )
				return false;
		
		return true;
	}
	
	void parse_rand_parameters( void )
	{
		//... parse random number options
		for( int i=0; i<=100; ++i )
		{
			char seedstr[128];
			std::string tempstr;
			sprintf(seedstr,"seed[%d]",i);
			if( pcf_->containsKey( "random", seedstr ) )
				tempstr = pcf_->getValue<std::string>( "random", seedstr );
			else
				tempstr = std::string("-2");
			
			if( is_number( tempstr ) )
			{	
				long ltemp;
				pcf_->convert( tempstr, ltemp );
				rngfnames_.push_back( "" );
				rngseeds_.push_back( ltemp );
			}else{
				rngfnames_.push_back( tempstr );
				rngseeds_.push_back(-1);
				std::cout << " - Random numbers for level " << std::setw(3) << i << " will be read from file.\n";
			}
			
		}
		
		//.. determine for which levels random seeds/random number files are given
		levelmin_seed_ = -1;
		for( unsigned ilevel = 0; ilevel < rngseeds_.size(); ++ilevel )
		{	
			if( levelmin_seed_ < 0 && (rngfnames_[ilevel].size() > 0 || rngseeds_[ilevel] > 0) )
				levelmin_seed_ = ilevel;
		}
		
	}
	
	void correct_avg( int icoarse, int ifine )
	{
		int shift[3], levelmin_poisson;
		shift[0] = pcf_->getValue<int>("setup","shift_x");
		shift[1] = pcf_->getValue<int>("setup","shift_y");
		shift[2] = pcf_->getValue<int>("setup","shift_z");
		
		levelmin_poisson = pcf_->getValue<unsigned>("setup","levelmin");
		
		int lfacc = 1<<(icoarse-levelmin_poisson);
		
		char fncoarse[128], fnfine[128];
		sprintf(fncoarse,"wnoise_%04d.bin",icoarse);
		sprintf(fnfine,"wnoise_%04d.bin",ifine);
		
		std::ifstream 
			iffine( fnfine, std::ios::binary ), 
			ifcoarse( fncoarse, std::ios::binary );
		
		int nc[3], i0c[3], nf[3], i0f[3];
		if( icoarse != levelmin_ )
		{
			nc[0]  = 2*prefh_->size(icoarse, 0);
			nc[1]  = 2*prefh_->size(icoarse, 1);
			nc[2]  = 2*prefh_->size(icoarse, 2);
			i0c[0] = prefh_->offset_abs(icoarse, 0) - lfacc*shift[0] - nc[0]/4;
			i0c[1] = prefh_->offset_abs(icoarse, 1) - lfacc*shift[1] - nc[1]/4;
			i0c[2] = prefh_->offset_abs(icoarse, 2) - lfacc*shift[2] - nc[2]/4;
			
		}
		else
		{
			nc[0]  = prefh_->size(icoarse, 0);
			nc[1]  = prefh_->size(icoarse, 1);
			nc[2]  = prefh_->size(icoarse, 2);
			i0c[0] = - lfacc*shift[0];
			i0c[1] = - lfacc*shift[1];
			i0c[2] = - lfacc*shift[2];
		}
		nf[0]  = 2*prefh_->size(ifine, 0);
		nf[1]  = 2*prefh_->size(ifine, 1);
		nf[2]  = 2*prefh_->size(ifine, 2);
		i0f[0] = prefh_->offset_abs(ifine, 0) - 2*lfacc*shift[0] - nf[0]/4;
		i0f[1] = prefh_->offset_abs(ifine, 1) - 2*lfacc*shift[1] - nf[1]/4;
		i0f[2] = prefh_->offset_abs(ifine, 2) - 2*lfacc*shift[2] - nf[2]/4;
		
		//.................................
		
		int nxc,nyc,nzc,nxf,nyf,nzf;
		iffine.read( reinterpret_cast<char*> (&nxf), sizeof(unsigned) );
		iffine.read( reinterpret_cast<char*> (&nyf), sizeof(unsigned) );
		iffine.read( reinterpret_cast<char*> (&nzf), sizeof(unsigned) );
		
		ifcoarse.read( reinterpret_cast<char*> (&nxc), sizeof(unsigned) );
		ifcoarse.read( reinterpret_cast<char*> (&nyc), sizeof(unsigned) );
		ifcoarse.read( reinterpret_cast<char*> (&nzc), sizeof(unsigned) );
		
		if( nxf!=nf[0] || nyf!=nf[1] || nzf!=nf[2] || nxc!=nc[0] || nyc!=nc[1] || nzc!=nc[2] )
			throw std::runtime_error("White noise file mismatch. This should not happen. Notify a developer!");
		
		int nxd(nxf/2),nyd(nyf/2),nzd(nzf/2);
		std::vector<T> deg_rand( nxd*nyd*nzd, 0.0 );
		double fac = 1.0/sqrt(8.0);
		for( int i=0; i<nxf; i+=2 )
		{	
			std::vector<T> fine_rand( 2*nyf*nzf, 0.0 );
			iffine.read( reinterpret_cast<char*> (&fine_rand[0]), 2*nyf*nzf*sizeof(T) );
			
			for( int j=0; j<nyf; j+=2 )
				for( int k=0; k<nzf; k+=2 )
				{
					unsigned qc = ((i/2)*nyd+(j/2))*nzd+(k/2);
					unsigned qf[8];
					qf[0] = (0*nyf+j+0)*nzf+k+0;
					qf[1] = (0*nyf+j+0)*nzf+k+1;
					qf[2] = (0*nyf+j+1)*nzf+k+0;
					qf[3] = (0*nyf+j+1)*nzf+k+1;
					qf[4] = (1*nyf+j+0)*nzf+k+0;
					qf[5] = (1*nyf+j+0)*nzf+k+1;
					qf[6] = (1*nyf+j+1)*nzf+k+0;
					qf[7] = (1*nyf+j+1)*nzf+k+1;
					
					for( int q=0; q<8; ++q )
						deg_rand[qc] += fac*fine_rand[qf[q]];
					
				}
		}
		
		//... now deg_rand holds the oct-averaged fine field, store this in the coarse field
		std::vector<T> coarse_rand(nxc*nyc*nzc,0.0);
		ifcoarse.read( reinterpret_cast<char*> (&coarse_rand[0]), nxc*nyc*nzc*sizeof(T) );
		
		int di,dj,dk;
		
		di = i0f[0]/2-i0c[0];
		dj = i0f[1]/2-i0c[1];
		dk = i0f[2]/2-i0c[2];
		
		for( int i=0; i<nxd; i++ )
			for( int j=0; j<nyd; j++ )
				for( int k=0; k<nxd; k++ )
				{
					unsigned qc = ((i+di)*nyc+(j+dj))*nzc+(k+dk);
					unsigned qcd = (i*nyd+j)*nzd+k;
					
					coarse_rand[qc] = deg_rand[qcd];
				}
		
		deg_rand.clear();
		
		ifcoarse.close();
		std::ofstream ofcoarse( fncoarse, std::ios::binary|std::ios::trunc );
		ofcoarse.write( reinterpret_cast<char*> (&nxc), sizeof(unsigned) );
		ofcoarse.write( reinterpret_cast<char*> (&nyc), sizeof(unsigned) );
		ofcoarse.write( reinterpret_cast<char*> (&nzc), sizeof(unsigned) );
		ofcoarse.write( reinterpret_cast<char*> (&coarse_rand[0]), nxc*nyc*nzc*sizeof(T) );
		ofcoarse.close();
	}
	
	void store_rnd( int ilevel, rng* prng )
	{
		int shift[3], levelmin_poisson;
		shift[0] = pcf_->getValue<int>("setup","shift_x");
		shift[1] = pcf_->getValue<int>("setup","shift_y");
		shift[2] = pcf_->getValue<int>("setup","shift_z");
		
		levelmin_poisson = pcf_->getValue<unsigned>("setup","levelmin");
		
		int lfac = 1<<(ilevel-levelmin_poisson);
		
		
		std::vector<T> data;
		if( ilevel == levelmin_ )
		{
			int N = 1<<levelmin_;
			int i0,j0,k0;
			
			i0 = -lfac*shift[0];
			j0 = -lfac*shift[1];
			k0 = -lfac*shift[2];
			
			char fname[128];
			sprintf(fname,"wnoise_%04d.bin",ilevel);

			std::ofstream ofs(fname,std::ios::binary|std::ios::trunc);
			
			ofs.write( reinterpret_cast<char*> (&N), sizeof(unsigned) );
			ofs.write( reinterpret_cast<char*> (&N), sizeof(unsigned) );
			ofs.write( reinterpret_cast<char*> (&N), sizeof(unsigned) );
			
			data.assign( N*N, 0.0 );
			for( int i=0; i<N; ++i )
			{	
				for( int j=0; j<N; ++j )
					for( int k=0; k<N; ++k )
						data[j*N+k] = (*prng)(i+i0,j+j0,k+k0);
				
				ofs.write(reinterpret_cast<char*> (&data[0]), N*N*sizeof(T) );
			}
			ofs.close();
		}
		else
		{
			int nx,ny,nz;
			int i0,j0,k0;
			
			nx = 2*prefh_->size(ilevel, 0);
			ny = 2*prefh_->size(ilevel, 1);
			nz = 2*prefh_->size(ilevel, 2);
			i0 = prefh_->offset_abs(ilevel, 0) - lfac*shift[0] - nx/4;
			j0 = prefh_->offset_abs(ilevel, 1) - lfac*shift[1] - nx/4;
			k0 = prefh_->offset_abs(ilevel, 2) - lfac*shift[2] - nx/4;
			
			char fname[128];
			sprintf(fname,"wnoise_%04d.bin",ilevel);
			
			std::ofstream ofs(fname,std::ios::binary|std::ios::trunc);
			
			ofs.write( reinterpret_cast<char*> (&nx), sizeof(unsigned) );
			ofs.write( reinterpret_cast<char*> (&ny), sizeof(unsigned) );
			ofs.write( reinterpret_cast<char*> (&nz), sizeof(unsigned) );
			
			data.assign( ny*nz, 0.0 );
			for( int i=0; i<nx; ++i )
			{	
				for( int j=0; j<ny; ++j )
					for( int k=0; k<nz; ++k )
						data[j*ny+k] = (*prng)(i+i0,j+j0,k+k0);
				
				ofs.write(reinterpret_cast<char*> (&data[0]), ny*nz*sizeof(T) );
			}
			ofs.close();
			
		}

	}
	
	void compute_random_numbers( void )
	{
		bool kavg = pcf_->getValueSafe<bool>("random","kaveraging",true);
		
		std::vector< rng* > randc(std::max(levelmax_+1,levelmin_seed_),(rng*)NULL);
		
		//--- FILL ALL WHITE NOISE ARRAYS FOR WHICH WE NEED THE FULL FIELD ---//
		
		//... seeds are given for a level coarser than levelmin
		if( levelmin_seed_ < levelmin_ )
		{
			if( rngfnames_[levelmin_seed_].size() > 0 )
				randc[levelmin_seed_] 
					= new rng( 1<<levelmin_seed_, rngfnames_[levelmin_seed_] );
			else
				randc[levelmin_seed_]
					= new rng( 1<<levelmin_seed_, ran_cube_size_, rngseeds_[levelmin_seed_], true );
			
			for( int i=levelmin_seed_+1; i<=levelmin_; ++i )
			{
#warning add possibility to read noise from file also here!
				randc[i] = new rng( *randc[i-1], ran_cube_size_, rngseeds_[i], kavg );
				delete randc[i-1];
				randc[i-1] = NULL;
			}
		}
		
		//... seeds are given for a level finer than levelmin, obtain by averaging
		if( levelmin_seed_ > levelmin_ )
		{
			randc[levelmin_seed_] = new rng( 1<<levelmin_seed_, ran_cube_size_, rngseeds_[levelmin_seed_], true );//, x0, lx );
			
			for( int ilevel = levelmin_seed_-1; ilevel >= (int)levelmin_; --ilevel ){
				if( rngseeds_[ilevel-levelmin_] > 0 )
					std::cerr << " - Warning: random seed for level " << ilevel << " will be ignored.\n"
					<< "            consistency requires that it is obtained by restriction from level " << levelmin_seed_ << std::endl;
				
			
				
				if( levelmin_ == levelmax_ )
					randc[ilevel] = new rng( *randc[ilevel+1], kavg );
				else
					randc[ilevel] = new rng( *randc[ilevel+1], false );
					
				if( ilevel+1 > levelmax_ )
				{
					delete randc[ilevel+1];
					randc[ilevel+1] = NULL;
				}
			}
			
		}
		
		//--- GENERATE AND STORE ALL LEVELS, INCLUDING REFINEMENTS ---//
		
		//... levelmin
		if( randc[levelmin_] == NULL )
			if( rngfnames_[levelmin_].size() > 0 )
				randc[levelmin_] = new rng( 1<<levelmin_, rngfnames_[levelmin_] );
			else
				randc[levelmin_] = new rng( 1<<levelmin_, ran_cube_size_, rngseeds_[levelmin_], true );
		
		store_rnd( levelmin_, randc[levelmin_] );
		
		
		
		//... refinement levels
		for( int ilevel=levelmin_+1; ilevel<=levelmax_; ++ilevel )
		{
			int lx[3], x0[3];
			int shift[3], levelmin_poisson;
			shift[0] = pcf_->getValue<int>("setup","shift_x");
			shift[1] = pcf_->getValue<int>("setup","shift_y");
			shift[2] = pcf_->getValue<int>("setup","shift_z");
			
			levelmin_poisson = pcf_->getValue<unsigned>("setup","levelmin");
			
			int lfac = 1<<(ilevel-levelmin_poisson);
			
			lx[0] = 2*prefh_->size(ilevel, 0);
			lx[1] = 2*prefh_->size(ilevel, 1);
			lx[2] = 2*prefh_->size(ilevel, 2);
			x0[0] = prefh_->offset_abs(ilevel, 0) - lfac*shift[0] - lx[0]/4;
			x0[1] = prefh_->offset_abs(ilevel, 1) - lfac*shift[1] - lx[1]/4;
			x0[2] = prefh_->offset_abs(ilevel, 2) - lfac*shift[2] - lx[2]/4;
			
			if( randc[ilevel] == NULL )
			{
				randc[ilevel] = new rng( *randc[ilevel-1], ran_cube_size_, rngseeds_[ilevel], kavg, x0, lx );
				store_rnd( ilevel, randc[ilevel] );
				delete randc[ilevel-1];
				randc[ilevel-1] = NULL;
			}
			else
			{
				store_rnd( ilevel, randc[ilevel] );
				delete randc[ilevel-1];
				randc[ilevel-1] = NULL;				
			}
		}
		
		delete randc[levelmax_];
		randc[levelmax_] = NULL;
		
		
#if 1
		//... make sure that the coarse grid contains oct averages where it overlaps with a fine grid
		for( int ilevel=levelmax_; ilevel>levelmin_; --ilevel )
			correct_avg( ilevel-1, ilevel );
#endif
		
		//.. we do not have random numbers for a coarse level, generate them
		/*if( levelmax_rand_ >= (int)levelmin_ )
		{
			std::cerr << "lmaxread >= (int)levelmin\n";
			randc[levelmax_rand_] = new rng( (unsigned)pow(2,levelmax_rand_), rngfnames_[levelmax_rand_] );
			for( int ilevel = levelmax_rand_-1; ilevel >= (int)levelmin_; --ilevel )
				randc[ilevel] = new rng( *randc[ilevel+1] );
		}*/
	}

public:
	random_number_generator( config_file& cf, refinement_hierarchy& refh )
	: pcf_( &cf ), prefh_( &refh ), constraints( cf )
	{
		levelmin_ = prefh_->levelmin();
		levelmax_ = prefh_->levelmax();
		
		ran_cube_size_	= pcf_->getValueSafe<unsigned>("random","cubesize",DEF_RAN_CUBE_SIZE);
		
		parse_rand_parameters();
		
		compute_random_numbers();
	}
	
	~random_number_generator()
	{  }
	
	template< typename array >
	void load( array& A, int ilevel )
	{
		char fname[128];
		sprintf(fname,"wnoise_%04d.bin",ilevel);
		
		std::ifstream ifs( fname, std::ios::binary );
		if( !ifs.good() )
			throw std::runtime_error("A white noise file was not found. This is an internal inconsistency. Inform a developer!");
		
		int nx,ny,nz;
		ifs.read( reinterpret_cast<char*> (&nx), sizeof(int) );
		ifs.read( reinterpret_cast<char*> (&ny), sizeof(int) );
		ifs.read( reinterpret_cast<char*> (&nz), sizeof(int) );
		
		if( nx!=A.size(0) || ny!=A.size(1) || nz!=A.size(2) )
			throw std::runtime_error("White noise file is not aligned with array. This is an internal inconsistency. Inform a developer!");
		
		for( int i=0; i<nx; ++i )
		{
			std::vector<T> slice( ny*nz, 0.0 );
			ifs.read( reinterpret_cast<char*> ( &slice[0] ), nx*nz*sizeof(T) );
			
			for( int j=0; j<ny; ++j )
				for( int k=0; k<nz; ++k )
					A(i,j,k) = slice[j*nz+k];
		}		
		ifs.close();
		
	}
};


#endif //__RANDOM_HH

