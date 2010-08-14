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
	
	//! vector of 3D meshes (the random number cbues) with random numbers
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
		unsigned i0cube[3], ncube[3];
		
		i0cube[0] = (int)((double)i0[0]/cubesize_);
		i0cube[1] = (int)((double)i0[1]/cubesize_);
		i0cube[2] = (int)((double)i0[2]/cubesize_);
		
		ncube[0] = (int)((double)(i0[0]+n[0])/cubesize_+1.0)-i0cube[0];
		ncube[1] = (int)((double)(i0[1]+n[1])/cubesize_+1.0)-i0cube[1];
		ncube[2] = (int)((double)(i0[2]+n[2])/cubesize_+1.0)-i0cube[2];
		
		double mean = 0.0;
		
		#pragma omp parallel for reduction(+:mean)
		for( int i=i0cube[0]; i<(int)ncube[0]; ++i )
			for( int j=i0cube[1]; j<(int)ncube[1]; ++j )
				for( int k=i0cube[2]; k<(int)ncube[2]; ++k )
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
			unsigned ixoff2=51,iyoff2=51,izoff2=51;
			unsigned nx2=52, ny2=52, nz2=52;
			
			unsigned ixoff1=34,iyoff1=34,izoff1=34;
			unsigned nx1=120,ny1=120,nz1=120;
			
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
	
	
public:
	
	//! constructor
	random_numbers( unsigned res, unsigned cubesize, long baseseed, int *x0, int *lx )
	: res_( res ), cubesize_( cubesize ), ncubes_( 1 ), baseseed_( baseseed )
	{
		std::cout << " - Generating random numbers (1) with seed " << baseseed << std::endl;
		initialize();
		fill_subvolume( x0, lx );
	}
	
	//! constructor
	random_numbers( unsigned res, unsigned cubesize, long baseseed, bool zeromean=true )
	: res_( res ), cubesize_( cubesize ), ncubes_( 1 ), baseseed_( baseseed )
	{
		std::cout << " - Generating random numbers (2) with seed " << baseseed << std::endl;
		
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
						
#ifdef RAND_DEBUG
						(*rnums_[0])(kk,jj,ii) = 0.0;
#else
						(*rnums_[0])(kk,jj,ii) = -in_float[q++];
#endif
						
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
#ifdef RAND_DEBUG						
						(*rnums_[0])(kk,jj,ii) = 0.0;
#else
						(*rnums_[0])(kk,jj,ii) = -in_double[q++];
#endif
						
					}
				ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
				if( blksz != nx*ny*sizeof(double) )
					throw std::runtime_error("corrupt random number file");
				
			}
		}
		
		double mean, var;
		mean = sum/count;
		var = sum2/count-mean*mean;
#ifdef RAND_DEBUG
		(*rnums_[0])(cubesize_/2-1,cubesize_/2-1,cubesize_/2-1) = 1.0;
#endif
		
		std::cout << " - Random numbers in file have \n"
				  << "     mean = " << mean << " and var = " << var << std::endl;
		
	}
	
	//... create constrained new field
	random_numbers( random_numbers<T>& rc, unsigned cubesize, long baseseed, bool zeromean=true )
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

	//... copy construct by averaging down
	explicit random_numbers( /*const*/ random_numbers <T>& rc )
	{
		//if( res > rc.m_res || (res/rc.m_res)%2 != 0 )
		//			throw std::runtime_error("Invalid restriction in random number container copy constructor.");
		
		double sum = 0.0, sum2 = 0.0;
		unsigned count = 0;
		
		//... initialize properties of container
		if( rc.rnums_.size() == 1 )
		{
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


#endif //__RANDOM_HH

