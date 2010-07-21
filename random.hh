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
		std::cerr << " - Generating random numbers w/ sample cube size of " << cubesize_ << std::endl;
		
		ncubes_ = std::max((int)((double)res_/cubesize_),1);
		if( res_ < cubesize_ )
		{	
			ncubes_ = 1;
			cubesize_ = res_;
		}
		
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
		
		return sum/(ncubes_*ncubes_*ncubes_);
	}

public:
	template< class C >
	double fill_all( C& dat )
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
		initialize();
		fill_subvolume( x0, lx );
	}
	
	//! constructor
	random_numbers( unsigned res, unsigned cubesize, long baseseed, bool zeromean=true )
	: res_( res ), cubesize_( cubesize ), ncubes_( 1 ), baseseed_( baseseed )
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
		
		ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
		
		//... read data ...//
		
		ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
		if( blksz == nx*ny*sizeof(float) )
			vartype = 4;
		else if( blksz == nx*ny*sizeof(double) )
			vartype = 8;
		else
			throw std::runtime_error("corrupt random number file");
		
			
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
				in_float.assign(nx*ny,0.0f);
				ifs.read( (char*)&in_float[0], nx*ny*sizeof(float) );
				
				for( int jj=0,q=0; jj<(int)ny; ++jj )
					for( int kk=0; kk<(int)nx; ++kk ){
						sum += in_float[q];
						sum2 += in_float[q]*in_float[q];
						++count;
						
#ifdef RAND_DEBUG
						(*rnums_[0])(ii,jj,kk) = 0.0;
#else
						(*rnums_[0])(ii,jj,kk) = in_float[q++];
#endif
						
					}
				
			}
		}
		else if( vartype == 8 )
		{
			for( int ii=0; ii<(int)nz; ++ii )
			{
				in_double.assign(nx*ny,0.0f);				
				ifs.read( (char*)&in_double[0], nx*ny*sizeof(double) );
				
				for( int jj=0,q=0; jj<(int)ny; ++jj )
					for( int kk=0; kk<(int)nx; ++kk )
					{
						sum += in_double[q];
						sum2 += in_double[q]*in_double[q];
						++count;
#ifdef RAND_DEBUG						
						(*rnums_[0])(ii,jj,kk) = 0.0;
#else
						(*rnums_[0])(ii,jj,kk) = in_double[q++];
#endif
						
					}
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

