/*
 
 random.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
*/

//... for testing purposes.............
//#define DEGRADE_RAND1
//#define DEGRADE_RAND2
//.....................................

#ifndef __RANDOM_HH
#define __RANDOM_HH

#define DEF_RAN_CUBE_SIZE	32

#include <fstream>
#include <algorithm>
#include <omp.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "general.hh"
#include "mesh.hh"
#include "mg_operators.hh"
#include "constraints.hh"


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
	double fill_cube( int i, int j, int k);
	
	//! subtract a constant from an entire cube
	void subtract_from_cube( int i, int j, int k, double val );
	
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
		
		size_t icube = (i*ncubes_+j)*ncubes_+k;
		
		for( int ii=0; ii<(int)cubesize_; ++ii )
			for( int jj=0; jj<(int)cubesize_; ++jj )
				for( int kk=0; kk<(int)cubesize_; ++kk )
					dat(offi+ii,offj+jj,offk+kk) = (*rnums_[icube])(ii,jj,kk);
	}
	
	//! free the memory associated with a subcube
	void free_cube( int i, int j, int k );
	
	//! initialize member variables and allocate memory
	void initialize( void );
	
	//! fill a cubic subvolume of the full grid with random numbers
	double fill_subvolume( int *i0, int *n );
	
	//! fill an entire grid with random numbers
	double fill_all( void );
	
	//! fill an external array instead of the internal field
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
	
	//! write the number of allocated random number cubes to stdout
	void print_allocated( void );
	
public:
	
	//! constructor
	random_numbers( unsigned res, unsigned cubesize, long baseseed, int *x0, int *lx );	
	
	//! constructor for constrained fine field
	random_numbers( random_numbers<T>& rc, unsigned cubesize, long baseseed, 
				    bool kspace=false, int *x0_=NULL, int *lx_=NULL, bool zeromean=true );
	
	//! constructor
	random_numbers( unsigned res, unsigned cubesize, long baseseed, bool zeromean=true );
	
	
	//! constructor to read white noise from file
	random_numbers( unsigned res, std::string randfname );
	

	//! copy constructor for averaged field (not copying) hence explicit!
	explicit random_numbers( /*const*/ random_numbers <T>& rc, bool kdegrade = true );
	
	//! destructor
	~random_numbers()
	{
		for( unsigned i=0; i<rnums_.size(); ++i )
			if( rnums_[i] != NULL )
				delete rnums_[i];
		rnums_.clear();
	}
	
	//! access a random number, this allocates a cube and fills it with consistent random numbers
	inline T& operator()( int i, int j, int k, bool fillrand=true )
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
			rnums_[ icube ] = new Meshvar<T>( cubesize_, 0, 0, 0 );

			if( fillrand )
				fill_cube(ic, jc, kc);
		}
		
		//... determine cell in cube
		is = (i - ic * cubesize_ + cubesize_) % cubesize_;
		js = (j - jc * cubesize_ + cubesize_) % cubesize_;
		ks = (k - kc * cubesize_ + cubesize_) % cubesize_;
		
		return (*rnums_[ icube ])(is,js,ks);
	}
	
	//! free all cubes
	void free_all_mem( void )
	{
		for( unsigned i=0; i<rnums_.size(); ++i )
			if( rnums_[i] != NULL )
			{
				delete rnums_[i];	
				rnums_[i] = NULL;
			}
	}
	
	
};


/*!
 * @brief encapsulates all things for multi-scale white noise generation
 */
template< typename rng, typename T >
class random_number_generator
{
protected:
	config_file						* pcf_;
	refinement_hierarchy			* prefh_;
	constraint_set					constraints;
	
	int								levelmin_, 
									levelmax_, 
									levelmin_seed_;
	std::vector<long>				rngseeds_;
	std::vector<std::string>		rngfnames_;
	
	bool							disk_cached_;
	std::vector< std::vector<T>* >	mem_cache_;
	
	unsigned						ran_cube_size_;
	

protected:
	
	//! checks if the specified string is numeric
	bool is_number(const std::string& s);
	
	//! parses the random number parameters in the conf file
	void parse_rand_parameters( void );
	
	//! correct coarse grid averages for the change in small scale when using Fourier interpolation
	void correct_avg( int icoarse, int ifine );
	
	//! the main driver routine for multi-scale white noise generation
	void compute_random_numbers( void );
	
	//! store the white noise fields in memory or on disk
	void store_rnd( int ilevel, rng* prng );
	

public:
	
	//! constructor
	random_number_generator( config_file& cf, refinement_hierarchy& refh, transfer_function *ptf = NULL );	
	
	//! destructor
	~random_number_generator();
	
	//! load random numbers to a new array
	template< typename array >
	void load( array& A, int ilevel )
	{
		if( disk_cached_ )
		{
			char fname[128];
			sprintf(fname,"wnoise_%04d.bin",ilevel);
			
			LOGUSER("Loading white noise from file \'%s\'...",fname);
			
			std::ifstream ifs( fname, std::ios::binary );
			if( !ifs.good() )
			{	
				LOGERR("White noise file \'%s\'was not found.",fname);
				throw std::runtime_error("A white noise file was not found. This is an internal inconsistency. Inform a developer!");
				
			}
			
			int nx,ny,nz;
			ifs.read( reinterpret_cast<char*> (&nx), sizeof(int) );
			ifs.read( reinterpret_cast<char*> (&ny), sizeof(int) );
			ifs.read( reinterpret_cast<char*> (&nz), sizeof(int) );
			
			if( nx!=(int)A.size(0) || ny!=(int)A.size(1) || nz!=(int)A.size(2) )
			{	
				LOGERR("White noise file is not aligned with array. File: [%d,%d,%d]. Mem: [%d,%d,%d].",nx,ny,nz,A.size(0),A.size(1),A.size(2));
				throw std::runtime_error("White noise file is not aligned with array. This is an internal inconsistency. Inform a developer!");
			}
			
			for( int i=0; i<nx; ++i )
			{
				std::vector<T> slice( ny*nz, 0.0 );
				ifs.read( reinterpret_cast<char*> ( &slice[0] ), ny*nz*sizeof(T) );
				
				#pragma omp parallel for
				for( int j=0; j<ny; ++j )
					for( int k=0; k<nz; ++k )
						A(i,j,k) = slice[j*nz+k];
				
			}		
			ifs.close();	
		}
		else
		{
			LOGUSER("Copying white noise from memory cache...");
			
			if( mem_cache_[ilevel-levelmin_] == NULL )
				LOGERR("Tried to access mem-cached random numbers for level %d. But these are not available!\n",ilevel);
			
			int nx( A.size(0) ), ny( A.size(1) ), nz( A.size(2) );
			
			if ( (size_t)nx*(size_t)ny*(size_t)nz != mem_cache_[ilevel-levelmin_]->size() )
			{
				LOGERR("White noise file is not aligned with array. File: [%d,%d,%d]. Mem: [%d,%d,%d].",nx,ny,nz,A.size(0),A.size(1),A.size(2));
				throw std::runtime_error("White noise file is not aligned with array. This is an internal inconsistency. Inform a developer!");
			}
			
			#pragma omp parallel for
			for( int i=0; i<nx; ++i )
				for( int j=0; j<ny; ++j )
					for( int k=0; k<nz; ++k )
						A(i,j,k) = (*mem_cache_[ilevel-levelmin_])[((size_t)i*ny+(size_t)j)*nz+(size_t)k];
			
			std::vector<T>().swap( *mem_cache_[ilevel-levelmin_] );
			delete mem_cache_[ilevel-levelmin_];
			mem_cache_[ilevel-levelmin_] = NULL;
			
		}

		
	}
};

typedef random_numbers<real_t> rand_nums;
typedef random_number_generator< rand_nums,real_t> rand_gen;


#endif //__RANDOM_HH

