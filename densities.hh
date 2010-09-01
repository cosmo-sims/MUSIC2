/*
 
 densities.hh - This file is part of MUSIC -
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

#ifndef __DENSITIES_HH
#define __DENSITIES_HH


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

#include <assert.h>

#include "config_file.hh"
#include "random.hh"
#include "cosmology.hh"
#include "transfer_function.hh"


void GenerateDensityHierarchy(	config_file& cf, transfer_function *ptf, tf_type type, 
							  refinement_hierarchy& refh, grid_hierarchy& delta, bool bdeconvolve, bool smooth );

void GenerateDensityUnigrid( config_file& cf, transfer_function *ptf, tf_type type, 
							refinement_hierarchy& refh, grid_hierarchy& delta, bool kspace, bool deconvolve, bool smooth );

void normalize_density( grid_hierarchy& delta );


/*!
 * @class DensityGrid
 * @brief provides infrastructure for computing the initial density field
 *
 * This class provides access and data management member functions that
 * are used when computing the initial density field by convolution with
 * transfer functions.
 */
template< typename real_t >
class DensityGrid
{
public:
	
	int nx_;	//!< number of grid cells in x-direction
	int ny_;	//!< number of grid cells in y-direction
	int nz_;	//!< number of grid cells in z-direction
	int nzp_;	//!< number of cells in memory (z-dir), used for Nyquist padding
	
	//! the actual data container in the form of a 1D array
	std::vector< real_t > data_;
	
	//! constructor
	/*! constructs an instance given the dimensions of the density field
	 * @param nx the number of cells in x
	 * @param ny the number of cells in y
	 * @param nz the number of cells in z
	 */
	DensityGrid( unsigned nx, unsigned ny, unsigned nz )
	: nx_(nx), ny_(ny), nz_(nz)
	{
		data_.assign(nx_*ny_*2*(nz_/2+1),0.0);
		nzp_ = 2*(nz_/2+1);
	}
	
	//! copy constructor
	explicit DensityGrid( const DensityGrid<real_t> & g )
	: nx_(g.nx_), ny_(g.ny_), nz_(g.nz_), nzp_(g.nzp_)
	{
		data_ = g.data_;
	}
	
	//!destructor
	~DensityGrid()
	{	}
	
	//! clears the density object
	/*! sets all dimensions to zero and frees the memory
	 */
	void clear( void )
	{
		nx_ = ny_ = nz_ = nzp_ = 0;
		data_.clear();
		std::vector<real_t>().swap(data_);
	}
	
	//! query the 3D array sizes of the density object
	/*! returns the size of the 3D density object along a specified dimension
	 * @param i the dimension for which size is to be returned
	 * @returns array size along dimension i
	 */
	int size( int i )
	{
		if(i==0) return nx_;
		if(i==1) return ny_;
		return ny_;
	}
	
	//! zeroes the density object
	/*! sets all values to 0.0
	 */
	void zero( void )
	{
		data_.assign(data_.size(),0.0);
	}
	
	//! assigns the contents of another DensityGrid to this
	DensityGrid& operator=( const DensityGrid<real_t>& g )
	{
		nx_ = g.nx_;
		ny_ = g.ny_;
		nz_ = g.nz_;
		nzp_= g.nzp_;
		data_ = g.data_;
		
		return *this;
	}
	
	//! 3D index based data access operator
	inline real_t& operator()( int i, int j, int k )
	{	return data_[(i*ny_+j)*nzp_+k]; 	}
	
	//! 3D index based const data access operator 
	inline const real_t& operator()( int i, int j, int k ) const
	{	return data_[(i*ny_+j)*nzp_+k]; 	}
	
	//! recover the pointer to the 1D data array
	inline real_t * get_data_ptr( void )
	{	return &data_[0];	}
	
	
	//! fills the density field with random number values
	/*! given a pointer to a random_numbers object, fills the field with random values
	 *  @param prc pointer to a random_numbers object
	 *  @param variance the variance of the random numbers (the values returned by prc are multiplied by this)
	 *  @param i0 x-offset (shift) in cells of the density field with respect to the random number field
	 *  @param j0 y-offset (shift) in cells of the density field with respect to the random number field
	 *  @param k0 z-offset (shift) in cells of the density field with respect to the random number field
	 *  @param setzero boolean, if true, the global mean will be subtracted
	 */
	void fill_rand( /*const*/ random_numbers<real_t>* prc, real_t variance, int i0, int j0, int k0, bool setzero=false )
	{
		double sum = 0.0;
		
		#pragma omp parallel for reduction(+:sum)
		for( int i=0; i<nx_; ++i )
			for( int j=0; j<ny_; ++j )
				for( int k=0; k<nz_; ++k )
				{
					(*this)(i,j,k) = (*prc)(i0+i,j0+j,k0+k) * variance;
					sum += (*this)(i,j,k);
				}
		
		sum /= nx_*ny_*nz_;
		
		if( setzero )
		{
			#pragma omp parallel for
			for( int i=0; i<nx_; ++i )
				for( int j=0; j<ny_; ++j )
					for( int k=0; k<nz_; ++k )
						(*this)(i,j,k) -= sum;
						
		}
					
	}
	
	//! copies the data from another field with 3D index-based access operator
	template< class array3 >
	void copy( array3& v )
	{
		for( int ix=0; ix<(int)nx_; ++ix )
			for( int iy=0; iy<(int)ny_; ++iy )	
				for( int iz=0; iz<(int)nz_; ++iz )
					v(ix,iy,iz) = (*this)(ix,iy,iz);
						
	}

	//! adds the data from another field with 3D index-based access operator
	template< class array3 >
	void copy_add( array3& v )
	{
		for( int ix=0; ix<(int)nx_; ++ix )
			for( int iy=0; iy<(int)ny_; ++iy )	
				for( int iz=0; iz<(int)nz_; ++iz )
					v(ix,iy,iz) += (*this)(ix,iy,iz);
	}
	
	//! subtract the mean of the field 
	/*! subtracting the total mean implies that a restriction does not change the mass
	 */
	double subtract_mean( void )
	{
		double sum = 0.0;
		unsigned count;
		
		for( int i=0; i<nx_; i++ )
			for( int j=0; j<ny_; j++ )
				for( int k=0; k<nz_; k++ )
				{
					sum += (*this)(i,j,k);
					count++;
				}
		sum /= count;
		
					
					
		for( int i=0; i<nx_; i++ )
			for( int j=0; j<ny_; j++ )
				for( int k=0; k<nz_; k++ )
					(*this)(i,j,k)			-= sum;
		
		return sum;
		
	}
	
	//! subtract the mean of each oct of the field 
	/*! subtracting the mean of each oct implies that a restriction of the field to the coarser level
	 *  is zero everywhere (needed for Hoffman-Ribak constraints)
	 */
	void subtract_oct_mean( void )
	{
		for( int i=0; i<nx_; i+=2 )
			for( int j=0; j<ny_; j+=2 )
				for( int k=0; k<nz_; k+=2 )
				{
					real_t avg = 
						0.125*((*this)(i,j,k)+(*this)(i+1,j,k)+(*this)(i,j+1,k)+(*this)(i,j,k+1)+
							   (*this)(i+1,j+1,k)+(*this)(i+1,j,k+1)+(*this)(i,j+1,k+1)+(*this)(i+1,j+1,k+1));
						
					(*this)(i,j,k)			-= avg;
					(*this)(i+1,j,k)		-= avg;
					(*this)(i,j+1,k)		-= avg;
					(*this)(i,j,k+1)		-= avg;
					(*this)(i+1,j+1,k)		-= avg;
					(*this)(i+1,j,k+1)		-= avg;
					(*this)(i,j+1,k+1)		-= avg;
					(*this)(i+1,j+1,k+1)	-= avg;
					
				}
		
	}
	
	//! replaces the value of all cells in an oct with the average value of the oct
	void set_to_oct_mean( void )
	{
		for( int i=0; i<nx_; i+=2 )
			for( int j=0; j<ny_; j+=2 )
				for( int k=0; k<nz_; k+=2 )
				{
					real_t avg = 
					0.125*((*this)(i,j,k)+(*this)(i+1,j,k)+(*this)(i,j+1,k)+(*this)(i,j,k+1)+
						   (*this)(i+1,j+1,k)+(*this)(i+1,j,k+1)+(*this)(i,j+1,k+1)+(*this)(i+1,j+1,k+1));
					
					(*this)(i,j,k)			= avg;
					(*this)(i+1,j,k)		= avg;
					(*this)(i,j+1,k)		= avg;
					(*this)(i,j,k+1)		= avg;
					(*this)(i+1,j+1,k)		= avg;
					(*this)(i+1,j,k+1)		= avg;
					(*this)(i,j+1,k+1)		= avg;
					(*this)(i+1,j+1,k+1)	= avg;
					
				}
		
	}
	
	//! sets the field to zero inside of a rectangular box
	void zero_subgrid( unsigned oxsub, unsigned oysub, unsigned ozsub, unsigned lxsub, unsigned lysub, unsigned lzsub )
	{
		//... correct offsets for padding (not needed for top grid)
		for( int ix=oxsub; ix<(int)(oxsub+lxsub); ++ix )
			for( int iy=oysub; iy<(int)(oysub+lysub); ++iy )	
				for( int iz=ozsub; iz<(int)(ozsub+lzsub); ++iz )
					(*this)(ix,iy,iz) = 0.0;
		
	}
	
	
	//! sets the field to zero inside of a rectangular box including the boundary of the box (used for isolated FFT)
	void zero_padded_subgrid( unsigned oxsub, unsigned oysub, unsigned ozsub, unsigned lxsub, unsigned lysub, unsigned lzsub )
	{
		//... correct offsets for padding (not needed for top grid)
		oxsub -= lxsub/2;
		oysub -= lysub/2;
		ozsub -= lzsub/2;
		lxsub *= 2;
		lysub *= 2;
		lzsub *= 2;
		for( int ix=oxsub; ix<(int)(oxsub+lxsub); ++ix )
			for( int iy=oysub; iy<(int)(oysub+lysub); ++iy )	
				for( int iz=ozsub; iz<(int)(ozsub+lzsub); ++iz )
					(*this)(ix,iy,iz) = 0.0;
					
	}
	
	//! sets the field to zero outside of a rectangular box including the boundary of the box (used for isolated FFT)
	void zero_but_padded_subgrid( int oxsub, int oysub, int ozsub, int lxsub, int lysub, int lzsub )
	{
		//... correct offsets for padding (not needed for top grid)
		oxsub -= lxsub/2;
		oysub -= lysub/2;
		ozsub -= lzsub/2;
		lxsub *= 2;
		lysub *= 2;
		lzsub *= 2;
		
		for( int ix=0; ix<nx_; ++ix )
			for( int iy=0; iy<ny_; ++iy )
				for( int iz=0; iz<nz_; ++iz )		
				{
					if( (ix<oxsub||ix>=oxsub+lxsub)
					   ||  (iy<oysub||iy>=oysub+lysub)
					   ||  (iz<ozsub||iz>=ozsub+lzsub) )
					(*this)(ix,iy,iz) = 0.0;	
					   
				}
	}
	
	//! sets the field to zero outside of a rectangular box
	void zero_but_subgrid( int oxsub, int oysub, int ozsub, int lxsub, int lysub, int lzsub )
	{
		for( int ix=0; ix<nx_; ++ix )
			for( int iy=0; iy<ny_; ++iy )
				for( int iz=0; iz<nz_; ++iz )		
				{
					if( (ix<oxsub||ix>=oxsub+lxsub)
					   ||  (iy<oysub||iy>=oysub+lysub)
					   ||  (iz<ozsub||iz>=ozsub+lzsub) )
						(*this)(ix,iy,iz) = 0.0;	
					
				}
	}
	
	//! sets the field to zero except on the boundary of a rectangular box
	void zero_but_subgrid_bnd( unsigned oxsub, unsigned oysub, unsigned ozsub, unsigned lxsub, unsigned lysub, unsigned lzsub )
	{
		//... correct offsets for padding (not needed for top grid)
		
		// zero the subgrid
		for( int ix=oxsub; ix<(int)(oxsub+lxsub); ++ix )
			for( int iy=oysub; iy<(int)(oysub+lysub); ++iy )	
				for( int iz=ozsub; iz<(int)(ozsub+lzsub); ++iz )
					(*this)(ix,iy,iz) = 0.0;
		
		oxsub -= lxsub/2;
		oysub -= lysub/2;
		ozsub -= lzsub/2;
		lxsub *= 2;
		lysub *= 2;
		lzsub *= 2;
		
		// zero the outside, except the boundary
		//#pragma parallel nowait
		{
			
		for( int ix=0; ix<oxsub; ++ix )
			for( int iy=0; iy<ny_; ++iy )
				for( int iz=0; iz<nz_; ++iz )
					(*this)(ix,iy,iz) = 0.0;
			
		for( int ix=oxsub+lxsub; ix<nx_; ++ix )
			for( int iy=0; iy<ny_; ++iy )
				for( int iz=0; iz<nz_; ++iz )
					(*this)(ix,iy,iz) = 0.0;
		
		for( int ix=0; ix<nx_; ++ix )
			for( int iy=0; iy<oysub; ++iy )
				for( int iz=0; iz<nz_; ++iz )
					(*this)(ix,iy,iz) = 0.0;
		for( int ix=0; ix<nx_; ++ix )
			for( int iy=oysub+lysub; iy<ny_; ++iy )
				for( int iz=0; iz<nz_; ++iz )
					(*this)(ix,iy,iz) = 0.0;
		
		for( int ix=0; ix<nx_; ++ix )
			for( int iy=0; iy<ny_; ++iy )
				for( int iz=0; iz<ozsub; ++iz )
					(*this)(ix,iy,iz) = 0.0;
		for( int ix=0; ix<nx_; ++ix )
			for( int iy=0; iy<ny_; ++iy )
				for( int iz=ozsub+lzsub; iz<nz_; ++iz )
					(*this)(ix,iy,iz) = 0.0;
		
		}
		
	}
	
};

template< typename real_t >
class PaddedDensitySubGrid : public DensityGrid<real_t>
{
public:
	using DensityGrid<real_t>::nx_;
	using DensityGrid<real_t>::ny_;
	using DensityGrid<real_t>::nz_;
	using DensityGrid<real_t>::data_;	
	
	using DensityGrid<real_t>::fill_rand;
	using DensityGrid<real_t>::get_data_ptr;
	
	int ox_, oy_, oz_;
public:
	
	PaddedDensitySubGrid( int ox, int oy, int oz, unsigned nx, unsigned ny, unsigned nz)
	: DensityGrid<real_t>(2*nx,2*ny,2*nz), ox_(ox), oy_(oy), oz_(oz)
	{

		//if( 2*ox-(int)nx/4 < 0 || 2*oy-(int)ny/4 < 0 || 2*oz-(int)nz/4 < 0 )
		//	throw std::runtime_error("subgrid extends across top grid");
		
		//.. the size in top grid cells is nx/2,ny/2,nz/2, so padding starts at ox-nx/4...
		//.. loop over relevant part of top grid and copy down grid values
		/*for( int ix=ox-nx/4,ixs=0; ix<ox+3*nx/4; ix++,ixs+=2 )
		 for( int iy=oy-ny/4,iys=0; iy<oy+3*ny/4; iy++,iys+=2 )
		 for( int iz=oz-nz/4,izs=0; iz<oz+3*nz/4; iz++,izs+=2 )
		 {
		 (*this)(ixs,iys,izs)		=
		 (*this)(ixs+1,iys,izs)		=
		 (*this)(ixs,iys+1,izs)		=
		 (*this)(ixs,iys,izs+1)		=
		 (*this)(ixs+1,iys+1,izs)	=
		 (*this)(ixs+1,iys,izs+1)	=
		 (*this)(ixs,iys+1,izs+1)	=
		 (*this)(ixs+1,iys+1,izs+1)	=	top(ix,iy,iz);
		 }*/
	}
	
	PaddedDensitySubGrid( const PaddedDensitySubGrid<real_t>& o )
	: DensityGrid<real_t>( o ), ox_(o.ox_), oy_(o.oy_), oz_(o.oz_)
	{ }
	
	void zero_padded_subgrid( unsigned oxsub, unsigned oysub, unsigned ozsub, unsigned lxsub, unsigned lysub, unsigned lzsub )
	{
		//... correct offsets for padding (not needed for top grid)
		
		oxsub += nx_/4;
		oysub += ny_/4;
		ozsub += nz_/4;
		
		//... correct offsets for padding (not needed for top grid)
		for( int ix=oxsub; ix<(int)(oxsub+lxsub); ++ix )
			for( int iy=oysub; iy<(int)(oysub+lysub); ++iy )	
				for( int iz=ozsub; iz<(int)(ozsub+lzsub); ++iz )
					(*this)(ix,iy,iz) = 0.0;
		
	}
	
	void zero_subgrid( unsigned oxsub, unsigned oysub, unsigned ozsub, unsigned lxsub, unsigned lysub, unsigned lzsub )
	{
		//... correct offsets for padding (not needed for top grid)
		
		oxsub += nx_/4;
		oysub += ny_/4;
		ozsub += nz_/4;
		
		for( int ix=oxsub; ix<(int)(oxsub+lxsub); ++ix )
			for( int iy=oysub; iy<(int)(oysub+lysub); ++iy )	
				for( int iz=ozsub; iz<(int)(ozsub+lzsub); ++iz )
					(*this)(ix,iy,iz) = 0.0;
		
	}
	
	void zero_but_subgrid_bnd( unsigned oxsub, unsigned oysub, unsigned ozsub, unsigned lxsub, unsigned lysub, unsigned lzsub )
	{
		//... correct offsets for padding (not needed for top grid)
		
		oxsub += nx_/4;
		oysub += ny_/4;
		ozsub += nz_/4;
		
		//lxsub += nx_/4;
		//lysub += nx_/4;
		//lzsub += nx_/4;
		
		// zero the subgrid
		for( int ix=oxsub; ix<(int)(oxsub+lxsub); ++ix )
			for( int iy=oysub; iy<(int)(oysub+lysub); ++iy )	
				for( int iz=ozsub; iz<(int)(ozsub+lzsub); ++iz )
					(*this)(ix,iy,iz) = 0.0;
		
		oxsub -= lxsub/2;
		oysub -= lysub/2;
		ozsub -= lzsub/2;
		lxsub *= 2;
		lysub *= 2;
		lzsub *= 2;
		
		// zero the outside, except the boundary
		//#pragma parallel block nowait
		{
			
			for( int ix=0; ix<oxsub; ++ix )
				for( int iy=0; iy<ny_; ++iy )
					for( int iz=0; iz<nz_; ++iz )
						(*this)(ix,iy,iz) = 0.0;
			
			std::cerr << oxsub+lxsub << " -> " << nx_ << std::endl;
			for( int ix=oxsub+lxsub; ix<nx_; ++ix )
				for( int iy=0; iy<ny_; ++iy )
					for( int iz=0; iz<nz_; ++iz )
						(*this)(ix,iy,iz) = 0.0;
			
			for( int ix=0; ix<nx_; ++ix )
				for( int iy=0; iy<oysub; ++iy )
					for( int iz=0; iz<nz_; ++iz )
						(*this)(ix,iy,iz) = 0.0;
			for( int ix=0; ix<nx_; ++ix )
				for( int iy=oysub+lysub; iy<ny_; ++iy )
					for( int iz=0; iz<nz_; ++iz )
						(*this)(ix,iy,iz) = 0.0;
			
			for( int ix=0; ix<nx_; ++ix )
				for( int iy=0; iy<ny_; ++iy )
					for( int iz=0; iz<ozsub; ++iz )
						(*this)(ix,iy,iz) = 0.0;
			for( int ix=0; ix<nx_; ++ix )
				for( int iy=0; iy<ny_; ++iy )
					for( int iz=ozsub+lzsub; iz<nz_; ++iz )
						(*this)(ix,iy,iz) = 0.0;
			
		}
	}		
	
	
	template< class array3 >
	void upload_bnd_add( array3& top )
	{
		int ox=ox_-nx_/8, oy=oy_-ny_/8, oz=oz_-nz_/8;
		int lx=ox+nx_/2, ly=oy+ny_/2, lz=oz+nz_/2;
		ox = std::max(ox,0);
		oy = std::max(oy,0);
		oz = std::max(oz,0);
		lx = std::min(lx,(int)top.size(0));
		ly = std::min(ly,(int)top.size(1));
		lz = std::min(lz,(int)top.size(2));
		
		
		if( ox < 0 || oy < 0 || oz < 0 ){
			
			std::cerr << "offset = " << ox << ", " << oy << ", " << oz << std::endl;
			std::cerr << "nx     = " << nx_ << ", " << ny_ << ", " << nz_ << std::endl;
			
			throw std::runtime_error("Subgrid extends beyond top grid. Increase levelmin or reduce size of refinement region!");
		}
		for( int ix=0,ixu=ox;ix<nx_&&ixu<lx;ix+=2,++ixu )
			for( int iy=0,iyu=oy;iy<ny_&&iyu<ly;iy+=2,++iyu )
				for( int iz=0,izu=oz;iz<nz_&&izu<lz;iz+=2,++izu )
				{
					if( (ix<nx_/4||ix>=3*nx_/4)
					   || (iy<ny_/4||iy>=3*ny_/4)
					   || (iz<nz_/4||iz>=3*nz_/4) )
					{
						top(ixu,iyu,izu) += 0.125* ((*this)(ix,iy,iz)+(*this)(ix+1,iy,iz)+(*this)(ix,iy+1,iz)+(*this)(ix,iy,iz+1)
											+(*this)(ix+1,iy+1,iz)+(*this)(ix+1,iy,iz+1)+(*this)(ix,iy+1,iz+1)+(*this)(ix+1,iy+1,iz+1) );

						
					}
				}
	}
	
	void constrain( const DensityGrid<real_t>& top  )
	{
		int ox=ox_-nx_/8, oy=oy_-ny_/8, oz=oz_-nz_/8;
		
		if( ox < 0 || oy < 0 || oz < 0 ){
		
			std::cerr << "offset = " << ox << ", " << oy << ", " << oz << std::endl;
			std::cerr << "nx     = " << nx_ << ", " << ny_ << ", " << nz_ << std::endl;
			
			throw std::runtime_error("Subgrid extends beyond top grid. Increase levelmin or reduce size of refinement region!");
		}
			
		//...top grid is not padded
		
		//.. fine grid will be multiplied by sqrt(8) later to account for increase in variance
		//.. so, we need to divide the coarse grid value by it
		
		real_t fac = 1.0/sqrt(8);
		for( int ixs=0,ix=ox;ixs<nx_;ixs+=2,++ix )
			for( int iys=0,iy=oy;iys<ny_;iys+=2,++iy )
				for( int izs=0,iz=oz;izs<nz_;izs+=2,++iz )
				{
					real_t mean = 0.125 * 
						((*this)(ixs,iys,izs)+
						 (*this)(ixs+1,iys,izs)+
						 (*this)(ixs,iys+1,izs)+
						 (*this)(ixs,iys,izs+1)+
						 (*this)(ixs+1,iys+1,izs)+
						 (*this)(ixs+1,iys,izs+1)+
						 (*this)(ixs,iys+1,izs+1)+
						 (*this)(ixs+1,iys+1,izs+1));
					 
					
					double aa = top(ix,iy,iz)*fac - mean;
					
					(*this)(ixs,iys,izs)		+=	aa;
					(*this)(ixs+1,iys,izs)		+=	aa;
					(*this)(ixs,iys+1,izs)		+=	aa;
					(*this)(ixs,iys,izs+1)		+=	aa;
					(*this)(ixs+1,iys+1,izs)	+=	aa;
					(*this)(ixs+1,iys,izs+1)	+=	aa;
					(*this)(ixs,iys+1,izs+1)	+=	aa;
					(*this)(ixs+1,iys+1,izs+1)	+=	aa;
				}
		
#if 0		
		unsigned nx=nx_/2, ny=ny_/2, nz=nz_/2;
		
		int ox=ox_, oy=oy_, oz=oz_;
		
		if( ox-nx/4 < 0 || oy-ny/4 < 0 || oz-nz/4 < 0 )
			throw("subgrid extends across top grid");
		
		//.. the size in top grid cells is nx/2,ny/2,nz/2, so padding starts at ox-nx/4...
		//.. loop over relevant part of top grid and copy down grid values
		for( int ix=ox-nx/4,ixs=0; ix<ox+3*nx/4; ix++,ixs+=2 )
			for( int iy=oy-ny/4,iys=0; iy<oy+3*ny/4; iy++,iys+=2 )
				for( int iz=oz-nz/4,izs=0; iz<oz+3*nz/4; iz++,izs+=2 )
				{
					//... removed -r, as the mean will be taken out right at the beginning
					(*this)(ixs,iys,izs)		+=	top(ix,iy,iz);
					(*this)(ixs+1,iys,izs)		+=	top(ix,iy,iz);
					(*this)(ixs,iys+1,izs)		+=	top(ix,iy,iz);
					(*this)(ixs,iys,izs+1)		+=	top(ix,iy,iz);
					(*this)(ixs+1,iys+1,izs)	+=	top(ix,iy,iz);
					(*this)(ixs+1,iys,izs+1)	+=	top(ix,iy,iz);
					(*this)(ixs,iys+1,izs+1)	+=	top(ix,iy,iz);
					(*this)(ixs+1,iys+1,izs+1)	+=	top(ix,iy,iz);
				}
#endif
	}
	
	void constrain( const PaddedDensitySubGrid<real_t>& top )
	{
		///... top grid is padded too...
		//.. ox is shifted by nx_/4, padded overlap starts at -nx_/8
		int ox=ox_+top.nx_/4-nx_/8, oy=oy_+top.ny_/4-ny_/8, oz=oz_+top.nz_/4-nz_/8;
		if( ox < 0 || oy < 0 || oz < 0 )
			throw std::runtime_error("subgrid extends beyond top grid");
		double fac = 1.0/sqrt(8.0);
		for( int ixs=0,ix=ox;ixs<nx_;ixs+=2,++ix )
			for( int iys=0,iy=oy;iys<ny_;iys+=2,++iy )
				for( int izs=0,iz=oz;izs<nz_;izs+=2,++iz )
				{
					
					double mean = 0.125 * 
					((*this)(ixs,iys,izs)+
					 (*this)(ixs+1,iys,izs)+
					 (*this)(ixs,iys+1,izs)+
					 (*this)(ixs,iys,izs+1)+
					 (*this)(ixs+1,iys+1,izs)+
					 (*this)(ixs+1,iys,izs+1)+
					 (*this)(ixs,iys+1,izs+1)+
					 (*this)(ixs+1,iys+1,izs+1));
					
					
					double aa = top(ix,iy,iz)*fac - mean;
					
					(*this)(ixs,iys,izs)		+=	aa;
					(*this)(ixs+1,iys,izs)		+=	aa;
					(*this)(ixs,iys+1,izs)		+=	aa;
					(*this)(ixs,iys,izs+1)		+=	aa;
					(*this)(ixs+1,iys+1,izs)	+=	aa;
					(*this)(ixs+1,iys,izs+1)	+=	aa;
					(*this)(ixs,iys+1,izs+1)	+=	aa;
					(*this)(ixs+1,iys+1,izs+1)	+=	aa;
				}

	}
	
	
	template< class array3 >
	void upload_bnd( unsigned nbnd, array3& v )
	{
		for( int ix=nx_/4-2*nbnd,ixu=-nbnd; ix<3*nx_/4+2*nbnd; ix+=2,++ixu)
			for( int iy=ny_/4-2*nbnd,iyu=-nbnd; iy<3*ny_/4+2*nbnd; iy+=2,++iyu )	
				for( int iz=nz_/4-2*nbnd,izu=-nbnd; iz<3*nz_/4+2*nbnd; iz+=2,++izu )
				{
					if( (ix<nx_/4||ix>=3*nx_/4)
					   || (iy<ny_/4||iy>=3*ny_/4)
					   || (iz<nz_/4||iz>=3*nz_/4) )
					   {
						   v(ixu+ox_,iyu+oy_,izu+oz_) 
						   = 1.0;
					   }
				}
	}
	
	template< class array3 >
	void copy_unpad( array3& v )
	{
		for( int ix=nx_/4,ixu=0; ix<3*nx_/4; ++ix,++ixu )
			for( int iy=ny_/4,iyu=0; iy<3*ny_/4; ++iy,++iyu )	
				for( int iz=nz_/4,izu=0; iz<3*nz_/4; ++iz,++izu )
					v(ixu,iyu,izu) = (*this)(ix,iy,iz);
	}
	
	template< class array3 >
	void copy_add_unpad( array3& v )
	{
		for( int ix=nx_/4,ixu=0; ix<3*nx_/4; ++ix,++ixu )
			for( int iy=ny_/4,iyu=0; iy<3*ny_/4; ++iy,++iyu )	
				for( int iz=nz_/4,izu=0; iz<3*nz_/4; ++iz,++izu )
					v(ixu,iyu,izu) += (*this)(ix,iy,iz);
	}
	
	template< class array3 >
	void copy_subtract_unpad( array3& v )
	{
		for( int ix=nx_/4,ixu=0; ix<3*nx_/4; ++ix,++ixu )
			for( int iy=ny_/4,iyu=0; iy<3*ny_/4; ++iy,++iyu )	
				for( int iz=nz_/4,izu=0; iz<3*nz_/4; ++iz,++izu )
					v(ixu,iyu,izu) -= (*this)(ix,iy,iz);
	}
	
	
	double oct_mean( int i, int j, int k )
	{
		
		return 0.125*((*this)(i,j,k)+(*this)(i+1,j,k)+(*this)(i,j+1,k)+(*this)(i,j,k+1)
					  +(*this)(i+1,j+1,k)+(*this)(i+1,j,k+1)+(*this)(i,j+1,k+1)+(*this)(i+1,j+1,k+1));
	}
	
	void set_oct_mean( int i, int j, int k, double val )
	{
		(*this)(i,j,k)			= val;
		(*this)(i+1,j,k)		= val;
		(*this)(i,j+1,k)		= val;
		(*this)(i,j,k+1)		= val;
		(*this)(i+1,j+1,k)		= val;
		(*this)(i+1,j,k+1)		= val;
		(*this)(i,j+1,k+1)		= val;
		(*this)(i+1,j+1,k+1)	= val;
	}
	
	void add_oct_val( int i, int j, int k, double val )
	{
		(*this)(i,j,k)		+= val;
		(*this)(i+1,j,k)		+= val;
		(*this)(i,j+1,k)		+= val;
		(*this)(i,j,k+1)		+= val;
		(*this)(i+1,j+1,k)	+= val;
		(*this)(i+1,j,k+1)	+= val;
		(*this)(i,j+1,k+1)	+= val;
		(*this)(i+1,j+1,k+1)	+= val;
	}
	
	void subtract_boundary_oct_mean( void )
	{
		#pragma omp parallel for
		for( int ix=0; ix<nx_/4-1; ix+=2 )
		{	
			for( int iy=0; iy<ny_; iy+=2 )
				for( int iz=0; iz<nz_; iz+=2 )
				{
					add_oct_val( ix,iy,iz, -oct_mean(ix,iy,iz) );
					add_oct_val( ix+3*nx_/4,iy,iz, -oct_mean(ix+3*nx_/4,iy,iz) );
				}
		}
		
		#pragma omp parallel for
		for( int ix=0; ix<nx_; ix+=2 )
			for( int iy=0; iy<ny_/4-1; iy+=2 )
			{
				for( int iz=0; iz<nz_; iz+=2 )
				{
					add_oct_val( ix,iy,iz, -oct_mean(ix,iy,iz) );
					add_oct_val( ix,iy+3*ny_/4,iz, -oct_mean(ix,iy+3*ny_/4,iz) );
				}
			}
		
		#pragma omp parallel for
		for( int ix=0; ix<nx_; ix+=2 )
			for( int iy=0; iy<ny_; iy+=2 )
				for( int iz=0; iz<nz_/4-1; iz+=2 )
				{
					add_oct_val( ix,iy,iz, -oct_mean(ix,iy,iz) );
					add_oct_val( ix,iy,iz+3*nz_/4, -oct_mean(ix,iy,iz+3*nz_/4) );
				}
	}
	
	void zero_boundary( void )
	{
		for( int ix=0; ix<nx_/4; ++ix )
			for( int iy=0; iy<ny_; ++iy )
				for( int iz=0; iz<nz_; ++iz )
			{
				(*this)(ix,iy,iz) = 0.0;
				(*this)(ix+3*nx_/4,iy,iz) = 0.0;
			}
		
		for( int ix=0; ix<nx_; ++ix )
			for( int iy=0; iy<ny_/4; ++iy )
				for( int iz=0; iz<nz_; ++iz )
				{
					(*this)(ix,iy,iz) = 0.0;
					(*this)(ix,iy+3*ny_/4,iz) = 0.0;
				}

		
		for( int ix=0; ix<nx_; ++ix )
			for( int iy=0; iy<ny_; ++iy )
				for( int iz=0; iz<nz_/4; ++iz )
				{
					(*this)(ix,iy,iz) = 0.0;
					(*this)(ix,iy,iz+3*nz_/4) = 0.0;
				}
		
	}
	

};

template< typename M >
inline void enforce_coarse_mean( M& v, M& V )
{
	double outersum =0.0, innersum = 0.0;
	unsigned 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);

	unsigned innercount = 0, outercount = 0;
	for( unsigned i=0; i<V.size(0); ++i )
		for( unsigned j=0; j<V.size(1); ++j )
			for( unsigned k=0; k<V.size(2); ++k )
			{
				if( (i >= ox && i < ox+nx) && 
				    (j >= oy && j < oy+ny) &&
				    (k >= oz && k < oz+nz ))
				{	
					innersum += V(i,j,k);
					++innercount;
				}
				else 
				{
					outersum += V(i,j,k);
					++outercount;
				}

			}
	innersum /= innercount;
	outersum /= outercount;
	
	std::cerr << "***\n";
	
	double dcorr = innersum * innercount / outercount;
	
	for( unsigned i=0; i<V.size(0); ++i )
		for( unsigned j=0; j<V.size(1); ++j )
			for( unsigned k=0; k<V.size(2); ++k )
			{
				if( !((i >= ox && i < ox+nx) && 
				   (j >= oy && j < oy+ny) &&
				   (k >= oz && k < oz+nz )))
					V(i,j,k) -= outersum + dcorr;//innersum;
			}
	
}

template< typename M >
inline void enforce_mean( M& v, M& V )
{

	int 
	nx = v.size(0)/2, 
	ny = v.size(1)/2, 
	nz = v.size(2)/2,
	ox = v.offset(0),
	oy = v.offset(1),
	oz = v.offset(2);

	
	double coarsemean = 0.0, finemean = 0.0;
	unsigned count = 0;
	#pragma omp parallel for reduction(+:coarsemean,finemean,count)
	for( int i=0; i<nx; ++i )
	{
		int i2 = 2*i;
		for( int j=0,j2=0; j<ny; ++j,j2+=2 )
			for( int k=0,k2=0; k<nz; ++k,k2+=2 )
			{
				
				finemean += 0.125 * ( v(i2+1,j2,k2) + v(i2,j2+1,k2) + v(i2,j2,k2+1) + v(i2+1,j2+1,k2) +
											v(i2+1,j2,k2+1) + v(i2+1,j2+1,k2+1) + v(i2,j2+1,k2+1) + v(i2,j2,k2) );
				
				coarsemean += V(i+ox,j+oy,k+oz);
				++count;
			}
	}
	
	coarsemean /= count;
	finemean /= count;
	
	double dmean = coarsemean-finemean;
	dmean = dmean/sqrt(2.0);
	
	std::cerr << " - enforce_mean correction : fine = " << finemean << ", coarse = " << coarsemean << ", diff = " << dmean << std::endl;
	
#pragma omp parallel for reduction(+:coarsemean,finemean)
	for( int i=0; i<nx; ++i )
	{
		int i2 = 2*i;
		for( int j=0,j2=0; j<ny; ++j,j2+=2 )
			for( int k=0,k2=0; k<nz; ++k,k2+=2 )
			{
				
				v(i2+1,j2,k2)		+= dmean;
				v(i2,j2+1,k2)		+= dmean; 
				v(i2,j2,k2+1)		+= dmean;
				v(i2+1,j2+1,k2)		+= dmean;
				v(i2+1,j2,k2+1)		+= dmean;
				v(i2+1,j2+1,k2+1)	+= dmean;
				v(i2,j2+1,k2+1)		+= dmean;
				v(i2,j2,k2)			+= dmean;
			}
	}
	
}
#endif

