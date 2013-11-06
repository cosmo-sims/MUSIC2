/*
 
 mesh.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
*/

#ifndef __MESH_HH
#define __MESH_HH

#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>

#include <math.h>

#include "config_file.hh"
#include "log.hh"


#include "region_generator.hh"

class refinement_mask
{
protected:
    std::vector<short> mask_;
    size_t nx_, ny_, nz_;
    
public:
    
    refinement_mask( void )
    : nx_( 0 ), ny_ ( 0 ), nz_( 0 )
    { }
    
    refinement_mask( size_t nx, size_t ny, size_t nz, short value = 0. )
    : nx_( nx ), ny_( ny ), nz_( nz )
    {
        mask_.assign( nx_*ny_*nz_, value );
    }
    
    refinement_mask( const refinement_mask& r )
    {
        nx_ = r.nx_;
        ny_ = r.ny_;
        nz_ = r.nz_;
        mask_ = r.mask_;
    }
    
    refinement_mask& operator=( const refinement_mask& r )
    {
        nx_ = r.nx_;
        ny_ = r.ny_;
        nz_ = r.nz_;
        mask_ = r.mask_;
        
        return *this;
    }
    
    void init( size_t nx, size_t ny, size_t nz, short value = 0. )
    {
        nx_ = nx;
        ny_ = ny;
        nz_ = nz;
        mask_.assign( nx_*ny_*nz_, value );
    }
    
    const short& operator()( size_t i, size_t j, size_t k ) const
    {
        return mask_[ (i*ny_+j)*nz_+k ];
    }
    
    short& operator()( size_t i, size_t j, size_t k )
    {
        return mask_[ (i*ny_+j)*nz_+k ];
    }
    
    size_t count_flagged( void )
    {
        size_t count = 0;
        for( size_t i=0; i<mask_.size(); ++i )
            if( mask_[i] )
                ++count;
        return count;
    }
    
    size_t count_notflagged( void )
    {
        size_t count = 0;
        for( size_t i=0; i<mask_.size(); ++i )
            if( !mask_[i] )
                ++count;
        return count;
    }
};

//! base class for all things that have rectangular mesh structure
template<typename T>
class Meshvar{
public:
	typedef T real_t;
	
	size_t 
		m_nx,	//!< x-extent of the rectangular mesh
		m_ny,	//!< y-extent of the rectangular mesh
		m_nz;	//!< z-extent of the rectangular mesh
	
	int 
		m_offx, //!< x-offset of the grid (just as a helper, not used inside the class)
		m_offy, //!< y-offset of the grid (just as a helper, not used inside the class)
		m_offz;	//!< z-offset of the grid (just as a helper, not used inside the class)
	
	real_t * m_pdata; //!< pointer to the dynamic data array
	
	//! constructor for cubic mesh
	explicit Meshvar( size_t n, int offx, int offy, int offz )
	: m_nx( n ), m_ny( n ), m_nz( n ), m_offx( offx ), m_offy( offy ), m_offz( offz )
	{
		m_pdata = new real_t[m_nx*m_ny*m_nz];
	}
	
	//! constructor for rectangular mesh
	Meshvar( size_t nx, size_t ny, size_t nz, int offx, int offy, int offz )
	: m_nx( nx ), m_ny( ny ), m_nz( nz ), m_offx( offx ), m_offy( offy ), m_offz( offz )
	{
		m_pdata = new real_t[m_nx*m_ny*m_nz];
	}
	
	//! variant copy constructor with optional copying of the actual data
	Meshvar( const Meshvar<real_t>& m, bool copy_over=true )
	{
		m_nx = m.m_nx;
		m_ny = m.m_ny;
		m_nz = m.m_nz;
		
		m_offx = m.m_offx;
		m_offy = m.m_offy;
		m_offz = m.m_offz;
		
		m_pdata = new real_t[m_nx*m_ny*m_nz];
		
		if( copy_over )
			for( size_t i=0; i<m_nx*m_ny*m_nz; ++i )
				m_pdata[i] = m.m_pdata[i];
	}
	
	//! standard copy constructor
	explicit Meshvar( const Meshvar<real_t>& m )
	{
		m_nx = m.m_nx;
		m_ny = m.m_ny;
		m_nz = m.m_nz;
		
		m_offx = m.m_offx;
		m_offy = m.m_offy;
		m_offz = m.m_offz;
		
		m_pdata = new real_t[m_nx*m_ny*m_nz];
		
		for( size_t i=0; i<m_nx*m_ny*m_nz; ++i )
			m_pdata[i] = m.m_pdata[i];
	}
	
	//! destructor
	~Meshvar()
	{
		if( m_pdata != NULL )
			delete[] m_pdata;
	}
	
	//! deallocate the data, but keep the structure
	inline void deallocate( void )
	{
		if( m_pdata != NULL )
			delete[] m_pdata;
		m_pdata = NULL;
	}
	
	//! get extent of the mesh along a specified dimension (const)
	inline size_t size( unsigned dim ) const
	{
		if( dim == 0 ) return m_nx;
		if( dim == 1 ) return m_ny;
		return m_nz;
	}
	
	//! get offset of the mesh along a specified dimension  (const)
	inline int offset( unsigned dim ) const
	{
		if( dim == 0 ) return m_offx;
		if( dim == 1 ) return m_offy;
		return m_offz;
	}
	
	//! get extent of the mesh along a specified dimension
	inline int& offset( unsigned dim )
	{
		if( dim == 0 ) return m_offx;
		if( dim == 1 ) return m_offy;
		return m_offz;
	}
	
	//! set all the data to zero values
	void zero( void )
	{
		for( size_t i=0; i<m_nx*m_ny*m_nz; ++i )
			m_pdata[i] = 0.0;
	}
	
	//! direct array random acces to the data block
	inline real_t * operator[]( const size_t i )
	{	return &m_pdata[i];	}

	//! direct array random acces to the data block (const)
	inline const real_t * operator[]( const size_t i ) const
	{	return &m_pdata[i];	}
	
	//! 3D random access to the data block via index 3-tuples
	inline real_t& operator()(const int ix, const int iy, const int iz )
	{
#ifdef DEBUG
        if( ix<0||ix>=(int)m_nx||iy<0||iy>=(int)m_ny||iz<0||iz>=(int)m_nz)
            LOGERR("Array index (%d,%d,%d) out of bounds",ix,iy,iz);
#endif
        
        return m_pdata[ ((size_t)ix*m_ny+(size_t)iy)*m_nz + (size_t)iz ];
    }
	
	//! 3D random access to the data block via index 3-tuples (const)
	inline const real_t& operator()(const int ix, const int iy, const int iz ) const
	{
#ifdef DEBUG
        if( ix<0||ix>=(int)m_nx||iy<0||iy>=(int)m_ny||iz<0||iz>=(int)m_nz)
            LOGERR("Array index (%d,%d,%d) out of bounds",ix,iy,iz);
#endif
        
        return m_pdata[ ((size_t)ix*m_ny+(size_t)iy)*m_nz + (size_t)iz ];
    }
	
	//! direct multiplication of the whole data block with a number
	Meshvar<real_t>& operator*=( real_t x )
	{
		for( size_t i=0; i<m_nx*m_ny*m_nz; ++i )
			m_pdata[i] *= x;
		return *this;
	}
	
	//! direct addition of a number to the whole data block
	Meshvar<real_t>& operator+=( real_t x )
	{
		for( size_t i=0; i<m_nx*m_ny*m_nz; ++i )
			m_pdata[i] += x;
		return *this;
	}

	//! direct element-wise division of the whole data block by a number
	Meshvar<real_t>& operator/=( real_t x )
	{
		for( size_t i=0; i<m_nx*m_ny*m_nz; ++i )
			m_pdata[i] /= x;
		return *this;
	}
	
	
	//! direct subtraction of a number from the whole data block
	Meshvar<real_t>& operator-=( real_t x )
	{
		for( size_t i=0; i<m_nx*m_ny*m_nz; ++i )
			m_pdata[i] -= x;
		return *this;
	}
	
	//! direct element-wise multiplication with another compatible mesh
	Meshvar<real_t>& operator*=( const Meshvar<real_t>& v )
	{
		if( v.m_nx*v.m_ny*v.m_nz != m_nx*m_ny*m_nz )
		{
            LOGERR("Meshvar::operator*= : attempt to operate on incompatible data");
            throw std::runtime_error("Meshvar::operator*= : attempt to operate on incompatible data");
		}
		for( size_t i=0; i<m_nx*m_ny*m_nz; ++i )
			m_pdata[i] *= v.m_pdata[i];
		
		return *this;
	}
	
	//! direct element-wise division with another compatible mesh
	Meshvar<real_t>& operator/=( const Meshvar<real_t>& v )
	{
		if( v.m_nx*v.m_ny*v.m_nz != m_nx*m_ny*m_nz )
		{
            LOGERR("Meshvar::operator/= : attempt to operate on incompatible data");
            throw std::runtime_error("Meshvar::operator/= : attempt to operate on incompatible data");
		}
        
		for( size_t i=0; i<m_nx*m_ny*m_nz; ++i )
			m_pdata[i] /= v.m_pdata[i];
		
		return *this;
	}
	
	//! direct element-wise addition of another compatible mesh
	Meshvar<real_t>& operator+=( const Meshvar<real_t>& v )
	{
		if( v.m_nx*v.m_ny*v.m_nz != m_nx*m_ny*m_nz )
		{
            LOGERR("Meshvar::operator+= : attempt to operate on incompatible data");
            throw std::runtime_error("Meshvar::operator+= : attempt to operate on incompatible data");
		}
		for( size_t i=0; i<m_nx*m_ny*m_nz; ++i )
			m_pdata[i] += v.m_pdata[i];
		
		return *this;
	}
	
	//! direct element-wise subtraction of another compatible mesh
	Meshvar<real_t>& operator-=( const Meshvar<real_t>& v )
	{
		if( v.m_nx*v.m_ny*v.m_nz != m_nx*m_ny*m_nz )
		{
            LOGERR("Meshvar::operator-= : attempt to operate on incompatible data");
            throw std::runtime_error("Meshvar::operator-= : attempt to operate on incompatible data");
		}
		for( size_t i=0; i<m_nx*m_ny*m_nz; ++i )
			m_pdata[i] -= v.m_pdata[i];
		
		return *this;
	}
	
	//! assignment operator for rectangular meshes
	Meshvar<real_t>& operator=( const Meshvar<real_t>& m )
	{
		m_nx = m.m_nx;
		m_ny = m.m_ny;
		m_nz = m.m_nz;
		
		m_offx = m.m_offx;
		m_offy = m.m_offy;
		m_offz = m.m_offz;
		
		if( m_pdata != NULL )
			delete m_pdata;
		
		m_pdata = new real_t[m_nx*m_ny*m_nz];
		
		for( size_t i=0; i<m_nx*m_ny*m_nz; ++i )
			m_pdata[i] = m.m_pdata[i];
		
		return *this;
	}
	
	real_t* get_ptr( void )
	{	return m_pdata;		}
};

//! MeshvarBnd derived class adding boundary ghost cell functionality
template< typename T >
class MeshvarBnd : public Meshvar< T >{
	using Meshvar<T>::m_nx;
	using Meshvar<T>::m_ny;
	using Meshvar<T>::m_nz;
	using Meshvar<T>::m_pdata;
	
	
	
public:
	typedef T real_t;
	
	//! number of boundary (ghost) cells
	int m_nbnd;		
	
	//! most general constructor
	MeshvarBnd( int nbnd, size_t nx, size_t ny, size_t nz, size_t xoff, size_t yoff, size_t zoff )
	: Meshvar<real_t>( nx+2*nbnd, ny+2*nbnd, nz+2*nbnd, xoff, yoff, zoff ), m_nbnd( nbnd )
	{ 	}

	//! zero-offset constructor
	MeshvarBnd( size_t nbnd, size_t nx, size_t ny, size_t nz )
	: Meshvar<real_t>( nx+2*nbnd, ny+2*nbnd, nz+2*nbnd, 0, 0, 0 ), m_nbnd( nbnd )
	{ 	}
	
	//! constructor for cubic meshes
	MeshvarBnd( size_t nbnd, size_t n, size_t xoff, size_t yoff, size_t zoff )
	: Meshvar<real_t>( n+2*nbnd, xoff, yoff, zoff ), m_nbnd( nbnd )
	{ 	}
	
	//! constructor for cubic meshes with zero offset
	MeshvarBnd( size_t nbnd, size_t n )
	: Meshvar<real_t>( n+2*nbnd, 0, 0, 0 ), m_nbnd( nbnd )
	{ 	}
	
	//! modified copy constructor, allows to avoid copying actual data
	MeshvarBnd( const MeshvarBnd<real_t>& v, bool copyover )
	: Meshvar<real_t>( v, copyover ), m_nbnd( v.m_nbnd )
	{   }
	
	//! copy constructor
	explicit MeshvarBnd( const MeshvarBnd<real_t>& v )
	: Meshvar<real_t>( v, true ), m_nbnd( v.m_nbnd )
	{   }
	
	//! get extent of the mesh along a specified dimension
	inline size_t size( unsigned dim=0 ) const
	{
		if( dim == 0 ) return m_nx-2*m_nbnd;
		if( dim == 1 ) return m_ny-2*m_nbnd;
		return m_nz-2*m_nbnd;
	}

	//! 3D random access to the data block via index 3-tuples
	inline real_t& operator()(const int ix, const int iy, const int iz )
	{
		size_t iix(ix+m_nbnd), iiy(iy+m_nbnd), iiz(iz+m_nbnd);
		return m_pdata[ (iix*m_ny+iiy)*m_nz + iiz ];
	}
	
	//! 3D random access to the data block via index 3-tuples (const)
	inline const real_t& operator()(const int ix, const int iy, const int iz ) const
	{
		size_t iix(ix+m_nbnd), iiy(iy+m_nbnd), iiz(iz+m_nbnd);
		return m_pdata[ (iix*m_ny+iiy)*m_nz + iiz ];
	}
	
	//! assignment operator for rectangular meshes with ghost zones
	MeshvarBnd<real_t>& operator=( const MeshvarBnd<real_t>& m )
	{
		if( this->m_nx != m.m_nx || this->m_ny != m.m_ny || this->m_nz != m.m_nz )
		{
			this->m_nx = m.m_nx;
			this->m_ny = m.m_ny;
			this->m_nz = m.m_nz;
			
			if( m_pdata != NULL )
				delete[] m_pdata;
			
			m_pdata = new real_t[m_nx*m_ny*m_nz];			
		}
		
		for( size_t i=0; i<m_nx*m_ny*m_nz; ++i )
			this->m_pdata[i] = m.m_pdata[i];
		
		return *this;
	}

	//! sets the value of all ghost zones to zero
	void zero_bnd( void )
	{
		
		int nx,ny,nz;
		nx = this->size(0);
		ny = this->size(1);
		nz = this->size(2);
		
		
		for( int j=-m_nbnd; j<ny+m_nbnd; ++j )
			for( int k=-m_nbnd; k<nz+m_nbnd; ++k ){
				for( int i=-m_nbnd;i<0;++i )
				{
					(*this)(i,j,k)  = 0.0;
					(*this)(nx-1-i,j,k) = 0.0;	
				}
				
			}
		
		for( int i=-m_nbnd; i<nx+m_nbnd; ++i )
			for( int k=-m_nbnd; k<nz+m_nbnd; ++k ){
				for( int j=-m_nbnd;j<0;++j )
				{
					(*this)(i,j,k) = 0.0;
					(*this)(i,ny-j-1,k) = 0.0;
				}
				
			}
		
		for( int i=-m_nbnd; i<nx+m_nbnd; ++i )
			for( int j=-m_nbnd; j<ny+m_nbnd; ++j ){
				for( int k=-m_nbnd;k<0;++k )
				{
					(*this)(i,j,k) = 0.0;
					(*this)(i,j,nz-k-1) = 0.0;	
				}
				
			}
	}
	
	//! outputs the data, for debugging only, not practical for large datasets
	void print( void ) const
	{
		int nbnd = m_nbnd;
		
		std::cout << "size is [" << this->size(0) << ", " << this->size(1) << ", " << this->size(2) << "]\n";
		std::cout << "ghost region has length of " << nbnd << std::endl;
		
		std::cout.precision(3);
		for(int i=-nbnd; i<(int)this->size(0)+nbnd; ++i )
		{
			std::cout << "ix = " << i << ": \n";
			
			for (int j=-nbnd; j<(int)this->size(1)+nbnd; ++j) {
				for (int k=-nbnd; k<(int)this->size(2)+nbnd; ++k) {
					if( i<0||i>=this->size(0)||j<0||j>=this->size(1)||k<0||k>=this->size(2))
						std::cout << "[" << std::setw(6) << (*this)(i,j,k) << "] ";
					else
						std::cout << std::setw(8) <<  (*this)(i,j,k) << " ";
				}
				std::cout << std::endl;
			}
			
			std::cout << std::endl;
			
		}
	}
};



//! class that subsumes a nested grid collection
template< typename T >
class GridHierarchy
{
public:
	
	//! number of ghost cells on boundary
	size_t m_nbnd;
	
	//! highest level without adaptive refinement
	unsigned m_levelmin;
	
	//! vector of pointers to the underlying rectangular mesh data for each level
	std::vector< MeshvarBnd<T>* > m_pgrids;
	
	std::vector<int> 
		m_xoffabs,		//!< vector of x-offsets of a level mesh relative to the coarser level
		m_yoffabs,		//!< vector of x-offsets of a level mesh relative to the coarser level
		m_zoffabs;		//!< vector of x-offsets of a level mesh relative to the coarser level
    
    std::vector< refinement_mask* > m_ref_masks;
    bool bhave_refmask;
	
protected:
	
	//! check whether a given grid has identical hierarchy, dimensions to this 
	bool is_consistent( const GridHierarchy<T>& gh )
	{
		if( gh.levelmax()!=levelmax() )
			return false;
		
		if( gh.levelmin()!=levelmin() )
			return false;
		
		for( unsigned i=levelmin(); i<=levelmax(); ++i )
			for( int j=0; j<3; ++j )
			{
				if( size(i,j) != gh.size(i,j) )
					return false;
				if( offset(i,j) != gh.offset(i,j) )
					return false;
			}
		
		return true;
	}
	
public:
	
	
	//! return a pointer to the MeshvarBnd object representing data for one level
	MeshvarBnd<T> *get_grid( unsigned ilevel )
	{	

		if( ilevel >= m_pgrids.size() )
		{
			LOGERR("Attempt to access level %d but maxlevel = %d", ilevel, m_pgrids.size()-1);
			throw std::runtime_error("Fatal: attempt to access non-existent grid");
		}
		return m_pgrids[ilevel];  
	}

	//! return a pointer to the MeshvarBnd object representing data for one level (const)	
	const MeshvarBnd<T> *get_grid( unsigned ilevel ) const
	{	
		if( ilevel >= m_pgrids.size() )
		{
            LOGERR("Attempt to access level %d but maxlevel = %d", ilevel, m_pgrids.size()-1 );
			throw std::runtime_error("Fatal: attempt to access non-existent grid");
		}

		return m_pgrids[ilevel];  
	}
	
	
	//! constructor for a collection of rectangular grids representing a multi-level hierarchy
	/*! creates an empty hierarchy, levelmin is initially zero, no grids are stored
	 * @param nbnd number of ghost zones added at the boundary
	 */
	explicit GridHierarchy( size_t nbnd )
	: m_nbnd( nbnd ), m_levelmin( 0 ), bhave_refmask( false )
	{
		m_pgrids.clear();
	}
	
	//! copy constructor
	explicit GridHierarchy( const GridHierarchy<T> & gh )
	{
		for( unsigned i=0; i<=gh.levelmax(); ++i )
			m_pgrids.push_back( new MeshvarBnd<T>( *gh.get_grid(i) ) );
		
		m_nbnd = gh.m_nbnd;
		m_levelmin = gh.m_levelmin;
		
		m_xoffabs = gh.m_xoffabs;
		m_yoffabs = gh.m_yoffabs;
		m_zoffabs = gh.m_zoffabs;
        
        //ref_mask   = gh.ref_mask;
        bhave_refmask = gh.bhave_refmask;
        
        if( bhave_refmask )
        {
            for( size_t i=0; i<gh.m_ref_masks.size(); ++i )
                m_ref_masks.push_back( new refinement_mask( *(gh.m_ref_masks[i]) ) );
        }
	}
	
	//! destructor
	~GridHierarchy()
	{
		this->deallocate();
	}
	
	//! free all memory occupied by the grid hierarchy
	void deallocate()
	{
		for( unsigned i=0; i<m_pgrids.size(); ++i )
			delete m_pgrids[i];
		m_pgrids.clear();
		std::vector< MeshvarBnd<T>* >().swap( m_pgrids );
		
		m_xoffabs.clear();
		m_yoffabs.clear();
		m_zoffabs.clear();
		m_levelmin = 0;
        
        for( size_t i=0; i<m_ref_masks.size(); ++i )
            delete m_ref_masks[i];
        m_ref_masks.clear();
	}
    
    void add_refinement_mask( const double *shift )
    {
        bhave_refmask = false;
        
        //! generate a mask
        if( m_levelmin != levelmax() )
        {
            for( int ilevel = (int)levelmax()-1; ilevel >= (int)levelmin(); --ilevel )
            {
                double xq[3], dx = 1.0/(1ul<<ilevel); //1.0/(1ul<<(levelmax()-1));
                
                m_ref_masks[ilevel]->init( size(ilevel,0), size(ilevel,1), size(ilevel,2), 0.0 );
                
                for( size_t i=0; i<size(ilevel,0); i++ )
                {
                    xq[0] = (offset_abs(ilevel,0) + i)*dx + 0.5*dx + shift[0];
                    for( size_t j=0; j<size(ilevel,1); ++j )
                    {
                        xq[1] = (offset_abs(ilevel,1) + j)*dx + 0.5*dx + shift[1];
                        for( size_t k=0; k<size(ilevel,2); ++k )
                        {
                            xq[2] = (offset_abs(ilevel,2) + k)*dx + 0.5*dx + shift[2];
                            
                            if( the_region_generator->query_point( xq, ilevel ) )
                                (*m_ref_masks[ilevel])(i,j,k) = 1.0;//true;
                            else
                                (*m_ref_masks[ilevel])(i,j,k) = 0.0;
                        }
                    }
                }
            }
            
            bhave_refmask = true;
            
            for( int ilevel = (int)levelmax()-1; ilevel >= (int)levelmin(); --ilevel )
            {
                for( size_t i=0; i<size(ilevel,0); i++ )
                {
                    for( size_t j=0; j<size(ilevel,1); ++j )
                    {
                        for( size_t k=0; k<size(ilevel,2); ++k )
                        {
                            int off[] = { offset(ilevel+1,0), offset(ilevel+1,1), offset(ilevel+1,2) };
                            
                            bool fine_is_flagged = false;
                            
                            /*fine_is_flagged |= (*m_ref_masks[ilevel+1])(off[0]+2*i+0,off[1]+2*j+0,off[2]+2*k+0);
                            fine_is_flagged |= (*m_ref_masks[ilevel+1])(off[0]+2*i+0,off[1]+2*j+0,off[2]+2*k+1);
                            fine_is_flagged |= (*m_ref_masks[ilevel+1])(off[0]+2*i+0,off[1]+2*j+1,off[2]+2*k+0);
                            fine_is_flagged |= (*m_ref_masks[ilevel+1])(off[0]+2*i+0,off[1]+2*j+1,off[2]+2*k+1);
                            fine_is_flagged |= (*m_ref_masks[ilevel+1])(off[0]+2*i+1,off[1]+2*j+0,off[2]+2*k+0);
                            fine_is_flagged |= (*m_ref_masks[ilevel+1])(off[0]+2*i+1,off[1]+2*j+0,off[2]+2*k+1);
                            fine_is_flagged |= (*m_ref_masks[ilevel+1])(off[0]+2*i+1,off[1]+2*j+1,off[2]+2*k+0);
                            fine_is_flagged |= (*m_ref_masks[ilevel+1])(off[0]+2*i+1,off[1]+2*j+1,off[2]+2*k+1);*/
                            
                            int ifine[] = {
                                2*(int)i-2*off[0],
                                2*(int)j-2*off[1],
                                2*(int)k-2*off[2],
                            };
                            
                            if( ifine[0]>=0 && ifine[0] < size(ilevel+1,0) &&
                                ifine[1]>=0 && ifine[1] < size(ilevel+1,1) &&
                                ifine[2]>=0 && ifine[2] < size(ilevel+1,2) )
                            {
                                fine_is_flagged |= (*m_ref_masks[ilevel+1])(ifine[0]+0,ifine[1]+0,ifine[2]+0);
                                fine_is_flagged |= (*m_ref_masks[ilevel+1])(ifine[0]+0,ifine[1]+0,ifine[2]+1);
                                fine_is_flagged |= (*m_ref_masks[ilevel+1])(ifine[0]+0,ifine[1]+1,ifine[2]+0);
                                fine_is_flagged |= (*m_ref_masks[ilevel+1])(ifine[0]+0,ifine[1]+1,ifine[2]+1);
                                fine_is_flagged |= (*m_ref_masks[ilevel+1])(ifine[0]+1,ifine[1]+0,ifine[2]+0);
                                fine_is_flagged |= (*m_ref_masks[ilevel+1])(ifine[0]+1,ifine[1]+0,ifine[2]+1);
                                fine_is_flagged |= (*m_ref_masks[ilevel+1])(ifine[0]+1,ifine[1]+1,ifine[2]+0);
                                fine_is_flagged |= (*m_ref_masks[ilevel+1])(ifine[0]+1,ifine[1]+1,ifine[2]+1);
                                
                            }
                            
                            if( fine_is_flagged )
                                (*m_ref_masks[ilevel])(i,j,k) = 0.0;
                        }
                    }
                }
            }
            
        }
    }
	
	
	//! get offset of a grid at specified refinement level
	/*! the offset describes the shift of a refinement grid with respect to its coarser parent grid
	 *  @param ilevel the level for which the offset is to be determined
	 *  @param idim the dimension along which the offset is to be determined
	 *  @return integer value denoting the offset in units of coarse grid cells
	 *  @sa offset_abs
	 */
	int offset( int ilevel, int idim ) const
	{
		return m_pgrids[ilevel]->offset(idim);
	}
	
	//! get size of a grid at specified refinement level
	/*! the size describes the number of cells along one dimension of a grid
	 *  @param ilevel the level for which the size is to be determined
	 *  @param idim the dimension along which the size is to be determined
	 *  @return integer value denoting the size of refinement grid at level ilevel along dimension idim
	 */
	size_t size( int ilevel, int idim ) const
	{
		return m_pgrids[ilevel]->size(idim);
	}
	
	
	//! get the absolute offset of a grid at specified refinement level
	/*! the absolute offset describes the shift of a refinement grid with respect to the simulation volume
	 *  @param ilevel the level for which the offset is to be determined
	 *  @param idim the dimension along which the offset is to be determined
	 *  @return integer value denoting the absolute offset in units of fine grid cells
	 *  @sa offset
	 */
	int offset_abs( int ilevel, int idim ) const
	{
		if( idim == 0 ) return m_xoffabs[ilevel];
		if( idim == 1 ) return m_yoffabs[ilevel];
		return m_zoffabs[ilevel];
	}
	
	
	//! get the coordinate posisition of a grid cell 
	/*! returns the position of a grid cell at specified level relative to the simulation volume
	 *  @param ilevel the refinement level of the grid cell
	 *  @param i the x-index of the cell in the level grid
	 *  @param j the y-index of the cell in the level grid
	 *  @param k the z-index of the cell in the level grid
	 *  @param ppos pointer to a double[3] array to which the coordinates are written
	 *  @return none
	 */
	void cell_pos( unsigned ilevel, int i, int j, int k, double* ppos ) const
	{
		double h = 1.0/(1<<ilevel);//, htop = h*2.0;
		ppos[0] = h*((double)offset_abs(ilevel,0)+(double)i+0.5);
		ppos[1] = h*((double)offset_abs(ilevel,1)+(double)j+0.5);
		ppos[2] = h*((double)offset_abs(ilevel,2)+(double)k+0.5);
		
		if( ppos[0] >= 1.0 || ppos[1] >= 1.0 || ppos[2] >= 1.0 )
			std::cerr << " - Cell seems outside domain! : (" << ppos[0] << ", " << ppos[1] << ", " << ppos[2] << "\n";
	}
    
    
    //! get the bounding box of a grid in code units
    /*! returns the bounding box of a grid at specified level in code units
     *  @param ilevel the refinement level of the grid
     *  @param left pointer to a double[3] array to which the left corner is written
     *  @param right pointer to a double[3] array to which the right corner is written
     *  @return none
     */
    void grid_bbox( unsigned ilevel, double *left, double *right ) const
    {
        double h = 1.0/(1<<ilevel);//, htop = h*2.0;
		left[0] = h*((double)offset_abs(ilevel,0));
		left[1] = h*((double)offset_abs(ilevel,1));
		left[2] = h*((double)offset_abs(ilevel,2));
        
        right[0] = left[0] + h*((double)size(ilevel,0));
        right[1] = left[1] + h*((double)size(ilevel,1));
        right[2] = left[2] + h*((double)size(ilevel,2));
    }
	
	//! checks whether a given grid cell is refined
	/*! a grid cell counts as refined if it is divided into 8 cells at the next higher level
	 *  @param ilevel the refinement level of the grid cell
	 *  @param i the x-index of the cell in the level grid
	 *  @param j the y-index of the cell in the level grid
	 *  @param k the z-index of the cell in the level grid
	 *  @return true if cell is refined, false otherwise
	 */
	bool is_refined( unsigned ilevel, int i, int j, int k ) const
	{
        //if( ilevel == levelmax() )
        //    return false;

        if( bhave_refmask )
            return !(*m_ref_masks[ilevel-1])(offset(ilevel,0)+i/2,offset(ilevel,1)+j/2,offset(ilevel,2)+k/2);

        if( ilevel == levelmax()-1 && bhave_refmask )
            return (*m_ref_masks[ilevel])(i,j,k);
		
		if( i < offset(ilevel+1,0) || i >= offset(ilevel+1, 0)+(int)size(ilevel+1,0)/2 ||
		    j < offset(ilevel+1,1) || j >= offset(ilevel+1, 1)+(int)size(ilevel+1,1)/2 ||
		    k < offset(ilevel+1,2) || k >= offset(ilevel+1, 2)+(int)size(ilevel+1,2)/2 )
			return false;
		
		return true;
	}
	
	//! sets the values of all grids on all levels to zero
	void zero( void )
	{
		for( unsigned i=0; i<m_pgrids.size(); ++i )
			m_pgrids[i]->zero();
	}
	
	
	//! count the number of cells that are not further refined (=leafs)
	/*! for allocation purposes it is useful to query the number of cells to be expected
	 *  @param lmin the minimum refinement level to consider
	 *  @param lmax the maximum refinement level to consider
	 *  @return the integer number of cells between lmin and lmax that are not further refined
	 */
	size_t count_leaf_cells( unsigned lmin, unsigned lmax ) const
	{
		size_t npcount = 0;
		
		for( int ilevel=lmax; ilevel>=(int)lmin; --ilevel )
			for( unsigned i=0; i<get_grid(ilevel)->size(0); ++i )
				for( unsigned j=0; j<get_grid(ilevel)->size(1); ++j )
					for( unsigned k=0; k<get_grid(ilevel)->size(2); ++k )
						if( ! is_refined(ilevel,i,j,k) )
							++npcount;
		
		return npcount;
	}
	
	//! count the number of cells that are not further refined (=leafs)
	/*! for allocation purposes it is useful to query the number of cells to be expected
	 *  @return the integer number of cells in the hierarchy that are not further refined
	 */
	size_t count_leaf_cells( void ) const
	{
		return count_leaf_cells( levelmin(), levelmax() );
	}
	
	//! creates a hierarchy of coextensive grids, refined by factors of 2
	/*! creates a hierarchy of lmax grids, each extending over the whole simulation volume with
	 *  grid length 2^n for level 0<=n<=lmax
	 *  @param lmax the maximum refinement level to be added (sets the resolution to 2^lmax for each dim)
	 *  @return none
	 */
	void create_base_hierarchy( unsigned lmax )
	{
		size_t n=1;
		
		this->deallocate();
		
		m_pgrids.clear();
		
		m_xoffabs.clear();
		m_yoffabs.clear();
		m_zoffabs.clear();
		
		for( unsigned i=0; i<= lmax; ++i )
		{
			//std::cout << "....adding level " << i << " (" << n << ", " << n << ", " << n << ")" << std::endl;
			m_pgrids.push_back( new MeshvarBnd<T>( m_nbnd, n, n, n, 0, 0, 0 ) );
			m_pgrids[i]->zero();
     		n *= 2;
			
			m_xoffabs.push_back( 0 );
			m_yoffabs.push_back( 0 );
			m_zoffabs.push_back( 0 );
		}
		
		m_levelmin = lmax;
        
        for( unsigned i=0; i<= lmax; ++i )
            m_ref_masks.push_back( new refinement_mask(size(i,0),size(i,1),size(i,2),(short)(i!=lmax)) );
	}
	
	//! multiply entire grid hierarchy by a constant
	GridHierarchy<T>& operator*=( T x )
	{
		for( unsigned i=0; i<m_pgrids.size(); ++i )
			(*m_pgrids[i]) *= x;
		return *this;
	}
	
	//! divide entire grid hierarchy by a constant
	GridHierarchy<T>& operator/=( T x )
	{
		for( unsigned i=0; i<m_pgrids.size(); ++i )
			(*m_pgrids[i]) /= x;
		return *this;
	}
	
	//! add a constant to the entire grid hierarchy
	GridHierarchy<T>& operator+=( T x )
	{
		for( unsigned i=0; i<m_pgrids.size(); ++i )
			(*m_pgrids[i]) += x;
		return *this;
	}
	
	//! subtract a constant from the entire grid hierarchy
	GridHierarchy<T>& operator-=( T x )
	{
		for( unsigned i=0; i<m_pgrids.size(); ++i )
			(*m_pgrids[i]) -= x;
		return *this;
	}
	
	//! multiply (element-wise) two grid hierarchies
	GridHierarchy<T>& operator*=( const GridHierarchy<T>& gh )
	{
		if( !is_consistent(gh) )
		{
            LOGERR("GridHierarchy::operator*= : attempt to operate on incompatible data");
            throw std::runtime_error("GridHierarchy::operator*= : attempt to operate on incompatible data");
		}
		for( unsigned i=0; i<m_pgrids.size(); ++i )
			(*m_pgrids[i]) *= *gh.get_grid(i);
		return *this;
	}
	
	//! divide (element-wise) two grid hierarchies
	GridHierarchy<T>& operator/=( const GridHierarchy<T>& gh )
	{
		if( !is_consistent(gh) )
		{
            LOGERR("GridHierarchy::operator/= : attempt to operate on incompatible data");
            throw std::runtime_error("GridHierarchy::operator/= : attempt to operate on incompatible data");
		}
		for( unsigned i=0; i<m_pgrids.size(); ++i )
			(*m_pgrids[i]) /= *gh.get_grid(i);
		return *this;
	}
	
	//! add (element-wise) two grid hierarchies
	GridHierarchy<T>& operator+=( const GridHierarchy<T>& gh )
	{
		if( !is_consistent(gh) )
			throw std::runtime_error("GridHierarchy::operator+= : attempt to operate on incompatible data");
		
		for( unsigned i=0; i<m_pgrids.size(); ++i )
			(*m_pgrids[i]) += *gh.get_grid(i);
		return *this;
	}
	
	//! subtract (element-wise) two grid hierarchies
	GridHierarchy<T>& operator-=( const GridHierarchy<T>& gh )
	{
		if( !is_consistent(gh) )
		{
            LOGERR("GridHierarchy::operator-= : attempt to operate on incompatible data");
            throw std::runtime_error("GridHierarchy::operator-= : attempt to operate on incompatible data");
		}
		for( unsigned i=0; i<m_pgrids.size(); ++i )
			(*m_pgrids[i]) -= *gh.get_grid(i);
		return *this;
	}
	
	//! assign (element-wise) two grid hierarchies
	GridHierarchy<T>& operator=( const GridHierarchy<T>& gh )
	{
        bhave_refmask = gh.bhave_refmask;

        if( bhave_refmask )
        {
            for( unsigned i=0; i<=gh.levelmax(); ++i )
                m_ref_masks.push_back( new refinement_mask( *(gh.m_ref_masks[i]) ) );
        }
      
		if( !is_consistent(gh) )
		{
			for( unsigned i=0; i<m_pgrids.size(); ++i )
				delete m_pgrids[i];
			m_pgrids.clear();
			
			for( unsigned i=0; i<=gh.levelmax(); ++i )
				m_pgrids.push_back( new MeshvarBnd<T>( *gh.get_grid(i) ) );
			m_levelmin = gh.levelmin();
			m_nbnd = gh.m_nbnd;
			
			m_xoffabs = gh.m_xoffabs;
			m_yoffabs = gh.m_yoffabs;
			m_zoffabs = gh.m_zoffabs;
			
			
			return *this;
		}//throw std::runtime_error("GridHierarchy::operator= : attempt to operate on incompatible data");
		
		for( unsigned i=0; i<m_pgrids.size(); ++i )
			(*m_pgrids[i]) = *gh.get_grid(i);
		return *this;
	}
	
	/*
	//! assignment operator
	GridHierarchy& operator=( const GridHierarchy<T>& gh )
	{
		for( unsigned i=0; i<m_pgrids.size(); ++i )
			delete m_pgrids[i];
		m_pgrids.clear();
		
		for( unsigned i=0; i<=gh.levelmax(); ++i )
			m_pgrids.push_back( new MeshvarBnd<T>( *gh.get_grid(i) ) );
		m_levelmin = gh.levelmin();
		m_nbnd = gh.m_nbnd;
		
		m_xoffabs = gh.m_xoffabs;
		m_yoffabs = gh.m_yoffabs;
		m_zoffabs = gh.m_zoffabs;
		
		
		return *this;
	}
	*/
	
	/*! add a new refinement patch to the so-far finest level
	 * @param xoff x-offset in units of the coarser grid (finest grid before adding new patch)
	 * @param yoff y-offset in units of the coarser grid (finest grid before adding new patch)
	 * @param zoff z-offset in units of the coarser grid (finest grid before adding new patch)
	 * @param nx x-extent in fine grid cells
	 * @param ny y-extent in fine grid cells
	 * @param nz z-extent in fine grid cells
	 */
	void add_patch( unsigned xoff, unsigned yoff, unsigned zoff, unsigned nx, unsigned ny, unsigned nz )
	{
		m_pgrids.push_back( new MeshvarBnd<T>( m_nbnd, nx, ny, nz, xoff, yoff, zoff ) );
		m_pgrids.back()->zero();
		
		//.. add absolute offsets (in units of current level grid cells)
		m_xoffabs.push_back( 2*(m_xoffabs.back() + xoff) );
		m_yoffabs.push_back( 2*(m_yoffabs.back() + yoff) );
		m_zoffabs.push_back( 2*(m_zoffabs.back() + zoff) );
        
        m_ref_masks.push_back( new refinement_mask(nx,ny,nz,0.0) );
	}
	
	/*! cut a refinement patch to the specified size
	 * @param ilevel grid level on which to perform the size adjustment
	 * @param xoff new x-offset in units of the coarser grid (finest grid before adding new patch)
	 * @param yoff new y-offset in units of the coarser grid (finest grid before adding new patch)
	 * @param zoff new z-offset in units of the coarser grid (finest grid before adding new patch)
	 * @param nx new x-extent in fine grid cells
	 * @param ny new y-extent in fine grid cells
	 * @param nz new z-extent in fine grid cells
	 */
	void cut_patch( unsigned ilevel, unsigned xoff, unsigned yoff, unsigned zoff, unsigned nx, unsigned ny, unsigned nz)
	{
		unsigned dx,dy,dz,dxtop,dytop,dztop;
		
		dx = xoff-m_xoffabs[ilevel];
		dy = yoff-m_yoffabs[ilevel];
		dz = zoff-m_zoffabs[ilevel];

		assert( dx%2==0 && dy%2==0 && dz%2==0 );
		
		dxtop = m_pgrids[ilevel]->offset(0)+dx/2;
		dytop = m_pgrids[ilevel]->offset(1)+dy/2;
		dztop = m_pgrids[ilevel]->offset(2)+dz/2;
		
		MeshvarBnd<T> *mnew = new MeshvarBnd<T>( m_nbnd, nx, ny, nz, dxtop, dytop, dztop );
		
		//... copy data
		for( unsigned i=0; i<nx; ++i )
			for( unsigned j=0; j<ny; ++j )
				for( unsigned k=0; k<nz; ++k )
					(*mnew)(i,j,k) = (*m_pgrids[ilevel])(i+dx,j+dy,k+dz);

		//... replace in hierarchy
		delete m_pgrids[ilevel];
		m_pgrids[ilevel] = mnew;
		
		//... update offsets
		m_xoffabs[ilevel] += dx;
		m_yoffabs[ilevel] += dy;
		m_zoffabs[ilevel] += dz;
		
		if( ilevel < levelmax() )
		{
			m_pgrids[ilevel+1]->offset(0) -= dx;
			m_pgrids[ilevel+1]->offset(1) -= dy;
			m_pgrids[ilevel+1]->offset(2) -= dz;
		}
		
		find_new_levelmin();
	}

  void cut_patch_enforce_top_density( unsigned ilevel, unsigned xoff, unsigned yoff, unsigned zoff, unsigned nx, unsigned ny, unsigned nz)
  {
    unsigned dx,dy,dz,dxtop,dytop,dztop;
		
    dx = xoff-m_xoffabs[ilevel];
    dy = yoff-m_yoffabs[ilevel];
    dz = zoff-m_zoffabs[ilevel];
    
    assert( dx%2==0 && dy%2==0 && dz%2==0 );
    
    dxtop = m_pgrids[ilevel]->offset(0)+dx/2;
    dytop = m_pgrids[ilevel]->offset(1)+dy/2;
    dztop = m_pgrids[ilevel]->offset(2)+dz/2;
    
    MeshvarBnd<T> *mnew = new MeshvarBnd<T>( m_nbnd, nx, ny, nz, dxtop, dytop, dztop );
    
    double coarsesum = 0.0, finesum = 0.0;
    size_t coarsecount = 0, finecount = 0;

    //... copy data
    for( unsigned i=0; i<nx; ++i )
      for( unsigned j=0; j<ny; ++j )
	for( unsigned k=0; k<nz; ++k ){
	  (*mnew)(i,j,k) = (*m_pgrids[ilevel])(i+dx,j+dy,k+dz);
	  finesum += (*mnew)(i,j,k);
	  finecount++;
	}

    //... replace in hierarchy
    delete m_pgrids[ilevel];
    m_pgrids[ilevel] = mnew;
    
    //... update offsets
    m_xoffabs[ilevel] += dx;
    m_yoffabs[ilevel] += dy;
    m_zoffabs[ilevel] += dz;
    
    if( ilevel < levelmax() )
      {
	m_pgrids[ilevel+1]->offset(0) -= dx;
	m_pgrids[ilevel+1]->offset(1) -= dy;
	m_pgrids[ilevel+1]->offset(2) -= dz;
      }

    //... enforce top mean density over same patch
    if( ilevel > levelmin() )
      {
	int ox = m_pgrids[ilevel]->offset(0);
	int oy = m_pgrids[ilevel]->offset(1);
	int oz = m_pgrids[ilevel]->offset(2);

	for( unsigned i=0; i<nx/2; ++i )
	  for( unsigned j=0; j<ny/2; ++j )
	    for( unsigned k=0; k<nz/2; ++k ){
	      coarsesum += (*m_pgrids[ilevel-1])(i+ox,j+oy,k+oz);
	      coarsecount++;
	    }

	coarsesum /= (double)coarsecount;
	finesum /= (double)finecount;

	for( unsigned i=0; i<nx; ++i )
	  for( unsigned j=0; j<ny; ++j )
	    for( unsigned k=0; k<nz; ++k )
	      (*mnew)(i,j,k) += (coarsesum-finesum);

	LOGINFO("level %d : corrected patch overlap mean density by %f",ilevel,coarsesum-finesum);
      }

    
    find_new_levelmin();
  }
  
	/*! determine level for which grid extends over entire domain */
	void find_new_levelmin( void )
	{
		for( unsigned i=0; i<=levelmax(); ++i )
		{
			unsigned n = 1<<i;
			if(	m_pgrids[i]->size(0) == n &&
			   	m_pgrids[i]->size(1) == n &&
			   	m_pgrids[i]->size(2) == n )
			{
				m_levelmin=i;
			}
		}
	}
	
	//! return maximum level in refinement hierarchy
	unsigned levelmax( void ) const
	{
		return m_pgrids.size()-1;
	}
	
	//! return minimum level in refinement hierarchy (the level which extends over the entire domain)
	unsigned levelmin( void ) const
	{
		return m_levelmin;
	}
	
};



//! class that computes the refinement structure given parameters
class refinement_hierarchy
{
	std::vector<double>
		x0_,	//!< x-coordinates of grid origins (in [0..1[)
		y0_,	//!< y-coordinates of grid origins (in [0..1[)
		z0_,	//!< z-coordinates of grid origins (in [0..1[)
		xl_,	//!< x-extent of grids (in [0..1[)
		yl_,	//!< y-extent of grids (in [0..1[)
		zl_;	//!< z-extent of grids (in [0..1[)
	
	std::vector<unsigned> 
		ox_,	//!< relative x-coordinates of grid origins (in coarser grid cells)
		oy_,	//!< relative y-coordinates of grid origins (in coarser grid cells)
		oz_,	//!< relative z-coordinates of grid origins (in coarser grid cells)
		oax_,	//!< absolute x-coordinates of grid origins (in fine grid cells)
		oay_,	//!< absolute y-coordinates of grid origins (in fine grid cells)
		oaz_,	//!< absolute z-coordinates of grid origins (in fine grid cells)
		nx_,	//!< x-extent of grids (in fine grid cells)
		ny_,	//!< y-extent of grids (in fine grid cells)
		nz_;	//!< z-extent of grids (in fine grid cells)
	
	unsigned 
		levelmin_,		//!< minimum grid level for Poisson solver
		levelmax_,		//!< maximum grid level for all operations
		levelmin_tf_,	//!< minimum grid level for density calculation
		padding_,		//!< padding in number of coarse cells between refinement levels
		blocking_factor_;
	
	
	config_file& cf_;	//!< reference to config_file
	
	bool align_top_,	//!< bool whether to align all grids with coarsest grid cells
	     equal_extent_; //!< bool whether the simulation code requires squared refinement regions (e.g. RAMSES)
	
	double 
		x0ref_[3],	//!< coordinates of refinement region origin (in [0..1[)
		lxref_[3];	//!< extent of refinement region (int [0..1[)

        size_t lnref_[3];
  bool   bhave_nref;
	
	int xshift_[3];	//!< shift of refinement region in coarse cells (in order to center it in the domain)
    double rshift_[3];
	
public:
	
	//! copy constructor
	refinement_hierarchy( const refinement_hierarchy& rh )
	: cf_( rh.cf_ )
	{
		*this = rh;
	}
	
	//! constructor from a config_file holding information about the desired refinement geometry
	explicit refinement_hierarchy( config_file& cf )
	: cf_( cf )
	{
	  //... query the parameter data we need
	  levelmin_	= cf_.getValue<unsigned>("setup","levelmin");
	  levelmax_	= cf_.getValue<unsigned>("setup","levelmax");
	  levelmin_tf_= cf_.getValueSafe<unsigned>("setup","levelmin_TF",levelmin_);
	  align_top_	= cf_.getValueSafe<bool>("setup","align_top",false);
	  equal_extent_ = cf_.getValueSafe<bool>("setup","force_equal_extent",false);
	  blocking_factor_= cf.getValueSafe<unsigned>( "setup", "blocking_factor",0);
		
	  bool bnoshift = cf_.getValueSafe<bool>("setup","no_shift",false);
	  bool force_shift = cf_.getValueSafe<bool>("setup","force_shift",false);
	  
        
	  //... call the region generator
	  if( levelmin_ != levelmax_ )
	    {
	      
	      double x1ref[3];
	      the_region_generator->get_AABB(x0ref_,x1ref,levelmax_);
	      for( int i=0; i<3; ++i )
		lxref_[i] = x1ref[i]-x0ref_[i];
	      bhave_nref = false;
	      
	      std::string region_type = cf.getValueSafe<std::string>("setup","region","box");
	      
	      LOGINFO("refinement region is \'%s\', w/ bounding box\n        left = [%f,%f,%f]\n       right = [%f,%f,%f]",
		      region_type.c_str(),x0ref_[0],x0ref_[1],x0ref_[2],x1ref[0],x1ref[1],x1ref[2]);
	      
	      
	      bhave_nref = the_region_generator->is_grid_dim_forced( lnref_ );
	    }
	  
	  
	  // if not doing any refinement levels, set extent to full box
	  if( levelmin_ == levelmax_ )
	    {
	      x0ref_[0] = 0.0;
	      x0ref_[1] = 0.0;
	      x0ref_[2] = 0.0;
	      
	      lxref_[0] = 1.0;
	      lxref_[1] = 1.0;
	      lxref_[2] = 1.0;
	    }
		
	  unsigned ncoarse = 1<<levelmin_;
	  
	  //... determine shift
	  
	  double xc[3];
	  xc[0] = fmod(x0ref_[0]+0.5*lxref_[0],1.0);
	  xc[1] = fmod(x0ref_[1]+0.5*lxref_[1],1.0);
	  xc[2] = fmod(x0ref_[2]+0.5*lxref_[2],1.0);
	  
	  if( (levelmin_ != levelmax_) && (!bnoshift || force_shift) )
	    {
	      xshift_[0] = (int)((0.5-xc[0])*ncoarse);
	      xshift_[1] = (int)((0.5-xc[1])*ncoarse);
	      xshift_[2] = (int)((0.5-xc[2])*ncoarse);
	    }else{
	    xshift_[0] = 0;
	    xshift_[1] = 0;
	    xshift_[2] = 0;
	  }
	  
	  char strtmp[32];
	  sprintf( strtmp, "%d", xshift_[0] );	cf_.insertValue( "setup", "shift_x", strtmp );
	  sprintf( strtmp, "%d", xshift_[1] );	cf_.insertValue( "setup", "shift_y", strtmp );
	  sprintf( strtmp, "%d", xshift_[2] );	cf_.insertValue( "setup", "shift_z", strtmp );
	  
	  rshift_[0] = -(double)xshift_[0]/ncoarse;
	  rshift_[1] = -(double)xshift_[1]/ncoarse;
	  rshift_[2] = -(double)xshift_[2]/ncoarse;
	  
	  x0ref_[0] += (double)xshift_[0]/ncoarse;
	  x0ref_[1] += (double)xshift_[1]/ncoarse;
	  x0ref_[2] += (double)xshift_[2]/ncoarse;
	  
	  
	  //... initialize arrays 
	  x0_.assign(levelmax_+1,0.0);	xl_.assign(levelmax_+1,1.0);
	  y0_.assign(levelmax_+1,0.0);	yl_.assign(levelmax_+1,1.0);
	  z0_.assign(levelmax_+1,0.0);	zl_.assign(levelmax_+1,1.0);
	  ox_.assign(levelmax_+1,0);	nx_.assign(levelmax_+1,0);
	  oy_.assign(levelmax_+1,0);	ny_.assign(levelmax_+1,0);
	  oz_.assign(levelmax_+1,0);	nz_.assign(levelmax_+1,0);
	  
	  oax_.assign(levelmax_+1,0);
	  oay_.assign(levelmax_+1,0);
	  oaz_.assign(levelmax_+1,0);
	  
		
	  nx_[levelmin_] = ncoarse;
	  ny_[levelmin_] = ncoarse;
	  nz_[levelmin_] = ncoarse;
	  
	  // set up base hierarchy sizes
	  for( unsigned ilevel=0; ilevel <=levelmin_; ++ilevel )
	    {
	      unsigned n = 1<<ilevel;
	      
	      xl_[ilevel] = yl_[ilevel] = zl_[ilevel] = 1.0;
	      nx_[ilevel] = ny_[ilevel] = nz_[ilevel] = n;
	    }
	  
	  // if no refinement, we can exit here
	  if( levelmax_ == levelmin_ )
            return;
	  
		
	  //... determine the position of the refinement region on the finest grid
	  int il,jl,kl,ir,jr,kr;
	  int nresmax = 1<<levelmax_;
	  
	  il = (int)(x0ref_[0] * nresmax);
	  jl = (int)(x0ref_[1] * nresmax);
	  kl = (int)(x0ref_[2] * nresmax);
	  ir = (int)((x0ref_[0]+lxref_[0]) * nresmax );//+ 1.0);
	  jr = (int)((x0ref_[1]+lxref_[1]) * nresmax );//+ 1.0);
	  kr = (int)((x0ref_[2]+lxref_[2]) * nresmax );//+ 1.0);
	  
	  //... align with coarser grids ...
	  if( align_top_ )
	    {
	      //... require alignment with top grid
	      unsigned nref = 1<<(levelmax_-levelmin_+1);
	      
	      if( bhave_nref )
		{
		  if( lnref_[0] % (1ul<<(levelmax_-levelmin_)) != 0 ||
		      lnref_[1] % (1ul<<(levelmax_-levelmin_)) != 0 ||
		      lnref_[2] % (1ul<<(levelmax_-levelmin_)) != 0 )
		    {
		      LOGERR("specified ref_dims and align_top=yes but cannot be aligned with coarse grid!");
		      throw std::runtime_error("specified ref_dims and align_top=yes but cannot be aligned with coarse grid!");
		    }
		}
	      
	      
	      il = (int)((double)il/nref)*nref;
	      jl = (int)((double)jl/nref)*nref;
	      kl = (int)((double)kl/nref)*nref;
	      
	      int irr = (int)((double)ir/nref)*nref;
	      int jrr = (int)((double)jr/nref)*nref;
	      int krr = (int)((double)kr/nref)*nref;
	      
	      if( irr < ir )
		ir = (int)((double)ir/nref + 1.0)*nref;
	      else
		ir = irr;
	      
	      if( jrr < jr )
		jr = (int)((double)jr/nref + 1.0)*nref;
	      else
		jr = jrr;
	      
	      if( krr < kr )
		kr = (int)((double)kr/nref + 1.0)*nref;
	      else
		kr = krr;
	      
	      
	    }else{
	    //... require alignment with coarser grid
	    il -= il%2; jl -= jl%2; kl -= kl%2;
	    ir += ir%2; jr += jr%2; kr += kr%2; 
	  }
	  
	  // if doing unigrid, set region to whole box
	  if( levelmin_ == levelmax_ )
	    {
	      il = jl = kl = 0;
	      ir = jr = kr = nresmax-1;
	    }
	  if( bhave_nref )
	    {
	      ir = il+lnref_[0];
	      jr = jl+lnref_[1];
	      kr = kl+lnref_[2];
	      
	    }
	  
	  //... make sure bounding box lies in domain
	  il = (il+nresmax)%nresmax; ir = (ir+nresmax)%nresmax;
	  jl = (jl+nresmax)%nresmax; jr = (jr+nresmax)%nresmax;
	  kl = (kl+nresmax)%nresmax; kr = (kr+nresmax)%nresmax;
	  
	  if( il>=ir || jl>=jr || kl>=kr )
	    {
	      LOGERR("Internal refinement bounding box error: [%d,%d]x[%d,%d]x[%d,%d]",il,ir,jl,jr,kl,kr);
	      throw std::runtime_error("refinement_hierarchy: Internal refinement bounding box error 1");
	    }
	  //... determine offsets
	  if( levelmin_ != levelmax_ )
	    {
	      
	      oax_[levelmax_] = (il+nresmax)%nresmax;
	      oay_[levelmax_] = (jl+nresmax)%nresmax;
	      oaz_[levelmax_] = (kl+nresmax)%nresmax;
	      nx_[levelmax_]  = ir-il;	
	      ny_[levelmax_]  = jr-jl;	
	      nz_[levelmax_]  = kr-kl;
	      
	      if( equal_extent_ )
		{
		  
		  if( bhave_nref && (lnref_[0]!=lnref_[1]||lnref_[0]!=lnref_[2]) )
		    {
		      LOGERR("Specified equal_extent=yes conflicting with ref_dims which are not equal.");
		      throw std::runtime_error("Specified equal_extent=yes conflicting with ref_dims which are not equal.");
		    }
		  size_t ilevel = levelmax_;
		  size_t nmax = std::max( nx_[ilevel], std::max( ny_[ilevel], nz_[ilevel] ) );				
		  int dx = (int)((double)(nmax-nx_[ilevel])*0.5);
		  int dy = (int)((double)(nmax-ny_[ilevel])*0.5);
		  int dz = (int)((double)(nmax-nz_[ilevel])*0.5);
		  
		  oax_[ilevel] -= dx;
		  oay_[ilevel] -= dy;
		  oaz_[ilevel] -= dz;
		  nx_[ilevel] = nmax;
		  ny_[ilevel] = nmax;
		  nz_[ilevel] = nmax;
		  
		  il = oax_[ilevel];
		  jl = oay_[ilevel];
		  kl = oaz_[ilevel];
		  ir = il + nmax;
		  jr = jl + nmax;
		  kr = kl + nmax;
		}
	    }
	  
	  padding_	= cf_.getValueSafe<unsigned>("setup","padding", 8);
	  
	  //... determine position of coarser grids
	  for( unsigned ilevel=levelmax_-1; ilevel> levelmin_; --ilevel )
	    {
	      il = (int)((double)il * 0.5 - padding_);
	      jl = (int)((double)jl * 0.5 - padding_);
	      kl = (int)((double)kl * 0.5 - padding_);
	      
	      ir = (int)((double)ir * 0.5 + padding_);
	      jr = (int)((double)jr * 0.5 + padding_);
	      kr = (int)((double)kr * 0.5 + padding_);
	      
	      //... align with coarser grids ...
	      if( align_top_ )
		{
		  //... require alignment with top grid
		  unsigned nref = 1<<(ilevel-levelmin_);
		  
		  il = (int)((double)il/nref)*nref;
		  jl = (int)((double)jl/nref)*nref;
		  kl = (int)((double)kl/nref)*nref;
		  
		  ir = (int)((double)ir/nref+1.0)*nref;
		  jr = (int)((double)jr/nref+1.0)*nref;
		  kr = (int)((double)kr/nref+1.0)*nref;
		  
		}
	      else
		{
		  //... require alignment with coarser grid
		  il -= il%2; jl -= jl%2; kl -= kl%2;
		  ir += ir%2; jr += jr%2; kr += kr%2; 
		}

			
	      if( il>=ir || jl>=jr || kl>=kr || il < 0 || jl < 0 || kl < 0)
		{
		  LOGERR("Internal refinement bounding box error: [%d,%d]x[%d,%d]x[%d,%d], level=%d",il,ir,jl,jr,kl,kr,ilevel);
		  throw std::runtime_error("refinement_hierarchy: Internal refinement bounding box error 2");
		}
	      oax_[ilevel] = il;
	      oay_[ilevel] = jl;
	      oaz_[ilevel] = kl;
	      nx_[ilevel]  = ir-il;
	      ny_[ilevel]  = jr-jl;
	      nz_[ilevel]  = kr-kl;
	      
	      if (blocking_factor_)
		{
		  nx_[ilevel] += nx_[ilevel] % blocking_factor_;
		  ny_[ilevel] += ny_[ilevel] % blocking_factor_;
		  nz_[ilevel] += nz_[ilevel] % blocking_factor_;
		}
	      
	      if( equal_extent_ )
		{
		  size_t nmax = std::max( nx_[ilevel], std::max( ny_[ilevel], nz_[ilevel] ) );				
		  int dx = (int)((double)(nmax-nx_[ilevel])*0.5);
		  int dy = (int)((double)(nmax-ny_[ilevel])*0.5);
		  int dz = (int)((double)(nmax-nz_[ilevel])*0.5);
				
		  oax_[ilevel] -= dx;
		  oay_[ilevel] -= dy;
		  oaz_[ilevel] -= dz;
		  nx_[ilevel] = nmax;
		  ny_[ilevel] = nmax;
		  nz_[ilevel] = nmax;
		  
		  il = oax_[ilevel];
		  jl = oay_[ilevel];
		  kl = oaz_[ilevel];
		  ir = il + nmax;
		  jr = jl + nmax;
		  kr = kl + nmax;
		}
	      
	    }
	  
	  //... determine relative offsets between grids
	  for( unsigned ilevel=levelmax_; ilevel>levelmin_; --ilevel )
	    {
	      ox_[ilevel] = (oax_[ilevel]/2 - oax_[ilevel-1]);
	      oy_[ilevel] = (oay_[ilevel]/2 - oay_[ilevel-1]);
	      oz_[ilevel] = (oaz_[ilevel]/2 - oaz_[ilevel-1]);
	    }
	  
	  //... do a forward sweep to ensure that absolute offsets are also correct now
	  for( unsigned ilevel=levelmin_+1; ilevel<=levelmax_; ++ilevel )
	    {
	      oax_[ilevel] = 2*oax_[ilevel-1]+2*ox_[ilevel];
	      oay_[ilevel] = 2*oay_[ilevel-1]+2*oy_[ilevel];
	      oaz_[ilevel] = 2*oaz_[ilevel-1]+2*oz_[ilevel];
	    }
	  
	  for( unsigned ilevel=levelmin_+1; ilevel<=levelmax_; ++ilevel )
	    {
	      double h = 1.0/(1ul<<ilevel);
	      
	      x0_[ilevel] = h*(double)oax_[ilevel];
	      y0_[ilevel] = h*(double)oay_[ilevel];
	      z0_[ilevel] = h*(double)oaz_[ilevel];
	      
	      xl_[ilevel] = h*(double)nx_[ilevel];
	      yl_[ilevel] = h*(double)ny_[ilevel];
	      zl_[ilevel] = h*(double)nz_[ilevel];
	    }
	  
	  
	  
	  // do a consistency check that largest subgrid in zoom is not larger than half the box size
	  for( unsigned ilevel=levelmin_+1; ilevel<=levelmax_; ++ilevel )
	    {	
	      if( nx_[ilevel] > (1ul<<(ilevel-1)) ||
		  ny_[ilevel] > (1ul<<(ilevel-1)) ||
		  nz_[ilevel] > (1ul<<(ilevel-1)) )
		{
		  LOGERR("On level %d, subgrid is larger than half the box. This is not allowed!",ilevel);
		  throw std::runtime_error("Fatal: Subgrid larger than half boxin zoom.");
		}
	    }
	  
	  // update the region generator with what has been actually created
	  double left[3] =  { x0_[levelmax_]+rshift_[0], y0_[levelmax_]+rshift_[1], z0_[levelmax_]+rshift_[2] };
	  double right[3] = { left[0]+xl_[levelmax_], left[1]+yl_[levelmax_], left[2]+zl_[levelmax_] };
	  the_region_generator->update_AABB( left, right );
	}
  
  //! asignment operator
  refinement_hierarchy& operator=( const refinement_hierarchy& o )
  {
    levelmin_ = o.levelmin_;
    levelmax_ = o.levelmax_;
    padding_  = o.padding_;
    cf_ = o.cf_;
    align_top_ = o.align_top_;
    for( int i=0; i<3; ++i )
      {
	x0ref_[i] = o.x0ref_[i];
	lxref_[i] = o.lxref_[i];
	xshift_[i] = o.xshift_[i];
	rshift_[i] = o.rshift_[i];
      }
    
    x0_ = o.x0_; y0_ = o.y0_; z0_ = o.z0_;
    xl_ = o.xl_; yl_ = o.yl_; zl_ = o.zl_;
    ox_ = o.ox_; oy_ = o.oy_; oz_ = o.oz_;
    oax_= o.oax_; oay_ = o.oay_; oaz_ = o.oaz_;
    nx_ = o.nx_; ny_=o.ny_; nz_=o.nz_;
    
    return *this;
  }
	
	/*! cut a grid level to the specified size
	 * @param ilevel grid level on which to perform the size adjustment
	 * @param nx new x-extent in fine grid cells
	 * @param ny new y-extent in fine grid cells
	 * @param nz new z-extent in fine grid cells
	 * @param oax new x-offset in units fine grid units
	 * @param oay new y-offset in units fine grid units
	 * @param oaz new z-offset in units fine grid units
	
	 */	
	void adjust_level( unsigned ilevel, int nx, int ny, int nz, int oax, int oay, int oaz )
	{
		double h = 1.0/(1<<ilevel);
		
		int dx,dy,dz;
		
		dx = oax_[ilevel] - oax;
		dy = oay_[ilevel] - oay;
		dz = oaz_[ilevel] - oaz;
		
		ox_[ilevel] -= dx/2;
		oy_[ilevel] -= dy/2;
		oz_[ilevel] -= dz/2;
		
		oax_[ilevel] = oax; 
		oay_[ilevel] = oay; 
		oaz_[ilevel] = oaz;
		
		nx_[ilevel] = nx;
		ny_[ilevel] = ny;
		nz_[ilevel] = nz;
		
		x0_[ilevel] = h*oax;
		y0_[ilevel] = h*oay;
		z0_[ilevel] = h*oaz;
		
		xl_[ilevel] = h*nx;
		yl_[ilevel] = h*ny;
		zl_[ilevel] = h*nz;
		
		if( ilevel < levelmax_ )
		{
			ox_[ilevel+1] += dx;
			oy_[ilevel+1] += dy;
			oz_[ilevel+1] += dz;
		}
		
		find_new_levelmin();
		
	}
	
	/*! determine level for which grid extends over entire domain */
	void find_new_levelmin( bool print=false )
	{
		unsigned old_levelmin( levelmin_ );
		
		for( unsigned i=0; i<=levelmax(); ++i )
		{
			unsigned n = 1<<i;
			
			if(	oax_[i]==0 && oay_[i]==0 && oaz_[i]==0
			   && nx_[i]==n && ny_[i]==n && nz_[i]==n )
			{
				levelmin_=i;
			}
		}
		
		if( (old_levelmin != levelmin_) && print)
			LOGINFO("refinement_hierarchy: set new levelmin to %d", levelmin_ );
	}
	
	//! get absolute grid offset for a specified level along a specified dimension (in fine grid units)
	unsigned offset_abs( unsigned ilevel, int dim ) const
	{
		if( dim==0 )
			return oax_.at(ilevel);
		if( dim==1 )
			return oay_.at(ilevel);
		return oaz_.at(ilevel);
	}
	
	//! get relative grid offset for a specified level along a specified dimension (in coarser grid units)
	int offset( unsigned ilevel, int dim ) const
	{
		if( dim==0 )
			return ox_.at(ilevel);
		if( dim==1 )
			return oy_.at(ilevel);
		return oz_.at(ilevel);
	}
	
	//! get grid size for a specified level along a specified dimension
	size_t size( unsigned ilevel, int dim ) const
	{
		if( dim==0 )
			return nx_.at(ilevel);
		if( dim==1 )
			return ny_.at(ilevel);
		return nz_.at(ilevel);
	}
	
	//! get minimum grid level (the level for which the grid covers the entire domain)
	unsigned levelmin( void ) const
	{	return levelmin_;	}
	
	//! get maximum grid level
	unsigned levelmax( void ) const
	{	return levelmax_;	}
	
	//! get the total shift of the coordinate system in units of coarse cells
	int get_shift( int idim ) const
	{	return xshift_[idim];  }
	
    //! get the total shift of the coordinate system in box coordinates
    const double* get_coord_shift( void ) const
    {   return rshift_; }
    
	//! write refinement hierarchy to stdout
	void output( void ) const
	{
		std::cout << "-------------------------------------------------------------\n";
		
		if( xshift_[0]!=0||xshift_[1]!=0||xshift_[2]!=0 )
			std::cout << " - Domain will be shifted by (" << xshift_[0] << ", " << xshift_[1] << ", " << xshift_[2] << ")\n" << std::endl;

		std::cout << " - Grid structure:\n";
		
		for( unsigned ilevel=levelmin_; ilevel<=levelmax_; ++ilevel )
		{
			std::cout 
			<< "     Level " << std::setw(3) << ilevel << " :   offset = (" << std::setw(5) << ox_[ilevel] << ", " << std::setw(5) << oy_[ilevel] << ", " << std::setw(5) << oz_[ilevel] << ")\n"
			<< "               offset_abs = (" << std::setw(5) << oax_[ilevel] << ", " << std::setw(5) << oay_[ilevel] << ", " << std::setw(5) << oaz_[ilevel] << ")\n"
			<< "                   size   = (" << std::setw(5) << nx_[ilevel] << ", " << std::setw(5) << ny_[ilevel] << ", " << std::setw(5) << nz_[ilevel] << ")\n";
		}
		std::cout << "-------------------------------------------------------------\n";
	}
	
	void output_log( void ) const
	{
		LOGUSER("   Domain shifted by      (%5d,%5d,%5d)",xshift_[0],xshift_[1],xshift_[2]);		
		for( unsigned ilevel=levelmin_; ilevel<=levelmax_; ++ilevel )
		{	
			LOGUSER("   Level %3d :   offset = (%5d,%5d,%5d)",ilevel,ox_[ilevel],oy_[ilevel],oz_[ilevel]);
			LOGUSER("                   size = (%5d,%5d,%5d)",nx_[ilevel],ny_[ilevel],nz_[ilevel]);
		}
	}
    
    
};


typedef GridHierarchy<real_t> grid_hierarchy;
typedef MeshvarBnd<real_t> meshvar_bnd;
typedef Meshvar<real_t> meshvar;


#endif

