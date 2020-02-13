/*
 *  schemes.hh
 *  GravitySolver
 *
 *  Created by Oliver Hahn on 2/1/10.
 *  Copyright 2010 KIPAC/Stanford University. All rights reserved.
 *
 */

#ifndef __SCHEME_HH
#define __SCHEME_HH

#include <vector>
#include <stdexcept>

#include "solver.hh"

//... abstract implementation of the Poisson/Force scheme
template< class L, class G, typename real_t=double >
class scheme
{
public:
	typedef L laplacian;
	typedef G gradient;
	
	laplacian m_laplacian;
	gradient m_gradient;
	
	template< class C >
	inline real_t grad_x( const C&c, const int i, const int j, const int k )
	{ return m_gradient.apply_x( c,i,j,k ); }
	
	template< class C >
	inline real_t grad_y( const C&c, const int i, const int j, const int k )
	{ return m_gradient.apply_y( c,i,j,k ); }
	
	template< class C >
	inline real_t grad_z( const C&c, const int i, const int j, const int k )
	{ return m_gradient.apply_z( c,i,j,k ); }
	
	template< class C >
	inline real_t L_apply( const C&c, const int i, const int j, const int k ) 
	{ return m_laplacian.apply( c,i,j,k ); }
	
	template< class C >
	inline real_t L_rhs( const C&c, const int i, const int j, const int k ) 
	{ return m_laplacian.rhs( c,i,j,k ); }
	
	inline real_t ccoeff( void )
	{ return m_laplacian.ccoeff(); }
	
};


template< int nextent, typename T >
class gradient
{
	typedef T real_t;
	std::vector<real_t> m_stencil;
	const unsigned nl;
public:
	
	gradient()
	: nl( 2*nextent+1 )
	{ 
		m_stencil.assign(nl*nl*nl,(real_t)0.0);
	}
	
	real_t& operator()(int i)
	{ return m_stencil[i+nextent]; }
	
	const real_t& operator()(int i) const
	{ return m_stencil[i+nextent]; }
	
	template< class C >
	inline void apply( const C& c, C& f, int dir )
	{
		f = c;
		
		int nx=c.size(0), ny=c.size(1), nz=c.size(2);		
		double hx = 1.0/(nx+1.0), hy = 1.0/(ny+1.0), hz = 1.0/(nz+1.0);
		
		f.zero();
		
		if( dir == 0 )
			for( int i=0; i<nx; ++i )
				for( int j=0; j<ny; ++j )
					for( int k=0; k<nz; ++k )
						for( int ii = -nextent; ii<=nextent; ++ii )
							f(i,j,k) += (*this)(ii) * c(i+ii,j,k)/hx;
		else if( dir == 1 )
			for( int i=0; i<nx; ++i )
				for( int j=0; j<ny; ++j )
					for( int k=0; k<nz; ++k )
						for( int jj = -nextent; jj<=nextent; ++jj )
							f(i,j,k) += (*this)(jj) * c(i,j+jj,k)/hy;
		else if( dir == 2 )
			for( int i=0; i<nx; ++i )
				for( int j=0; j<ny; ++j )
					for( int k=0; k<nz; ++k )
						for( int kk = -nextent; kk<=nextent; ++kk )
							f(i,j,k) += (*this)(kk) * c(i,j,k+kk)/hz;
		
	}
};

template< int nextent, typename real_t >
class base_stencil
{
protected:
	std::vector<real_t> m_stencil;
	const unsigned nl;
public:
	bool m_modsource;
	
public:
	base_stencil( bool amodsource = false )
	: nl( 2*nextent+1 ), m_modsource( amodsource )
	{
		m_stencil.assign(nl*nl*nl,(real_t)0.0);
	}
	
	real_t& operator()(int i, int j, int k)
	{ return m_stencil[((i+nextent)*nl+(j+nextent))*nl+(k+nextent)]; }
	
	const real_t& operator()(unsigned i, unsigned j, unsigned k) const
	{ return m_stencil[((i+nextent)*nl+(j+nextent))*nl+(k+nextent)]; }
	
	template< class C >
	inline real_t rhs( const C& c, const int i, const int j, const int k )
	{
		real_t sum = this->apply( c, i, j, k );
		sum -= (*this)(0,0,0) * c(i,j,k);
		return sum;
	}
	
	inline real_t ccoeff( void )
	{
		return (*this)(0,0,0);
	}
	
	
	template< class C >
	inline real_t apply( const C& c, const int i, const int j, const int k )
	{
		real_t sum = 0.0;
		
		for( int ii=-nextent; ii<=nextent; ++ii )
			for( int jj=-nextent; jj<=nextent; ++jj )
				for( int kk=-nextent; kk<=nextent; ++kk )
					sum += (*this)(ii,jj,kk) * c(i+ii,j+jj,k+kk);
		
		return sum;
	}
	
	template< class C >
	inline real_t modsource( const C& c, const int i, const int j, const int k )
	{
		return 0.0;
	}
	
};


/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/


//... Implementation of the Gradient schemes............................................


template< typename real_t >
class deriv_2P : public gradient<1,real_t>
{
	
public:
	deriv_2P( void )
	{
		(*this)( 0 ) =  0.0;
		(*this)(-1 ) = -0.5;
		(*this)(+1 ) = +0.5;		
	}
	
	
};

//... Implementation of the Laplacian schemes..........................................


template< typename real_t >
class stencil_7P : public base_stencil<1,real_t>
{
	
public:
	stencil_7P( void )
	{
		(*this)( 0, 0, 0) = -6.0;
		(*this)(-1, 0, 0) = +1.0;
		(*this)(+1, 0, 0) = +1.0;
		(*this)( 0,-1, 0) = +1.0;
		(*this)( 0,+1, 0) = +1.0;
		(*this)( 0, 0,-1) = +1.0;
		(*this)( 0, 0,+1) = +1.0;
	}
	
	template< class C >
	inline real_t apply( const C& c, const int i, const int j, const int k ) const
	{
		return c(i-1,j,k)+c(i+1,j,k)+c(i,j-1,k)+c(i,j+1,k)+c(i,j,k-1)+c(i,j,k+1)-6.0*c(i,j,k);
	}
	
	template< class C >
	inline real_t rhs( const C& c, const int i, const int j, const int k ) const
	{
		return c(i-1,j,k)+c(i+1,j,k)+c(i,j-1,k)+c(i,j+1,k)+c(i,j,k-1)+c(i,j,k+1);
	}
	
	inline real_t ccoeff( void )
	{
		return -6.0;
	}
};


template< typename real_t >
class stencil_13P : public base_stencil<2,real_t>
{
	
public:
	stencil_13P( void )
	{
		(*this)( 0, 0, 0) = -90.0/12.;
		
		(*this)(-1, 0, 0) = 
		(*this)(+1, 0, 0) = 
		(*this)( 0,-1, 0) = 
		(*this)( 0,+1, 0) = 
		(*this)( 0, 0,-1) = 
		(*this)( 0, 0,+1) = 16./12.;
		
		(*this)(-2, 0, 0) = 
		(*this)(+2, 0, 0) = 
		(*this)( 0,-2, 0) = 
		(*this)( 0,+2, 0) = 
		(*this)( 0, 0,-2) = 
		(*this)( 0, 0,+2) = -1./12.;
	}
	
	template< class C >
	inline real_t apply( const C& c, const int i, const int j, const int k )
	{
		return 
			(-1.0*(c(i-2,j,k)+c(i+2,j,k)+c(i,j-2,k)+c(i,j+2,k)+c(i,j,k-2)+c(i,j,k+2))
			 +16.0*(c(i-1,j,k)+c(i+1,j,k)+c(i,j-1,k)+c(i,j+1,k)+c(i,j,k-1)+c(i,j,k+1))
			 -90.0*c(i,j,k))/12.0;
	}
	
	template< class C >
	inline real_t rhs( const C& c, const int i, const int j, const int k )
	{
		return 
			(-1.0*(c(i-2,j,k)+c(i+2,j,k)+c(i,j-2,k)+c(i,j+2,k)+c(i,j,k-2)+c(i,j,k+2))
			 +16.0*(c(i-1,j,k)+c(i+1,j,k)+c(i,j-1,k)+c(i,j+1,k)+c(i,j,k-1)+c(i,j,k+1)))/12.0;
	}
	
	inline real_t ccoeff( void )
	{
		return -90.0/12.0;
	}
};

#endif


