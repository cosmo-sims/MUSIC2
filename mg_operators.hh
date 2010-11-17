/*
 
 mg_operators.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
*/

#ifndef __MG_OPERATORS_HH
#define __MG_OPERATORS_HH

//! injection operator based on Catmull-Rom splines with arbitrary refinement factor
class mg_cubic_mult
{
protected:
	template< typename real_t >
	inline double CUBE( real_t x ) const
	{ return x*x*x; }
	
	template< typename T >
	inline T SQR( T x ) const
	{ return x*x; }
	
public:
	template<typename M >
	inline double interp_cubic( double dx, double dy, double dz, M& V, int ox, int oy, int oz ) const
	{
		register int	i, j, k;
		double           u[4], v[4], w[4];
		double           r[4], q[4];
		double           vox = 0;
		
		int sx=0,sy=0,sz=0;
		if( dx > 0.0 ) sx=1;
		if( dy > 0.0 ) sy=1;
		if( dz > 0.0 ) sz=1;
		
		
		if( dx > 0 )
		{	
			w[0] = -0.5*CUBE(dx)+SQR(dx)-0.5*dx;
			w[1] = 1.5*CUBE(dx)-2.5*SQR(dx)+1.0;
			w[2] = -1.5*CUBE(dx)+2.0*SQR(dx)+0.5*dx;
			w[3] = 0.5*CUBE(dx)-0.5*SQR(dx);
		} else {
			w[0] = -0.5*CUBE(dx)-0.5*SQR(dx);
			w[1] = 1.5*CUBE(dx)+2.0*SQR(dx)-0.5*dx;
			w[2] = -1.5*CUBE(dx)-2.5*SQR(dx)+1.0;
			w[3] = 0.5*CUBE(dx)+SQR(dx)+0.5*dx;
		}
		
		if( dy > 0 )
		{
			v[0] = -0.5*CUBE(dy)+SQR(dy)-0.5*dy;
			v[1] = 1.5*CUBE(dy)-2.5*SQR(dy)+1.0;
			v[2] = -1.5*CUBE(dy)+2.0*SQR(dy)+0.5*dy;
			v[3] = 0.5*CUBE(dy)-0.5*SQR(dy);
		} else {
			v[0] = -0.5*CUBE(dy)-0.5*SQR(dy);
			v[1] = 1.5*CUBE(dy)+2.0*SQR(dy)-0.5*dy;
			v[2] = -1.5*CUBE(dy)-2.5*SQR(dy)+1.0;
			v[3] = 0.5*CUBE(dy)+SQR(dy)+0.5*dy;
		}
		
		if( dz > 0 )
		{
			u[0] = -0.5*CUBE(dz)+SQR(dz)-0.5*dz;
			u[1] = 1.5*CUBE(dz)-2.5*SQR(dz)+1.0;
			u[2] = -1.5*CUBE(dz)+2.0*SQR(dz)+0.5*dz;
			u[3] = 0.5*CUBE(dz)-0.5*SQR(dz);
		} else {
			u[0] = -0.5*CUBE(dz)-0.5*SQR(dz);
			u[1] = 1.5*CUBE(dz)+2.0*SQR(dz)-0.5*dz;
			u[2] = -1.5*CUBE(dz)-2.5*SQR(dz)+1.0;
			u[3] = 0.5*CUBE(dz)+SQR(dz)+0.5*dz;
		}
		
		for (k = 0; k < 4; k++)
		{
			q[k] = 0;
			for (j = 0; j < 4; j++)
			{
				r[j] = 0;
				for (i = 0; i < 4; i++)
				{
					r[j] += u[i] * V(ox+k-2+sx,oy+j-2+sy,oz+i-2+sz);
				}
				q[k] += v[j] * r[j];
			}
			vox += w[k] * q[k];
		}
		
		return vox;
		
	}
	
	template< typename m1, typename m2 >
	inline void prolong( m1& V, m2& v ) const
	{
		//int N = 1<<R;
		int Nl= 0;//, Nr=N/2;
		double dx = 0.5;//1.0/(double)N;
		int nx = v.size(0), ny = v.size(1), nz = v.size(2);
		double finemean = 0.0, coarsemean = 0.0;
		unsigned finecount = 0, coarsecount = 0;
		
		int oxc = V.offset(0), oyc = V.offset(1), ozc = V.offset(2);
		int oxf = v.offset(0), oyf = v.offset(1), ozf = v.offset(2);
		
		
		
		#pragma omp parallel for reduction(+:finemean,finecount,coarsemean,coarsecount)
		for( int i=0; i<nx; i+=2 )
			for( int j=0; j<ny; j+=2 )
				for( int k=0; k<nz; k+=2 )
				{
					//determine the coarse cell index
					
					int iic,jjc,kkc;
					
					iic = oxf+(int)(dx*(double)i);
					jjc = oyf+(int)(dx*(double)j);
					kkc = ozf+(int)(dx*(double)k);
					
					double xc,yc,zc,xf,yf,zf;
					
					
					
					double coarseval = V(iic,jjc,kkc);
					double finesum   = 0.0;
					for( int ii=0; ii<2; ++ii )
						for( int jj=0; jj<2; ++jj )
							for( int kk=0; kk<2; ++kk )
							{
								xf = oxf+(int)(dx*(double)(i+ii))-0.5*dx;
								yf = oyf+(int)(dx*(double)(j+jj))-0.5*dx;
								zf = ozf+(int)(dx*(double)(k+kk))-0.5*dx;
								
								xc = (double)(iic);
								yc = (double)(jjc);
								zc = (double)(kkc);
								
								double val = interp_cubic(xf-xc,yf-yc,zf-zc,V,iic,jjc,kkc);
								
								finesum += val;
								v(i+ii,j+jj,k+kk) = val;
							}
					
					finesum /= 8.0;
					
					/*double dv = coarseval-finesum;
					
					for( int ii=0; ii<2; ++ii )
						for( int jj=0; jj<2; ++jj )
							for( int kk=0; kk<2; ++kk )
								v(i+ii,j+jj,k+kk) += dv;
					 */
				}
		
		
	}
	
	template< typename m1, typename m2 >
	inline void prolong( m1& V, m2& v, int oxc, int oyc, int ozc, int oxf, int oyf, int ozf, int R ) const
	{
		int N = 1<<R;
		int Nl= -N/2+1;//, Nr=N/2;
		double dx = 1.0/(double)N;
		int nx = v.size(0), ny = v.size(1), nz = v.size(2);
		double finemean = 0.0, coarsemean = 0.0;
		unsigned finecount = 0, coarsecount = 0;
		
		#pragma omp parallel for reduction(+:finemean,finecount,coarsemean,coarsecount)
		for( int i=0; i<nx; ++i )
			for( int j=0; j<ny; ++j )
				for( int k=0; k<nz; ++k )
				{
					//determine the coarse cell index
					int iic,jjc,kkc;
					
					iic = (int)(dx*(double)(oxf+i)-oxc);
					jjc = (int)(dx*(double)(oyf+j)-oyc);
					kkc = (int)(dx*(double)(ozf+k)-ozc);
					
					double xc,yc,zc,xf,yf,zf;
					
					xc = iic+0.5;
					yc = jjc+0.5;
					zc = kkc+0.5;
					
					xf = (dx*(double)(oxf+i+Nl)-oxc)-0.5*dx + 0.5;
					yf = (dx*(double)(oyf+j+Nl)-oyc)-0.5*dx + 0.5;
					zf = (dx*(double)(ozf+k+Nl)-ozc)-0.5*dx + 0.5;
					
					
					
					double val = interp_cubic(xf-xc,yf-yc,zf-zc,V,iic,jjc,kkc);
					v(i,j,k) = val;
					
					finemean += val; finecount++;
					coarsemean += V(iic,jjc,kkc); coarsecount++;
				}
		
		
		coarsemean /= coarsecount;
		finemean /= finecount;
		
		
		//... subtract the mean difference caused by interpolation
		/*double dmean = coarsemean-finemean;
		
		#pragma omp parallel for
		for( int i=0; i<nx; ++i )
			for( int j=0; j<ny; ++j )
				for( int k=0; k<nz; ++k )
					v(i,j,k) += dmean;*/
	}
	
	template< typename m1, typename m2 >
	inline void prolong_add( m1& V, m2& v, int oxc, int oyc, int ozc, int oxf, int oyf, int ozf, int R ) const
	{
		int N = 1<<R;
		int Nl= -N/2+1;//, Nr=N/2;
		double dx = 1.0/(double)N;
		
		int nx = v.size(0), ny = v.size(1), nz = v.size(2);
		
		double finemean = 0.0, coarsemean = 0.0;
		unsigned finecount = 0, coarsecount = 0;
		
		#pragma omp parallel for reduction(+:finemean,finecount,coarsemean,coarsecount)
		for( int i=0; i<nx; ++i )
			for( int j=0; j<ny; ++j )
				for( int k=0; k<nz; ++k )
				{
					//determine the coarse cell index
					int iic,jjc,kkc;
					
					iic = (int)(dx*(double)(oxf+i)-oxc);
					jjc = (int)(dx*(double)(oyf+j)-oyc);
					kkc = (int)(dx*(double)(ozf+k)-ozc);
					
					double xc,yc,zc,xf,yf,zf;
					
					xc = iic+0.5;
					yc = jjc+0.5;
					zc = kkc+0.5;
					
					xf = (dx*(double)(oxf+i+Nl)-oxc)-0.5*dx + 0.5;
					yf = (dx*(double)(oyf+j+Nl)-oyc)-0.5*dx + 0.5;
					zf = (dx*(double)(ozf+k+Nl)-ozc)-0.5*dx + 0.5;
					
					double val = interp_cubic(xf-xc,yf-yc,zf-zc,V,iic,jjc,kkc); 
					
					v(i,j,k) += val;
					
					finemean += val; finecount++;
					coarsemean += V(iic,jjc,kkc); coarsecount++;
				}
		
		coarsemean /= coarsecount;
		finemean /= finecount;
		
		/*double dmean = coarsemean-finemean;
		
		//... subtract the mean difference caused by interpolation
		#pragma omp parallel for
		for( int i=0; i<nx; ++i )
			for( int j=0; j<ny; ++j )
				for( int k=0; k<nz; ++k )
					v(i,j,k) += dmean;
		*/	
	}
};


//! injection operator based on Catmull-Rom splines
class mg_cubic
{
protected:
	template< typename real_t >
	inline double CUBE( real_t x ) const
	{ return x*x*x; }
	
	template< typename T >
	inline T SQR( T x ) const
	{ return x*x; }

public:
	template< int sx, int sy, int sz, typename M >
	inline double interp_cubic( int x, int y, int z, M& V, double s=1.0 ) const
	{
		int				i, j, k;
		//double           dx, dy, dz;
		double           u[4], v[4], w[4];
		double           r[4], q[4];
		double           vox = 0;
		
		//dx = 0.5*((double)sz - 0.5)*s;
		//dy = 0.5*((double)sy - 0.5)*s;
		//dz = 0.5*((double)sx - 0.5)*s;
				
		if( sz == 1 )
		{	
			u[0] = -4.5/64.0;
			u[1] = 55.5/64.0;
			u[2] = 14.5/64.0;
			u[3] = -1.5/64.0;
		} else {
			u[0] = -1.5/64.0;
			u[1] = 14.5/64.0;
			u[2] = 55.5/64.0;
			u[3] = -4.5/64.0;
		}
		
		if( sy == 1 )
		{
			v[0] = -4.5/64.0;
			v[1] = 55.5/64.0;
			v[2] = 14.5/64.0;
			v[3] = -1.5/64.0;
		} else{
			v[0] = -1.5/64.0;
			v[1] = 14.5/64.0;
			v[2] = 55.5/64.0;
			v[3] = -4.5/64.0;
		}
		
		if( sx == 1 )
		{
			w[0] = -4.5/64.0;
			w[1] = 55.5/64.0;
			w[2] = 14.5/64.0;
			w[3] = -1.5/64.0;
		} else {
			w[0] = -1.5/64.0;
			w[1] = 14.5/64.0;
			w[2] = 55.5/64.0;
			w[3] = -4.5/64.0;
		}


		for (k = 0; k < 4; k++)
		{
			q[k] = 0;
			for (j = 0; j < 4; j++)
			{
				r[j] = 0;
				for (i = 0; i < 4; i++)
				{
					r[j] += u[i] * V(x+k-2+sx,y+j-2+sy,z+i-2+sz);
				}
				q[k] += v[j] * r[j];
			}
			vox += w[k] * q[k];
		}
		
		
		return vox;
		
	}
	
public:
	
	
	//... restricts v to V
	template< typename m1, typename m2 >
	inline void restrict( m1& v, m2& V )
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
		
		#pragma omp parallel for
		for( int i=0; i<nx; ++i )
		{
			int i2 = 2*i;
			for( int j=0,j2=0; j<ny; ++j,j2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 )
					V(ox+i,oy+j,oz+k) = 0.125 * ( interp_cubic<1,0,0>(i2,j2,k2,v)
												 +interp_cubic<0,1,0>(i2,j2,k2,v)
												 +interp_cubic<0,0,1>(i2,j2,k2,v)
												 +interp_cubic<1,1,0>(i2,j2,k2,v)
												 +interp_cubic<0,1,1>(i2,j2,k2,v)
												 +interp_cubic<1,0,1>(i2,j2,k2,v)
												 +interp_cubic<1,1,1>(i2,j2,k2,v)
												 +interp_cubic<0,0,0>(i2,j2,k2,v));
		} 
												 		
		
	}
	
	//... restricts v to V
	template< typename m1, typename m2 >
	inline void restrict_add( m1& v, m2& V )
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
		#pragma omp parallel for
		for( int i=0; i<nx; ++i )
		{
			int i2 = 2*i;
			for( int j=0,j2=0; j<ny; ++j,j2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 )
					V(ox+i,oy+j,oz+k) += 0.125 * ( interp_cubic<1,0,0>(i2+1,j2,k2,v)
												 +interp_cubic<0,1,0>(i2,j2+1,k2,v)
												 +interp_cubic<0,0,1>(i2,j2,k2+1,v)
												 +interp_cubic<1,1,0>(i2+1,j2+1,k2,v)
												 +interp_cubic<0,1,1>(i2,j2+1,k2+1,v)
												 +interp_cubic<1,0,1>(i2+1,j2,k2+1,v)
												 +interp_cubic<1,1,1>(i2+1,j2+1,k2+1,v)
												 +interp_cubic<0,0,0>(i2,j2,k2,v));
		} 
		
		
	}
	
	
	//.. straight restriction on boundary
	template< typename m1, typename m2 >
	inline void restrict_bnd( const m1& v, m2& V ) const
	{	
		unsigned 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
		
		//... boundary points
		for( int j=0,j2=0; j<ny; ++j,j2+=2 )
			for( int k=0,k2=0; k<nz; ++k,k2+=2 ){
				V(ox-1,j+oy,k+oz)  = 0.125*(v(-1,j2,k2)+v(-1,j2+1,k2)+v(-1,j2,k2+1)+v(-1,j2+1,k2+1)+
											v(-2,j2,k2)+v(-2,j2+1,k2)+v(-2,j2,k2+1)+v(-2,j2+1,k2+1));
				V(ox+nx,j+oy,k+oz) = 0.125*(v(2*nx,j2,k2)+v(2*nx,j2+1,k2)+v(2*nx,j2,k2+1)+v(2*nx,j2+1,k2+1)+
											v(2*nx+1,j2,k2)+v(2*nx+1,j2+1,k2)+v(2*nx+1,j2,k2+1)+v(2*nx+1,j2+1,k2+1));
			}
		
		for( int i=0,i2=0; i<nx; ++i,i2+=2 )
			for( int k=0,k2=0; k<nz; ++k,k2+=2 ){
				V(i+ox,oy-1,k+oz)  = 0.125*(v(i2,-1,k2)+v(i2+1,-1,k2)+v(i2,-1,k2+1)+v(i2+1,-1,k2+1)+
										    v(i2,-2,k2)+v(i2+1,-2,k2)+v(i2,-2,k2+1)+v(i2+1,-2,k2+1) );
				V(i+ox,oy+ny,k+oz) = 0.125*(v(i2,2*ny,k2)+v(i2+1,2*ny,k2)+v(i2,2*ny,k2+1)+v(i2+1,2*ny,k2+1)+
											v(i2,2*ny+1,k2)+v(i2+1,2*ny+1,k2)+v(i2,2*ny+1,k2+1)+v(i2+1,2*ny+1,k2+1));
			}
		
		for( int i=0,i2=0; i<nx; ++i,i2+=2 )
			for( int j=0,j2=0; j<ny; ++j,j2+=2 ){
				V(i+ox,j+oy,oz-1)  = 0.125*(v(i2,j2,-1)+v(i2+1,j2,-1)+v(i2,j2+1,-1)+v(i2+1,j2+1,-1)+
											v(i2,j2,-2)+v(i2+1,j2,-2)+v(i2,j2+1,-2)+v(i2+1,j2+1,-2));
				V(i+ox,j+oy,oz+nz) = 0.125*(v(i2,j2,2*nz)+v(i2+1,j2,2*nz)+v(i2,j2+1,2*nz)+v(i2+1,j2+1,2*nz)+
											v(i2,j2,2*nz+1)+v(i2+1,j2,2*nz+1)+v(i2,j2+1,2*nz+1)+v(i2+1,j2+1,2*nz+1));
			}
		
		
	}
	
	template< typename m1, typename m2 >
	inline void restrict_bnd_add( const m1& v, m2& V ) const
	{
		
		unsigned 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
		
		//... boundary points
		for( int j=0,j2=0; j<ny; ++j,j2+=2 )
			for( int k=0,k2=0; k<nz; ++k,k2+=2 ){
				V(ox-1,j+oy,k+oz)  += 0.125*(v(-1,j2,k2)+v(-1,j2+1,k2)+v(-1,j2,k2+1)+v(-1,j2+1,k2+1)+
									 		 v(-2,j2,k2)+v(-2,j2+1,k2)+v(-2,j2,k2+1)+v(-2,j2+1,k2+1));
				V(ox+nx,j+oy,k+oz) += 0.125*(v(2*nx,j2,k2)+v(2*nx,j2+1,k2)+v(2*nx,j2,k2+1)+v(2*nx,j2+1,k2+1)+
											 v(2*nx+1,j2,k2)+v(2*nx+1,j2+1,k2)+v(2*nx+1,j2,k2+1)+v(2*nx+1,j2+1,k2+1));
			}
		
		for( int i=0,i2=0; i<nx; ++i,i2+=2 )
			for( int k=0,k2=0; k<nz; ++k,k2+=2 ){
				V(i+ox,oy-1,k+oz)  += 0.125*(v(i2,-1,k2)+v(i2+1,-1,k2)+v(i2,-1,k2+1)+v(i2+1,-1,k2+1)+
										     v(i2,-2,k2)+v(i2+1,-2,k2)+v(i2,-2,k2+1)+v(i2+1,-2,k2+1) );
				V(i+ox,oy+ny,k+oz) += 0.125*(v(i2,2*ny,k2)+v(i2+1,2*ny,k2)+v(i2,2*ny,k2+1)+v(i2+1,2*ny,k2+1)+
											 v(i2,2*ny+1,k2)+v(i2+1,2*ny+1,k2)+v(i2,2*ny+1,k2+1)+v(i2+1,2*ny+1,k2+1));
			}
		
		for( int i=0,i2=0; i<nx; ++i,i2+=2 )
			for( int j=0,j2=0; j<ny; ++j,j2+=2 ){
				V(i+ox,j+oy,oz-1)  += 0.125*(v(i2,j2,-1)+v(i2+1,j2,-1)+v(i2,j2+1,-1)+v(i2+1,j2+1,-1)+
											 v(i2,j2,-2)+v(i2+1,j2,-2)+v(i2,j2+1,-2)+v(i2+1,j2+1,-2));
				V(i+ox,j+oy,oz+nz) += 0.125*(v(i2,j2,2*nz)+v(i2+1,j2,2*nz)+v(i2,j2+1,2*nz)+v(i2+1,j2+1,2*nz)+
											 v(i2,j2,2*nz+1)+v(i2+1,j2,2*nz+1)+v(i2,j2+1,2*nz+1)+v(i2+1,j2+1,2*nz+1));
			}
		
		
	}
	
	
	

			
	template< typename m1, typename m2 >
	inline void prolong( m1& V, m2& v ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
		#pragma omp parallel for
		for( int i=0; i<nx; ++i )
		{
			int i2(2*i);
			for( int j=0,j2=0; j<ny; ++j,j2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 )
				{
					double sum = 0.0;
					sum+= (v(i2+1,j2,k2)		= interp_cubic<1,0,0>(i+ox,j+oy,k+oz,V));
					sum+= (v(i2,j2+1,k2)		= interp_cubic<0,1,0>(i+ox,j+oy,k+oz,V));
					sum+= (v(i2,j2,k2+1)		= interp_cubic<0,0,1>(i+ox,j+oy,k+oz,V));
					sum+= (v(i2+1,j2+1,k2)		= interp_cubic<1,1,0>(i+ox,j+oy,k+oz,V));
					sum+= (v(i2+1,j2,k2+1)		= interp_cubic<1,0,1>(i+ox,j+oy,k+oz,V));
					sum+= (v(i2+1,j2+1,k2+1)	= interp_cubic<1,1,1>(i+ox,j+oy,k+oz,V));
					sum+= (v(i2,j2+1,k2+1)		= interp_cubic<0,1,1>(i+ox,j+oy,k+oz,V));
					sum+= (v(i2,j2,k2)			= interp_cubic<0,0,0>(i+ox,j+oy,k+oz,V));
					sum *= 0.125;
					
					double dd = V(i+ox,j+oy,k+oz)-sum;
					v(i2+1,j2,k2)		+= dd;
					v(i2,j2+1,k2)		+= dd;
					v(i2,j2,k2+1)		+= dd;
					v(i2+1,j2+1,k2)		+= dd;
					v(i2+1,j2,k2+1)		+= dd;
					v(i2+1,j2+1,k2+1)	+= dd;
					v(i2,j2+1,k2+1)		+= dd;
					v(i2,j2,k2)			+= dd;
				}
		}
	}

	template< typename m1, typename m2 >
	inline void prolong_add( m1& V, m2& v ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
#pragma omp parallel for
		for( int i=0; i<nx; ++i )
		{
			int i2(2*i);
			for( int j=0,j2=0; j<ny; ++j,j2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 )
				{
					double val, sum = 0.0;
					v(i2+1,j2,k2)		+= (val=interp_cubic<1,0,0>(i+ox,j+oy,k+oz,V)); sum +=val;
					v(i2,j2+1,k2)		+= (val=interp_cubic<0,1,0>(i+ox,j+oy,k+oz,V)); sum +=val;
					v(i2,j2,k2+1)		+= (val=interp_cubic<0,0,1>(i+ox,j+oy,k+oz,V)); sum +=val;
					v(i2+1,j2+1,k2)		+= (val=interp_cubic<1,1,0>(i+ox,j+oy,k+oz,V)); sum +=val;
					v(i2+1,j2,k2+1)		+= (val=interp_cubic<1,0,1>(i+ox,j+oy,k+oz,V)); sum +=val;
					v(i2+1,j2+1,k2+1)	+= (val=interp_cubic<1,1,1>(i+ox,j+oy,k+oz,V)); sum +=val;
					v(i2,j2+1,k2+1)		+= (val=interp_cubic<0,1,1>(i+ox,j+oy,k+oz,V)); sum +=val;
					v(i2,j2,k2)			+= (val=interp_cubic<0,0,0>(i+ox,j+oy,k+oz,V)); sum +=val;
					
					sum *= 0.125;
					
					double dd = V(i+ox,j+oy,k+oz)-sum;
					v(i2+1,j2,k2)		+= dd;
					v(i2,j2+1,k2)		+= dd;
					v(i2,j2,k2+1)		+= dd;
					v(i2+1,j2+1,k2)		+= dd;
					v(i2+1,j2,k2+1)		+= dd;
					v(i2+1,j2+1,k2+1)	+= dd;
					v(i2,j2+1,k2+1)		+= dd;
					v(i2,j2,k2)			+= dd;
				}
		}
	}
	
	template< typename m1, typename m2 >
	inline void prolong_bnd( m1& V, m2& v )
	{
		int 
			nx = v.size(0)/2, 
			ny = v.size(1)/2, 
			nz = v.size(2)/2,
			ox = v.offset(0),
			oy = v.offset(1),
			oz = v.offset(2);
			
		int nbnd = V.m_nbnd;
		int nbndtop = nbnd/2;
		
		
		#pragma omp parallel for
		for( int i=-nbndtop; i<nx+nbndtop; ++i )
		{	
			int i2(2*i);
			for( int j=-nbndtop,j2=-2*nbndtop; j<ny+nbndtop; ++j,j2+=2 )
				for( int k=-nbndtop,k2=-2*nbndtop; k<nz+nbndtop; ++k,k2+=2 )
				{
					
					if( i>=0&&i<nx&&j>=0&&j<ny&&k>=0&&k<nz )
						continue;
					
					v(i2+1,j2,k2)		= interp_cubic<1,0,0>(i+ox,j+oy,k+oz,V);
					v(i2,j2+1,k2)		= interp_cubic<0,1,0>(i+ox,j+oy,k+oz,V);
					v(i2,j2,k2+1)		= interp_cubic<0,0,1>(i+ox,j+oy,k+oz,V);
					v(i2+1,j2+1,k2)		= interp_cubic<1,1,0>(i+ox,j+oy,k+oz,V);
					v(i2+1,j2,k2+1)		= interp_cubic<1,0,1>(i+ox,j+oy,k+oz,V);
					v(i2+1,j2+1,k2+1)	= interp_cubic<1,1,1>(i+ox,j+oy,k+oz,V);
					v(i2,j2+1,k2+1)		= interp_cubic<0,1,1>(i+ox,j+oy,k+oz,V);
					v(i2,j2,k2)			= interp_cubic<0,0,0>(i+ox,j+oy,k+oz,V);
				}
		}
	}
	
	template< typename m1, typename m2 >
	inline void prolong_add_bnd( m1& V, m2& v )
	{
		unsigned 
		nx = V.size(0), 
		ny = V.size(1), 
		nz = V.size(2);
		
		#pragma omp parallel for
		for( int i=-1; i<=nx; ++i )
		{	
			int i2(2*i);
			for( int j=-1,j2=-2; j<=ny; ++j,j2+=2 )
				for( int k=-1,k2=-2; k<=nz; ++k,k2+=2 )
				{
					
					if( i>=0&&i<nx&&j>=0&&j<ny&&k>=0&&k<nz )
						continue;
					
					v(i2+1,j2,k2)		+= interp_cubic<1,0,0>(i,j,k,V);
					v(i2,j2+1,k2)		+= interp_cubic<0,1,0>(i,j,k,V);
					v(i2,j2,k2+1)		+= interp_cubic<0,0,1>(i,j,k,V);
					v(i2+1,j2+1,k2)		+= interp_cubic<1,1,0>(i,j,k,V);
					v(i2+1,j2,k2+1)		+= interp_cubic<1,0,1>(i,j,k,V);
					v(i2+1,j2+1,k2+1)	+= interp_cubic<1,1,1>(i,j,k,V);
					v(i2,j2+1,k2+1)		+= interp_cubic<0,1,1>(i,j,k,V);
					v(i2,j2,k2)			+= interp_cubic<0,0,0>(i,j,k,V);
				}
		}
	}
	
	
};


/*
class mg_lin
{
public:
	
	template< typename M >
	inline double restr_lin( int x, int y, int z, M& V ) const
	{
		double w[4] = { 1.0/4.0, 3.0/4.0, 3.0/4.0, 1.0/4.0 };
		double vox = 0.0;
		int i,j,k;
		
		for( i=0; i<4; ++i )
			for( j=0; j<4; ++j )
				for( k=0; k<4; ++k )
					vox += w[i]*w[j]*w[k] * V(x+i-1,y+j-1,z+k-1);
		return vox;
	}
	
	template< int sx, int sy, int sz, typename M >
	inline double interp_lin( int x, int y, int z, M& V, double s=1.0 ) const
	{
		double           u[2], v[2], w[2];
		double           vox = 0;
		int				 i,j,k;
		
		
		if( sx==0 ){
			u[0] = 1.0/4.0;
			u[1] = 3.0/4.0;
		}else{
			u[0] = 3.0/4.0;
			u[1] = 1.0/4.0;
		}
		if( sy==0 ){
			v[0] = 1.0/4.0;
			v[1] = 3.0/4.0;
		}else{
			v[0] = 3.0/4.0;
			v[1] = 1.0/4.0;
		}
		if( sz==0 ){
			w[0] = 1.0/4.0;
			w[1] = 3.0/4.0;
		}else{
			w[0] = 3.0/4.0;
			w[1] = 1.0/4.0;
		}
		
		for( i=0; i<2; ++i )
			for( j=0; j<2; ++j )
				for( k=0; k<2; ++k )
					vox += u[i]*v[j]*w[k] * V(x+i+sx-1,y+j+sy-1,z+k+sz-1);
		
		return vox;
	}
	
	
	template< typename m1, typename m2 >
	inline void prolong( const m1& V, m2& v ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
		#pragma omp parallel for
		for( int i=0; i<nx; ++i )
		{
			int i2(2*i);
			for( int j=0,j2=0; j<ny; ++j,j2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 )
				{
					v(i2+1,j2,k2)		= interp_lin<1,0,0>(i+ox,j+oy,k+oz,V);
					v(i2,j2+1,k2)		= interp_lin<0,1,0>(i+ox,j+oy,k+oz,V);
					v(i2,j2,k2+1)		= interp_lin<0,0,1>(i+ox,j+oy,k+oz,V);
					v(i2+1,j2+1,k2)		= interp_lin<1,1,0>(i+ox,j+oy,k+oz,V);
					v(i2+1,j2,k2+1)		= interp_lin<1,0,1>(i+ox,j+oy,k+oz,V);
					v(i2+1,j2+1,k2+1)	= interp_lin<1,1,1>(i+ox,j+oy,k+oz,V);
					v(i2,j2+1,k2+1)		= interp_lin<0,1,1>(i+ox,j+oy,k+oz,V);
					v(i2,j2,k2)			= interp_lin<0,0,0>(i+ox,j+oy,k+oz,V);
				}
		}
	}
	
	template< typename m1, typename m2 >
	inline void prolong_add( const m1& V, m2& v ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
#pragma omp parallel for
		for( int i=0; i<nx; ++i )
		{
			int i2(2*i);
			for( int j=0,j2=0; j<ny; ++j,j2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 )
				{
					v(i2+1,j2,k2)		+= interp_lin<1,0,0>(i+ox,j+oy,k+oz,V);
					v(i2,j2+1,k2)		+= interp_lin<0,1,0>(i+ox,j+oy,k+oz,V);
					v(i2,j2,k2+1)		+= interp_lin<0,0,1>(i+ox,j+oy,k+oz,V);
					v(i2+1,j2+1,k2)		+= interp_lin<1,1,0>(i+ox,j+oy,k+oz,V);
					v(i2+1,j2,k2+1)		+= interp_lin<1,0,1>(i+ox,j+oy,k+oz,V);
					v(i2+1,j2+1,k2+1)	+= interp_lin<1,1,1>(i+ox,j+oy,k+oz,V);
					v(i2,j2+1,k2+1)		+= interp_lin<0,1,1>(i+ox,j+oy,k+oz,V);
					v(i2,j2,k2)			+= interp_lin<0,0,0>(i+ox,j+oy,k+oz,V);
				}
		}
	}
	
	
	template< typename m1, typename m2 >
	inline void restrict( const m1& v, m2& V ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
	#pragma omp parallel for
		for( int i=0; i<nx; ++i )
		{
			int i2 = 2*i;
			for( int j=0,j2=0; j<ny; ++j,j2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 )
					V(i+ox,j+oy,k+oz) = 0.125 * restr_lin(i2, j2, k2, v);
		}		
		
	}	
};
*/

//! linear grid injection/restriction operator
//template< typename T >
class mg_linear
{
public:
	//typedef T real_t;
	
	
	//inline void prolong_bnd( const MeshvarBnd<real_t>& V, MeshvarBnd<real_t>& v ) const
	template< typename m1, typename m2 >
	inline void prolong_bnd( const m1& V, m2& v ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
#pragma omp parallel 
		{
			
			{
				for(int q=0;q<2; ++q)
				{
					int i=nx;
					if(q==0) i=-1;
					
					int i2 = 2*i;
					for( int j=0,j2=0; j<ny; ++j,j2+=2 )
						for( int k=0,k2=0; k<nz; ++k,k2+=2 )
						{
							int ii=i+ox,jj=j+oy,kk=k+oz;
							v(i2,j2,k2)			= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii-1,jj-1,kk-1))/64.;
							v(i2+1,j2,k2)		= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii+1,jj-1,kk-1))/64.;
							v(i2,j2+1,k2)		= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj+1,kk)+V(ii-1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii-1,jj+1,kk-1))/64.;
							v(i2,j2,k2+1)		= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii-1,jj-1,kk+1))/64.;
							
							v(i2+1,j2+1,k2)		= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii+1,jj+1,kk-1))/64.;
							v(i2+1,j2,k2+1)		= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii+1,jj-1,kk+1))/64.;
							v(i2,j2+1,k2+1)		= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii-1,jj+1,kk+1))/64.;
							v(i2+1,j2+1,k2+1)	= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii+1,jj+1,kk+1))/64.;
						}
				}
			}
			
			{
				for(int q=0;q<2; ++q)
				{
					int j=ny;
					if(q==0) j=-1;
					
					int j2 = 2*j;
					for( int i=0,i2=0; i<nx; ++i,i2+=2 )
						for( int k=0,k2=0; k<nz; ++k,k2+=2 )
						{
							int ii=i+ox,jj=j+oy,kk=k+oz;
							v(i2,j2,k2)			= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii-1,jj-1,kk-1))/64.;
							v(i2+1,j2,k2)		= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii+1,jj-1,kk-1))/64.;
							v(i2,j2+1,k2)		= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj+1,kk)+V(ii-1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii-1,jj+1,kk-1))/64.;
							v(i2,j2,k2+1)		= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii-1,jj-1,kk+1))/64.;
							
							v(i2+1,j2+1,k2)		= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii+1,jj+1,kk-1))/64.;
							v(i2+1,j2,k2+1)		= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii+1,jj-1,kk+1))/64.;
							v(i2,j2+1,k2+1)		= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii-1,jj+1,kk+1))/64.;
							v(i2+1,j2+1,k2+1)	= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii+1,jj+1,kk+1))/64.;
						}
				}
			}
			
			{
				for(int q=0;q<2; ++q)
				{
					int k = nz;
					if(q==0) k=-1;
					
					int k2 = 2*k;
					for( int i=0,i2=0; i<nx; ++i,i2+=2 )
						for( int j=0,j2=0; j<nx; ++i,i2+=2 )
						{
							int ii=i+ox,jj=j+oy,kk=k+oz;
							v(i2,j2,k2)			= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii-1,jj-1,kk-1))/64.;
							v(i2+1,j2,k2)		= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii+1,jj-1,kk-1))/64.;
							v(i2,j2+1,k2)		= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj+1,kk)+V(ii-1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii-1,jj+1,kk-1))/64.;
							v(i2,j2,k2+1)		= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii-1,jj-1,kk+1))/64.;
							
							v(i2+1,j2+1,k2)		= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii+1,jj+1,kk-1))/64.;
							v(i2+1,j2,k2+1)		= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii+1,jj-1,kk+1))/64.;
							v(i2,j2+1,k2+1)		= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii-1,jj+1,kk+1))/64.;
							v(i2+1,j2+1,k2+1)	= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii+1,jj+1,kk+1))/64.;
							
							
						}
				}
			}
			
			
#pragma omp barrier
		}
	}
	
	template< typename m1, typename m2 >
	inline void prolong_add_bnd( const m1& V, m2& v ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
#pragma omp parallel 
		{
			
			{
				for(int q=0;q<2; ++q)
				{
					int i=nx;
					if(q==0) i=-1;
					
					int i2 = 2*i;
					for( int j=0,j2=0; j<ny; ++j,j2+=2 )
						for( int k=0,k2=0; k<nz; ++k,k2+=2 )
						{
							int ii=i+ox,jj=j+oy,kk=k+oz;
							v(i2,j2,k2)			+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii-1,jj-1,kk-1))/64.;
							v(i2+1,j2,k2)		+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii+1,jj-1,kk-1))/64.;
							v(i2,j2+1,k2)		+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj+1,kk)+V(ii-1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii-1,jj+1,kk-1))/64.;
							v(i2,j2,k2+1)		+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii-1,jj-1,kk+1))/64.;
							
							v(i2+1,j2+1,k2)		+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii+1,jj+1,kk-1))/64.;
							v(i2+1,j2,k2+1)		+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii+1,jj-1,kk+1))/64.;
							v(i2,j2+1,k2+1)		+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii-1,jj+1,kk+1))/64.;
							v(i2+1,j2+1,k2+1)	+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii+1,jj+1,kk+1))/64.;
						}
				}
			}
			
			{
				for(int q=0;q<2; ++q)
				{
					int j=ny;
					if(q==0) j=-1;
					
					int j2 = 2*j;
					for( int i=0,i2=0; i<nx; ++i,i2+=2 )
						for( int k=0,k2=0; k<nz; ++k,k2+=2 )
						{
							int ii=i+ox,jj=j+oy,kk=k+oz;
							v(i2,j2,k2)			+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii-1,jj-1,kk-1))/64.;
							v(i2+1,j2,k2)		+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii+1,jj-1,kk-1))/64.;
							v(i2,j2+1,k2)		+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj+1,kk)+V(ii-1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii-1,jj+1,kk-1))/64.;
							v(i2,j2,k2+1)		+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii-1,jj-1,kk+1))/64.;
							
							v(i2+1,j2+1,k2)		+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii+1,jj+1,kk-1))/64.;
							v(i2+1,j2,k2+1)		+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii+1,jj-1,kk+1))/64.;
							v(i2,j2+1,k2+1)		+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii-1,jj+1,kk+1))/64.;
							v(i2+1,j2+1,k2+1)	+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii+1,jj+1,kk+1))/64.;
						}
				}
			}
			
			{
				for(int q=0;q<2; ++q)
				{
					int k = nz;
					if(q==0) k=-1;
					
					int k2 = 2*k;
					for( int i=0,i2=0; i<nx; ++i,i2+=2 )
						for( int j=0,j2=0; j<nx; ++i,i2+=2 )
						{
							int ii=i+ox,jj=j+oy,kk=k+oz;
							v(i2,j2,k2)			+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii-1,jj-1,kk-1))/64.;
							v(i2+1,j2,k2)		+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii+1,jj-1,kk-1))/64.;
							v(i2,j2+1,k2)		+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj+1,kk)+V(ii-1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii-1,jj+1,kk-1))/64.;
							v(i2,j2,k2+1)		+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii-1,jj-1,kk+1))/64.;
							
							v(i2+1,j2+1,k2)		+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii+1,jj+1,kk-1))/64.;
							v(i2+1,j2,k2+1)		+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii+1,jj-1,kk+1))/64.;
							v(i2,j2+1,k2+1)		+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii-1,jj+1,kk+1))/64.;
							v(i2+1,j2+1,k2+1)	+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii+1,jj+1,kk+1))/64.;
							
							
						}
				}
			}
			
			
#pragma omp barrier
		}
	}
	
	template< typename m1, typename m2 >
	inline void prolong( const m1& V, m2& v ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
#pragma omp parallel for
		for( int i=0; i<nx; ++i ){
			int i2 = 2*i;
			for( int j=0,j2=0; j<ny; ++j,j2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 )
				{
					int ii=i+ox,jj=j+oy,kk=k+oz;
					v(i2,j2,k2)			= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii-1,jj-1,kk-1))/64.;
					v(i2+1,j2,k2)		= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii+1,jj-1,kk-1))/64.;
					v(i2,j2+1,k2)		= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj+1,kk)+V(ii-1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii-1,jj+1,kk-1))/64.;
					v(i2,j2,k2+1)		= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii-1,jj-1,kk+1))/64.;
					
					v(i2+1,j2+1,k2)		= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii+1,jj+1,kk-1))/64.;
					v(i2+1,j2,k2+1)		= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii+1,jj-1,kk+1))/64.;
					v(i2,j2+1,k2+1)		= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii-1,jj+1,kk+1))/64.;
					v(i2+1,j2+1,k2+1)	= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii+1,jj+1,kk+1))/64.;
					
					
				}
		}
	}
	
	template< typename m1, typename m2 >
	inline void prolong_add( const m1& V, m2& v ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
#pragma omp parallel for
		for( int i=0; i<nx; ++i ){
			int i2 = 2*i;
			for( int j=0,j2=0; j<ny; ++j,j2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 )
				{
					int ii=i+ox,jj=j+oy,kk=k+oz;
					v(i2,j2,k2)			+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii-1,jj-1,kk-1))/64.;
					v(i2+1,j2,k2)		+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk-1)+V(ii,jj-1,kk-1))+V(ii+1,jj-1,kk-1))/64.;
					v(i2,j2+1,k2)		+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii-1,jj+1,kk)+V(ii-1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii-1,jj+1,kk-1))/64.;
					v(i2,j2,k2+1)		+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj-1,kk)+V(ii-1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii-1,jj-1,kk+1))/64.;
					
					v(i2+1,j2+1,k2)		+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk-1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk-1)+V(ii,jj+1,kk-1))+V(ii+1,jj+1,kk-1))/64.;
					v(i2+1,j2,k2+1)		+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj-1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj-1,kk)+V(ii+1,jj,kk+1)+V(ii,jj-1,kk+1))+V(ii+1,jj-1,kk+1))/64.;
					v(i2,j2+1,k2+1)		+= (27.*V(ii,jj,kk)+9.*(V(ii-1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii-1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii-1,jj+1,kk+1))/64.;
					v(i2+1,j2+1,k2+1)	+= (27.*V(ii,jj,kk)+9.*(V(ii+1,jj,kk)+V(ii,jj+1,kk)+V(ii,jj,kk+1))+3.*(V(ii+1,jj+1,kk)+V(ii+1,jj,kk+1)+V(ii,jj+1,kk+1))+V(ii+1,jj+1,kk+1))/64.;
					
					
				}
		}
	}
	
};

//! zero order grid injection/restriction operator
class mg_straight
{
public:
	
	template< typename m1, typename m2 >
	inline void restrict_bnd( const m1& v, m2& V ) const
	{	
		unsigned 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
		
		//... boundary points
		for( int j=0,j2=0; j<ny; ++j,j2+=2 )
			for( int k=0,k2=0; k<nz; ++k,k2+=2 ){
				V(ox-1,j+oy,k+oz)  = 0.25*(v(-1,j2,k2)+v(-1,j2+1,k2)+v(-1,j2,k2+1)+v(-1,j2+1,k2+1));
				V(ox+nx,j+oy,k+oz) = 0.25*(v(2*nx,j2,k2)+v(2*nx,j2+1,k2)+v(2*nx,j2,k2+1)+v(2*nx,j2+1,k2+1));
			}
		
		for( int i=0,i2=0; i<nx; ++i,i2+=2 )
			for( int k=0,k2=0; k<nz; ++k,k2+=2 ){
				V(i+ox,oy-1,k+oz)  = 0.25*(v(i2,-1,k2)+v(i2+1,-1,k2)+v(i2,-1,k2+1)+v(i2+1,-1,k2+1));
				V(i+ox,oy+ny,k+oz) = 0.25*(v(i2,2*ny,k2)+v(i2+1,2*ny,k2)+v(i2,2*ny,k2+1)+v(i2+1,2*ny,k2+1));
			}
		
		for( int i=0,i2=0; i<nx; ++i,i2+=2 )
			for( int j=0,j2=0; j<ny; ++j,j2+=2 ){
				V(i+ox,j+oy,oz-1)  = 0.25*(v(i2,j2,-1)+v(i2+1,j2,-1)+v(i2,j2+1,-1)+v(i2+1,j2+1,-1));
				V(i+ox,j+oy,oz+nz) = 0.25*(v(i2,j2,2*nz)+v(i2+1,j2,2*nz)+v(i2,j2+1,2*nz)+v(i2+1,j2+1,2*nz));
			}
		
		
	}
	
	
	//... restricts v to V
	template< typename m1, typename m2 >
	inline void restrict( const m1& v, m2& V ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
#pragma omp parallel for
		for( int i=0; i<nx; ++i )
		{
			int i2 = 2*i;
			for( int j=0,j2=0; j<ny; ++j,j2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 )
					V(i+ox,j+oy,k+oz) = 0.125 * ( v(i2+1,j2,k2) + v(i2,j2+1,k2) + v(i2,j2,k2+1) + v(i2+1,j2+1,k2) +
												 v(i2+1,j2,k2+1) + v(i2+1,j2+1,k2+1) + v(i2,j2+1,k2+1) + v(i2,j2,k2) );
		}
	}	
	
	template< typename m1, typename m2 >
	inline void restrict_add( const m1& v, m2& V ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
#pragma omp parallel for
		for( int i=0; i<nx; ++i )
		{
			int i2 = 2*i;
			for( int j=0,j2=0; j<ny; ++j,j2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 )
					V(i+ox,j+oy,k+oz) += 0.125 * ( v(i2+1,j2,k2) + v(i2,j2+1,k2) + v(i2,j2,k2+1) + v(i2+1,j2+1,k2) +
												 v(i2+1,j2,k2+1) + v(i2+1,j2+1,k2+1) + v(i2,j2+1,k2+1) + v(i2,j2,k2) );
		}
	}	
	
	template< typename m1, typename m2 >
	inline void restrict_bnd_add( const m1& v, m2& V ) const
	{	
		unsigned 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
		
		//... boundary points
		for( int j=0,j2=0; j<ny; ++j,j2+=2 )
			for( int k=0,k2=0; k<nz; ++k,k2+=2 ){
				V(ox-1,j+oy,k+oz)  += 0.25*(v(-1,j2,k2)+v(-1,j2+1,k2)+v(-1,j2,k2+1)+v(-1,j2+1,k2+1));
				V(ox+nx,j+oy,k+oz) += 0.25*(v(2*nx,j2,k2)+v(2*nx,j2+1,k2)+v(2*nx,j2,k2+1)+v(2*nx,j2+1,k2+1));
			}
		
		for( int i=0,i2=0; i<nx; ++i,i2+=2 )
			for( int k=0,k2=0; k<nz; ++k,k2+=2 ){
				V(i+ox,oy-1,k+oz)  += 0.25*(v(i2,-1,k2)+v(i2+1,-1,k2)+v(i2,-1,k2+1)+v(i2+1,-1,k2+1));
				V(i+ox,oy+ny,k+oz) += 0.25*(v(i2,2*ny,k2)+v(i2+1,2*ny,k2)+v(i2,2*ny,k2+1)+v(i2+1,2*ny,k2+1));
			}
		
		for( int i=0,i2=0; i<nx; ++i,i2+=2 )
			for( int j=0,j2=0; j<ny; ++j,j2+=2 ){
				V(i+ox,j+oy,oz-1)  += 0.25*(v(i2,j2,-1)+v(i2+1,j2,-1)+v(i2,j2+1,-1)+v(i2+1,j2+1,-1));
				V(i+ox,j+oy,oz+nz) += 0.25*(v(i2,j2,2*nz)+v(i2+1,j2,2*nz)+v(i2,j2+1,2*nz)+v(i2+1,j2+1,2*nz));
			}
		
		
	}
	
	//... prolongs V to v
	template< typename m1, typename m2 >
	inline void prolong( const m1& V, m2& v ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
#pragma omp parallel 
		{
			
#pragma omp for nowait
			for( int i=0; i<nx; ++i ){
				int i2 = 2*i;
				for( int j=0,j2=0; j<ny; ++j,j2+=2 )
					for( int k=0,k2=0; k<nz; ++k,k2+=2 )
					{
						v(i2+1,j2,k2)		= 
						v(i2,j2+1,k2)		= 
						v(i2,j2,k2+1)		= 
						v(i2+1,j2+1,k2)		= 
						v(i2+1,j2,k2+1)		= 
						v(i2+1,j2+1,k2+1)	= 
						v(i2,j2+1,k2+1)		= 
						v(i2,j2,k2)			= V(i+ox,j+oy,k+oz);
					}
			}
			
			
#pragma omp barrier
		}
		
	}
	
	//... prolongs V to v on grid boundary cells
	template< typename m1, typename m2 >
	inline void prolong_bnd( const m1& V, m2& v ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
#pragma omp parallel 
		{
			
			
			
			for( int j=0,j2=0; j<ny; ++j,j2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 ){
					v(-1,j2,k2) = 
					v(-1,j2+1,k2)=
					v(-1,j2,k2+1)=
					v(-1,j2+1,k2+1) = V(ox-1,j+oy,k+oz);
					
					v(2*nx,j2,k2)=
					v(2*nx,j2+1,k2)=
					v(2*nx,j2,k2+1)=
					v(2*nx,j2+1,k2+1)=V(ox+nx,j+oy,k+oz);
				}
			
			for( int i=0,i2=0; i<nx; ++i,i2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 ){
					v(i2,-1,k2)=
					v(i2+1,-1,k2)=
					v(i2,-1,k2+1)=
					v(i2+1,-1,k2+1)=V(i+ox,oy-1,k+oz);
					v(i2,2*ny,k2)=
					v(i2+1,2*ny,k2)=
					v(i2,2*ny,k2+1)=
					v(i2+1,2*ny,k2+1)=V(i+ox,oy+ny,k+oz);
				}
			
			for( int i=0,i2=0; i<nx; ++i,i2+=2 )
				for( int j=0,j2=0; j<ny; ++j,j2+=2 ){
					v(i2,j2,-1)=
					v(i2+1,j2,-1)=
					v(i2,j2+1,-1)=
					v(i2+1,j2+1,-1)=V(i+ox,j+oy,oz-1);
					
					v(i2,j2,2*nz)=
					v(i2+1,j2,2*nz)=
					v(i2,j2+1,2*nz)=
					v(i2+1,j2+1,2*nz)=V(i+ox,j+oy,oz+nz);
				}
			
#pragma omp barrier
		}
		
	}
	
	
	//... prolongs V to v
	template< typename m1, typename m2 >
	inline void prolong_add_bnd( const m1& V, m2& v ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
#pragma omp parallel 
		{
			
			
			
			for( int j=0,j2=0; j<ny; ++j,j2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 ){
					v(-1,j2,k2) +=V(ox-1,j+oy,k+oz);
					v(-1,j2+1,k2)+=V(ox-1,j+oy,k+oz);
					v(-1,j2,k2+1)+=V(ox-1,j+oy,k+oz);
					v(-1,j2+1,k2+1)+=V(ox-1,j+oy,k+oz);
					
					v(2*nx,j2,k2)+=V(ox+nx,j+oy,k+oz);
					v(2*nx,j2+1,k2)+=V(ox+nx,j+oy,k+oz);
					v(2*nx,j2,k2+1)+=V(ox+nx,j+oy,k+oz);
					v(2*nx,j2+1,k2+1)+=V(ox+nx,j+oy,k+oz);
				}
			
			for( int i=0,i2=0; i<nx; ++i,i2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 ){
					v(i2,-1,k2)+=V(i+ox,oy-1,k+oz);
					v(i2+1,-1,k2)+=V(i+ox,oy-1,k+oz);
					v(i2,-1,k2+1)+=V(i+ox,oy-1,k+oz);
					v(i2+1,-1,k2+1)+=V(i+ox,oy-1,k+oz);
					v(i2,2*ny,k2)+=V(i+ox,oy+ny,k+oz);
					v(i2+1,2*ny,k2)+=V(i+ox,oy+ny,k+oz);
					v(i2,2*ny,k2+1)+=V(i+ox,oy+ny,k+oz);
					v(i2+1,2*ny,k2+1)+=V(i+ox,oy+ny,k+oz);
				}
			
			for( int i=0,i2=0; i<nx; ++i,i2+=2 )
				for( int j=0,j2=0; j<ny; ++j,j2+=2 ){
					v(i2,j2,-1)+=V(i+ox,j+oy,oz-1);
					v(i2+1,j2,-1)+=V(i+ox,j+oy,oz-1);
					v(i2,j2+1,-1)+=V(i+ox,j+oy,oz-1);
					v(i2+1,j2+1,-1)+=V(i+ox,j+oy,oz-1);
					
					v(i2,j2,2*nz)+=V(i+ox,j+oy,oz+nz);
					v(i2+1,j2,2*nz)+=V(i+ox,j+oy,oz+nz);
					v(i2,j2+1,2*nz)+=V(i+ox,j+oy,oz+nz);
					v(i2+1,j2+1,2*nz)+=V(i+ox,j+oy,oz+nz);
				}
			
#pragma omp barrier
		}
		
	}
	
	//! prolongs V and adds it to v
	template< typename m1, typename m2 >
	inline void prolong_add( const m1& V, m2& v ) const
	{
		int 
		nx = v.size(0)/2, 
		ny = v.size(1)/2, 
		nz = v.size(2)/2,
		ox = v.offset(0),
		oy = v.offset(1),
		oz = v.offset(2);
		
		#pragma omp parallel for
		for( int i=0; i<nx; ++i ){
			int i2 = 2*i;
			for( int j=0,j2=0; j<ny; ++j,j2+=2 )
				for( int k=0,k2=0; k<nz; ++k,k2+=2 )
				{
					v(i2+1,j2,k2)		+= V(i+ox,j+oy,k+oz);
					v(i2,j2+1,k2)		+= V(i+ox,j+oy,k+oz);
					v(i2,j2,k2+1)		+= V(i+ox,j+oy,k+oz);
					v(i2+1,j2+1,k2)		+= V(i+ox,j+oy,k+oz);
					v(i2+1,j2,k2+1)		+= V(i+ox,j+oy,k+oz);
					v(i2+1,j2+1,k2+1)	+= V(i+ox,j+oy,k+oz);
					v(i2,j2+1,k2+1)		+= V(i+ox,j+oy,k+oz);
					v(i2,j2,k2)			+= V(i+ox,j+oy,k+oz);
				}
		}
	}
	
};


#endif //__MG_OPERATORS_HH

