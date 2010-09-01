/*
 *  mg_interp.hh
 *  FROLIC_mg
 *
 *  Created by Oliver Hahn on 5/27/10.
 *  Copyright 2010 KIPAC/Stanford University. All rights reserved.
 *
 */

#ifndef __MG_INTERP_HH
#define __MG_INTERP_HH

#include "mg_operators.hh"


struct coarse_fine_interpolation
{
	template< class G >
	void interp_coarse_fine( unsigned ilevel, G& coarse, G& fine )
	{ }
};


//! general 2nd order polynomial interpolation
inline double interp2( double x1, double x2, double x3, double f1, double f2, double f3, double x )
{
	double a,b,c;	
	a = (x1 * f3 - x3 * f1 - x2 * f3 - x1 * f2 + x2 * f1 + x3 * f2) / (x1 * x3 * x3 - x2 * x3 * x3 + x2 * x1 * x1 - x3 * x1 * x1 + x3 * x2 * x2 - x1 * x2 * x2);
	b = -(x1 * x1 * f3 - x1 * x1 * f2 - f1 * x3 * x3 + f2 * x3 * x3 - x2 * x2 * f3 + f1 * x2 * x2) / (x1 - x2) / (x1 * x2 - x1 * x3 + x3 * x3 - x2 * x3);
	c = (x1 * x1 * x2 * f3 - x1 * x1 * x3 * f2 - x2 * x2 * x1 * f3 + f2 * x1 * x3 * x3 + x2 * x2 * x3 * f1 - f1 * x2 * x3 * x3) / (x1 - x2) / (x1 * x2 - x1 * x3 + x3 * x3 - x2 * x3);
	
	return a*x*x+b*x+c;
}


//! optimized 4th order polynomial interpolation across a left boundary for i-1 values
inline double interp4left( double f0, double f1, double f2, double f3, double f4 )
{
	//return -4.0/231.0*f0+4.0/7.0*f1+5.0/7.0*f2-1.0/3.0*f3+5./77.0*f4; 
	return 1.0/13.0*f0-21./55.*f1+7.0/9.0*f2+8./15.*f3-8./1287*f4;
}

//! optimized 4th order polynomial interpolation across a right boundary for i+1 values
inline double interp4right( double f0, double f1, double f2, double f3, double f4 )
{
	return interp4left(f4,f3,f2,f1,f0);
}

//! optimized 4th order polynomial interpolation across a left boundary for i-2 values
inline double interp4lleft( double f0, double f1, double f2, double f3, double f4 )
{
	//return 16./231.*f0+48.0/35.0*f1-6.0/7.0*f2-8.0/15.0*f3-9./77.*f4;
	return -15./91.*f0+8./11.*f1+-10./9.*f2+32./21.*f3+32./1287.*f4;
}

//! optimized 4th order polynomial interpolation across a right boundary for i+2 values
inline double interp4rright( double f0, double f1, double f2, double f3, double f4 )
{
	return interp4lleft(f4,f3,f2,f1,f0);
}

//! general 4th order polynomial interpolation
inline double interp4( double fm2, double fm1, double f0, double fp1, double fp2, double x )
{
	double x2 = x*x, x3=x2*x, x4=x3*x;
	double a,b,c,d,e;
	
	a= 1.0/24.0*fm2-1.0/6.0*fm1+0.25*f0-1.0/6.0*fp1+1.0/24.0*fp2;
	b=-1.0/12.0*fm2+1.0/6.0*fm1-1.0/6.0*fp1+1.0/12.0*fp2;
	c=-1.0/24.0*fm2+2.0/3.0*fm1-5.0/4.0*f0+2.0/3.0*fp1-1.0/24.0*fp2;
	d= 1.0/12.0*fm2-2.0/3.0*fm1+2.0/3.0*fp1-1.0/12.0*fp2;
	e=f0;
	
	return a*x4+b*x3+c*x2+d*x+e;
}

//! general 4th order polynomial interpolation
inline double interp4( double* f, double x )
{
	double x2 = x*x, x3=x2*x, x4=x3*x;
	double a,b,c,d,e;
	
	a= 1.0/24.0*f[0]-1.0/6.0*f[1]+0.25*f[2]-1.0/6.0*f[3]+1.0/24.0*f[4];
	b=-1.0/12.0*f[0]+1.0/6.0*f[1]-1.0/6.0*f[3]+1.0/12.0*f[4];
	c=-1.0/24.0*f[0]+2.0/3.0*f[1]-5.0/4.0*f[2]+2.0/3.0*f[3]-1.0/24.0*f[4];
	d= 1.0/12.0*f[0]-2.0/3.0*f[1]+2.0/3.0*f[3]-1.0/12.0*f[4];
	e=f[2];
	
	return a*x4+b*x3+c*x2+d*x+e;
}

//! general 6th order polynomial interpolation
inline double interp6( double *f, double x )
{
	double x2 = x*x, x3=x2*x, x4=x3*x, x5=x4*x, x6=x5*x;
	double a,b,c,d,e,g,h;
	
	a=(f[0]-6.*f[1]+15.*f[2]-20.*f[3]+15.*f[4]-6.*f[5]+f[6])/720.;
	b=(-3.*f[0]+12.*f[1]-15.*f[2]+15*f[4]-12.*f[5]+3.*f[6])/720.;
	c=(-5.*f[0]+60.*f[1]-195.*f[2]+280.*f[3]-195.*f[4]+60.*f[5]-5.*f[6])/720.;
	d=(15.*f[0]-120.*f[1]+195.*f[2]-195.*f[4]+120.*f[5]-15.*f[6])/720.;
	e=(4.*f[0]-54.*f[1]+540.*f[2]-980.*f[3]+540.*f[4]-54.*f[5]+4.*f[6])/720.;
	g=(-12.*f[0]+108.*f[1]-540.*f[2]+540.*f[4]-108.*f[5]+12.*f[6])/720.;
	h=f[3];
	
	return a*x6+b*x5+c*x4+d*x3+e*x2+g*x+h;
	
}

//! general 6th order polynomial interpolation
inline double interp6( double f0, double f1, double f2, double f3, double f4, double f5, double f6, double x )
{
	double x2 = x*x, x3=x2*x, x4=x3*x, x5=x4*x, x6=x5*x;
	double a,b,c,d,e,g,h;
	double f[7]={
		f0,f1,f2,f3,f4,f5,f6
	};
	
	a=(f[0]-6.*f[1]+15.*f[2]-20.*f[3]+15.*f[4]-6.*f[5]+f[6])/720.;
	b=(-3.*f[0]+12.*f[1]-15.*f[2]+15*f[4]-12.*f[5]+3.*f[6])/720.;
	c=(-5.*f[0]+60.*f[1]-195.*f[2]+280.*f[3]-195.*f[4]+60.*f[5]-5.*f[6])/720.;
	d=(15.*f[0]-120.*f[1]+195.*f[2]-195.*f[4]+120.*f[5]-15.*f[6])/720.;
	e=(4.*f[0]-54.*f[1]+540.*f[2]-980.*f[3]+540.*f[4]-54.*f[5]+4.*f[6])/720.;
	g=(-12.*f[0]+108.*f[1]-540.*f[2]+540.*f[4]-108.*f[5]+12.*f[6])/720.;
	h=f[3];
	
	return a*x6+b*x5+c*x4+d*x3+e*x2+g*x+h;
	
}

//! general 2nd order polynomial interpolation
inline double interp2( double fleft, double fcenter, double fright, double x )
{
	double a,b,c;
	a = 0.5*(fleft+fright)-fcenter;
	b = 0.5*(fright-fleft);
	c = fcenter;
	
	return a*x*x+b*x+c;
}

//! optimized 2nd order polynomial interpolation for i-1 values
inline double interp2left( double fleft, double fcenter, double fright )
{
	double a,b,c;
	a = (6.0*fright-10.0*fcenter+4.0*fleft)/15.0;
	b = (-4.0*fleft+9.0*fright-5.0*fcenter)/15.0;
	c = fcenter;
	
	return a-b+c;
}

//! optimized 2nd order polynomial interpolation for i+1 values
inline double interp2right( double fleft, double fcenter, double fright )
{
	double a,b,c;
	a = (6.0*fleft-10.0*fcenter+4.0*fright)/15.0;
	b = (4.0*fright-9.0*fleft+5.0*fcenter)/15.0;
	c = fcenter;
	
	return a+b+c;
}

//! optimized 2nd order polynomial interpolation for i-2 values
inline double interp2lleft( double fleft, double fcenter, double fright )
{
	double a,b,c;
	a = (-10.0*fcenter+4.0*fleft+6.0*fright)/15.0;
	b = (-12.0*fleft+15.0*fcenter-3.0*fright)/15.0;
	c = (-3.0*fright+10.0*fcenter+8.0*fleft)/15.0;
	
	return a-b+c;
}

//! optimized 2nd order polynomial interpolation for i+2 values
inline double interp2rright( double fleft, double fcenter, double fright )
{
	double a,b,c;
	a = (-10.0*fcenter+4.0*fleft+6.0*fright)/15.0;
	b = (-12.0*fleft+15.0*fcenter-3.0*fright)/15.0;
	c = (-3.0*fright+10.0*fcenter+8.0*fleft)/15.0;
	
	return a-b+c;
}

//! optimized 6th order polynomial interpolation for i-1 values
inline double interp6left( double f0, double f1, double f2, double f3, double f4, double f5, double f6 )
{
	return 4./2431.*f0-24./1001.*f1+4./7.*f2+60./77.*f3-6./13.*f4+12./77.*f5-5./221.*f6;
}

//! optimized 6th order polynomial interpolation for i-2 values
inline double interp6lleft( double f0, double f1, double f2, double f3, double f4, double f5, double f6 )
{
	return -12./2431.*f0+40./429.*f1+4./3.*f2-10./11.*f3+28./39.*f4-3./11.*f5+28./663.*f6;
}

//! optimized 6th order polynomial interpolation for i-3 values
inline double interp6llleft( double f0, double f1, double f2, double f3, double f4, double f5, double f6 )
{
	return -36./2431.*f0+600./1001.*f1+20./21.*f2-100./77.*f3+15./13.*f4-36./77.*f5+50./663.*f6;
}

//! optimized 6th order polynomial interpolation for i+1 values
inline double interp6right( double f0, double f1, double f2, double f3, double f4, double f5, double f6 )
{
	return interp6left(f6,f5,f4,f3,f2,f1,f0);
}

//! optimized 6th order polynomial interpolation for i+2 values
inline double interp6rright( double f0, double f1, double f2, double f3, double f4, double f5, double f6 )
{
	return interp6lleft(f6,f5,f4,f3,f2,f1,f0);
}

//! optimized 6th order polynomial interpolation for i+3 values
inline double interp6rrright( double f0, double f1, double f2, double f3, double f4, double f5, double f6 )
{
	return interp6llleft(f6,f5,f4,f3,f2,f1,f0);
}




#include "fd_schemes.hh"

//! tri-cubic interpolation for non-conservative injection
struct cubic_interp
: public coarse_fine_interpolation
{
	template< class G >
	void interp_coarse_fine( unsigned ilevel, G& coarse, G& fine )
	{
		mg_cubic().prolong_bnd( coarse, fine );
		
		
		//return;
		//......................................................
		
		G *u    = &fine;
		G *utop = &coarse;
		
		int
		xoff = u->offset(0),
		yoff = u->offset(1),
		zoff = u->offset(2);
		
		//... don't do anything if we are not an additional refinement region
		if( xoff == 0 && yoff == 0 && zoff == 0 )
			return;
		
		int
		nx = u->size(0), 
		ny = u->size(1), 
		nz = u->size(2);
		
#pragma omp parallel for schedule(dynamic)
		for( int ix=-1; ix<=nx; ++ix )
			for( int iy=-1; iy<=ny; ++iy )
				for( int iz=-1; iz<=nz; ++iz )
				{
					bool xbnd=(ix==-1||ix==nx),ybnd=(iy==-1||iy==ny),zbnd=(iz==-1||iz==nz);
					
					//if(ix==-1||ix==nx||iy==-1||iy==ny||iz==-1||iz==nz)
					if( xbnd || ybnd || zbnd )
						//if( xbnd ^ ybnd ^ zbnd )
					{
						
						//... only deal with proper ghostzones
						if( (xbnd&&ybnd) || (xbnd&&zbnd) || (ybnd&&zbnd) || (xbnd&&ybnd&&zbnd))
							continue;
						
						int ixtop = (int)(0.5*(double)(ix))+xoff;
						int iytop = (int)(0.5*(double)(iy))+yoff;
						int iztop = (int)(0.5*(double)(iz))+zoff;
						
						if( ix==-1 ) ixtop=xoff-1;
						if( iy==-1 ) iytop=yoff-1;
						if( iz==-1 ) iztop=zoff-1;

						double coarse_flux, fine_flux, dflux;
				
						//double fac = 0.125*27.0/24.0, fac2 = -0.125*1.0/24.0;//0.125;//24.0/26.0*0.125*0.5;
						double fac = 0.5*0.125*27.0/24.0, fac2 = -0.5*0.125*1.0/24.0;
						//double fac = 0.125*15.0/12.0, fac2 = -0.125*1.0/12.0;
						
						if(xbnd)
						{
							if( ix==-1 )
							{
								
								fine_flux = 0.0;
								fine_flux += Laplace_flux_O4<double>().apply_x(-1,*u,ix+1,iy,iz);
								fine_flux += Laplace_flux_O4<double>().apply_x(-1,*u,ix+1,iy+1,iz);
								fine_flux += Laplace_flux_O4<double>().apply_x(-1,*u,ix+1,iy,iz+1);
								fine_flux += Laplace_flux_O4<double>().apply_x(-1,*u,ix+1,iy+1,iz+1);
								
								coarse_flux = Laplace_flux_O4<double>().apply_x(-1,*utop,ixtop+1,iytop,iztop)/2.0;
								fine_flux /= 4.0;
							
								dflux = coarse_flux - fine_flux;
								
								for(int j=0;j<2;++j)
									for( int k=0;k<2;++k)
									{
										(*u)(ix,iy+j,iz+k)   += fac*dflux;
										(*u)(ix-1,iy+j,iz+k) += fac2*dflux;
									}
								//(*u)(ix,iy+j,iz+k) += 24.0/27.0*dflux;
							}
							else if( ix==nx )									
							{
								
								fine_flux = 0.0;
								fine_flux += Laplace_flux_O4<double>().apply_x(+1,*u,ix,iy,iz);
								fine_flux += Laplace_flux_O4<double>().apply_x(+1,*u,ix,iy+1,iz);
								fine_flux += Laplace_flux_O4<double>().apply_x(+1,*u,ix,iy,iz+1);
								fine_flux += Laplace_flux_O4<double>().apply_x(+1,*u,ix,iy+1,iz+1);
								
								coarse_flux = Laplace_flux_O4<double>().apply_x(+1,*utop,ixtop,iytop,iztop)/2.0;
								fine_flux /= 4.0;
								
								dflux = coarse_flux - fine_flux;
								
								for(int j=0;j<2;++j)
									for( int k=0;k<2;++k)
									{
										(*u)(ix,iy+j,iz+k)   += fac*dflux;
										(*u)(ix+1,iy+j,iz+k) += fac2*dflux;
									}
								//(*u)(ix,iy+j,iz+k) += 24.0/27.0*dflux;
							}
						}
						
						if(ybnd)
						{
							if( iy==-1 )
							{
								
								fine_flux = 0.0;
								fine_flux += Laplace_flux_O4<double>().apply_y(-1,*u,ix,iy+1,iz);
								fine_flux += Laplace_flux_O4<double>().apply_y(-1,*u,ix+1,iy+1,iz);
								fine_flux += Laplace_flux_O4<double>().apply_y(-1,*u,ix,iy+1,iz+1);
								fine_flux += Laplace_flux_O4<double>().apply_y(-1,*u,ix+1,iy+1,iz+1);
								
								coarse_flux = Laplace_flux_O4<double>().apply_y(-1,*utop,ixtop,iytop+1,iztop)/2.0;
								fine_flux /= 4.0;
								
								dflux = coarse_flux - fine_flux;
								
								for(int i=0;i<2;++i)
									for( int k=0;k<2;++k)
									{
										(*u)(ix+i,iy,iz+k)   += fac*dflux;
										(*u)(ix+i,iy-1,iz+k) += fac2*dflux;
									}
							}
							else if( iy==ny )									
							{
								
								fine_flux = 0.0;
								fine_flux += Laplace_flux_O4<double>().apply_y(+1,*u,ix,iy,iz);
								fine_flux += Laplace_flux_O4<double>().apply_y(+1,*u,ix+1,iy,iz);
								fine_flux += Laplace_flux_O4<double>().apply_y(+1,*u,ix,iy,iz+1);
								fine_flux += Laplace_flux_O4<double>().apply_y(+1,*u,ix+1,iy,iz+1);
								
								coarse_flux = Laplace_flux_O4<double>().apply_y(+1,*utop,ixtop,iytop,iztop)/2.0;
								fine_flux /= 4.0;
								
								dflux = coarse_flux - fine_flux;
								
								for(int i=0;i<2;++i)
									for( int k=0;k<2;++k)
									{
										(*u)(ix+i,iy,iz+k)   += fac*dflux;
										(*u)(ix+i,iy+1,iz+k) += fac2*dflux;
									}
							}
						}
						
						if(zbnd)
						{
							if( iz==-1 )
							{
								
								fine_flux = 0.0;
								fine_flux += Laplace_flux_O4<double>().apply_z(-1,*u,ix,iy,iz+1);
								fine_flux += Laplace_flux_O4<double>().apply_z(-1,*u,ix+1,iy,iz+1);
								fine_flux += Laplace_flux_O4<double>().apply_z(-1,*u,ix,iy+1,iz+1);
								fine_flux += Laplace_flux_O4<double>().apply_z(-1,*u,ix+1,iy+1,iz+1);
								
								coarse_flux = Laplace_flux_O4<double>().apply_z(-1,*utop,ixtop,iytop,iztop+1)/2.0;
								fine_flux /= 4.0;
								
								dflux = coarse_flux - fine_flux;
								
								for(int i=0;i<2;++i)
									for( int j=0;j<2;++j)
									{
										(*u)(ix+i,iy+j,iz)   += fac*dflux;
										(*u)(ix+i,iy+j,iz-1) += fac2*dflux;
									}
							}
							else if( iz==nz )									
							{
								
								fine_flux = 0.0;
								fine_flux += Laplace_flux_O4<double>().apply_z(+1,*u,ix,iy,iz);
								fine_flux += Laplace_flux_O4<double>().apply_z(+1,*u,ix+1,iy,iz);
								fine_flux += Laplace_flux_O4<double>().apply_z(+1,*u,ix,iy+1,iz);
								fine_flux += Laplace_flux_O4<double>().apply_z(+1,*u,ix+1,iy+1,iz);
								
								coarse_flux = Laplace_flux_O4<double>().apply_z(+1,*utop,ixtop,iytop,iztop)/2.0;
								fine_flux /= 4.0;
								
								dflux = coarse_flux - fine_flux;
								
								for(int i=0;i<2;++i)
									for( int j=0;j<2;++j)
									{
										(*u)(ix+i,iy+j,iz)   += fac*dflux;
										(*u)(ix+i,iy+j,iz+1) += fac2*dflux;
									}
							}
						}
					}
				}
		
	}
};


//! 3rd order flux-corrected coarse-fine interpolation 
struct interp_O3_fluxcorr
: public coarse_fine_interpolation
{
	template< class G >
	void interp_coarse_fine( unsigned ilevel, G& coarse, G& fine )
	{
		
		//... use cubic interpolation to get all values on boundary
		mg_cubic().prolong_bnd( coarse, fine );

		//... use flux corrected quadratic interpolation for the
		//... the boundary overlapped by the Laplace stencil
		G *u    = &fine;
		G *utop = &coarse;
		
		int
			xoff = u->offset(0),
			yoff = u->offset(1),
			zoff = u->offset(2);
		
		//... don't do anything if we are not an additional refinement region
		if( xoff == 0 && yoff == 0 && zoff == 0 )
			return;
		
		int
			nx = u->size(0), 
			ny = u->size(1), 
			nz = u->size(2);
		
		//... set boundary condition for fine grid
		#pragma omp parallel for schedule(dynamic)
		for( int ix=-1; ix<=nx; ++ix )
			for( int iy=-1; iy<=ny; ++iy )
				for( int iz=-1; iz<=nz; ++iz )
				{
					bool xbnd=(ix==-1||ix==nx),ybnd=(iy==-1||iy==ny),zbnd=(iz==-1||iz==nz);
					
					if( xbnd || ybnd || zbnd )
					{
						
						//... only deal with proper ghostzones
						if( (xbnd&&ybnd) || (xbnd&&zbnd) || (ybnd&&zbnd) || (xbnd&&ybnd&&zbnd))
							continue;
						
						int ixtop = (int)(0.5*(double)(ix))+xoff;
						int iytop = (int)(0.5*(double)(iy))+yoff;
						int iztop = (int)(0.5*(double)(iz))+zoff;
						
						if( ix==-1 ) ixtop=xoff-1;
						if( iy==-1 ) iytop=yoff-1;
						if( iz==-1 ) iztop=zoff-1;
						
						double ustar1, ustar2, ustar3, uhat;			
						double fac = 0.5;//0.25;
						double flux;;
						if( ix == -1 && iy%2==0 && iz%2==0 )
						{
							flux = 0.0;
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{		
									ustar1 = interp2( (*utop)(ixtop,iytop-1,iztop-1),(*utop)(ixtop,iytop,iztop-1),(*utop)(ixtop,iytop+1,iztop-1), fac*((double)j-0.5) );
									ustar2 = interp2( (*utop)(ixtop,iytop-1,iztop),(*utop)(ixtop,iytop,iztop),(*utop)(ixtop,iytop+1,iztop), fac*((double)j-0.5) );
									ustar3 = interp2( (*utop)(ixtop,iytop-1,iztop+1),(*utop)(ixtop,iytop,iztop+1),(*utop)(ixtop,iytop+1,iztop+1), fac*((double)j-0.5) );
									
									uhat   = interp2( ustar1, ustar2, ustar3, fac*((double)k-0.5) );
									
									(*u)(ix,iy+j,iz+k) = interp2left( uhat, (*u)(ix+1,iy+j,iz+k), (*u)(ix+2,iy+j,iz+k) );
									
									flux += ((*u)(ix+1,iy+j,iz+k)-(*u)(ix,iy+j,iz+k));
								}
							
							flux /= 4.0;
							
							double dflux = ((*utop)(ixtop+1,iytop,iztop)-(*utop)(ixtop,iytop,iztop))/2.0 - flux;
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
									(*u)(ix,iy+j,iz+k) -= dflux;
						}
						// right boundary
						if( ix == nx && iy%2==0 && iz%2==0 )
						{
							flux = 0.0;
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{		
									ustar1 = interp2( (*utop)(ixtop,iytop-1,iztop-1),(*utop)(ixtop,iytop,iztop-1),(*utop)(ixtop,iytop+1,iztop-1), fac*((double)j-0.5) );
									ustar2 = interp2( (*utop)(ixtop,iytop-1,iztop),(*utop)(ixtop,iytop,iztop),(*utop)(ixtop,iytop+1,iztop), fac*((double)j-0.5) );
									ustar3 = interp2( (*utop)(ixtop,iytop-1,iztop+1),(*utop)(ixtop,iytop,iztop+1),(*utop)(ixtop,iytop+1,iztop+1), fac*((double)j-0.5) );
									
									uhat   = interp2( -1.0, 0.0, 1.0, ustar1, ustar2, ustar3, fac*((double)k-0.5) );
									
									(*u)(ix,iy+j,iz+k) = interp2right( (*u)(ix-2,iy+j,iz+k), (*u)(ix-1,iy+j,iz+k), uhat );
									flux += ((*u)(ix,iy+j,iz+k)-(*u)(ix-1,iy+j,iz+k));
								}
							flux /= 4.0;
							double dflux = ((*utop)(ixtop,iytop,iztop)-(*utop)(ixtop-1,iytop,iztop))/2.0 - flux;
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
									(*u)(ix,iy+j,iz+k) += dflux;
						}

						// bottom boundary
						if( iy == -1 && ix%2==0 && iz%2==0 )
						{
							flux = 0.0;
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{
									ustar1 = interp2( (*utop)(ixtop-1,iytop,iztop-1),(*utop)(ixtop,iytop,iztop-1),(*utop)(ixtop+1,iytop,iztop-1), fac*(j-0.5) );
									ustar2 = interp2( (*utop)(ixtop-1,iytop,iztop),(*utop)(ixtop,iytop,iztop),(*utop)(ixtop+1,iytop,iztop), fac*(j-0.5) );
									ustar3 = interp2( (*utop)(ixtop-1,iytop,iztop+1),(*utop)(ixtop,iytop,iztop+1),(*utop)(ixtop+1,iytop,iztop+1), fac*(j-0.5) );
									
									uhat   = interp2( -1.0, 0.0, 1.0, ustar1, ustar2, ustar3, fac*((double)k-0.5) );
									
									(*u)(ix+j,iy,iz+k) = interp2left( uhat, (*u)(ix+j,iy+1,iz+k), (*u)(ix+j,iy+2,iz+k) );
									
									flux += ((*u)(ix+j,iy+1,iz+k)-(*u)(ix+j,iy,iz+k));
								}
							flux /= 4.0;
							double dflux = ((*utop)(ixtop,iytop+1,iztop)-(*utop)(ixtop,iytop,iztop))/2.0 - flux;
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
									(*u)(ix+j,iy,iz+k) -= dflux;
						}
						// top boundary
						if( iy == ny && ix%2==0 && iz%2==0 )
						{
							flux = 0.0;
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{		
									ustar1 = interp2( (*utop)(ixtop-1,iytop,iztop-1),(*utop)(ixtop,iytop,iztop-1),(*utop)(ixtop+1,iytop,iztop-1), fac*(j-0.5) );
									ustar2 = interp2( (*utop)(ixtop-1,iytop,iztop),(*utop)(ixtop,iytop,iztop),(*utop)(ixtop+1,iytop,iztop), fac*(j-0.5) );
									ustar3 = interp2( (*utop)(ixtop-1,iytop,iztop+1),(*utop)(ixtop,iytop,iztop+1),(*utop)(ixtop+1,iytop,iztop+1), fac*(j-0.5) );
									
									uhat   = interp2( -1.0, 0.0, 1.0, ustar1, ustar2, ustar3, fac*((double)k-0.5) );
									
									(*u)(ix+j,iy,iz+k) = interp2right( (*u)(ix+j,iy-2,iz+k), (*u)(ix+j,iy-1,iz+k), uhat  );
									
									flux += ((*u)(ix+j,iy,iz+k)-(*u)(ix+j,iy-1,iz+k));
								}
							flux /= 4.0;
							double dflux = ((*utop)(ixtop,iytop,iztop)-(*utop)(ixtop,iytop-1,iztop))/2.0 - flux;
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
									(*u)(ix+j,iy,iz+k) += dflux;
						}

						// front boundary
						if( iz == -1 && ix%2==0 && iy%2==0 )
						{
							flux = 0.0;
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{		
									ustar1 = interp2( (*utop)(ixtop-1,iytop-1,iztop),(*utop)(ixtop,iytop-1,iztop),(*utop)(ixtop+1,iytop-1,iztop), fac*(j-0.5) );
									ustar2 = interp2( (*utop)(ixtop-1,iytop,iztop),(*utop)(ixtop,iytop,iztop),(*utop)(ixtop+1,iytop,iztop), fac*(j-0.5) );
									ustar3 = interp2( (*utop)(ixtop-1,iytop+1,iztop),(*utop)(ixtop,iytop+1,iztop),(*utop)(ixtop+1,iytop+1,iztop), fac*(j-0.5) );
									
									uhat   = interp2( -1.0, 0.0, 1.0, ustar1, ustar2, ustar3, fac*((double)k-0.5) );
									
									(*u)(ix+j,iy+k,iz) = interp2left( uhat, (*u)(ix+j,iy+k,iz+1), (*u)(ix+j,iy+k,iz+2) );
									
									flux += ((*u)(ix+j,iy+k,iz+1)-(*u)(ix+j,iy+k,iz));
								}
							flux /= 4.0;
							double dflux = ((*utop)(ixtop,iytop,iztop+1)-(*utop)(ixtop,iytop,iztop))/2.0 - flux;
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
									(*u)(ix+j,iy+k,iz) -= dflux;
						}

						// back boundary
						if( iz == nz && ix%2==0 && iy%2==0 )
						{
							flux = 0.0;
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{		
									ustar1 = interp2( (*utop)(ixtop-1,iytop-1,iztop),(*utop)(ixtop,iytop-1,iztop),(*utop)(ixtop+1,iytop-1,iztop), fac*(j-0.5) );
									ustar2 = interp2( (*utop)(ixtop-1,iytop,iztop),(*utop)(ixtop,iytop,iztop),(*utop)(ixtop+1,iytop,iztop), fac*(j-0.5) );
									ustar3 = interp2( (*utop)(ixtop-1,iytop+1,iztop),(*utop)(ixtop,iytop+1,iztop),(*utop)(ixtop+1,iytop+1,iztop), fac*(j-0.5) );
									
									uhat   = interp2( -1.0, 0.0, 1.0, ustar1, ustar2, ustar3, fac*((double)k-0.5) );
									
									(*u)(ix+j,iy+k,iz) = interp2right( (*u)(ix+j,iy+k,iz-2), (*u)(ix+j,iy+k,iz-1), uhat );
									
									flux += ((*u)(ix+j,iy+k,iz)-(*u)(ix+j,iy+k,iz-1));
								}
							flux /= 4.0;
							double dflux = ((*utop)(ixtop,iytop,iztop)-(*utop)(ixtop,iytop,iztop-1))/2.0 - flux;
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
									(*u)(ix+j,iy+k,iz) += dflux;
						}
						
					}
				}
		
	}
};

//! 5th order flux-corrected coarse-fine interpolation
struct interp_O5_fluxcorr
: public coarse_fine_interpolation
{
	
		
	template< class G >
	void interp_coarse_fine( unsigned ilevel, G& coarse, G& fine )
	{
		
		//... use cubic interpolation to get all values on boundary
		mg_cubic().prolong_bnd( coarse, fine );
		
		//... use flux corrected quadratic interpolation for the
		//... the boundary overlapped by the Laplace stencil
		G *u    = &fine;
		G *utop = &coarse;
		
		
		int
			xoff = u->offset(0),
			yoff = u->offset(1),
			zoff = u->offset(2);
			
		//... don't do anything if we are not an additional refinement region
		if( xoff == 0 && yoff == 0 && zoff == 0 )
			return;
		
		int
		nx = u->size(0), 
		ny = u->size(1), 
		nz = u->size(2);
		
		//... set boundary condition for fine grid
		#pragma omp parallel for schedule(dynamic)
		for( int ix=-1; ix<=nx; ++ix )
			for( int iy=-1; iy<=ny; ++iy )
				for( int iz=-1; iz<=nz; ++iz )
				{
					bool xbnd=(ix==-1||ix==nx),ybnd=(iy==-1||iy==ny),zbnd=(iz==-1||iz==nz);
					bool bnd=xbnd|ybnd|zbnd;
					
					if( bnd )
					{
						
						int ixtop = (int)(0.5*(double)(ix))+xoff;
						int iytop = (int)(0.5*(double)(iy))+yoff;
						int iztop = (int)(0.5*(double)(iz))+zoff;
						
						if( ix==-1 ) ixtop=xoff-1;
						if( iy==-1 ) iytop=yoff-1;
						if( iz==-1 ) iztop=zoff-1;
						
						double ustar[5], uhat[2];			
						double fac = 0.5;
						
						double coarse_flux, fine_flux, dflux;
						
						double ffac = 12./14.;
									
						// left boundary
						if( ix == -1 && iy%2==0 && iz%2==0 )
						{
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{		
									for( int p=0; p<2; ++p )
									{
										for( int q=-2;q<=2;++q )
											ustar[q+2] = interp4( (*utop)(ixtop+p-1,iytop-2,iztop+q), (*utop)(ixtop+p-1,iytop-1,iztop+q), 
																(*utop)(ixtop+p-1,iytop,iztop+q),   (*utop)(ixtop+p-1,iytop+1,iztop+q), 
																(*utop)(ixtop+p-1,iytop+2,iztop+q), fac*((double)j-0.5) );
										uhat[p] = interp4( ustar, fac*((double)k-0.5) );//-1.5 );
									}
									
									(*u)(ix,iy+j,iz+k)   = interp4left( uhat[0], uhat[1], (*u)(ix+1,iy+j,iz+k), 
																	   (*u)(ix+2,iy+j,iz+k), (*u)(ix+3,iy+j,iz+k) );
									(*u)(ix-1,iy+j,iz+k) = interp4lleft( uhat[0], uhat[1], (*u)(ix+1,iy+j,iz+k), 
																		(*u)(ix+2,iy+j,iz+k), (*u)(ix+3,iy+j,iz+k) );
								}
							
							fine_flux = 0.0;
							fine_flux += Laplace_flux_O4<double>().apply_x(-1,*u,ix+1,iy,iz);
							fine_flux += Laplace_flux_O4<double>().apply_x(-1,*u,ix+1,iy+1,iz);
							fine_flux += Laplace_flux_O4<double>().apply_x(-1,*u,ix+1,iy,iz+1);
							fine_flux += Laplace_flux_O4<double>().apply_x(-1,*u,ix+1,iy+1,iz+1);
							fine_flux /= 4.0;
							
							coarse_flux = Laplace_flux_O4<double>().apply_x(-1,*utop,ixtop+1,iytop,iztop)/2.0;
							
							dflux = coarse_flux - fine_flux;
							
							for(int j=0;j<2;++j)
								for( int k=0;k<2;++k)
								{
									(*u)(ix,iy+j,iz+k)   += ffac*dflux;
									(*u)(ix-1,iy+j,iz+k) += ffac*dflux;
								}
							
							
						}
						// right boundary
						if( ix == nx && iy%2==0 && iz%2==0 )
						{
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{		
									for( int p=0; p<2; ++p )
									{
										for( int q=-2;q<=2;++q )
											ustar[q+2] = interp4( (*utop)(ixtop+p,iytop-2,iztop+q), (*utop)(ixtop+p,iytop-1,iztop+q), 
																 (*utop)(ixtop+p,iytop,iztop+q),   (*utop)(ixtop+p,iytop+1,iztop+q), 
																 (*utop)(ixtop+p,iytop+2,iztop+q), fac*((double)j-0.5) );
										uhat[p] = interp4( ustar, fac*((double)k-0.5));//-1.5 );
									}
									
									(*u)(ix,iy+j,iz+k)   = interp4right( (*u)(ix-3,iy+j,iz+k), (*u)(ix-2,iy+j,iz+k), 
																		(*u)(ix-1,iy+j,iz+k), uhat[0], uhat[1] );
									(*u)(ix+1,iy+j,iz+k) = interp4rright( (*u)(ix-3,iy+j,iz+k), (*u)(ix-2,iy+j,iz+k), 
																		 (*u)(ix-1,iy+j,iz+k), uhat[0], uhat[1] );
								}
							
							fine_flux = 0.0;
							fine_flux += Laplace_flux_O4<double>().apply_x(+1,*u,ix,iy,iz);
							fine_flux += Laplace_flux_O4<double>().apply_x(+1,*u,ix,iy+1,iz);
							fine_flux += Laplace_flux_O4<double>().apply_x(+1,*u,ix,iy,iz+1);
							fine_flux += Laplace_flux_O4<double>().apply_x(+1,*u,ix,iy+1,iz+1);
							
							coarse_flux = Laplace_flux_O4<double>().apply_x(+1,*utop,ixtop,iytop,iztop)/2.0;
							fine_flux /= 4.0;
							
							dflux = coarse_flux - fine_flux;
							
							for(int j=0;j<2;++j)
								for( int k=0;k<2;++k)
								{
									(*u)(ix,iy+j,iz+k)   += ffac*dflux;
									(*u)(ix+1,iy+j,iz+k) += ffac*dflux;
								}
							
						}
						// bottom boundary
						if( iy == -1 && ix%2==0 && iz%2==0 )
						{
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{	
									for( int p=0; p<2; ++p )
									{
										for( int q=-2;q<=2;++q )
											ustar[q+2] = interp4( (*utop)(ixtop-2,iytop+p-1,iztop+q), (*utop)(ixtop-1,iytop+p-1,iztop+q), 
																 (*utop)(ixtop,iytop+p-1,iztop+q),   (*utop)(ixtop+1,iytop+p-1,iztop+q), 
																 (*utop)(ixtop+2,iytop+p-1,iztop+q), fac*((double)j-0.5) );
										uhat[p] = interp4( ustar, fac*((double)k-0.5));//-1.5 );
									}
									
									(*u)(ix+j,iy,iz+k)   = interp4left( uhat[0], uhat[1], (*u)(ix+j,iy+1,iz+k), 
																	   (*u)(ix+j,iy+2,iz+k), (*u)(ix+j,iy+3,iz+k) );									
									(*u)(ix+j,iy-1,iz+k) = interp4lleft( uhat[0], uhat[1], (*u)(ix+j,iy+1,iz+k), 
																		(*u)(ix+j,iy+2,iz+k), (*u)(ix+j,iy+3,iz+k) );
								}
							
							fine_flux = 0.0;
							fine_flux += Laplace_flux_O4<double>().apply_y(-1,*u,ix,iy+1,iz);
							fine_flux += Laplace_flux_O4<double>().apply_y(-1,*u,ix+1,iy+1,iz);
							fine_flux += Laplace_flux_O4<double>().apply_y(-1,*u,ix,iy+1,iz+1);
							fine_flux += Laplace_flux_O4<double>().apply_y(-1,*u,ix+1,iy+1,iz+1);
							
							coarse_flux = Laplace_flux_O4<double>().apply_y(-1,*utop,ixtop,iytop+1,iztop)/2.0;
							fine_flux /= 4.0;
							
							dflux = coarse_flux - fine_flux;
							
							for(int i=0;i<2;++i)
								for( int k=0;k<2;++k)
								{
									(*u)(ix+i,iy,iz+k)   += ffac*dflux;
									(*u)(ix+i,iy-1,iz+k) += ffac*dflux;
								}
							
						}
						// top boundary
						if( iy == ny && ix%2==0 && iz%2==0 )
						{
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{	
									for( int p=0; p<2; ++p )
									{
										for( int q=-2;q<=2;++q )
											ustar[q+2] = interp4( (*utop)(ixtop-2,iytop+p,iztop+q), (*utop)(ixtop-1,iytop+p,iztop+q), 
																 (*utop)(ixtop,iytop+p,iztop+q),   (*utop)(ixtop+1,iytop+p,iztop+q), 
																 (*utop)(ixtop+2,iytop+p,iztop+q), fac*((double)j-0.5) );
										uhat[p] = interp4( ustar, fac*((double)k-0.5));//+1.5 );
									}
									
									(*u)(ix+j,iy,iz+k)   = interp4right( (*u)(ix+j,iy-3,iz+k), (*u)(ix+j,iy-2,iz+k), 
																		(*u)(ix+j,iy-1,iz+k), uhat[0], uhat[1] );									
									(*u)(ix+j,iy+1,iz+k) = interp4rright( (*u)(ix+j,iy-3,iz+k), (*u)(ix+j,iy-2,iz+k), 
																		 (*u)(ix+j,iy-1,iz+k), uhat[0], uhat[1] );									
								}
							
							fine_flux = 0.0;
							fine_flux += Laplace_flux_O4<double>().apply_y(+1,*u,ix,iy,iz);
							fine_flux += Laplace_flux_O4<double>().apply_y(+1,*u,ix+1,iy,iz);
							fine_flux += Laplace_flux_O4<double>().apply_y(+1,*u,ix,iy,iz+1);
							fine_flux += Laplace_flux_O4<double>().apply_y(+1,*u,ix+1,iy,iz+1);
							
							coarse_flux = Laplace_flux_O4<double>().apply_y(+1,*utop,ixtop,iytop,iztop)/2.0;
							fine_flux /= 4.0;
							
							dflux = coarse_flux - fine_flux;
							
							for(int i=0;i<2;++i)
								for( int k=0;k<2;++k)
								{
									(*u)(ix+i,iy,iz+k)   += ffac*dflux;
									(*u)(ix+i,iy+1,iz+k) += ffac*dflux;
								}
							
							
						}
						// front boundary
						if( iz == -1 && ix%2==0 && iy%2==0 )
						{
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{		
									for( int p=0; p<2; ++p )
									{
										for( int q=-2;q<=2;++q )
											ustar[q+2] = interp4( (*utop)(ixtop-2,iytop+q,iztop+p-1), (*utop)(ixtop-1,iytop+q,iztop+p-1), 
																 (*utop)(ixtop,iytop+q,iztop+p-1),   (*utop)(ixtop+1,iytop+q,iztop+p-1), 
																 (*utop)(ixtop+2,iytop+q,iztop+p-1), fac*((double)j-0.5) );
										uhat[p] = interp4( ustar, fac*((double)k-0.5));//-1.5 );
									}
									
									(*u)(ix+j,iy+k,iz)   = interp4left( uhat[0], uhat[1], (*u)(ix+j,iy+k,iz+1), 
																	   (*u)(ix+j,iy+k,iz+2), (*u)(ix+j,iy+k,iz+3) );									
									(*u)(ix+j,iy+k,iz-1) = interp4lleft( uhat[0], uhat[1], (*u)(ix+j,iy+k,iz+1), 
																		(*u)(ix+j,iy+k,iz+2), (*u)(ix+j,iy+k,iz+3) );
								}

							
							fine_flux = 0.0;
							fine_flux += Laplace_flux_O4<double>().apply_z(-1,*u,ix,iy,iz+1);
							fine_flux += Laplace_flux_O4<double>().apply_z(-1,*u,ix+1,iy,iz+1);
							fine_flux += Laplace_flux_O4<double>().apply_z(-1,*u,ix,iy+1,iz+1);
							fine_flux += Laplace_flux_O4<double>().apply_z(-1,*u,ix+1,iy+1,iz+1);
							
							coarse_flux = Laplace_flux_O4<double>().apply_z(-1,*utop,ixtop,iytop,iztop+1)/2.0;
							fine_flux /= 4.0;
							
							dflux = coarse_flux - fine_flux;
							
							for(int i=0;i<2;++i)
								for( int j=0;j<2;++j)
								{
									(*u)(ix+i,iy+j,iz)   += ffac*dflux;
									(*u)(ix+i,iy+j,iz-1) += ffac*dflux;
								}
							
						}
						// back boundary
						if( iz == nz && ix%2==0 && iy%2==0 )
						{
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{		
									for( int p=0; p<2; ++p )
									{
										for( int q=-2;q<=2;++q )
											ustar[q+2] = interp4( (*utop)(ixtop-2,iytop+q,iztop+p), (*utop)(ixtop-1,iytop+q,iztop+p), 
																 (*utop)(ixtop,iytop+q,iztop+p),   (*utop)(ixtop+1,iytop+q,iztop+p), 
																 (*utop)(ixtop+2,iytop+q,iztop+p), fac*((double)j-0.5) );
										uhat[p] = interp4( ustar, fac*((double)k-0.5));//+1.5 );
									}
									
									(*u)(ix+j,iy+k,iz)   = interp4right( (*u)(ix+j,iy+k,iz-3), (*u)(ix+j,iy+k,iz-2), 
																		(*u)(ix+j,iy+k,iz-1), uhat[0], uhat[1] );									
									(*u)(ix+j,iy+k,iz+1) = interp4rright((*u)(ix+j,iy+k,iz-3), (*u)(ix+j,iy+k,iz-2), 
																		 (*u)(ix+j,iy+k,iz-1), uhat[0], uhat[1] );
								}
							
							fine_flux = 0.0;
							fine_flux += Laplace_flux_O4<double>().apply_z(+1,*u,ix,iy,iz);
							fine_flux += Laplace_flux_O4<double>().apply_z(+1,*u,ix+1,iy,iz);
							fine_flux += Laplace_flux_O4<double>().apply_z(+1,*u,ix,iy+1,iz);
							fine_flux += Laplace_flux_O4<double>().apply_z(+1,*u,ix+1,iy+1,iz);
							
							coarse_flux = Laplace_flux_O4<double>().apply_z(+1,*utop,ixtop,iytop,iztop)/2.0;
							fine_flux /= 4.0;
							
							dflux = coarse_flux - fine_flux;
							
							for(int i=0;i<2;++i)
								for( int j=0;j<2;++j)
								{
									(*u)(ix+i,iy+j,iz)   += ffac*dflux;
									(*u)(ix+i,iy+j,iz+1) += ffac*dflux;
								}
						}
					}
				}
	}
};

//! 7th order flux-corrected coarse-fine interpolation
struct interp_O7_fluxcorr
: public coarse_fine_interpolation
{
	
	
	template< class G >
	void interp_coarse_fine( unsigned ilevel, G& coarse, G& fine )
	{
		
		//... use cubic interpolation to get all values on boundary
		mg_cubic().prolong_bnd( coarse, fine );
		
		//... use flux corrected quadratic interpolation for the
		//... the boundary overlapped by the Laplace stencil
		G *u    = &fine;
		G *utop = &coarse;
		
		
		int
			xoff = u->offset(0),
			yoff = u->offset(1),
			zoff = u->offset(2);
			
		//... don't do anything if we are not an additional refinement region
		if( xoff == 0 && yoff == 0 && zoff == 0 )
			return;
		
		int
		nx = u->size(0), 
		ny = u->size(1), 
		nz = u->size(2);
		
		//... set boundary condition for fine grid
#pragma omp parallel for schedule(dynamic)
		for( int ix=-1; ix<=nx; ++ix )
			for( int iy=-1; iy<=ny; ++iy )
				for( int iz=-1; iz<=nz; ++iz )
				{
					bool xbnd=(ix==-1||ix==nx),ybnd=(iy==-1||iy==ny),zbnd=(iz==-1||iz==nz);
					bool bnd=xbnd|ybnd|zbnd;
					
					if( bnd )
					{
						
						int ixtop = (int)(0.5*(double)(ix))+xoff;
						int iytop = (int)(0.5*(double)(iy))+yoff;
						int iztop = (int)(0.5*(double)(iz))+zoff;
						
						if( ix==-1 ) ixtop=xoff-1;
						if( iy==-1 ) iytop=yoff-1;
						if( iz==-1 ) iztop=zoff-1;
						
						double ustar[7], uhat[3];			
						double fac = 0.5;
						
						double coarse_flux, fine_flux, dflux;
						
						double ffac = 180./222.;
						
						// left boundary
						if( ix == -1 && iy%2==0 && iz%2==0 )
						{
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{		
									for( int p=0; p<3; ++p )
									{
										for( int q=-3;q<=3;++q )
											ustar[q+3] = interp6( (*utop)(ixtop+p-2,iytop-3,iztop+q), (*utop)(ixtop+p-2,iytop-2,iztop+q), 
																 (*utop)(ixtop+p-2,iytop-1,iztop+q), (*utop)(ixtop+p-2,iytop,iztop+q),   
																 (*utop)(ixtop+p-2,iytop+1,iztop+q), (*utop)(ixtop+p-2,iytop+2,iztop+q), 
																 (*utop)(ixtop+p-2,iytop+3,iztop+q), fac*((double)j-0.5) );
										uhat[p] = interp6( ustar, fac*((double)k-0.5));//-1.5 );
									}
									
									(*u)(ix,iy+j,iz+k)   = interp6left( uhat[0], uhat[1], uhat[2], (*u)(ix+1,iy+j,iz+k), 
																	   (*u)(ix+2,iy+j,iz+k), (*u)(ix+3,iy+j,iz+k), (*u)(ix+4,iy+j,iz+k) );
									(*u)(ix-1,iy+j,iz+k) = interp6lleft( uhat[0], uhat[1], uhat[2], (*u)(ix+1,iy+j,iz+k), 
																		(*u)(ix+2,iy+j,iz+k), (*u)(ix+3,iy+j,iz+k),(*u)(ix+4,iy+j,iz+k) );
									(*u)(ix-2,iy+j,iz+k) = interp6llleft( uhat[0], uhat[1], uhat[2], (*u)(ix+1,iy+j,iz+k), 
																		(*u)(ix+2,iy+j,iz+k), (*u)(ix+3,iy+j,iz+k),(*u)(ix+4,iy+j,iz+k) );
								}
							
							fine_flux = 0.0;
							fine_flux += Laplace_flux_O6<double>().apply_x(-1,*u,ix+1,iy,iz);
							fine_flux += Laplace_flux_O6<double>().apply_x(-1,*u,ix+1,iy+1,iz);
							fine_flux += Laplace_flux_O6<double>().apply_x(-1,*u,ix+1,iy,iz+1);
							fine_flux += Laplace_flux_O6<double>().apply_x(-1,*u,ix+1,iy+1,iz+1);
							fine_flux /= 4.0;
							
							coarse_flux = Laplace_flux_O6<double>().apply_x(-1,*utop,ixtop+1,iytop,iztop)/2.0;
							
							dflux = coarse_flux - fine_flux;
							
							for(int j=0;j<2;++j)
								for( int k=0;k<2;++k)
								{
									(*u)(ix,iy+j,iz+k)   += ffac*dflux;
									(*u)(ix-1,iy+j,iz+k) += ffac*dflux;
									(*u)(ix-2,iy+j,iz+k) += ffac*dflux;
								}
							
							
						}
						// right boundary
						if( ix == nx && iy%2==0 && iz%2==0 )
						{
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{		
									for( int p=0; p<3; ++p )
									{
										for( int q=-3;q<=3;++q )
											ustar[q+3] = interp6( (*utop)(ixtop+p,iytop-3,iztop+q), (*utop)(ixtop+p,iytop-2,iztop+q), 
																 (*utop)(ixtop+p,iytop-1,iztop+q), (*utop)(ixtop+p,iytop,iztop+q),
																 (*utop)(ixtop+p,iytop+1,iztop+q), (*utop)(ixtop+p,iytop+2,iztop+q),
																 (*utop)(ixtop+p,iytop+3,iztop+q), fac*((double)j-0.5) );
										uhat[p] = interp6( ustar, fac*((double)k-0.5) );//-1.5 );
									}
									
									(*u)(ix,iy+j,iz+k)   = interp6right( (*u)(ix-4,iy+j,iz+k), (*u)(ix-3,iy+j,iz+k), (*u)(ix-2,iy+j,iz+k), 
																		(*u)(ix-1,iy+j,iz+k), uhat[0], uhat[1], uhat[2] );
									(*u)(ix+1,iy+j,iz+k)   = interp6rright( (*u)(ix-4,iy+j,iz+k), (*u)(ix-3,iy+j,iz+k), (*u)(ix-2,iy+j,iz+k), 
																		(*u)(ix-1,iy+j,iz+k), uhat[0], uhat[1], uhat[2] );
									(*u)(ix+2,iy+j,iz+k)   = interp6rrright( (*u)(ix-4,iy+j,iz+k), (*u)(ix-3,iy+j,iz+k), (*u)(ix-2,iy+j,iz+k), 
																		(*u)(ix-1,iy+j,iz+k), uhat[0], uhat[1], uhat[2] );

									
								}
							
							fine_flux = 0.0;
							fine_flux += Laplace_flux_O6<double>().apply_x(+1,*u,ix,iy,iz);
							fine_flux += Laplace_flux_O6<double>().apply_x(+1,*u,ix,iy+1,iz);
							fine_flux += Laplace_flux_O6<double>().apply_x(+1,*u,ix,iy,iz+1);
							fine_flux += Laplace_flux_O6<double>().apply_x(+1,*u,ix,iy+1,iz+1);
							
							coarse_flux = Laplace_flux_O6<double>().apply_x(+1,*utop,ixtop,iytop,iztop)/2.0;
							fine_flux /= 4.0;
							
							dflux = coarse_flux - fine_flux;
							
							for(int j=0;j<2;++j)
								for( int k=0;k<2;++k)
								{
									(*u)(ix,iy+j,iz+k)   += ffac*dflux;
									(*u)(ix+1,iy+j,iz+k) += ffac*dflux;
									(*u)(ix+2,iy+j,iz+k) += ffac*dflux;
								}
							
						}
						// bottom boundary
						if( iy == -1 && ix%2==0 && iz%2==0 )
						{
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{	
									for( int p=0; p<3; ++p )
									{
										for( int q=-3;q<=3;++q )
											ustar[q+3] = interp6( (*utop)(ixtop-3,iytop+p-2,iztop+q), (*utop)(ixtop-2,iytop+p-2,iztop+q),
																 (*utop)(ixtop-1,iytop+p-2,iztop+q), (*utop)(ixtop,iytop+p-2,iztop+q),   
																 (*utop)(ixtop+1,iytop+p-2,iztop+q), (*utop)(ixtop+2,iytop+p-2,iztop+q),
																 (*utop)(ixtop+3,iytop+p-2,iztop+q), fac*((double)j-0.5) );
										uhat[p] = interp6( ustar, fac*((double)k-0.5));//-1.5 );
									}
									
									(*u)(ix+j,iy,iz+k)   = interp6left( uhat[0], uhat[1], uhat[2], (*u)(ix+j,iy+1,iz+k), 
																	   (*u)(ix+j,iy+2,iz+k), (*u)(ix+j,iy+3,iz+k),(*u)(ix+j,iy+4,iz+k) );									
									(*u)(ix+j,iy-1,iz+k)   = interp6lleft( uhat[0], uhat[1], uhat[2], (*u)(ix+j,iy+1,iz+k), 
																	   (*u)(ix+j,iy+2,iz+k), (*u)(ix+j,iy+3,iz+k),(*u)(ix+j,iy+4,iz+k) );									
									(*u)(ix+j,iy-2,iz+k)   = interp6llleft( uhat[0], uhat[1], uhat[2], (*u)(ix+j,iy+1,iz+k), 
																	   (*u)(ix+j,iy+2,iz+k), (*u)(ix+j,iy+3,iz+k),(*u)(ix+j,iy+4,iz+k) );									

								}
							
							fine_flux = 0.0;
							fine_flux += Laplace_flux_O6<double>().apply_y(-1,*u,ix,iy+1,iz);
							fine_flux += Laplace_flux_O6<double>().apply_y(-1,*u,ix+1,iy+1,iz);
							fine_flux += Laplace_flux_O6<double>().apply_y(-1,*u,ix,iy+1,iz+1);
							fine_flux += Laplace_flux_O6<double>().apply_y(-1,*u,ix+1,iy+1,iz+1);
							
							coarse_flux = Laplace_flux_O6<double>().apply_y(-1,*utop,ixtop,iytop+1,iztop)/2.0;
							fine_flux /= 4.0;
							
							dflux = coarse_flux - fine_flux;
							
							for(int i=0;i<2;++i)
								for( int k=0;k<2;++k)
								{
									(*u)(ix+i,iy,iz+k)   += ffac*dflux;
									(*u)(ix+i,iy-1,iz+k) += ffac*dflux;
									(*u)(ix+i,iy-2,iz+k) += ffac*dflux;
								}
							
						}
						// top boundary
						if( iy == ny && ix%2==0 && iz%2==0 )
						{
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{	
									for( int p=0; p<3; ++p )
									{
										for( int q=-3;q<=3;++q )
											ustar[q+3] = interp6( (*utop)(ixtop-3,iytop+p,iztop+q),  (*utop)(ixtop-2,iytop+p,iztop+q), 
																 (*utop)(ixtop-1,iytop+p,iztop+q), (*utop)(ixtop,iytop+p,iztop+q), 
																 (*utop)(ixtop+1,iytop+p,iztop+q), (*utop)(ixtop+2,iytop+p,iztop+q),
																  (*utop)(ixtop+3,iytop+p,iztop+q), fac*((double)j-0.5) );
										uhat[p] = interp6( ustar, fac*((double)k-0.5));//+1.5 );
									}
									
									(*u)(ix+j,iy,iz+k)   = interp6right( (*u)(ix+j,iy-4,iz+k), (*u)(ix+j,iy-3,iz+k), (*u)(ix+j,iy-2,iz+k), 
																		(*u)(ix+j,iy-1,iz+k), uhat[0], uhat[1], uhat[2] );									
									(*u)(ix+j,iy+1,iz+k)   = interp6rright( (*u)(ix+j,iy-4,iz+k), (*u)(ix+j,iy-3,iz+k), (*u)(ix+j,iy-2,iz+k), 
																		(*u)(ix+j,iy-1,iz+k), uhat[0], uhat[1], uhat[2] );									
									(*u)(ix+j,iy+2,iz+k)   = interp6rrright( (*u)(ix+j,iy-4,iz+k), (*u)(ix+j,iy-3,iz+k), (*u)(ix+j,iy-2,iz+k), 
																		(*u)(ix+j,iy-1,iz+k), uhat[0], uhat[1], uhat[2] );									

								}
							
							fine_flux = 0.0;
							fine_flux += Laplace_flux_O6<double>().apply_y(+1,*u,ix,iy,iz);
							fine_flux += Laplace_flux_O6<double>().apply_y(+1,*u,ix+1,iy,iz);
							fine_flux += Laplace_flux_O6<double>().apply_y(+1,*u,ix,iy,iz+1);
							fine_flux += Laplace_flux_O6<double>().apply_y(+1,*u,ix+1,iy,iz+1);
							
							coarse_flux = Laplace_flux_O6<double>().apply_y(+1,*utop,ixtop,iytop,iztop)/2.0;
							fine_flux /= 4.0;
							
							dflux = coarse_flux - fine_flux;
							
							for(int i=0;i<2;++i)
								for( int k=0;k<2;++k)
								{
									(*u)(ix+i,iy,iz+k)   += ffac*dflux;
									(*u)(ix+i,iy+1,iz+k) += ffac*dflux;
									(*u)(ix+i,iy+2,iz+k) += ffac*dflux;
								}
							
							
						}
						// front boundary
						if( iz == -1 && ix%2==0 && iy%2==0 )
						{
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{		
									for( int p=0; p<3; ++p )
									{
										for( int q=-3;q<=3;++q )
											ustar[q+3] = interp6( (*utop)(ixtop-3,iytop+q,iztop+p-2), (*utop)(ixtop-2,iytop+q,iztop+p-2),
																 (*utop)(ixtop-1,iytop+q,iztop+p-2), (*utop)(ixtop,iytop+q,iztop+p-2), 
																 (*utop)(ixtop+1,iytop+q,iztop+p-2), (*utop)(ixtop+2,iytop+q,iztop+p-2),
																 (*utop)(ixtop+3,iytop+q,iztop+p-2), fac*((double)j-0.5) );
										uhat[p] = interp6( ustar, fac*((double)k-0.5));//-1.5 );
									}
									
									(*u)(ix+j,iy+k,iz)   = interp6left( uhat[0], uhat[1], uhat[2], (*u)(ix+j,iy+k,iz+1), 
																	   (*u)(ix+j,iy+k,iz+2), (*u)(ix+j,iy+k,iz+3),(*u)(ix+j,iy+k,iz+4) );									
									(*u)(ix+j,iy+k,iz-1)   = interp6lleft( uhat[0], uhat[1], uhat[2], (*u)(ix+j,iy+k,iz+1), 
																	   (*u)(ix+j,iy+k,iz+2), (*u)(ix+j,iy+k,iz+3), (*u)(ix+j,iy+k,iz+4) );									
									(*u)(ix+j,iy+k,iz-2)   = interp6llleft( uhat[0], uhat[1], uhat[2], (*u)(ix+j,iy+k,iz+1), 
																	   (*u)(ix+j,iy+k,iz+2), (*u)(ix+j,iy+k,iz+3), (*u)(ix+j,iy+k,iz+4) );									
								}
							
							
							fine_flux = 0.0;
							fine_flux += Laplace_flux_O6<double>().apply_z(-1,*u,ix,iy,iz+1);
							fine_flux += Laplace_flux_O6<double>().apply_z(-1,*u,ix+1,iy,iz+1);
							fine_flux += Laplace_flux_O6<double>().apply_z(-1,*u,ix,iy+1,iz+1);
							fine_flux += Laplace_flux_O6<double>().apply_z(-1,*u,ix+1,iy+1,iz+1);
							
							coarse_flux = Laplace_flux_O6<double>().apply_z(-1,*utop,ixtop,iytop,iztop+1)/2.0;
							fine_flux /= 4.0;
							
							dflux = coarse_flux - fine_flux;
							
							for(int i=0;i<2;++i)
								for( int j=0;j<2;++j)
								{
									(*u)(ix+i,iy+j,iz)   += ffac*dflux;
									(*u)(ix+i,iy+j,iz-1) += ffac*dflux;
									(*u)(ix+i,iy+j,iz-2) += ffac*dflux;
								}
							
						}
						// back boundary
						if( iz == nz && ix%2==0 && iy%2==0 )
						{
							for( int j=0;j<=1;j++)
								for( int k=0;k<=1;k++)
								{		
									for( int p=0; p<3; ++p )
									{
										for( int q=-3;q<=3;++q )
											ustar[q+3] = interp6( (*utop)(ixtop-3,iytop+q,iztop+p), (*utop)(ixtop-2,iytop+q,iztop+p), 
																 (*utop)(ixtop-1,iytop+q,iztop+p), (*utop)(ixtop,iytop+q,iztop+p),   
																 (*utop)(ixtop+1,iytop+q,iztop+p), (*utop)(ixtop+2,iytop+q,iztop+p), 
																 (*utop)(ixtop+3,iytop+q,iztop+p), fac*((double)j-0.5) );
										uhat[p] = interp6( ustar, fac*((double)k-0.5));//+1.5 );
									}
									
									(*u)(ix+j,iy+k,iz)   = interp6right( (*u)(ix+j,iy+k,iz-4), (*u)(ix+j,iy+k,iz-3), (*u)(ix+j,iy+k,iz-2), 
																		(*u)(ix+j,iy+k,iz-1), uhat[0], uhat[1], uhat[2] );
									(*u)(ix+j,iy+k,iz+1)   = interp6rright( (*u)(ix+j,iy+k,iz-4), (*u)(ix+j,iy+k,iz-3), (*u)(ix+j,iy+k,iz-2), 
																		(*u)(ix+j,iy+k,iz-1), uhat[0], uhat[1], uhat[2] );
									(*u)(ix+j,iy+k,iz+2)   = interp6rrright( (*u)(ix+j,iy+k,iz-4), (*u)(ix+j,iy+k,iz-3), (*u)(ix+j,iy+k,iz-2), 
																		(*u)(ix+j,iy+k,iz-1), uhat[0], uhat[1], uhat[2] );

								}
							
							fine_flux = 0.0;
							fine_flux += Laplace_flux_O6<double>().apply_z(+1,*u,ix,iy,iz);
							fine_flux += Laplace_flux_O6<double>().apply_z(+1,*u,ix+1,iy,iz);
							fine_flux += Laplace_flux_O6<double>().apply_z(+1,*u,ix,iy+1,iz);
							fine_flux += Laplace_flux_O6<double>().apply_z(+1,*u,ix+1,iy+1,iz);
							
							coarse_flux = Laplace_flux_O6<double>().apply_z(+1,*utop,ixtop,iytop,iztop)/2.0;
							fine_flux /= 4.0;
							
							dflux = coarse_flux - fine_flux;
							
							for(int i=0;i<2;++i)
								for( int j=0;j<2;++j)
								{
									(*u)(ix+i,iy+j,iz)   += ffac*dflux;
									(*u)(ix+i,iy+j,iz+1) += ffac*dflux;
									(*u)(ix+i,iy+j,iz+2) += ffac*dflux;
								}
						}
					}
				}
	}
};


#endif // __MG_INTERP_HH
