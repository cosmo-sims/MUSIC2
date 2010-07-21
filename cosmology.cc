/*
 
 cosmology.cc - This file is part of MUSIC -
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

#include "cosmology.hh"
#include "mesh.hh"
#include "mg_operators.hh"

#define ACC(i,j,k) ((*u.get_grid((ilevel)))((i),(j),(k)))
#define SQR(x)	((x)*(x))


void compute_LLA_density( const grid_hierarchy& u, grid_hierarchy& fnew, unsigned order )
{
	fnew = u;
	
	for( unsigned ilevel=u.levelmin(); ilevel<=u.levelmax(); ++ilevel )
	{
		double h = pow(2.0,ilevel), h2 = h*h, h2_4 = 0.25*h2;
		meshvar_bnd *pvar = fnew.get_grid(ilevel);
		
		
		if( order == 2 )
		{
			#pragma omp parallel for //reduction(+:sum_corr,sum,sum2)
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
					{
						double D[3][3];
						
						D[0][0] = (ACC(ix-1,iy,iz)-2.0*ACC(ix,iy,iz)+ACC(ix+1,iy,iz)) * h2;
						D[1][1] = (ACC(ix,iy-1,iz)-2.0*ACC(ix,iy,iz)+ACC(ix,iy+1,iz)) * h2;
						D[2][2] = (ACC(ix,iy,iz-1)-2.0*ACC(ix,iy,iz)+ACC(ix,iy,iz+1)) * h2;
											
						D[0][1] = D[1][0] = (ACC(ix-1,iy-1,iz)-ACC(ix-1,iy+1,iz)-ACC(ix+1,iy-1,iz)+ACC(ix+1,iy+1,iz))*h2_4;
						D[0][2] = D[2][0] = (ACC(ix-1,iy,iz-1)-ACC(ix-1,iy,iz+1)-ACC(ix+1,iy,iz-1)+ACC(ix+1,iy,iz+1))*h2_4;
						D[1][2] = D[2][1] = (ACC(ix,iy-1,iz-1)-ACC(ix,iy-1,iz+1)-ACC(ix,iy+1,iz-1)+ACC(ix,iy+1,iz+1))*h2_4;
						
						(*pvar)(ix,iy,iz) = D[0][0]+D[1][1]+D[2][2] +
							( D[0][0]*D[1][1] + D[0][0]*D[2][2] + D[1][1]*D[2][2] +
							  D[0][1]*D[1][0] + D[0][2]*D[2][0] + D[1][2]*D[2][1] +
							  D[0][0]*D[0][0] + D[1][1]*D[1][1] + D[2][2]*D[2][2] );
						
					}
		}
		else if ( order == 4 )
		{
			#pragma omp parallel for 
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
					{
						double D[3][3];
						
						D[0][0] = (-ACC(ix-2,iy,iz)+16.*ACC(ix-1,iy,iz)-30.0*ACC(ix,iy,iz)+16.*ACC(ix+1,iy,iz)-ACC(ix+2,iy,iz)) * h2/12.0;
						D[1][1] = (-ACC(ix,iy-2,iz)+16.*ACC(ix,iy-1,iz)-30.0*ACC(ix,iy,iz)+16.*ACC(ix,iy+1,iz)-ACC(ix,iy+2,iz)) * h2/12.0;
						D[2][2] = (-ACC(ix,iy,iz-2)+16.*ACC(ix,iy,iz-1)-30.0*ACC(ix,iy,iz)+16.*ACC(ix,iy,iz+1)-ACC(ix,iy,iz+2)) * h2/12.0;
						
						D[0][1] = D[1][0] = (ACC(ix-1,iy-1,iz)-ACC(ix-1,iy+1,iz)-ACC(ix+1,iy-1,iz)+ACC(ix+1,iy+1,iz))*h2_4;
						D[0][2] = D[2][0] = (ACC(ix-1,iy,iz-1)-ACC(ix-1,iy,iz+1)-ACC(ix+1,iy,iz-1)+ACC(ix+1,iy,iz+1))*h2_4;
						D[1][2] = D[2][1] = (ACC(ix,iy-1,iz-1)-ACC(ix,iy-1,iz+1)-ACC(ix,iy+1,iz-1)+ACC(ix,iy+1,iz+1))*h2_4;
						
						(*pvar)(ix,iy,iz) = D[0][0]+D[1][1]+D[2][2] +
						( D[0][0]*D[1][1] + D[0][0]*D[2][2] + D[1][1]*D[2][2] +
						 D[0][1]*D[1][0] + D[0][2]*D[2][0] + D[1][2]*D[2][1] +
						 D[0][0]*D[0][0] + D[1][1]*D[1][1] + D[2][2]*D[2][2] );
						
					}
		}
		else if ( order == 6 )
		{
			h2_4/=36.;
			h2/=180.;
			#pragma omp parallel for 
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
					{
						double D[3][3];
						
						D[0][0] = (2.*ACC(ix-3,iy,iz)-27.*ACC(ix-2,iy,iz)+270.*ACC(ix-1,iy,iz)-490.0*ACC(ix,iy,iz)+270.*ACC(ix+1,iy,iz)-27.*ACC(ix+2,iy,iz)+2.*ACC(ix+3,iy,iz)) * h2;
						D[1][1] = (2.*ACC(ix,iy-3,iz)-27.*ACC(ix,iy-2,iz)+270.*ACC(ix,iy-1,iz)-490.0*ACC(ix,iy,iz)+270.*ACC(ix,iy+1,iz)-27.*ACC(ix,iy+2,iz)+2.*ACC(ix,iy+3,iz)) * h2;
						D[2][2] = (2.*ACC(ix,iy,iz-3)-27.*ACC(ix,iy,iz-2)+270.*ACC(ix,iy,iz-1)-490.0*ACC(ix,iy,iz)+270.*ACC(ix,iy,iz+1)-27.*ACC(ix,iy,iz+2)+2.*ACC(ix,iy,iz+3)) * h2;
						
						//.. this is actually 8th order accurate
						D[0][1] = D[1][0] = (64.*(ACC(ix-1,iy-1,iz)-ACC(ix-1,iy+1,iz)-ACC(ix+1,iy-1,iz)+ACC(ix+1,iy+1,iz))
											 -8.*(ACC(ix-2,iy-1,iz)-ACC(ix+2,iy-1,iz)-ACC(ix-2,iy+1,iz)+ACC(ix+2,iy+1,iz)
												+ ACC(ix-1,iy-2,iz)-ACC(ix-1,iy+2,iz)-ACC(ix+1,iy-2,iz)+ACC(ix+1,iy+2,iz))
											 +1.*(ACC(ix-2,iy-2,iz)-ACC(ix-2,iy+2,iz)-ACC(ix+2,iy-2,iz)+ACC(ix+2,iy+2,iz)))*h2_4;
						D[0][2] = D[2][0] = (64.*(ACC(ix-1,iy,iz-1)-ACC(ix-1,iy,iz+1)-ACC(ix+1,iy,iz-1)+ACC(ix+1,iy,iz+1))
											 -8.*(ACC(ix-2,iy,iz-1)-ACC(ix+2,iy,iz-1)-ACC(ix-2,iy,iz+1)+ACC(ix+2,iy,iz+1)
												+ ACC(ix-1,iy,iz-2)-ACC(ix-1,iy,iz+2)-ACC(ix+1,iy,iz-2)+ACC(ix+1,iy,iz+2))
											 +1.*(ACC(ix-2,iy,iz-2)-ACC(ix-2,iy,iz+2)-ACC(ix+2,iy,iz-2)+ACC(ix+2,iy,iz+2)))*h2_4;
						D[1][2] = D[2][1] = (64.*(ACC(ix,iy-1,iz-1)-ACC(ix,iy-1,iz+1)-ACC(ix,iy+1,iz-1)+ACC(ix,iy+1,iz+1))
											 -8.*(ACC(ix,iy-2,iz-1)-ACC(ix,iy+2,iz-1)-ACC(ix,iy-2,iz+1)+ACC(ix,iy+2,iz+1)
												+ ACC(ix,iy-1,iz-2)-ACC(ix,iy-1,iz+2)-ACC(ix,iy+1,iz-2)+ACC(ix,iy+1,iz+2))
											 +1.*(ACC(ix,iy-2,iz-2)-ACC(ix,iy-2,iz+2)-ACC(ix,iy+2,iz-2)+ACC(ix,iy+2,iz+2)))*h2_4;
						
						(*pvar)(ix,iy,iz) = D[0][0]+D[1][1]+D[2][2] +
						( D[0][0]*D[1][1] + D[0][0]*D[2][2] + D[1][1]*D[2][2] +
						 D[0][1]*D[1][0] + D[0][2]*D[2][0] + D[1][2]*D[2][1] +
						 D[0][0]*D[0][0] + D[1][1]*D[1][1] + D[2][2]*D[2][2] );
						
					}
			//TODO: test sixth order
			
		}else
			throw std::runtime_error("compute_LLA_density : invalid operator order specified");

	}
	
}


void compute_Lu_density( const grid_hierarchy& u, grid_hierarchy& fnew )
{
	fnew = u;
	
	for( unsigned ilevel=u.levelmin(); ilevel<=u.levelmax(); ++ilevel )
	{
		double h = pow(2.0,ilevel), h2 = h*h;
		meshvar_bnd *pvar = fnew.get_grid(ilevel);
		
		#pragma omp parallel for
		for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
			for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
				for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
				{
					double D[3][3];
					
					D[0][0] = 1.0 + (ACC(ix-1,iy,iz)-2.0*ACC(ix,iy,iz)+ACC(ix+1,iy,iz)) * h2;
					D[1][1] = 1.0 + (ACC(ix,iy-1,iz)-2.0*ACC(ix,iy,iz)+ACC(ix,iy+1,iz)) * h2;
					D[2][2] = 1.0 + (ACC(ix,iy,iz-1)-2.0*ACC(ix,iy,iz)+ACC(ix,iy,iz+1)) * h2;
					
					(*pvar)(ix,iy,iz) = D[0][0]+D[1][1]+D[2][2] - 3.0;
					
				}
	}
	
}


void compute_2LPT_source( const grid_hierarchy& u, grid_hierarchy& fnew, unsigned order )
{
	fnew = u;
	
	for( unsigned ilevel=u.levelmin(); ilevel<=u.levelmax(); ++ilevel )
	{
		double h = pow(2.0,ilevel), h2 = h*h, h2_4 = 0.25*h2;
		meshvar_bnd *pvar = fnew.get_grid(ilevel);
		
		if ( order == 2 )
		{
			
			#pragma omp parallel for
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
					{
						double D[3][3];
						
						D[0][0] = (ACC(ix-1,iy,iz)-2.0*ACC(ix,iy,iz)+ACC(ix+1,iy,iz)) * h2;
						D[1][1] = (ACC(ix,iy-1,iz)-2.0*ACC(ix,iy,iz)+ACC(ix,iy+1,iz)) * h2;
						D[2][2] = (ACC(ix,iy,iz-1)-2.0*ACC(ix,iy,iz)+ACC(ix,iy,iz+1)) * h2;
						D[0][1] = D[1][0] = (ACC(ix-1,iy-1,iz)-ACC(ix-1,iy+1,iz)-ACC(ix+1,iy-1,iz)+ACC(ix+1,iy+1,iz))*h2_4;
						D[0][2] = D[2][0] = (ACC(ix-1,iy,iz-1)-ACC(ix-1,iy,iz+1)-ACC(ix+1,iy,iz-1)+ACC(ix+1,iy,iz+1))*h2_4;
						D[1][2] = D[2][1] = (ACC(ix,iy-1,iz-1)-ACC(ix,iy-1,iz+1)-ACC(ix,iy+1,iz-1)+ACC(ix,iy+1,iz+1))*h2_4;
						
						
						(*pvar)(ix,iy,iz) =   D[0][0]*D[1][1] - SQR( D[0][1] )
											+ D[0][0]*D[2][2] - SQR( D[0][2] )
											+ D[1][1]*D[2][2] - SQR( D[1][2] );
						
					}
		}
		else if ( order == 4 )
		{
			#pragma omp parallel for
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
					{
						double D[3][3];
						
						D[0][0] = (-ACC(ix-2,iy,iz)+16.*ACC(ix-1,iy,iz)-30.0*ACC(ix,iy,iz)+16.*ACC(ix+1,iy,iz)-ACC(ix+2,iy,iz)) * h2/12.0;
						D[1][1] = (-ACC(ix,iy-2,iz)+16.*ACC(ix,iy-1,iz)-30.0*ACC(ix,iy,iz)+16.*ACC(ix,iy+1,iz)-ACC(ix,iy+2,iz)) * h2/12.0;
						D[2][2] = (-ACC(ix,iy,iz-2)+16.*ACC(ix,iy,iz-1)-30.0*ACC(ix,iy,iz)+16.*ACC(ix,iy,iz+1)-ACC(ix,iy,iz+2)) * h2/12.0;
						D[0][1] = D[1][0] = (ACC(ix-1,iy-1,iz)-ACC(ix-1,iy+1,iz)-ACC(ix+1,iy-1,iz)+ACC(ix+1,iy+1,iz))*h2_4;
						D[0][2] = D[2][0] = (ACC(ix-1,iy,iz-1)-ACC(ix-1,iy,iz+1)-ACC(ix+1,iy,iz-1)+ACC(ix+1,iy,iz+1))*h2_4;
						D[1][2] = D[2][1] = (ACC(ix,iy-1,iz-1)-ACC(ix,iy-1,iz+1)-ACC(ix,iy+1,iz-1)+ACC(ix,iy+1,iz+1))*h2_4;
						
						
						(*pvar)(ix,iy,iz) =   D[0][0]*D[1][1] - SQR( D[0][1] )
						+ D[0][0]*D[2][2] - SQR( D[0][2] )
						+ D[1][1]*D[2][2] - SQR( D[1][2] );
						
					}
		}
		else if ( order == 6 )
		{
			h2_4/=36.;
			h2/=180.;
#pragma omp parallel for 
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
					{
						double D[3][3];
						
						D[0][0] = (2.*ACC(ix-3,iy,iz)-27.*ACC(ix-2,iy,iz)+270.*ACC(ix-1,iy,iz)-490.0*ACC(ix,iy,iz)+270.*ACC(ix+1,iy,iz)-27.*ACC(ix+2,iy,iz)+2.*ACC(ix+3,iy,iz)) * h2;
						D[1][1] = (2.*ACC(ix,iy-3,iz)-27.*ACC(ix,iy-2,iz)+270.*ACC(ix,iy-1,iz)-490.0*ACC(ix,iy,iz)+270.*ACC(ix,iy+1,iz)-27.*ACC(ix,iy+2,iz)+2.*ACC(ix,iy+3,iz)) * h2;
						D[2][2] = (2.*ACC(ix,iy,iz-3)-27.*ACC(ix,iy,iz-2)+270.*ACC(ix,iy,iz-1)-490.0*ACC(ix,iy,iz)+270.*ACC(ix,iy,iz+1)-27.*ACC(ix,iy,iz+2)+2.*ACC(ix,iy,iz+3)) * h2;
						
						//.. this is actually 8th order accurate
						D[0][1] = D[1][0] = (64.*(ACC(ix-1,iy-1,iz)-ACC(ix-1,iy+1,iz)-ACC(ix+1,iy-1,iz)+ACC(ix+1,iy+1,iz))
											 -8.*(ACC(ix-2,iy-1,iz)-ACC(ix+2,iy-1,iz)-ACC(ix-2,iy+1,iz)+ACC(ix+2,iy+1,iz)
												  + ACC(ix-1,iy-2,iz)-ACC(ix-1,iy+2,iz)-ACC(ix+1,iy-2,iz)+ACC(ix+1,iy+2,iz))
											 +1.*(ACC(ix-2,iy-2,iz)-ACC(ix-2,iy+2,iz)-ACC(ix+2,iy-2,iz)+ACC(ix+2,iy+2,iz)))*h2_4;
						D[0][2] = D[2][0] = (64.*(ACC(ix-1,iy,iz-1)-ACC(ix-1,iy,iz+1)-ACC(ix+1,iy,iz-1)+ACC(ix+1,iy,iz+1))
											 -8.*(ACC(ix-2,iy,iz-1)-ACC(ix+2,iy,iz-1)-ACC(ix-2,iy,iz+1)+ACC(ix+2,iy,iz+1)
												  + ACC(ix-1,iy,iz-2)-ACC(ix-1,iy,iz+2)-ACC(ix+1,iy,iz-2)+ACC(ix+1,iy,iz+2))
											 +1.*(ACC(ix-2,iy,iz-2)-ACC(ix-2,iy,iz+2)-ACC(ix+2,iy,iz-2)+ACC(ix+2,iy,iz+2)))*h2_4;
						D[1][2] = D[2][1] = (64.*(ACC(ix,iy-1,iz-1)-ACC(ix,iy-1,iz+1)-ACC(ix,iy+1,iz-1)+ACC(ix,iy+1,iz+1))
											 -8.*(ACC(ix,iy-2,iz-1)-ACC(ix,iy+2,iz-1)-ACC(ix,iy-2,iz+1)+ACC(ix,iy+2,iz+1)
												  + ACC(ix,iy-1,iz-2)-ACC(ix,iy-1,iz+2)-ACC(ix,iy+1,iz-2)+ACC(ix,iy+1,iz+2))
											 +1.*(ACC(ix,iy-2,iz-2)-ACC(ix,iy-2,iz+2)-ACC(ix,iy+2,iz-2)+ACC(ix,iy+2,iz+2)))*h2_4;
						
						(*pvar)(ix,iy,iz) =   D[0][0]*D[1][1] - SQR( D[0][1] )
											+ D[0][0]*D[2][2] - SQR( D[0][2] )
											+ D[1][1]*D[2][2] - SQR( D[1][2] );
						
					}
			//TODO: test sixth order
			
		}
		else
			throw std::runtime_error("compute_2LPT_source : invalid operator order specified");


	}
	
	
	//.. subtract global mean so the multi-grid poisson solver behaves well
	
	for( int i=fnew.levelmax(); i>(int)fnew.levelmin(); --i )
		mg_straight().restrict( (*fnew.get_grid(i)), (*fnew.get_grid(i-1)) );
	
	double sum = 0.0;
	int nx,ny,nz;
	
	nx = fnew.get_grid(fnew.levelmin())->size(0);
	ny = fnew.get_grid(fnew.levelmin())->size(1);
	nz = fnew.get_grid(fnew.levelmin())->size(2);
	
	for( int ix=0; ix<nx; ++ix )
		for( int iy=0; iy<ny; ++iy )
			for( int iz=0; iz<nz; ++iz )
				sum += (*fnew.get_grid(fnew.levelmin()))(ix,iy,iz);
	
	sum /= (nx*ny*nz);
	
	for( unsigned i=fnew.levelmin(); i<=fnew.levelmax(); ++i )
	{		
		nx = fnew.get_grid(i)->size(0);
		ny = fnew.get_grid(i)->size(1);
		nz = fnew.get_grid(i)->size(2);
		
		for( int ix=0; ix<nx; ++ix )
			for( int iy=0; iy<ny; ++iy )
				for( int iz=0; iz<nz; ++iz )
					(*fnew.get_grid(i))(ix,iy,iz) -= sum;
	}
	
}
#undef SQR
#undef ACC

