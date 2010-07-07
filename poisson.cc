/*
 
 poisson.cc - This file is part of MUSIC -
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

/****** ABSTRACT FACTORY PATTERN IMPLEMENTATION *******/

#include "poisson.hh"
std::map< std::string, poisson_plugin_creator *>& 
get_poisson_plugin_map()
{
	static std::map< std::string, poisson_plugin_creator* > poisson_plugin_map;
	return poisson_plugin_map;
}

void print_poisson_plugins()
{
	std::map< std::string, poisson_plugin_creator *>& m = get_poisson_plugin_map();
	
	std::map< std::string, poisson_plugin_creator *>::iterator it;
	it = m.begin();
	std::cout << " - Available poisson solver plug-ins:\n";
	while( it!=m.end() )
	{
		if( (*it).second )
			std::cout << "\t\'" << (*it).first << "\'\n";
		++it;
	}
	
	
}


/****** CALL IMPLEMENTATIONS OF POISSON SOLVER CLASSES ******/

#include "mg_solver.hh"
#include "fd_schemes.hh"

#ifdef SINGLE_PRECISION
typedef multigrid::solver< stencil_7P<float>, interp_O3_fluxcorr, mg_straight, float > poisson_solver_O2;
typedef multigrid::solver< stencil_13P<float>, interp_O5_fluxcorr, mg_straight, float > poisson_solver_O4;
typedef multigrid::solver< stencil_19P<float>, interp_O7_fluxcorr, mg_straight, float > poisson_solver_O6;
#else
typedef multigrid::solver< stencil_7P<double>, interp_O3_fluxcorr, mg_straight, double > poisson_solver_O2;
typedef multigrid::solver< stencil_13P<double>, interp_O5_fluxcorr, mg_straight, double > poisson_solver_O4;
typedef multigrid::solver< stencil_19P<double>, interp_O7_fluxcorr, mg_straight, double > poisson_solver_O6;
#endif


/**************************************************************************************/
/**************************************************************************************/
#pragma mark -


double multigrid_poisson_plugin::solve( grid_hierarchy& f, grid_hierarchy& u )
{
	unsigned verbosity = cf_.getValueSafe<unsigned>("setup","verbosity",2);
	
	if( verbosity > 0 )
		std::cout << " - Invoking multi-grid Poisson solver..." << std::endl;
	
	double acc = 1e-5, err;
	std::string ps_smoother_name;
	unsigned ps_presmooth, ps_postsmooth, order;
	
	acc                 = cf_.getValueSafe<double>("poisson","accuracy",acc);
	ps_presmooth		= cf_.getValueSafe<unsigned>("poisson","pre_smooth",3);
	ps_postsmooth		= cf_.getValueSafe<unsigned>("poisson","post_smooth",3);
	ps_smoother_name	= cf_.getValueSafe<std::string>("poisson","smoother","gs");
	order				= cf_.getValueSafe<unsigned>( "poisson", "laplace_order", 4 );
	
	multigrid::opt::smtype ps_smtype = multigrid::opt::sm_gauss_seidel;
	
	if ( ps_smoother_name == std::string("gs") )
		ps_smtype = multigrid::opt::sm_gauss_seidel;
	else if ( ps_smoother_name == std::string("jacobi") )
		ps_smtype = multigrid::opt::sm_jacobi;
	else if ( ps_smoother_name == std::string("sor") )
		ps_smtype = multigrid::opt::sm_sor;
	else
		std::cerr << " - Warning: unknown smoother \'" << ps_smoother_name << "\' for multigrid solver!\n"
		<< "            reverting to \'gs\' (Gauss-Seidel)" << std::endl;
	
	
	double tstart, tend;
	tstart = omp_get_wtime();
	
	
	//----- run Poisson solver -----//
	if( order == 2 )
	{
		poisson_solver_O2 ps( f, ps_smtype, ps_presmooth, ps_postsmooth );
		err = ps.solve( u, acc, true );	
	}
	else if( order == 4 )
	{
		poisson_solver_O4 ps( f, ps_smtype, ps_presmooth, ps_postsmooth );
		err = ps.solve( u, acc, true );	
	}
	else if( order == 6 )
	{
		poisson_solver_O6 ps( f, ps_smtype, ps_presmooth, ps_postsmooth );
		err = ps.solve( u, acc, true );	
	}	
	else
	{	
		throw std::runtime_error("Invalid order specified for Laplace operator");
	}
	
	//------------------------------//
	
	tend = omp_get_wtime();
	if( verbosity > 1 )
		std::cout << " - Poisson solver took " << tend-tstart << "s with " << omp_get_max_threads() << " threads." << std::endl;
	
	return err;
}

double multigrid_poisson_plugin::gradient( int dir, grid_hierarchy& u, grid_hierarchy& Du )
{
	Du = u;
	
	unsigned order = cf_.getValueSafe<unsigned>( "poisson", "grad_order", 4 );
	
	switch( order )
	{
		case 2: 
			implementation().gradient_O2( dir, u, Du );
			break;	
		case 4: 
			implementation().gradient_O4( dir, u, Du );
			break;
		case 6: 
			implementation().gradient_O6( dir, u, Du );
			break;
		default:
			throw std::runtime_error("Invalid order specified for gradient operator!");
	}
	
	return 0.0;
}

void multigrid_poisson_plugin::implementation::gradient_O2( int dir, grid_hierarchy& u, grid_hierarchy& Du )
{
	for( unsigned ilevel=u.levelmin(); ilevel<=u.levelmax(); ++ilevel )
	{
		double h = pow(2.0,ilevel);
		meshvar_bnd *pvar = Du.get_grid(ilevel);
		
		if( dir == 0 )
#pragma omp parallel for
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
						(*pvar)(ix,iy,iz) = 0.5*((*u.get_grid(ilevel))(ix+1,iy,iz)-(*u.get_grid(ilevel))(ix-1,iy,iz))*h;
		
		else if( dir == 1 )
#pragma omp parallel for
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
						(*pvar)(ix,iy,iz) = 0.5*((*u.get_grid(ilevel))(ix,iy+1,iz)-(*u.get_grid(ilevel))(ix,iy-1,iz))*h;
		
		else if( dir == 2 )
#pragma omp parallel for
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
						(*pvar)(ix,iy,iz) = 0.5*((*u.get_grid(ilevel))(ix,iy,iz+1)-(*u.get_grid(ilevel))(ix,iy,iz-1))*h;	
	}
}

void multigrid_poisson_plugin::implementation::gradient_O4( int dir, grid_hierarchy& u, grid_hierarchy& Du )
{
	for( unsigned ilevel=u.levelmin(); ilevel<=u.levelmax(); ++ilevel )
	{
		double h = pow(2.0,ilevel);
		meshvar_bnd *pvar = Du.get_grid(ilevel);
		
		h /= 12.0;
		
		if( dir == 0 )
#pragma omp parallel for
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
						(*pvar)(ix,iy,iz) = ((*u.get_grid(ilevel))(ix-2,iy,iz)
											-8.0*(*u.get_grid(ilevel))(ix-1,iy,iz)
											 +8.0*(*u.get_grid(ilevel))(ix+1,iy,iz)
											 -(*u.get_grid(ilevel))(ix+2,iy,iz))*h;
		
		else if( dir == 1 )
#pragma omp parallel for
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
						(*pvar)(ix,iy,iz) = ((*u.get_grid(ilevel))(ix,iy-2,iz)
											 -8.0*(*u.get_grid(ilevel))(ix,iy-1,iz)
											 +8.0*(*u.get_grid(ilevel))(ix,iy+1,iz)
											 -(*u.get_grid(ilevel))(ix,iy+2,iz))*h;
		
		else if( dir == 2 )
#pragma omp parallel for
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
						(*pvar)(ix,iy,iz) = ((*u.get_grid(ilevel))(ix,iy,iz-2)
											 -8.0*(*u.get_grid(ilevel))(ix,iy,iz-1)
											 +8.0*(*u.get_grid(ilevel))(ix,iy,iz+1)
											 -(*u.get_grid(ilevel))(ix,iy,iz+2))*h;
	}		
	
}

void multigrid_poisson_plugin::implementation::gradient_O6( int dir, grid_hierarchy& u, grid_hierarchy& Du )
{
	for( unsigned ilevel=u.levelmin(); ilevel<=u.levelmax(); ++ilevel )
	{
		double h = pow(2.0,ilevel);
		meshvar_bnd *pvar = Du.get_grid(ilevel);
		
		h /= 60.;
		if( dir == 0 )
#pragma omp parallel for
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
						(*pvar)(ix,iy,iz) = 
						(-(*u.get_grid(ilevel))(ix-3,iy,iz)
						 +9.0*(*u.get_grid(ilevel))(ix-2,iy,iz)
						 -45.0*(*u.get_grid(ilevel))(ix-1,iy,iz)
						 +45.0*(*u.get_grid(ilevel))(ix+1,iy,iz)
						 -9.0*(*u.get_grid(ilevel))(ix+2,iy,iz)
						 +(*u.get_grid(ilevel))(ix+3,iy,iz))*h;
		
		else if( dir == 1 )
#pragma omp parallel for
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
						(*pvar)(ix,iy,iz) = 
						(-(*u.get_grid(ilevel))(ix,iy-3,iz)
						 +9.0*(*u.get_grid(ilevel))(ix,iy-2,iz)
						 -45.0*(*u.get_grid(ilevel))(ix,iy-1,iz)
						 +45.0*(*u.get_grid(ilevel))(ix,iy+1,iz)
						 -9.0*(*u.get_grid(ilevel))(ix,iy+2,iz)
						 +(*u.get_grid(ilevel))(ix,iy+3,iz))*h;
		
		else if( dir == 2 )
#pragma omp parallel for
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
						(*pvar)(ix,iy,iz) = 
						(-(*u.get_grid(ilevel))(ix,iy,iz-3)
						 +9.0*(*u.get_grid(ilevel))(ix,iy,iz-2)
						 -45.0*(*u.get_grid(ilevel))(ix,iy,iz-1)
						 +45.0*(*u.get_grid(ilevel))(ix,iy,iz+1)
						 -9.0*(*u.get_grid(ilevel))(ix,iy,iz+2)
						 +(*u.get_grid(ilevel))(ix,iy,iz+3))*h;
	}
		
}
	
/**************************************************************************************/
/**************************************************************************************/
#pragma mark -

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

double fft_poisson_plugin::solve( grid_hierarchy& f, grid_hierarchy& u )
{
	unsigned verbosity = cf_.getValueSafe<unsigned>("setup","verbosity",2);
	
	if( f.levelmin() != f.levelmax() )
		throw std::runtime_error("fft_poisson_plugin::solve : k-space method can only be used in unigrid mode (levelmin=levelmax)");
	
	if( verbosity > 0 )
		std::cout << " - Invoking unigrid FFT Poisson solver..." << std::endl;
	
	int nx,ny,nz,nzp;
	nx = f.get_grid(f.levelmax())->size(0);
	ny = f.get_grid(f.levelmax())->size(1);
	nz = f.get_grid(f.levelmax())->size(2);
	nzp = 2*(nz/2+1);
	
	
	//... copy data ..................................................
	fftw_real *data = new fftw_real[nx*ny*nzp];
	fftw_complex *cdata = reinterpret_cast<fftw_complex*> (data);
	
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )	
			for( int k=0; k<nz; ++k )
			{
				unsigned idx = (i*ny+j)*nzp+k;
				data[idx] = (*f.get_grid(f.levelmax()))(i,j,k);
			}
	
	//... perform FFT and Poisson solve................................
	rfftwnd_plan 
	plan = rfftw3d_create_plan( nx,ny,nz,
							   FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE),
	iplan = rfftw3d_create_plan( nx,ny,nz,
								FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
	
	
#ifndef SINGLETHREAD_FFTW		
	rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), plan, data, NULL );
#else
	rfftwnd_one_real_to_complex( plan, data, NULL );
#endif
	
	double boxlength = cf_.getValue<double>("setup","boxlength");
	double kfac = 2.0*M_PI/boxlength;
	double fac = -1.0/(nx*ny*nz)/boxlength;
	
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )	
			for( int k=0; k<nz/2+1; ++k )
			{
				int ii = i; if(ii>nx/2) ii-=nx;
				int jj = j; if(jj>ny/2) jj-=ny;
				double ki = (double)ii;
				double kj = (double)jj;
				double kk = (double)k;
				
				double kk2 = kfac*kfac*(ki*ki+kj*kj+kk*kk);
				
				unsigned idx = (i*ny+j)*nzp/2+k;
				double re = cdata[idx].re;
				double im = cdata[idx].im;
				
				cdata[idx].re = -re/kk2*fac;
				cdata[idx].im = -im/kk2*fac;
			}
	
	cdata[0].re = 0.0;
	cdata[0].im = 0.0;
	
#ifndef SINGLETHREAD_FFTW		
	rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), iplan, cdata, NULL );
#else
	rfftwnd_one_complex_to_real( iplan, cdata, NULL );
#endif
	
	rfftwnd_destroy_plan(plan);
	rfftwnd_destroy_plan(iplan);
	
	//... copy data ..........................................
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )	
			for( int k=0; k<nz; ++k )
			{
				unsigned idx = (i*ny+j)*nzp+k;
				(*u.get_grid(u.levelmax()))(i,j,k) = data[idx];
			}
	
	delete[] data;
	
	return 0.0;
}


double fft_poisson_plugin::gradient( int dir, grid_hierarchy& u, grid_hierarchy& Du )
{
	
	if( u.levelmin() != u.levelmax() )
		throw std::runtime_error("fft_poisson_plugin::gradient : k-space method can only be used in unigrid mode (levelmin=levelmax)");
	
	Du = u;
	int nx,ny,nz,nzp;
	nx = u.get_grid(u.levelmax())->size(0);
	ny = u.get_grid(u.levelmax())->size(1);
	nz = u.get_grid(u.levelmax())->size(2);
	nzp = 2*(nz/2+1);
	
	//... copy data ..................................................
	fftw_real *data = new fftw_real[nx*ny*nzp];
	fftw_complex *cdata = reinterpret_cast<fftw_complex*> (data);
	
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )	
			for( int k=0; k<nz; ++k )
			{
				unsigned idx = (i*ny+j)*nzp+k;
				data[idx] = (*u.get_grid(u.levelmax()))(i,j,k);
			}
	
	//... perform FFT and Poisson solve................................
	rfftwnd_plan 
	plan = rfftw3d_create_plan( nx,ny,nz,
							   FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE),
	iplan = rfftw3d_create_plan( nx,ny,nz,
								FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
	
	
#ifndef SINGLETHREAD_FFTW		
	rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), plan, data, NULL );
#else
	rfftwnd_one_real_to_complex( plan, data, NULL );
#endif
	
	double fac = -1.0/(nx*ny*nz);
	double boxlength = cf_.getValue<double>("setup","boxlength");
	double kfac = 2.0*M_PI/boxlength;
	
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )	
			for( int k=0; k<nz/2+1; ++k )
			{
				int ii = i; if(ii>nx/2) ii-=nx;
				int jj = j; if(jj>ny/2) jj-=ny;
				double ki = (double)ii;
				double kj = (double)jj;
				double kk = (double)k;
				
				double kdir;
				if( dir == 0 )
					kdir = kfac*ki;
				else if( dir == 1 )
					kdir = kfac*kj;
				else if( dir == 2 )
					kdir = kfac*kk;
				
				unsigned idx = (i*ny+j)*nzp/2+k;
				double re = cdata[idx].re;
				double im = cdata[idx].im;
				
				cdata[idx].re = fac*im*kdir;
				cdata[idx].im = -fac*re*kdir;
			}
	
	cdata[0].re = 0.0;
	cdata[0].im = 0.0;
	
#ifndef SINGLETHREAD_FFTW		
	rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), iplan, cdata, NULL );
#else
	rfftwnd_one_complex_to_real( iplan, cdata, NULL );
#endif
	
	rfftwnd_destroy_plan(plan);
	rfftwnd_destroy_plan(iplan);
	
	//... copy data ..........................................
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )	
			for( int k=0; k<nz; ++k )
			{
				unsigned idx = (i*ny+j)*nzp+k;
				(*Du.get_grid(u.levelmax()))(i,j,k) = data[idx];
			}
	
	delete[] data;
	return 0.0;
}


/**************************************************************************************/
/**************************************************************************************/
#pragma mark -



namespace{
	poisson_plugin_creator_concrete<multigrid_poisson_plugin> multigrid_poisson_creator("mg_poisson");
	poisson_plugin_creator_concrete<fft_poisson_plugin> fft_poisson_creator("fft_poisson");
}
