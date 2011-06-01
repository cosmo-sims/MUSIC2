/*
 
 poisson.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

/****** ABSTRACT FACTORY PATTERN IMPLEMENTATION *******/

#include "poisson.hh"
#include "Numerics.hh"

#if defined(FFTW3) && defined(SINGLE_PRECISION)
#define fftw_complex fftwf_complex
#endif


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
		
		LOGINFO("Poisson plug-in :: %s",std::string((*it).first).c_str());
		
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
	LOGUSER("Initializing multi-grid Poisson solver...");
	
	unsigned verbosity = cf_.getValueSafe<unsigned>("setup","verbosity",2);
	
	if( verbosity > 0 )
	{	
		std::cout << "-------------------------------------------------------------\n";
		std::cout << " - Invoking multi-grid Poisson solver..." << std::endl;
	}
	
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
	{	
		ps_smtype = multigrid::opt::sm_gauss_seidel;
		LOGUSER("Selected Gauss-Seidel multigrid smoother");
	}
	else if ( ps_smoother_name == std::string("jacobi") )
	{	
		ps_smtype = multigrid::opt::sm_jacobi;
		LOGUSER("Selected Jacobi multigrid smoother");	
	}
	else if ( ps_smoother_name == std::string("sor") )
	{	
		ps_smtype = multigrid::opt::sm_sor;
		LOGUSER("Selected SOR multigrid smoother");
	}
	else
	{	
		LOGWARN("Unknown multigrid smoother \'%s\' specified. Reverting to Gauss-Seidel.",ps_smoother_name.c_str());
		std::cerr << " - Warning: unknown smoother \'" << ps_smoother_name << "\' for multigrid solver!\n"
			<< "            reverting to \'gs\' (Gauss-Seidel)" << std::endl;
	}
		
	
	
	double tstart, tend;
	tstart = omp_get_wtime();
	
	
	//----- run Poisson solver -----//
	if( order == 2 )
	{
		LOGUSER("Running multigrid solver with 2nd order Laplacian...");
		poisson_solver_O2 ps( f, ps_smtype, ps_presmooth, ps_postsmooth );
		err = ps.solve( u, acc, true );	
	}
	else if( order == 4 )
	{
		LOGUSER("Running multigrid solver with 4th order Laplacian...");
		poisson_solver_O4 ps( f, ps_smtype, ps_presmooth, ps_postsmooth );
		err = ps.solve( u, acc, true );	
	}
	else if( order == 6 )
	{
		LOGUSER("Running multigrid solver with 6th order Laplacian..");
		poisson_solver_O6 ps( f, ps_smtype, ps_presmooth, ps_postsmooth );
		err = ps.solve( u, acc, true );	
	}	
	else
	{	
		LOGERR("Invalid order specified for Laplace operator");
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
			LOGERR("Invalid order %d specified for gradient operator", order);
			throw std::runtime_error("Invalid order specified for gradient operator!");
	}
	
	return 0.0;
}

double multigrid_poisson_plugin::gradient_add( int dir, grid_hierarchy& u, grid_hierarchy& Du )
{
	//Du = u;
	
	unsigned order = cf_.getValueSafe<unsigned>( "poisson", "grad_order", 4 );
	
	switch( order )
	{
		case 2: 
			implementation().gradient_add_O2( dir, u, Du );
			break;	
		case 4: 
			implementation().gradient_add_O4( dir, u, Du );
			break;
		case 6: 
			implementation().gradient_add_O6( dir, u, Du );
			break;
		default:
			LOGERR("Invalid order %d specified for gradient operator!",order);
			throw std::runtime_error("Invalid order specified for gradient operator!");
	}
	
	return 0.0;
}

void multigrid_poisson_plugin::implementation::gradient_O2( int dir, grid_hierarchy& u, grid_hierarchy& Du )
{
	LOGUSER("Computing a 2nd order finite difference gradient...");
	
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
	
	LOGUSER("Done computing a 2nd order finite difference gradient.");
}

void multigrid_poisson_plugin::implementation::gradient_add_O2( int dir, grid_hierarchy& u, grid_hierarchy& Du )
{
	LOGUSER("Computing a 2nd order finite difference gradient...");
	
	for( unsigned ilevel=u.levelmin(); ilevel<=u.levelmax(); ++ilevel )
	{
		double h = pow(2.0,ilevel);
		meshvar_bnd *pvar = Du.get_grid(ilevel);
		
		if( dir == 0 )
			#pragma omp parallel for
			for( int ix = 0; ix < (int)(*Du.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*Du.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*Du.get_grid(ilevel)).size(2); ++iz )
						(*pvar)(ix,iy,iz) += 0.5*((*u.get_grid(ilevel))(ix+1,iy,iz)-(*u.get_grid(ilevel))(ix-1,iy,iz))*h;
		
		else if( dir == 1 )
			#pragma omp parallel for
			for( int ix = 0; ix < (int)(*Du.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*Du.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*Du.get_grid(ilevel)).size(2); ++iz )
						(*pvar)(ix,iy,iz) += 0.5*((*u.get_grid(ilevel))(ix,iy+1,iz)-(*u.get_grid(ilevel))(ix,iy-1,iz))*h;
		
		else if( dir == 2 )
			#pragma omp parallel for
			for( int ix = 0; ix < (int)(*Du.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*Du.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*Du.get_grid(ilevel)).size(2); ++iz )
						(*pvar)(ix,iy,iz) += 0.5*((*u.get_grid(ilevel))(ix,iy,iz+1)-(*u.get_grid(ilevel))(ix,iy,iz-1))*h;	
	}
	
	LOGUSER("Done computing a 4th order finite difference gradient.");
}

void multigrid_poisson_plugin::implementation::gradient_O4( int dir, grid_hierarchy& u, grid_hierarchy& Du )
{
	LOGUSER("Computing a 4th order finite difference gradient...");
	
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
	
	LOGUSER("Done computing a 4th order finite difference gradient.");
}

void multigrid_poisson_plugin::implementation::gradient_add_O4( int dir, grid_hierarchy& u, grid_hierarchy& Du )
{
	LOGUSER("Computing a 4th order finite difference gradient...");
	
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
						(*pvar)(ix,iy,iz) += ((*u.get_grid(ilevel))(ix-2,iy,iz)
											 -8.0*(*u.get_grid(ilevel))(ix-1,iy,iz)
											 +8.0*(*u.get_grid(ilevel))(ix+1,iy,iz)
											 -(*u.get_grid(ilevel))(ix+2,iy,iz))*h;
		
		else if( dir == 1 )
#pragma omp parallel for
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
						(*pvar)(ix,iy,iz) += ((*u.get_grid(ilevel))(ix,iy-2,iz)
											 -8.0*(*u.get_grid(ilevel))(ix,iy-1,iz)
											 +8.0*(*u.get_grid(ilevel))(ix,iy+1,iz)
											 -(*u.get_grid(ilevel))(ix,iy+2,iz))*h;
		
		else if( dir == 2 )
#pragma omp parallel for
			for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
				for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
					for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
						(*pvar)(ix,iy,iz) += ((*u.get_grid(ilevel))(ix,iy,iz-2)
											 -8.0*(*u.get_grid(ilevel))(ix,iy,iz-1)
											 +8.0*(*u.get_grid(ilevel))(ix,iy,iz+1)
											 -(*u.get_grid(ilevel))(ix,iy,iz+2))*h;
	}		
	
	
	LOGUSER("Done computing a 4th order finite difference gradient.");
}


void multigrid_poisson_plugin::implementation::gradient_O6( int dir, grid_hierarchy& u, grid_hierarchy& Du )
{
	LOGUSER("Computing a 6th order finite difference gradient...");
	
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
		
	LOGUSER("Done computing a 6th order finite difference gradient.");
}
	

void multigrid_poisson_plugin::implementation::gradient_add_O6( int dir, grid_hierarchy& u, grid_hierarchy& Du )
{
	LOGUSER("Computing a 6th order finite difference gradient...");
	
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
						(*pvar)(ix,iy,iz) += 
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
						(*pvar)(ix,iy,iz) += 
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
						(*pvar)(ix,iy,iz) += 
						(-(*u.get_grid(ilevel))(ix,iy,iz-3)
						 +9.0*(*u.get_grid(ilevel))(ix,iy,iz-2)
						 -45.0*(*u.get_grid(ilevel))(ix,iy,iz-1)
						 +45.0*(*u.get_grid(ilevel))(ix,iy,iz+1)
						 -9.0*(*u.get_grid(ilevel))(ix,iy,iz+2)
						 +(*u.get_grid(ilevel))(ix,iy,iz+3))*h;
	}
	
	LOGUSER("Done computing a 6th order finite difference gradient.");
}


/**************************************************************************************/
/**************************************************************************************/
#pragma mark -
#include "general.hh"

double fft_poisson_plugin::solve( grid_hierarchy& f, grid_hierarchy& u )
{
	LOGUSER("Entering k-space Poisson solver...");
	
	unsigned verbosity = cf_.getValueSafe<unsigned>("setup","verbosity",2);
	
	if( f.levelmin() != f.levelmax() )
	{	
		LOGERR("Attempt to run k-space Poisson solver on non unigrid mesh.");
		throw std::runtime_error("fft_poisson_plugin::solve : k-space method can only be used in unigrid mode (levelmin=levelmax)");
	}
	
	if( verbosity > 0 )
	{	
		std::cout << "-------------------------------------------------------------\n";
		std::cout << " - Invoking unigrid FFT Poisson solver..." << std::endl;
	}
	
	int nx,ny,nz,nzp;
	nx = f.get_grid(f.levelmax())->size(0);
	ny = f.get_grid(f.levelmax())->size(1);
	nz = f.get_grid(f.levelmax())->size(2);
	nzp = 2*(nz/2+1);
	
	
	//... copy data ..................................................
	fftw_real *data = new fftw_real[nx*ny*nzp];
	fftw_complex *cdata = reinterpret_cast<fftw_complex*> (data);
	
	#pragma omp parallel for
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )	
			for( int k=0; k<nz; ++k )
			{
				size_t idx = ((size_t)i*ny+(size_t)j)*nzp+(size_t)k;
				data[idx] = (*f.get_grid(f.levelmax()))(i,j,k);
			}
	
	//... perform FFT and Poisson solve................................
	LOGUSER("Performing forward transform.");

#ifdef FFTW3
	#ifdef SINGLE_PRECISION
	fftwf_plan 
		plan  = fftwf_plan_dft_r2c_3d( nx, ny, nz, data, cdata, FFTW_ESTIMATE ),
		iplan = fftwf_plan_dft_c2r_3d( nx, ny, nz, cdata, data, FFTW_ESTIMATE );
	
	fftwf_execute(plan);
	#else
	fftw_plan 
	plan  = fftw_plan_dft_r2c_3d( nx, ny, nz, data, cdata, FFTW_ESTIMATE ),
	iplan = fftw_plan_dft_c2r_3d( nx, ny, nz, cdata, data, FFTW_ESTIMATE );
	
	fftw_execute(plan);
	#endif
	
#else
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
	
#endif
	double kfac = 2.0*M_PI;
	double fac = -1.0/(nx*ny*nz);
	
	#pragma omp parallel for
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
				
				size_t idx = ((size_t)i*ny+(size_t)j)*nzp/2+(size_t)k;
#ifdef FFTW3
				cdata[idx][0] *= -1.0/kk2*fac;
				cdata[idx][1] *= -1.0/kk2*fac;

#else
				cdata[idx].re *= -1.0/kk2*fac;
				cdata[idx].im *= -1.0/kk2*fac;
#endif			
			}

	LOGUSER("Performing backward transform.");
	
#ifdef FFTW3
	cdata[0][0] = 0.0;
	cdata[0][1] = 0.0;

	#ifdef SINGLE_PRECISION
	fftwf_execute(iplan);
	fftwf_destroy_plan(plan);
	fftwf_destroy_plan(iplan);
	#else
	fftw_execute(iplan);
	fftw_destroy_plan(plan);
	fftw_destroy_plan(iplan);
	#endif
#else
	cdata[0].re = 0.0;
	cdata[0].im = 0.0;
	
	#ifndef SINGLETHREAD_FFTW		
	rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), iplan, cdata, NULL );
	#else
	rfftwnd_one_complex_to_real( iplan, cdata, NULL );
	#endif
	
	rfftwnd_destroy_plan(plan);
	rfftwnd_destroy_plan(iplan);
#endif
	
	


	
	//... copy data ..........................................
	#pragma omp parallel for
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )	
			for( int k=0; k<nz; ++k )
			{
				size_t idx = ((size_t)i*ny+(size_t)j)*nzp+(size_t)k;
				(*u.get_grid(u.levelmax()))(i,j,k) = data[idx];
			}
	
	delete[] data;
	
	//... set boundary values ................................
	int nb = u.get_grid(u.levelmax())->m_nbnd;
	for( int iy=-nb; iy<ny+nb; ++iy )
		for( int iz=-nb; iz<nz+nb; ++iz )
		{
			int iiy( (iy+ny)%ny ), iiz( (iz+nz)%nz );
			
			for( int i=-nb; i<0; ++i )
			{
				(*u.get_grid(u.levelmax()))(i,iy,iz) = (*u.get_grid(u.levelmax()))(nx+i,iiy,iiz);
				(*u.get_grid(u.levelmax()))(nx-1-i,iy,iz) = (*u.get_grid(u.levelmax()))(-1-i,iiy,iiz);	
			}
			
		}
		
	for( int ix=-nb; ix<nx+nb; ++ix )
		for( int iz=-nb; iz<nz+nb; ++iz )
		{
			int iix( (ix+nx)%nx ), iiz( (iz+nz)%nz );
			
			for( int i=-nb; i<0; ++i )
			{
				(*u.get_grid(u.levelmax()))(ix,i,iz) = (*u.get_grid(u.levelmax()))(iix,ny+i,iiz);
				(*u.get_grid(u.levelmax()))(ix,ny-1-i,iz) = (*u.get_grid(u.levelmax()))(iix,-1-i,iiz);
			}
		}
		
	for( int ix=-nb; ix<nx+nb; ++ix )
		for( int iy=-nb; iy<ny+nb; ++iy )
		{
			int iix( (ix+nx)%nx ), iiy( (iy+ny)%ny );
			
			for( int i=-nb; i<0; ++i )
			{
				(*u.get_grid(u.levelmax()))(ix,iy,i) = (*u.get_grid(u.levelmax()))(iix,iiy,nz+i);
				(*u.get_grid(u.levelmax()))(ix,iy,nz-1-i) = (*u.get_grid(u.levelmax()))(iix,iiy,-1-i);
			}
		}
		
		

	
	
	LOGUSER("Done with k-space Poisson solver.");
	return 0.0;
}


double fft_poisson_plugin::gradient( int dir, grid_hierarchy& u, grid_hierarchy& Du )
{
	
	LOGUSER("Computing a gradient in k-space...\n");
	
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
	
	#pragma omp parallel for
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )	
			for( int k=0; k<nz; ++k )
			{
				size_t idx = ((size_t)i*ny+(size_t)j)*nzp+(size_t)k;
				data[idx] = (*u.get_grid(u.levelmax()))(i,j,k);
			}
	
	//... perform FFT and Poisson solve................................
	
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
	fftwf_plan 
		plan  = fftwf_plan_dft_r2c_3d(nx, ny, nz, data, cdata, FFTW_ESTIMATE),
		iplan = fftwf_plan_dft_c2r_3d(nx, ny, nz, cdata, data, FFTW_ESTIMATE);
	
	fftwf_execute(plan);
	#else	
	fftw_plan 
	plan  = fftw_plan_dft_r2c_3d(nx, ny, nz, data, cdata, FFTW_ESTIMATE),
	iplan = fftw_plan_dft_c2r_3d(nx, ny, nz, cdata, data, FFTW_ESTIMATE);
	
	fftw_execute(plan);
	#endif
#else
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
	
#endif
	
	double fac = -1.0/(nx*ny*nz);
	double kfac = 2.0*M_PI;
	
	
	
	bool do_glass = cf_.getValueSafe<bool>("output","glass",false);
	bool deconvolve_cic = do_glass;
	
	#pragma omp parallel for
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )	
			for( int k=0; k<nz/2+1; ++k )
			{
				size_t idx = (i*(size_t)ny+j)*(size_t)nzp/2+(size_t)k;
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
				else //if( dir == 2 )
					kdir = kfac*kk;
				
#ifdef FFTW3
				double re = cdata[idx][0];
				double im = cdata[idx][1];
				
				cdata[idx][0] = fac*im*kdir;
				cdata[idx][1] = -fac*re*kdir;	
				
				if( deconvolve_cic )
				{
					double dfx, dfy, dfz;
					dfx = M_PI*ki/(double)nx; dfx = (i!=0)? sin(dfx)/dfx : 1.0;
					dfy = M_PI*kj/(double)ny; dfy = (j!=0)? sin(dfy)/dfy : 1.0;
					dfz = M_PI*kk/(double)nz; dfz = (k!=0)? sin(dfz)/dfz : 1.0;
					
					dfx = 1.0/(dfx*dfy*dfz); dfx = dfx*dfx;
					cdata[idx][0] *= dfx;
					cdata[idx][1] *= dfx;
					
				}
#else
				double re = cdata[idx].re;
				double im = cdata[idx].im;
				
				cdata[idx].re = fac*im*kdir;
				cdata[idx].im = -fac*re*kdir;	
				
				if( deconvolve_cic )
				{
					double dfx, dfy, dfz;
					dfx = M_PI*ki/(double)nx; dfx = (i!=0)? sin(dfx)/dfx : 1.0;
					dfy = M_PI*kj/(double)ny; dfy = (j!=0)? sin(dfy)/dfy : 1.0;
					dfz = M_PI*kk/(double)nz; dfz = (k!=0)? sin(dfz)/dfz : 1.0;
					
					dfx = 1.0/(dfx*dfy*dfz); dfx = dfx*dfx;

					cdata[idx].re *= dfx;
					cdata[idx].im *= dfx;
				}
#endif			
				
				/*double ktot = sqrt(ii*ii+jj*jj+k*k);
				if( ktot >= nx/2 )//dir == 0 && i==nx/2 || dir == 1 && j==ny/2 || dir == 2 && k==nz/2 )
				{
					cdata[idx].re = 0.0;
					cdata[idx].im = 0.0;
				}*/
			}
	
	
#ifdef FFTW3
	cdata[0][0] = 0.0;
	cdata[0][1] = 0.0;
	
	#ifdef SINGLE_PRECISION
	fftwf_execute(iplan);
	fftwf_destroy_plan(plan);
	fftwf_destroy_plan(iplan);
	#else
	fftw_execute(iplan);
	fftw_destroy_plan(plan);
	fftw_destroy_plan(iplan);
	#endif

#else
	cdata[0].re = 0.0;
	cdata[0].im = 0.0;
	
	#ifndef SINGLETHREAD_FFTW		
	rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), iplan, cdata, NULL );
	#else
	rfftwnd_one_complex_to_real( iplan, cdata, NULL );
	#endif
	
	rfftwnd_destroy_plan(plan);
	rfftwnd_destroy_plan(iplan);
#endif
	
	//... copy data ..........................................
	double dmax = 0.0;
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )	
			for( int k=0; k<nz; ++k )
			{
				size_t idx = ((size_t)i*ny+(size_t)j)*nzp+(size_t)k;
				(*Du.get_grid(u.levelmax()))(i,j,k) = data[idx];
				if(fabs(data[idx])>dmax)
					dmax = fabs(data[idx]);
			}

	//std::cerr << " - component max. is " << dmax*nx << std::endl;
	
	delete[] data;
	
	LOGUSER("Done with k-space gradient.\n");
	
	return 0.0;
}

/**************************************************************************************/
/**************************************************************************************/
#pragma mark -


template<int order>
double poisson_hybrid_kernel( int idir, int i, int j, int k, int n )
{
	return 1.0;
}

template<>
inline double poisson_hybrid_kernel<2>(int idir, int i, int j, int k, int n )
{
	if(i==0&&j==0&&k==0)
		return 0.0;
	
	double 
	ki(M_PI*(double)i/(double)n), 
	kj(M_PI*(double)j/(double)n), 
	kk(M_PI*(double)k/(double)n), 
	kr(sqrt(ki*ki+kj*kj+kk*kk));
	
	double grad = 1.0, laplace = 1.0;
	
	if( idir==0 )
		grad = sin(ki);
	else if( idir==1 )
		grad = sin(kj);
	else 
		grad = sin(kk);
	
	laplace = 2.0*((-cos(ki)+1.0)+(-cos(kj)+1.0)+(-cos(kk)+1.0));
	
	double kgrad = 1.0;
	if( idir==0 )
		kgrad = ki;
	else if( idir ==1)
		kgrad = kj;
	else if( idir ==2)
		kgrad = kk;
	
	return kgrad/kr/kr-grad/laplace;
}

template<>
inline double poisson_hybrid_kernel<4>(int idir, int i, int j, int k, int n )
{
	
	if(i==0&&j==0&&k==0)
	return 0.0;

	double 
	ki(M_PI*(double)i/(double)n), 
	kj(M_PI*(double)j/(double)n), 
	kk(M_PI*(double)k/(double)n), 
	kr(sqrt(ki*ki+kj*kj+kk*kk));

	double grad = 1.0, laplace = 1.0;

	if( idir==0 )
	   grad = 0.166666666667*(-sin(2.*ki)+8.*sin(ki));
	else if( idir==1 )
	   grad = 0.166666666667*(-sin(2.*kj)+8.*sin(kj));
	else if( idir==2 )
	   grad = 0.166666666667*(-sin(2.*kk)+8.*sin(kk));

	laplace = 0.1666666667*((cos(2*ki)-16.*cos(ki)+15.)
						   +(cos(2*kj)-16.*cos(kj)+15.)
						   +(cos(2*kk)-16.*cos(kk)+15.));

	double kgrad = 1.0;
	if( idir==0 )
	kgrad = ki;
	else if( idir ==1)
	kgrad = kj;
	else if( idir ==2)
	kgrad = kk;

	return kgrad/kr/kr-grad/laplace;
}
	   
template<>
inline double poisson_hybrid_kernel<6>(int idir, int i, int j, int k, int n )
{
	double 
	ki(M_PI*(double)i/(double)n), 
	kj(M_PI*(double)j/(double)n), 
	kk(M_PI*(double)k/(double)n), 
	kr(sqrt(ki*ki+kj*kj+kk*kk));

	if(i==0&&j==0&&k==0)
		return 0.0;

	double grad = 1.0, laplace = 1.0;

	if( idir==0 )
		grad = 0.0333333333333*(sin(3.*ki)-9.*sin(2.*ki)+45.*sin(ki));
	else if( idir==1 )
		grad = 0.0333333333333*(sin(3.*kj)-9.*sin(2.*kj)+45.*sin(kj));
	else if( idir==2 )
		grad = 0.0333333333333*(sin(3.*kk)-9.*sin(2.*kk)+45.*sin(kk));

	laplace = 0.01111111111111*(
								(-2.*cos(3.0*ki)+27.*cos(2.*ki)-270.*cos(ki)+245.)
								+(-2.*cos(3.0*kj)+27.*cos(2.*kj)-270.*cos(kj)+245.)
								+(-2.*cos(3.0*kk)+27.*cos(2.*kk)-270.*cos(kk)+245.));

	double kgrad = 1.0;
	if( idir==0 )
		kgrad = ki;
	else if( idir ==1)
		kgrad = kj;
	else if( idir ==2)
		kgrad = kk;

//	if( i*i+j*j+k*k >= n*n )
//		kgrad = 0.0;
	
	return kgrad/kr/kr-grad/laplace;
}
	   
	   
template<int order>
void do_poisson_hybrid( fftw_real* data, int idir, int nxp, int nyp, int nzp, bool periodic )
{
	double fftnorm = 1.0/(nxp*nyp*nzp);
		
	fftnorm = 1.0/(nxp*nyp*nzp);
	
	fftw_complex	*cdata = reinterpret_cast<fftw_complex*>(data);
	
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
	fftwf_plan iplan, plan;
	plan  = fftwf_plan_dft_r2c_3d(nxp, nyp, nzp, data, cdata, FFTW_ESTIMATE);
	iplan = fftwf_plan_dft_c2r_3d(nxp, nyp, nzp, cdata, data, FFTW_ESTIMATE);
	fftwf_execute(plan);
	#else
	fftw_plan iplan, plan;
	plan  = fftw_plan_dft_r2c_3d(nxp, nyp, nzp, data, cdata, FFTW_ESTIMATE);
	iplan = fftw_plan_dft_c2r_3d(nxp, nyp, nzp, cdata, data, FFTW_ESTIMATE);
	fftw_execute(plan);
	#endif
#else
	rfftwnd_plan	iplan, plan;
	
	plan  = rfftw3d_create_plan( nxp, nyp, nzp,
								FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE);
	
	iplan = rfftw3d_create_plan( nxp, nyp, nzp, 
								FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
	
	#ifndef SINGLETHREAD_FFTW		
	rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), plan, data, NULL );
	#else
	rfftwnd_one_real_to_complex( plan, data, NULL );
	#endif
#endif
	
	double ksum = 0.0;
	size_t kcount = 0;
	
	#pragma omp parallel for reduction(+:ksum,kcount)
	for( int i=0; i<nxp; ++i )
		for( int j=0; j<nyp; ++j )
			for( int k=0; k<nzp/2+1; ++k )
			{
				unsigned ii = (i*nyp + j) * (nzp/2+1) + k;
#ifdef FFTW3
				if( k==0 || k==nzp/2 )
				{
					ksum  += cdata[ii][0];
					kcount++;
				}else{
					ksum  += 2.0*(cdata[ii][0]);
					kcount+=2;
				}
#else
				if( k==0 || k==nzp/2 )
				{
					ksum  += cdata[ii].re;
					kcount++;
				}else{
					ksum  += 2.0*(cdata[ii].re);
					kcount+=2;
				}
#endif
			}
	
	ksum /= kcount;
	kcount = 0;
	
	
	#pragma omp parallel for
	for( int i=0; i<nxp; ++i )
		for( int j=0; j<nyp; ++j )
			for( int k=0; k<nzp/2+1; ++k )
			{
				unsigned ii = (i*nyp + j) * (nzp/2+1) + k;
				
				int ki(i), kj(j);
				if( ki > nxp/2 ) ki-=nxp;
				if( kj > nyp/2 ) kj-=nyp;
				
				//... apply hybrid correction
				double dk = poisson_hybrid_kernel<order>(idir, ki, kj, k, nxp/2 );
				//cdata[ii].re -= ksum;
				
#ifdef FFTW3
				fftw_real re = cdata[ii][0], im = cdata[ii][1];
				
				cdata[ii][0] = -im*dk*fftnorm;
				cdata[ii][1] = re*dk*fftnorm;
#else
				fftw_real re = cdata[ii].re, im = cdata[ii].im;
				
				cdata[ii].re = -im*dk*fftnorm;
				cdata[ii].im = re*dk*fftnorm;
#endif
				//cdata[ii].re += ksum*fftnorm;
				
				//if( i==nxp/2||j==nyp/2||k==nzp/2 )
				//{
				//	cdata[ii].re = 0.0;
				//	cdata[ii].im = 0.0;
				//}
				
				
				
			}
#ifdef FFTW3
	cdata[0][0] = 0.0;
	cdata[0][1] = 0.0;
	
	#ifdef SINGLE_PRECISION
	fftwf_execute(iplan);
	fftwf_destroy_plan(plan);
	fftwf_destroy_plan(iplan);
	#else
	fftw_execute(iplan);
	fftw_destroy_plan(plan);
	fftw_destroy_plan(iplan);
	#endif	
#else
	cdata[0].re = 0.0;
	cdata[0].im = 0.0;
	
	#ifndef SINGLETHREAD_FFTW		
	rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), iplan, cdata, NULL);
	#else		
	rfftwnd_one_complex_to_real(iplan, cdata, NULL);
	#endif
	
	rfftwnd_destroy_plan(plan);
	rfftwnd_destroy_plan(iplan);
#endif
	
}
   
template< typename T >
void poisson_hybrid( T& f, int idir, int order, bool periodic )

{
	int nx=f.size(0), ny=f.size(1), nz=f.size(2), nxp, nyp, nzp;
	fftw_real		*data;
	int xo=0,yo=0,zo=0;
	int nmax = std::max(nx,std::max(ny,nz));
	
	LOGUSER("Entering hybrid Poisson solver...");
	
	if(!periodic)
	{
		nxp = 2*nmax;
		nyp = 2*nmax;
		nzp = 2*nmax;
		xo  = nmax/2;
		yo  = nmax/2;
		zo  = nmax/2;
	}
	else
	{
		nxp = nmax;
		nyp = nmax;
		nzp = nmax;
	}
	
	
	
	data		= new fftw_real[(size_t)nxp*(size_t)nyp*(size_t)(nzp+2)];
	
	if(idir==0)
		std::cout << "   - Performing hybrid Poisson step... (" << nxp <<  ", " << nyp << ", " << nzp << ")\n";
	
	//size_t N = (size_t)nxp*(size_t)nyp*2*((size_t)nzp/2+1);
	
	//#pragma omp parallel for
	//for( size_t i=0; i<N; ++i )
	//	data[i]=0.0;
	
	#pragma omp parallel for
	for( int i=0; i<nxp; ++i )
		for( int j=0; j<nyp; ++j )
			for( int k=0; k<=nzp; ++k )
			{
				size_t idx = ((size_t)i*(size_t)nxp+(size_t)j)*(size_t)(nzp+2)+(size_t)k;
				data[idx] = 0.0;
			}
	
	#pragma omp parallel for
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )
			for( int k=0; k<nz; ++k )
			{
				size_t idx = ((i+xo)*nyp + j+yo) * 2*(nzp/2+1) + k+zo;
				data[idx] = f(i,j,k);
			}
	
	switch (order) {
		case 2:
			do_poisson_hybrid<2>(data, idir, nxp, nyp, nzp, periodic);
			break;
		case 4:
			do_poisson_hybrid<4>(data, idir, nxp, nyp, nzp, periodic);
			break;
		case 6:
			do_poisson_hybrid<6>(data, idir, nxp, nyp, nzp, periodic);
			break;
		default:
			std::cerr << " - ERROR: invalid operator order specified in deconvolution.";
			LOGERR("Invalid operator order specified in deconvolution.");
			break;
	}
	
	LOGUSER("Copying hybrid correction factor...");
	
	#pragma omp parallel for
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )
			for( int k=0; k<nz; ++k )
			{
				size_t idx = (((size_t)i+xo)*nyp + (size_t)j+yo) * 2*(nzp/2+1) + (size_t)k+zo;	
				f(i,j,k) = data[idx];
			}
	
	delete[] data;

	LOGUSER("Done with hybrid Poisson solve.");
}
	   
	   
/**************************************************************************************/
/**************************************************************************************/
#pragma mark -

template void poisson_hybrid< MeshvarBnd<double> >( MeshvarBnd<double>& f, int idir, int order, bool periodic );
template void poisson_hybrid< MeshvarBnd<float> >( MeshvarBnd<float>& f, int idir, int order, bool periodic );

namespace{
	poisson_plugin_creator_concrete<multigrid_poisson_plugin> multigrid_poisson_creator("mg_poisson");
	poisson_plugin_creator_concrete<fft_poisson_plugin> fft_poisson_creator("fft_poisson");
}
