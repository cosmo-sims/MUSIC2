/*
 
 main.cc - This file is part of MUSIC -
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


#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>



#include "general.hh"
#include "defaults.hh"
#include "output.hh"

#include "config_file.hh"

#include "poisson.hh"
#include "mg_solver.hh"
#include "fd_schemes.hh"
#include "random.hh"
#include "densities.hh"

#include "convolution_kernel.hh"
#include "cosmology.hh"
#include "transfer_function.hh"

#define THE_CODE_NAME "music!"
#define THE_CODE_VERSION "0.7.3a"


namespace music
{

	struct framework
	{
		transfer_function *the_transfer_function;
		//poisson_solver *the_poisson_solver;
		config_file *the_config_file;
		refinement_hierarchy *the_refinement_hierarchy;
	};
	
}
 
transfer_function *TransferFunction_real::ptf_ = NULL;
//transfer_function *TransferFunction_real::ptf = NULL;
transfer_function *TransferFunction_k::ptf_ = NULL;

tf_type TransferFunction_k::type_;
tf_type TransferFunction_real::type_;


real_t TransferFunction_real::nspec_ = -1.0;
//real_t TransferFunction_real::nspec = -1.0;
real_t TransferFunction_k::nspec_ = -1.0;

void splash(void)
{
	
	std::cout 
	<< "\n    __    __     __  __     ______     __     ______      \n"
	<< "   /\\ \"-./  \\   /\\ \\/\\ \\   /\\  ___\\   /\\ \\   /\\  ___\\  \n"   
	<< "   \\ \\ \\-./\\ \\  \\ \\ \\_\\ \\  \\ \\___  \\  \\ \\ \\  \\ \\ \\____ \n"  
	<< "    \\ \\_\\ \\ \\_\\  \\ \\_____\\  \\/\\_____\\  \\ \\_\\  \\ \\_____\\ \n"
	<< "     \\/_/  \\/_/   \\/_____/   \\/_____/   \\/_/   \\/_____/ \n\n"
	<< "                            this is " << THE_CODE_NAME << " version " << THE_CODE_VERSION << "\n\n\n";
	
	
}

void modify_grid_for_TF( const refinement_hierarchy& rh_full, refinement_hierarchy& rh_TF, config_file& cf )
{
	unsigned lbase, lbaseTF, lmax, overlap;
	
	lbase				= cf.getValue<unsigned>( "setup", "levelmin" );
	lmax				= cf.getValue<unsigned>( "setup", "levelmax" );
	lbaseTF				= cf.getValueSafe<unsigned>( "setup", "levelmin_TF", lbase );
	overlap				= cf.getValueSafe<unsigned>( "setup", "overlap", 4 );
	
	rh_TF = rh_full;
	
	unsigned pad = overlap;
	
	for( unsigned i=lbase+1; i<=lmax; ++i )
	{
		int x0[3], lx[3], lxmax = 0;
		
		for( int j=0; j<3; ++j )
		{
			lx[j] = rh_TF.size(i,j)+2*pad;
			x0[j] = rh_TF.offset_abs(i,j)-pad;
			
			if( lx[j] > lxmax )
				lxmax = lx[j];
		}
		
		//... make sure that grids are divisible by 4 for convolution.
		lxmax += lxmax%4;
		
		for( int j=0; j<3; ++j )
		{
			double dl = 0.5*((double)(lxmax-lx[j]));
			int add_left = (int)ceil(dl);
			
			lx[j] = lxmax;
			x0[j] -= add_left;
			x0[j] += x0[j]%2;
		}
		
		rh_TF.adjust_level(i, lx[0], lx[1], lx[2], x0[0], x0[1], x0[2] );
	}
	
	if( lbaseTF > lbase )
	{
		std::cout << " - Will use levelmin = " << lbaseTF << " to compute density field...\n";
	
		for( unsigned i=lbase; i<=lbaseTF; ++i )
		{
			unsigned nfull = (unsigned)pow(2,i);
			rh_TF.adjust_level(i, nfull, nfull, nfull, 0, 0, 0);
		}
	}
	
	
		
}

void coarsen_density( const refinement_hierarchy& rh, grid_hierarchy& u )
{
	for( int i=rh.levelmax(); i>0; --i )
		mg_straight().restrict( *(u.get_grid(i)), *(u.get_grid(i-1)) );

	for( unsigned i=1; i<=rh.levelmax(); ++i )
	{
		if(	  rh.offset(i,0) != u.get_grid(i)->offset(0)
   		   || rh.offset(i,1) != u.get_grid(i)->offset(1)
		   || rh.offset(i,2) != u.get_grid(i)->offset(2)
		   || rh.size(i,0) != u.get_grid(i)->size(0)
		   || rh.size(i,1) != u.get_grid(i)->size(1)
		   || rh.size(i,2) != u.get_grid(i)->size(2) )
		{
			u.cut_patch(i, rh.offset_abs(i,0), rh.offset_abs(i,1), rh.offset_abs(i,2), 
						rh.size(i,0), rh.size(i,1), rh.size(i,2) );
		}
	}
	
	
}

void store_grid_structure( config_file& cf, const refinement_hierarchy& rh )
{
	char str1[128], str2[128];
	for( unsigned i=rh.levelmin(); i<=rh.levelmax(); ++i )
	{
		for( int j=0; j<3; ++j )
		{
			sprintf(str1,"offset(%d,%d)",i,j);	
			sprintf(str2,"%d",rh.offset(i,j));
			cf.insertValue("setup",str1,str2);

			sprintf(str1,"size(%d,%d)",i,j);	
			sprintf(str2,"%ld",rh.size(i,j));
			cf.insertValue("setup",str1,str2);
			
		}		
	}
}

void subtract_finest_mean( grid_hierarchy& u )
{
	std::cout << " - Subtracting component mean...\n";
	double sum = 0.0;
	for( int ix = 0; ix < (int)(*u.get_grid(u.levelmax())).size(0); ++ix )
		for( int iy = 0; iy < (int)(*u.get_grid(u.levelmax())).size(1); ++iy )
			for( int iz = 0; iz < (int)(*u.get_grid(u.levelmax())).size(2); ++iz )
				sum += 0.5*(*u.get_grid(u.levelmax()))(ix,iy,iz);
	
	sum /= (*u.get_grid(u.levelmax())).size(0)
			* (*u.get_grid(u.levelmax())).size(1)
			* (*u.get_grid(u.levelmax())).size(2);
	
	std::cout << "     component mean is " << sum << std::endl;
	
	for( unsigned ilevel=u.levelmin(); ilevel<=u.levelmax(); ++ilevel )
		#pragma omp parallel for
		for( int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix )
			for( int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy )
				for( int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz )
					(*u.get_grid(ilevel))(ix,iy,iz) -= sum;
	
}


/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/



int main (int argc, const char * argv[]) 
{
	const unsigned nbnd = 4;
	
	unsigned lbase, lmax, lbaseTF;
	double   err;
	
	cosmology cosmo;
	double boxlength, zstart;
	std::vector<long> rngseeds;
	std::vector<std::string> rngfnames;
	
	double x0[3], lx[3];
	unsigned npad = 8;
	
	splash();
	if( argc != 2 ){
		std::cout << " This version is compiled with the following plug-ins:\n";
		
		print_transfer_function_plugins();
		print_output_plugins();
		
		std::cerr << "\n In order to run, you need to specify a parameter file!\n\n";
		exit(0);
	}
	
	//... open log file
	char logfname[128];
	sprintf(logfname,"%s_log.txt",argv[1]);
	MUSIC::log::setOutput(logfname);
	time_t ltime=time(NULL);
	LOGINFO("Opening log file \'%s\'.",logfname);
	LOGUSER("Running %s, version %s",THE_CODE_NAME,THE_CODE_VERSION);
	LOGUSER("Running with a maximum of %d OpenMP threads", omp_get_max_threads() );
	LOGUSER("Log is for run started %s",asctime( localtime(&ltime) ));
	
#ifdef SINGLETHREAD_FFTW
	LOGUSER("Code was compiled for single-threaded FFTW");
#else
	LOGUSER("Code was compiled for multi-threaded FFTW");
#endif
	
#ifdef SINGLE_PRECISION
	LOGUSER("Code was compiled for single precision.");
#else
	LOGUSER("Code was compiled for double precision.");
#endif
	
	
	/******************************************************************************************************/
	/* read and interpret config file *********************************************************************/
	/******************************************************************************************************/
	config_file cf(argv[1]);
	std::string tfname,randfname,temp, outformat, outfname, poisson_solver_name;;
	bool shift_back(false), align_top(false), kspace(false), force_shift(false);
	float tf0,tf1,tf2;
	
	
	boxlength           = cf.getValue<double>( "setup", "boxlength" );
	lbase				= cf.getValue<unsigned>( "setup", "levelmin" );
	lmax				= cf.getValue<unsigned>( "setup", "levelmax" );
	lbaseTF				= cf.getValueSafe<unsigned>( "setup", "levelmin_TF", lbase );
	force_shift			= cf.getValueSafe<bool>("setup", "force_shift", force_shift );
	
	
	if( lbase == lmax && !force_shift )
		cf.insertValue("setup","no_shift","yes");
	
	if( lbaseTF < lbase )
	{
		std::cout << " - WARNING: levelminTF < levelmin. This is not good!\n"
				  << "            I will set levelminTF = levelmin.\n";
		
		LOGUSER("levelminTF < levelmin. set levelminTF = levelmin.");
		
		lbaseTF = lbase;
		cf.insertValue("setup","levelmin_TF",cf.getValue<std::string>("setup","levelmin"));
	}
	
	temp				= cf.getValue<std::string>( "setup", "ref_offset" );
	sscanf( temp.c_str(), "%g,%g,%g", &tf0, &tf1, &tf2 ); x0[0] = tf0; x0[1] = tf1; x0[2] = tf2;
	
	
	temp				= cf.getValue<std::string>( "setup", "ref_extent" );
	sscanf( temp.c_str(), "%g,%g,%g", &tf0, &tf1, &tf2 ); lx[0] = tf0; lx[1] = tf1; lx[2] = tf2;
	
	npad                = cf.getValue<unsigned>( "setup", "padding" );
	align_top			= cf.getValueSafe<bool>( "setup", "align_top", false );
	kspace				= cf.getValueSafe<bool>( "poisson", "kspace", false );
	
	if( kspace )
		poisson_solver_name = std::string("fft_poisson");
	else
		poisson_solver_name = std::string("mg_poisson");
	
	// TODO: move cosmology parameters reading to cosmo_calc
	zstart				= cf.getValue<double>( "setup", "zstart" );
	cosmo.astart		= 1.0/(1.0+zstart);
	cosmo.Omega_b		= cf.getValue<double>( "cosmology", "Omega_b" );
	cosmo.Omega_m		= cf.getValue<double>( "cosmology", "Omega_m" );
	cosmo.Omega_L		= cf.getValue<double>( "cosmology", "Omega_L" );
	cosmo.H0			= cf.getValue<double>( "cosmology", "H0" );
	cosmo.sigma8		= cf.getValue<double>( "cosmology", "sigma_8" );
	cosmo.nspect		= cf.getValue<double>( "cosmology", "nspec" );
	cosmo.WDMg_x		= cf.getValueSafe<double>( "cosmology", "WDMg_x", 1.5 );
	cosmo.WDMmass		= cf.getValueSafe<double>( "cosmology", "WDMmass", 0.0 );
	cosmo.dplus			= 0.0;
	cosmo.pnorm			= 0.0;
	cosmo.vfact			= 0.0;
	
	//cosmo.Gamma			= cf.getValueSafe<double>( "cosmology", "Gamma", -1.0 );
	
	/******************************************************************************************************/
	/******************************************************************************************************/
	
	shift_back			= cf.getValueSafe<bool>( "output", "shift_back", shift_back );
	outformat			= cf.getValue<std::string>( "output", "format" );
	outfname			= cf.getValue<std::string>( "output", "filename" );
	
	
	unsigned grad_order = cf.getValueSafe<unsigned> ( "poisson" , "grad_order", 4 );

	bool bdefd = cf.getValueSafe<bool> ( "poisson" , "fft_fine", true );
	
	//... if in unigrid mode, use k-space instead
	//if(bdefd&lbase==lmax)
	//kspace=true;
	
	//... switch off if using kspace anyway
	bdefd &= !kspace;
	
	/******************************************************************************************************/
	/******************************************************************************************************/
	/******************************************************************************************************/

#if not defined(SINGLETHREAD_FFTW)
#ifdef FFTW3
	fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());
#else
	fftw_threads_init();
#endif
#endif
	
	transfer_function_plugin *the_transfer_function_plugin
		= select_transfer_function_plugin( cf );
	
	CosmoCalc ccalc(cosmo,the_transfer_function_plugin);
	cosmo.pnorm	= ccalc.ComputePNorm( 2.0*M_PI/boxlength );
	cosmo.dplus	= ccalc.CalcGrowthFactor( cosmo.astart )/ccalc.CalcGrowthFactor( 1.0 );
	cosmo.vfact = ccalc.ComputeVFact( cosmo.astart );
	{
		char tmpstr[128];
		sprintf(tmpstr,"%.12g",cosmo.pnorm);
		cf.insertValue("cosmology","pnorm",tmpstr);
		sprintf(tmpstr,"%.12g",cosmo.dplus);
		cf.insertValue("cosmology","dplus",tmpstr);
		sprintf(tmpstr,"%.12g",cosmo.vfact);
		cf.insertValue("cosmology","vfact",tmpstr);
		
	}
	
	/******************************************************************************************************/
	/******************************************************************************************************/
	
	bool 
		do_baryons	= cf.getValue<bool>("setup","baryons"),
		do_2LPT		= cf.getValue<bool>("setup","use_2LPT"),
		do_LLA		= cf.getValue<bool>("setup","use_LLA"),
		do_CVM		= cf.getValueSafe<bool>("setup","center_velocities",false);
	

	//... determine the refinement hierarchy
	refinement_hierarchy rh_Poisson( cf );
	store_grid_structure(cf, rh_Poisson);
	rh_Poisson.output();
	
	refinement_hierarchy rh_TF( rh_Poisson );
	modify_grid_for_TF( rh_Poisson, rh_TF, cf );
	//rh_TF.output();

	LOGUSER("Grid structure for Poisson solver:");
	rh_Poisson.output_log();
	LOGUSER("Grid structure for density convolution:");
	rh_TF.output_log();
	
	
	if( !the_transfer_function_plugin->tf_is_distinct() && do_baryons )
		std::cout	<< " - WARNING: The selected transfer function does not support\n"
					<< "            distinct amplitudes for baryon and DM fields!\n"
					<< "            Perturbation amplitudes will be identical!" << std::endl;
	
	
	//... initialize the output plug-in
	output_plugin *the_output_plugin = select_output_plugin( cf );
	
	//... initialize the random numbers
	rand_gen rand( cf, rh_TF );
	
	//... initialize the Poisson solver
	poisson_plugin_creator *the_poisson_plugin_creator = get_poisson_plugin_map()[ poisson_solver_name ];
	poisson_plugin *the_poisson_solver = the_poisson_plugin_creator->create( cf );
	
	//... THIS IS THE MAIN DRIVER BRANCHING TREE RUNNING THE VARIOUS PARTS OF THE CODE
	bool bfatal = false;
	try{
		if( ! do_2LPT )
		{
			LOGUSER("Entering 1LPT branch");
			
			//... cdm density and displacements
			std::cout << "=============================================================\n";
			std::cout << "   COMPUTING DARK MATTER DISPLACEMENTS\n";
			std::cout << "-------------------------------------------------------------\n";
			LOGUSER("Computing dark matter displacements...");
			
			grid_hierarchy f( nbnd ), u(nbnd);
			
			GenerateDensityHierarchy(	cf, the_transfer_function_plugin, cdm , rh_TF, rand, f, true, false );
			coarsen_density(rh_Poisson, f);
			normalize_density(f);
			
			u = f;	u.zero();
			
			the_output_plugin->write_dm_mass(f);
			the_output_plugin->write_dm_density(f);
			
			err = the_poisson_solver->solve(f, u);
			
			if(!bdefd)
				f.deallocate();
			
			the_output_plugin->write_dm_potential(u);
			
			
			//... DM displacements
			{
				grid_hierarchy data_forIO(u);
				for( int icoord = 0; icoord < 3; ++icoord )
				{
					if( bdefd )
					{
						data_forIO.zero();
						*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
						poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order, data_forIO.levelmin()==data_forIO.levelmax());					
						*data_forIO.get_grid(data_forIO.levelmax()) /= 1<<f.levelmax();
						the_poisson_solver->gradient_add(icoord, u, data_forIO );
						
					}
					else
						//... displacement
						the_poisson_solver->gradient(icoord, u, data_forIO );

					the_output_plugin->write_dm_position(icoord, data_forIO );
				}
			}
			
			//... gas density
			if( do_baryons )
			{
				std::cout << "=============================================================\n";
				std::cout << "   COMPUTING BARYON DENSITY\n";
				std::cout << "-------------------------------------------------------------\n";
				LOGUSER("Computing baryon density...");
				
				GenerateDensityHierarchy(	cf, the_transfer_function_plugin, baryon , rh_TF, rand, f, false, true );
				coarsen_density(rh_Poisson, f);
				normalize_density(f);
				
				if( do_LLA )
				{
					u = f;	u.zero();
					err = the_poisson_solver->solve(f, u);
					compute_LLA_density( u, f,grad_order );
					normalize_density(f);
				}
				
				the_output_plugin->write_gas_density(f);
			}
			
			std::cout << "=============================================================\n";
			std::cout << "   COMPUTING VELOCITIES\n";
			std::cout << "-------------------------------------------------------------\n";
			LOGUSER("Computing velocitites...");
			
			//... velocities
			if( do_baryons )
			{
				GenerateDensityHierarchy(	cf, the_transfer_function_plugin, total , rh_TF, rand, f, true, false );
				coarsen_density(rh_Poisson, f);
				normalize_density(f);
				
				u = f;	u.zero();
					
				err = the_poisson_solver->solve(f, u);
				
				if(!bdefd)
					f.deallocate();
			}
			grid_hierarchy data_forIO(u);
			for( int icoord = 0; icoord < 3; ++icoord )
			{
				//... displacement
				if(bdefd)
				{
					data_forIO.zero();
					*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
					poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order, data_forIO.levelmin()==data_forIO.levelmax());					
					*data_forIO.get_grid(data_forIO.levelmax()) /= 1<<f.levelmax();
					the_poisson_solver->gradient_add(icoord, u, data_forIO );
				}
				else 
					the_poisson_solver->gradient(icoord, u, data_forIO );

				//... multiply to get velocity
				data_forIO *= cosmo.vfact;
				
				if(do_CVM)
					subtract_finest_mean(data_forIO);
				
				the_output_plugin->write_dm_velocity(icoord, data_forIO);
				if( do_baryons )
					the_output_plugin->write_gas_velocity(icoord, data_forIO);
			}
		
		}else {
			//.. use 2LPT ...
			LOGUSER("Entering 2LPT branch");
			
			grid_hierarchy f( nbnd ), u1(nbnd), u2LPT(nbnd), f2LPT( nbnd );
			
			std::cout << "=============================================================\n";
			std::cout << "   COMPUTING VELOCITIES\n";
			std::cout << "-------------------------------------------------------------\n";	
			LOGUSER("Computing velocities...");
			
			GenerateDensityHierarchy(	cf, the_transfer_function_plugin, total , rh_TF, rand, f, true, false );
			coarsen_density(rh_Poisson, f);
			normalize_density(f);
			
			u1 = f;	u1.zero();
			
			if(bdefd)
				f2LPT=f;

			//... compute 1LPT term
			err = the_poisson_solver->solve(f, u1);
			the_output_plugin->write_dm_potential(u1);
			
			//... compute 2LPT term
			u2LPT = f; u2LPT.zero();
			
			if( !kspace )
				compute_2LPT_source(u1, f2LPT, grad_order );
			else
				compute_2LPT_source_FFT(cf, u1, f2LPT);
			
				
			err = the_poisson_solver->solve(f2LPT, u2LPT);
			
			//... if doing the hybrid step, we need a combined source term
			if( bdefd )
			{
				f2LPT*=6.0/7.0;
				f+=f2LPT;
				
				if( do_baryons )
					f2LPT.deallocate();
			}
			
			//... add the 2LPT contribution
			u2LPT *= 6.0/7.0;
			u1 += u2LPT;
			
			if( do_baryons )
				u2LPT.deallocate();
			
			grid_hierarchy data_forIO(u1);
			for( int icoord = 0; icoord < 3; ++icoord )
			{
				if(bdefd)
				{
					data_forIO.zero();
					*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
					poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order, data_forIO.levelmin()==data_forIO.levelmax());					
					*data_forIO.get_grid(data_forIO.levelmax()) /= 1<<f.levelmax();
					the_poisson_solver->gradient_add(icoord, u1, data_forIO );
				}
				else 
					the_poisson_solver->gradient(icoord, u1, data_forIO );
				
				data_forIO *= cosmo.vfact;
				
				if( do_CVM )
					subtract_finest_mean(data_forIO);
					
				the_output_plugin->write_dm_velocity(icoord, data_forIO);					
				
				if( do_baryons )
					the_output_plugin->write_gas_velocity(icoord, data_forIO);				
			}
			data_forIO.deallocate();
			
			
			std::cout << "=============================================================\n";
			std::cout << "   COMPUTING DARK MATTER DISPLACEMENTS\n";
			std::cout << "-------------------------------------------------------------\n";
			LOGUSER("Computing dark matter displacements...");
			
			//... if baryons are enabled, the displacements have to be recomputed
			//... otherwise we can compute them directly from the velocities
			if( do_baryons )
			{
				GenerateDensityHierarchy(	cf, the_transfer_function_plugin, cdm , rh_TF, rand, f, true, false );
				coarsen_density(rh_Poisson, f);
				normalize_density(f);
				the_output_plugin->write_dm_density(f);
				the_output_plugin->write_dm_mass(f);
				u1 = f;	u1.zero();
				
				if(bdefd)
					f2LPT=f;
				
				//... compute 1LPT term
				err = the_poisson_solver->solve(f, u1);
				
				//... compute 2LPT term
				u2LPT = f; u2LPT.zero();
				
				if( !kspace )
					compute_2LPT_source(u1, f2LPT, grad_order );
				else
					compute_2LPT_source_FFT(cf, u1, f2LPT);
				
				err = the_poisson_solver->solve(f2LPT, u2LPT);
				
				if( bdefd )
				{
					f2LPT*=3.0/7.0;
					f+=f2LPT;
					f2LPT.deallocate();
				}
				
				u2LPT *= 3.0/7.0;
				u1 += u2LPT;
				u2LPT.deallocate();
			}else{
				//... reuse prior data
				f-=f2LPT;
				the_output_plugin->write_dm_density(f);
				the_output_plugin->write_dm_mass(f);
				f+=f2LPT;
				
				u2LPT *= 0.5;
				u1 -= u2LPT;
				u2LPT.deallocate();
				
				if(bdefd)
				{
					f2LPT *= 0.5;
					f-=f2LPT;
					f2LPT.deallocate();
				}
			}
						
			data_forIO = u1;
			
			for( int icoord = 0; icoord < 3; ++icoord )
			{
				//... displacement
				if(bdefd)
				{
					data_forIO.zero();
					*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
					poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order, data_forIO.levelmin()==data_forIO.levelmax());					
					*data_forIO.get_grid(data_forIO.levelmax()) /= 1<<f.levelmax();
					the_poisson_solver->gradient_add(icoord, u1, data_forIO );
				}
				else 
					the_poisson_solver->gradient(icoord, u1, data_forIO );
				
				the_output_plugin->write_dm_position(icoord, data_forIO );	
			}
			

			if( do_baryons )
			{	
				std::cout << "=============================================================\n";
				std::cout << "   COMPUTING BARYON DENSITY\n";
				std::cout << "-------------------------------------------------------------\n";
				LOGUSER("Computing baryon density...");
				
				GenerateDensityHierarchy(	cf, the_transfer_function_plugin, baryon , rh_TF, rand, f, false, true );
				coarsen_density(rh_Poisson, f);
				normalize_density(f);
				
				if( !do_LLA )
					the_output_plugin->write_gas_density(f);
				else 
				{	
					u1 = f;	u1.zero();
					
					//... compute 1LPT term
					err = the_poisson_solver->solve(f, u1);
					
					//... compute 2LPT term
					u2LPT = f; u2LPT.zero();
					
					if( !kspace )
						compute_2LPT_source(u1, f2LPT, grad_order );
					else
						compute_2LPT_source_FFT(cf, u1, f2LPT);
					
					err = the_poisson_solver->solve(f2LPT, u2LPT);
					u2LPT *= 3.0/7.0;
					u1 += u2LPT;
					u2LPT.deallocate();
					
					compute_LLA_density( u1, f, grad_order );
					normalize_density(f);
					the_output_plugin->write_gas_density(f);
				}
			}
		}
		
		
		
	}catch(std::runtime_error& excp){
		LOGERR("Fatal error occured. Code will exit.");
		std::cerr << " - " << excp.what() << std::endl;
		std::cerr << " - A fatal error occured. We need to exit...\n";
		bfatal = true;
	}

	std::cout << "=============================================================\n";
	
	//... clean up
	the_output_plugin->finalize();
	delete the_output_plugin;
	
	
	if( !bfatal )
	{	
		std::cout << " - Wrote output file \'" << outfname << "\'\n     using plugin \'" << outformat << "\'...\n";
		LOGUSER("Wrote output file \'%s\'.",outfname.c_str());
	}
	
	delete the_transfer_function_plugin;
	delete the_poisson_solver;
	
	/** we are done ! **/
	std::cout << " - Done!" << std::endl << std::endl;
	ltime=time(NULL);
	LOGUSER("Run finished succesfully on %s",asctime( localtime(&ltime) ));
	
	///*****************************************///
	
	/*std::string save_fname(std::string(argv[1])+std::string("_stats"));
	std::ofstream ofs(save_fname.c_str());
	time_t ltime=time(NULL);
	
	ofs << "Parameter dump for the run on " << asctime( localtime(&ltime) );
	ofs << "You ran " << THE_CODE_NAME << " version " << THE_CODE_VERSION << std::endl << std::endl;
	
	cf.dump( ofs );
*/
	
#ifdef FFTW3
	fftw_cleanup_threads();
#endif
	
	cf.log_dump();
	
	
	return 0;
}
