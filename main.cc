/*
 
 main.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
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
#define THE_CODE_VERSION "1.5"


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
 
//... declare static class members here
transfer_function *TransferFunction_real::ptf_ = NULL;
transfer_function *TransferFunction_k::ptf_ = NULL;
tf_type TransferFunction_k::type_;
tf_type TransferFunction_real::type_;
real_t TransferFunction_real::nspec_ = -1.0;
real_t TransferFunction_k::nspec_ = -1.0;


//... prototypes for routines used in main driver routine
void splash(void);
void modify_grid_for_TF( const refinement_hierarchy& rh_full, refinement_hierarchy& rh_TF, config_file& cf );
void print_hierarchy_stats( config_file& cf, const refinement_hierarchy& rh );
void store_grid_structure( config_file& cf, const refinement_hierarchy& rh );
double compute_finest_mean( grid_hierarchy& u );
double compute_finest_sigma( grid_hierarchy& u );


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

void print_hierarchy_stats( config_file& cf, const refinement_hierarchy& rh )
{
	double omegam = cf.getValue<double>("cosmology","Omega_m");
	double omegab = cf.getValue<double>("cosmology","Omega_b");
	bool bbaryons = cf.getValue<bool>("setup","baryons");
	double boxlength = cf.getValue<double>("setup","boxlength");
	
	unsigned levelmin = rh.levelmin();
	double dx = boxlength/(double)(1<<levelmin), dx3=dx*dx*dx;
	double rhom = 2.77519737e11; // h-1 M_o / (h-1 Mpc)**3
	double cmass, bmass(0.0), mtotgrid;
	if( bbaryons )
	{
		cmass = (omegam-omegab)*rhom*dx3;
		bmass = omegab*rhom*dx3;
	}else
		cmass = omegam*rhom*dx3;
	
	std::cout << "-------------------------------------------------------------\n";
	
	if( rh.get_shift(0)!=0||rh.get_shift(1)!=0||rh.get_shift(2)!=0 )
		std::cout << " - Domain will be shifted by (" << rh.get_shift(0) << ", " << rh.get_shift(1) << ", " << rh.get_shift(2) << ")\n" << std::endl;
	
	std::cout << " - Grid structure:\n";
	
	 
	
	for( unsigned ilevel=rh.levelmin(); ilevel<=rh.levelmax(); ++ilevel )
	{
		double rfac = 1.0/(1<<(ilevel-rh.levelmin())), rfac3=rfac*rfac*rfac;
		
		mtotgrid = omegam*rhom*dx3*rfac3*rh.size(ilevel, 0)*rh.size(ilevel, 1)*rh.size(ilevel, 2);
		std::cout 
		<< "     Level " << std::setw(3) << ilevel << " :   offset = (" << std::setw(5) << rh.offset(ilevel,0) << ", " << std::setw(5) << rh.offset(ilevel,1) << ", " << std::setw(5) << rh.offset(ilevel,2) << ")\n"
		<< "                     size = (" << std::setw(5) << rh.size(ilevel,0) << ", " << std::setw(5) << rh.size(ilevel,1) << ", " << std::setw(5) << rh.size(ilevel,2) << ")\n";
		
		if( ilevel == rh.levelmax() )
		{
			std::cout << "-------------------------------------------------------------\n";
			std::cout << " - Finest level :\n";

			if( dx*rfac > 0.1 )	
			  std::cout << "                   extent =  " << dx*rfac*rh.size(ilevel,0) << " x " << dx*rfac*rh.size(ilevel,1) << " x " << dx*rfac * rh.size(ilevel,2) << " h-3 Mpc**3\n";
			else if( dx*rfac > 1e-4 )
			  std::cout << "                   extent =  " << dx*rfac*1000.0*rh.size(ilevel,0) << " x " << dx*rfac*1000.0*rh.size(ilevel,1) << " x " << dx*rfac*1000.0*rh.size(ilevel,2) << " h-3 kpc**3\n";
			else
			  std::cout << "                   extent =  " << dx*rfac*1.e6*rh.size(ilevel,0) << " x " << dx*rfac*1.e6*rh.size(ilevel,1) << " x " << dx*rfac*1.e6 * rh.size(ilevel,2) << " h-3 pc**3\n";

			std::cout << "                 mtotgrid =  " << mtotgrid << " h-1 M_o\n";
			std::cout << "            particle mass =  " << cmass*rfac3 << " h-1 M_o\n";
			if( bbaryons )
				std::cout << "         baryon mass/cell =  " << bmass*rfac3 << " h-1 M_o\n";
			if( dx*rfac > 0.1 )
				std::cout << "                       dx =  " << dx*rfac << " h-1 Mpc\n";
			else if( dx*rfac > 1e-4 )
				std::cout << "                       dx =  " << dx*rfac*1000.0 << " h-1 kpc\n";
			else
				std::cout << "                       dx =  " << dx*rfac*1.e6 << " h-1 pc\n";
		}
		
	}
	std::cout << "-------------------------------------------------------------\n";
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

double compute_finest_mean( grid_hierarchy& u )
{

	double sum = 0.0;
    size_t count = 0;
	for( int ix = 0; ix < (int)(*u.get_grid(u.levelmax())).size(0); ++ix )
		for( int iy = 0; iy < (int)(*u.get_grid(u.levelmax())).size(1); ++iy )
			for( int iz = 0; iz < (int)(*u.get_grid(u.levelmax())).size(2); ++iz )
                if( ! u.is_refined(u.levelmax(),ix,iy,iz) )
                {
                    sum += (*u.get_grid(u.levelmax()))(ix,iy,iz);
                    ++count;
                }
	sum /= count;
	return sum;
	
}

double compute_finest_sigma( grid_hierarchy& u )
{
	double sum = 0.0, sum2 = 0.0;
	for( int ix = 0; ix < (int)(*u.get_grid(u.levelmax())).size(0); ++ix )
		for( int iy = 0; iy < (int)(*u.get_grid(u.levelmax())).size(1); ++iy )
			for( int iz = 0; iz < (int)(*u.get_grid(u.levelmax())).size(2); ++iz )
			{
				sum +=  (*u.get_grid(u.levelmax()))(ix,iy,iz);
				sum2 +=  (*u.get_grid(u.levelmax()))(ix,iy,iz)* (*u.get_grid(u.levelmax()))(ix,iy,iz);
			}

	size_t N = (size_t)(*u.get_grid(u.levelmax())).size(0)
		 * (size_t)(*u.get_grid(u.levelmax())).size(1)
		 * (size_t)(*u.get_grid(u.levelmax())).size(2);
	sum /= N;
	sum2 /= N;

	return sqrt(sum2-sum*sum);
}


/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/

region_generator_plugin *the_region_generator;
RNG_plugin *the_random_number_generator;

int main (int argc, const char * argv[]) 
{
	const unsigned nbnd = 4;
	
	unsigned lbase, lmax, lbaseTF;
	double   err = 1.0;
	
	//------------------------------------------------------------------------------
	//... parse command line options
	//------------------------------------------------------------------------------
	
	splash();
	if( argc != 2 ){
		std::cout << " This version is compiled with the following plug-ins:\n";
		
		print_region_generator_plugins();
		print_transfer_function_plugins();
		print_RNG_plugins();
		print_output_plugins();
		
		std::cerr << "\n In order to run, you need to specify a parameter file!\n\n";
		exit(0);
	}
	
	//------------------------------------------------------------------------------
	//... open log file
	//------------------------------------------------------------------------------

	char logfname[128];
	sprintf(logfname,"%s_log.txt",argv[1]);
	MUSIC::log::setOutput(logfname);
	time_t ltime=time(NULL);
	LOGINFO("Opening log file \'%s\'.",logfname);
	LOGUSER("Running %s, version %s",THE_CODE_NAME,THE_CODE_VERSION);
	LOGUSER("Log is for run started %s",asctime( localtime(&ltime) ));
	
#ifdef FFTW3
	LOGUSER("Code was compiled using FFTW version 3.x");
#else
	LOGUSER("Code was compiled using FFTW version 2.x");
#endif
	
#ifdef SINGLETHREAD_FFTW
	LOGUSER("Code was compiled for single-threaded FFTW");
#else
	LOGUSER("Code was compiled for multi-threaded FFTW");
	LOGUSER("Running with a maximum of %d OpenMP threads", omp_get_max_threads() );
#endif
	
#ifdef SINGLE_PRECISION
	LOGUSER("Code was compiled for single precision.");
#else
	LOGUSER("Code was compiled for double precision.");
#endif
	
	
	//------------------------------------------------------------------------------
	//... read and interpret config file
	//------------------------------------------------------------------------------
	config_file cf(argv[1]);
	std::string tfname,randfname,temp;
	bool force_shift(false);
	double boxlength;
	
	//------------------------------------------------------------------------------
	//... initialize some parameters about grid set-up
	//------------------------------------------------------------------------------
	
	boxlength           = cf.getValue<double>( "setup", "boxlength" );
	lbase				= cf.getValue<unsigned>( "setup", "levelmin" );
	lmax				= cf.getValue<unsigned>( "setup", "levelmax" );
	lbaseTF				= cf.getValueSafe<unsigned>( "setup", "levelmin_TF", lbase );
	
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
	
    // .. determine if spectral sampling should be used
    if( !cf.containsKey( "setup", "kspace_TF" ))
        cf.insertValue( "setup", "kspace_TF", "yes");
    
    bool bspectral_sampling = cf.getValue<bool>( "setup", "kspace_TF" );
    
	if( bspectral_sampling )
	  LOGINFO("Using k-space sampled transfer functions...");
	else
	  LOGINFO("Using real space sampled transfer functions...");
		
	//------------------------------------------------------------------------------
	//... initialize multithread FFTW
	//------------------------------------------------------------------------------
	
#if not defined(SINGLETHREAD_FFTW)
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
	fftwf_init_threads();
	fftwf_plan_with_nthreads(omp_get_max_threads());
	#else
	fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());
	#endif
#else
	fftw_threads_init();
#endif
#endif
	
	//------------------------------------------------------------------------------
	//... initialize cosmology
	//------------------------------------------------------------------------------
	bool 
		do_baryons	= cf.getValue<bool>("setup","baryons"),
		do_2LPT		= cf.getValueSafe<bool>("setup","use_2LPT",false),
		do_LLA		= cf.getValueSafe<bool>("setup","use_LLA",false),
		do_CVM		= cf.getValueSafe<bool>("setup","center_vel",false),
        do_CVMfrompos = false;
	
	transfer_function_plugin *the_transfer_function_plugin
		= select_transfer_function_plugin( cf );
    
	cosmology cosmo( cf );
	
	std::cout << " - starting at a=" << cosmo.astart << std::endl;
	
	CosmoCalc ccalc(cosmo,the_transfer_function_plugin);
	cosmo.pnorm	= ccalc.ComputePNorm( 2.0*M_PI/boxlength );
	cosmo.dplus	= ccalc.CalcGrowthFactor( cosmo.astart )/ccalc.CalcGrowthFactor( 1.0 );
	cosmo.vfact = ccalc.ComputeVFact( cosmo.astart );
	
	if( !the_transfer_function_plugin->tf_has_total0() )
	        cosmo.pnorm *= cosmo.dplus*cosmo.dplus;
	
	//... directly use the normalisation via a parameter rather than the calculated one
	cosmo.pnorm = cf.getValueSafe<double>("setup","force_pnorm",cosmo.pnorm);
	

	double vfac2lpt = 1.0;
	
	if( the_transfer_function_plugin->tf_velocity_units() && do_baryons )
	{
		vfac2lpt = cosmo.vfact; // if the velocities are in velocity units, we need to divide by vfact for the 2lPT term
		cosmo.vfact = 1.0;
	}

	
	//
	{
		char tmpstr[128];
		sprintf(tmpstr,"%.12g",cosmo.pnorm);
		cf.insertValue("cosmology","pnorm",tmpstr);
		sprintf(tmpstr,"%.12g",cosmo.dplus);
		cf.insertValue("cosmology","dplus",tmpstr);
		sprintf(tmpstr,"%.12g",cosmo.vfact);
		cf.insertValue("cosmology","vfact",tmpstr);
		
	}
	
    the_region_generator = select_region_generator_plugin( cf );
  
    the_random_number_generator = select_RNG_plugin( cf );
	//------------------------------------------------------------------------------
	//... determine run parameters
	//------------------------------------------------------------------------------
	
	if( !the_transfer_function_plugin->tf_is_distinct() && do_baryons )
		std::cout	<< " - WARNING: The selected transfer function does not support\n"
		<< "            distinct amplitudes for baryon and DM fields!\n"
		<< "            Perturbation amplitudes will be identical!" << std::endl;
	

	double kickfac = 0.0, kickfac_2LPT = 0.0;
    double velocity_restframe[3] = {0.0,0.0,0.0};
    
	if( do_CVM )
	{
		double ztarget = 0.0;//cf.getValueSafe<double>("setup","center_vel_zfinal",0.0);
		kickfac = ccalc.ComputeVelocityCompensation( cosmo.astart, 1./(1.+ztarget) )/cosmo.vfact;
		kickfac_2LPT = ccalc.ComputeVelocityCompensation_2LPT( cosmo.astart, 1./(1.+ztarget) )/cosmo.vfact;
		std::cout	<< " - Will center velocities for target redshift " << ztarget << std::endl;
		LOGUSER("Will center velocities for target redshift %f.",ztarget);
        
        double atarget = 1.0;
        
        if( cf.containsKey("setup", "center_vel_posfinal" ) )
        {
            double CVM_x0final[3], CVM_x0initial[3], dx[3];
            
            do_CVMfrompos = true;
            atarget = cf.getValueSafe("setup", "center_vel_afinal", atarget );
            
            std::string strtemp;
            strtemp            = cf.getValue<std::string>( "setup", "center_vel_posfinal" );
            sscanf( strtemp.c_str(), "%lf,%lf,%lf", &CVM_x0final[0], &CVM_x0final[1], &CVM_x0final[2] );
            
            kickfac = ccalc.CalcTFac( cosmo.astart, atarget );
            
            LOGINFO("Will center velocities for target position\n     (%f,%f,%f) at a=%f.",CVM_x0final[0],CVM_x0final[1],CVM_x0final[2],atarget);
            
            
            LOGINFO("kick factor is %f",kickfac);
            // initial center definition goes here:
            the_region_generator->get_center_unshifted( CVM_x0initial );
            
            LOGINFO("Will center velocities for initial position\n     (%f,%f,%f) at a=%f.",CVM_x0initial[0],CVM_x0initial[1],CVM_x0initial[2],cosmo.astart);
            
            dx[0] = (CVM_x0final[0] - CVM_x0initial[0]);
            dx[1] = (CVM_x0final[1] - CVM_x0initial[1]);
            dx[2] = (CVM_x0final[2] - CVM_x0initial[2]);
            
            velocity_restframe[0] = -dx[0] * kickfac;
            velocity_restframe[1] = -dx[1] * kickfac;
            velocity_restframe[2] = -dx[2] * kickfac;
            
            LOGINFO("Setting initial velocity rest frame to (%f,%f,%f)", velocity_restframe[0],velocity_restframe[1],velocity_restframe[2]);
            
            
        }else{
            LOGINFO("Will zero initial velocities, no target center specified.");
            
        }
        
	}
	
	
	
	//------------------------------------------------------------------------------
	//... determine the refinement hierarchy
	//------------------------------------------------------------------------------
	
    refinement_hierarchy rh_Poisson( cf );
	store_grid_structure(cf, rh_Poisson);
	//rh_Poisson.output();
	print_hierarchy_stats( cf, rh_Poisson );
	
	refinement_hierarchy rh_TF( rh_Poisson );
	modify_grid_for_TF( rh_Poisson, rh_TF, cf );
	//rh_TF.output();

	LOGUSER("Grid structure for Poisson solver:");
	rh_Poisson.output_log();
	LOGUSER("Grid structure for density convolution:");
	rh_TF.output_log();
	
	//------------------------------------------------------------------------------
	//... initialize the output plug-in
	//------------------------------------------------------------------------------
	std::string outformat, outfname;
	outformat			= cf.getValue<std::string>( "output", "format" );
	outfname			= cf.getValue<std::string>( "output", "filename" );
	output_plugin *the_output_plugin = select_output_plugin( cf );
	
	//------------------------------------------------------------------------------
	//... initialize the random numbers
	//------------------------------------------------------------------------------
	std::cout << "=============================================================\n";
	std::cout << "   GENERATING WHITE NOISE\n";
	std::cout << "-------------------------------------------------------------\n";
	LOGUSER("Computing white noise...");
	rand_gen rand( cf, rh_TF, the_transfer_function_plugin );
	
	//------------------------------------------------------------------------------
	//... initialize the Poisson solver
	//------------------------------------------------------------------------------
	bool bdefd	= cf.getValueSafe<bool> ( "poisson" , "fft_fine", true );
	bool bglass = cf.getValueSafe<bool>("output","glass", false);
	bool bsph	= cf.getValueSafe<bool>("setup","do_SPH",false) && do_baryons;
	bool bbshift= bsph && !bglass;
	
	bool kspace	= cf.getValueSafe<bool>( "poisson", "kspace", false );
        bool kspace2LPT = kspace;

	bool decic_DM = cf.getValueSafe<bool>( "output", "glass_cicdeconvolve", false );
	bool decic_baryons = cf.getValueSafe<bool>( "output", "glass_cicdeconvolve", false ) & bsph;

	//... if in unigrid mode, use k-space instead of hybrid
	if(bdefd && (lbase==lmax))
	{
		kspace=true;
		bdefd=false;
		kspace2LPT=false;
	}
	
	std::string poisson_solver_name;
	if( kspace )
		poisson_solver_name = std::string("fft_poisson");
	else
		poisson_solver_name = std::string("mg_poisson");
	
	unsigned grad_order = cf.getValueSafe<unsigned> ( "poisson" , "grad_order", 4 );
	
	
	
	
	//... switch off if using kspace anyway
	//bdefd &= !kspace;
	
	poisson_plugin_creator *the_poisson_plugin_creator = get_poisson_plugin_map()[ poisson_solver_name ];
	poisson_plugin *the_poisson_solver = the_poisson_plugin_creator->create( cf );
	
	//---------------------------------------------------------------------------------
	//... THIS IS THE MAIN DRIVER BRANCHING TREE RUNNING THE VARIOUS PARTS OF THE CODE
	//---------------------------------------------------------------------------------
	bool bfatal = false;
	try{
		if( ! do_2LPT )
		{
			LOGUSER("Entering 1LPT branch");
			
			//------------------------------------------------------------------------------
			//... cdm density and displacements
			//------------------------------------------------------------------------------
			std::cout << "=============================================================\n";
			std::cout << "   COMPUTING DARK MATTER DISPLACEMENTS\n";
			std::cout << "-------------------------------------------------------------\n";
			LOGUSER("Computing dark matter displacements...");
			
			grid_hierarchy f( nbnd );//, u(nbnd);
			tf_type my_tf_type = cdm;
			if( !do_baryons || !the_transfer_function_plugin->tf_is_distinct() )
				my_tf_type = total;
			
			
			GenerateDensityHierarchy(	cf, the_transfer_function_plugin, my_tf_type , rh_TF, rand, f, false, false );
			coarsen_density(rh_Poisson, f, bspectral_sampling);
            f.add_refinement_mask( rh_Poisson.get_coord_shift() );
            
			normalize_density(f);
			
			LOGUSER("Writing CDM data");
			the_output_plugin->write_dm_mass(f);
			the_output_plugin->write_dm_density(f);
			
			grid_hierarchy u( f );	u.zero();
			err = the_poisson_solver->solve(f, u);
			
			if(!bdefd)
				f.deallocate();	
			
			LOGUSER("Writing CDM potential");
			the_output_plugin->write_dm_potential(u);
			
			
			//------------------------------------------------------------------------------
			//... DM displacements
			//------------------------------------------------------------------------------
			{
				grid_hierarchy data_forIO(u);
				for( int icoord = 0; icoord < 3; ++icoord )
				{
					if( bdefd )
					{
						data_forIO.zero();
						*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
						poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order,
							       data_forIO.levelmin()==data_forIO.levelmax(), decic_DM );					
						*data_forIO.get_grid(data_forIO.levelmax()) /= 1<<f.levelmax();
						the_poisson_solver->gradient_add(icoord, u, data_forIO );
						
					}
					else
						//... displacement
						the_poisson_solver->gradient(icoord, u, data_forIO );

                    coarsen_density( rh_Poisson, data_forIO, false );
					LOGUSER("Writing CDM displacements");
					the_output_plugin->write_dm_position(icoord, data_forIO );
				}
				if( do_baryons )
					u.deallocate();
				data_forIO.deallocate();
			}
			
			
			//------------------------------------------------------------------------------
			//... gas density
			//------------------------------------------------------------------------------
			if( do_baryons )
			{
				std::cout << "=============================================================\n";
				std::cout << "   COMPUTING BARYON DENSITY\n";
				std::cout << "-------------------------------------------------------------\n";
				LOGUSER("Computing baryon density...");
				GenerateDensityHierarchy(	cf, the_transfer_function_plugin, baryon , rh_TF, rand, f, false, bbshift );
				coarsen_density(rh_Poisson, f, bspectral_sampling);
                f.add_refinement_mask( rh_Poisson.get_coord_shift() );
				normalize_density(f);
				
				if( !do_LLA )
				{	
					LOGUSER("Writing baryon density");
					the_output_plugin->write_gas_density(f);
				}
				
				if( bsph )
				{
					u = f;	u.zero();
					err = the_poisson_solver->solve(f, u);					
					
					if(!bdefd)
						f.deallocate();
					
					grid_hierarchy data_forIO(u);
					for( int icoord = 0; icoord < 3; ++icoord )
					{
						if( bdefd )
						{
							data_forIO.zero();
							*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
							poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order, 
								       data_forIO.levelmin()==data_forIO.levelmax(), decic_baryons);					
							*data_forIO.get_grid(data_forIO.levelmax()) /= 1<<f.levelmax();
							the_poisson_solver->gradient_add(icoord, u, data_forIO );
							
						}
						else
							//... displacement
							the_poisson_solver->gradient(icoord, u, data_forIO );
						
						coarsen_density( rh_Poisson, data_forIO, false );
                        LOGUSER("Writing baryon displacements");
						the_output_plugin->write_gas_position(icoord, data_forIO );
						
					}	
					u.deallocate();
					data_forIO.deallocate();
					if( bdefd )
						f.deallocate();
				}
				else if( do_LLA )
				{
					u = f;	u.zero();
					err = the_poisson_solver->solve(f, u);
					compute_LLA_density( u, f,grad_order );
					u.deallocate();
					normalize_density(f);
					LOGUSER("Writing baryon density");
					the_output_plugin->write_gas_density(f);
				}
				
				f.deallocate();
			}
			
			
			
			//------------------------------------------------------------------------------
			//... velocities
			//------------------------------------------------------------------------------
			if( (!the_transfer_function_plugin->tf_has_velocities() || !do_baryons) && !bsph )
			{	
				std::cout << "=============================================================\n";
				std::cout << "   COMPUTING VELOCITIES\n";
				std::cout << "-------------------------------------------------------------\n";
				LOGUSER("Computing velocitites...");
				
				if( do_baryons )
				{
					GenerateDensityHierarchy(	cf, the_transfer_function_plugin, total , rh_TF, rand, f, false, false );
					coarsen_density(rh_Poisson, f, bspectral_sampling);
					f.add_refinement_mask( rh_Poisson.get_coord_shift() );
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
						poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order, 
							       data_forIO.levelmin()==data_forIO.levelmax(), decic_baryons );					
						*data_forIO.get_grid(data_forIO.levelmax()) /= 1<<f.levelmax();
						the_poisson_solver->gradient_add(icoord, u, data_forIO );
					}
					else 
						the_poisson_solver->gradient(icoord, u, data_forIO );
					
					
					
					//... multiply to get velocity
					data_forIO *= cosmo.vfact;
					
					//... velocity kick to keep refined region centered?
					if(do_CVM)
					{	
					    
						double ukick = 0.0;
                        if( !do_CVMfrompos )
                            ukick = compute_finest_mean(data_forIO);
                        else{
                            ukick = compute_finest_mean(data_forIO);
                            LOGINFO("mean of %c-velocity of high-res particles is %f",'x'+icoord,ukick);
                            ukick = -velocity_restframe[icoord];
                        }
						data_forIO -= ukick;
					}
					double sigv = compute_finest_sigma( data_forIO );
					LOGINFO("sigma of %c-velocity of high-res particles is %f",'x'+icoord, sigv);

					coarsen_density( rh_Poisson, data_forIO, false );
					LOGUSER("Writing CDM velocities");
					the_output_plugin->write_dm_velocity(icoord, data_forIO);

					if( do_baryons )
					{	
						LOGUSER("Writing baryon velocities");
						the_output_plugin->write_gas_velocity(icoord, data_forIO);
					}
				
				}
				
				u.deallocate();
				data_forIO.deallocate();
				
			}
			else
			{
				LOGINFO("Computing separate velocities for CDM and baryons:");
				std::cout << "=============================================================\n";
				std::cout << "   COMPUTING DARK MATTER VELOCITIES\n";
				std::cout << "-------------------------------------------------------------\n";
				LOGUSER("Computing dark matter velocitites...");
				
				//... we do baryons and have velocity transfer functions, or we do SPH and not to shift
				//... do DM first
				GenerateDensityHierarchy(	cf, the_transfer_function_plugin, vcdm , rh_TF, rand, f, false, false );
				coarsen_density(rh_Poisson, f, bspectral_sampling);
				f.add_refinement_mask( rh_Poisson.get_coord_shift() );
                normalize_density(f);
				
				u = f;	u.zero();
				
				err = the_poisson_solver->solve(f, u);
				
				if(!bdefd)
					f.deallocate();
				
				double uref[3];
				
				grid_hierarchy data_forIO(u);
				for( int icoord = 0; icoord < 3; ++icoord )
				{
					//... displacement
					if(bdefd)
					{
						data_forIO.zero();
						*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
						poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order, 
							       data_forIO.levelmin()==data_forIO.levelmax(), decic_DM );					
						*data_forIO.get_grid(data_forIO.levelmax()) /= 1<<f.levelmax();
						the_poisson_solver->gradient_add(icoord, u, data_forIO );
					}
					else 
						the_poisson_solver->gradient(icoord, u, data_forIO );
					
					//... multiply to get velocity
					data_forIO *= cosmo.vfact;
					
					//... velocity kick to keep refined region centered?
					if(do_CVM)
					{
					    if( !do_CVMfrompos )
                            uref[icoord] = compute_finest_mean(data_forIO);
                        else{
                            double umean = compute_finest_mean(data_forIO);
                            LOGINFO("mean of %c-velocity of high-res particles is %f",'x'+icoord,umean);
                            
                            uref[icoord] = -velocity_restframe[icoord];
                        }
						data_forIO -= uref[icoord];
					}
				
					double sigv = compute_finest_sigma( data_forIO );
					LOGINFO("sigma of %c-velocity of high-res DM is %f",'x'+icoord, sigv);

					coarsen_density( rh_Poisson, data_forIO, false );
					LOGUSER("Writing CDM velocities");
					the_output_plugin->write_dm_velocity(icoord, data_forIO);
				}
				u.deallocate();
				data_forIO.deallocate();
				f.deallocate();
				
				
				std::cout << "=============================================================\n";
				std::cout << "   COMPUTING BARYON VELOCITIES\n";
				std::cout << "-------------------------------------------------------------\n";
				LOGUSER("Computing baryon velocitites...");
				//... do baryons
				GenerateDensityHierarchy(	cf, the_transfer_function_plugin, vbaryon , rh_TF, rand, f, false, bbshift );
				coarsen_density(rh_Poisson, f, bspectral_sampling);
				f.add_refinement_mask( rh_Poisson.get_coord_shift() );
                normalize_density(f);
				
				u = f;	u.zero();
				
				err = the_poisson_solver->solve(f, u);
				
				if(!bdefd)
					f.deallocate();
				
				data_forIO = u;
				for( int icoord = 0; icoord < 3; ++icoord )
				{
					//... displacement
					if(bdefd)
					{
						data_forIO.zero();
						*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
						poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order, 
							       data_forIO.levelmin()==data_forIO.levelmax(), decic_baryons );					
						*data_forIO.get_grid(data_forIO.levelmax()) /= 1<<f.levelmax();
						the_poisson_solver->gradient_add(icoord, u, data_forIO );
					}
					else 
						the_poisson_solver->gradient(icoord, u, data_forIO );
					
					//... multiply to get velocity
					data_forIO *= cosmo.vfact;
					
					//... velocity kick to keep refined region centered?
					if( do_CVM )
						data_forIO -= uref[icoord];
					
					double sigv = compute_finest_sigma( data_forIO );
                    LOGINFO("sigma of %c-velocity of high-res baryons is %f",'x'+icoord, sigv);
					
					coarsen_density( rh_Poisson, data_forIO, false );
					LOGUSER("Writing baryon velocities");
					the_output_plugin->write_gas_velocity(icoord, data_forIO);
				}
				u.deallocate();
				f.deallocate();
				data_forIO.deallocate();
			}
		/*********************************************************************************************/
		/*********************************************************************************************/
		/*** 2LPT ************************************************************************************/
		/*********************************************************************************************/
		}else {
			//.. use 2LPT ...
			LOGUSER("Entering 2LPT branch");
			
			grid_hierarchy f( nbnd ), u1(nbnd), u2LPT(nbnd), f2LPT( nbnd );
			
			
			
			tf_type my_tf_type = vcdm;
			bool dm_only = !do_baryons;
			if( !do_baryons || !the_transfer_function_plugin->tf_has_velocities() )
				my_tf_type = total;
			
			std::cout << "=============================================================\n";
			if( my_tf_type == total )
			{
				std::cout << "   COMPUTING VELOCITIES\n";
				LOGUSER("Computing velocities...");				
			}else{
				std::cout << "   COMPUTING DARK MATTER VELOCITIES\n";
				LOGUSER("Computing dark matter velocities...");	
			}
			std::cout << "-------------------------------------------------------------\n";	

			
			GenerateDensityHierarchy(	cf, the_transfer_function_plugin, my_tf_type , rh_TF, rand, f, false, false );
			coarsen_density(rh_Poisson, f, bspectral_sampling);
			f.add_refinement_mask( rh_Poisson.get_coord_shift() );
            normalize_density(f);
			
			if( dm_only )
			{
				the_output_plugin->write_dm_density(f);
				the_output_plugin->write_dm_mass(f);	
			}
			
			u1 = f;	u1.zero();
			
			//... compute 1LPT term
			err = the_poisson_solver->solve(f, u1);
			
			
			//... compute 2LPT term
			if(bdefd)
				f2LPT=f;
			else
				f.deallocate();
		
			LOGINFO("Computing 2LPT term....");
			if( !kspace2LPT )
				compute_2LPT_source(u1, f2LPT, grad_order );
			else{
				LOGUSER("computing term using FFT");
				compute_2LPT_source_FFT(cf, u1, f2LPT);
			}
			LOGINFO("Solving 2LPT Poisson equation");
			u2LPT = u1; u2LPT.zero();
			err = the_poisson_solver->solve(f2LPT, u2LPT);
			
			//... if doing the hybrid step, we need a combined source term
			if( bdefd )
			{
				f2LPT*=6.0/7.0/vfac2lpt;
				f+=f2LPT;
				
				if( !dm_only )
					f2LPT.deallocate();
			}
			
			//... add the 2LPT contribution
			u2LPT *= 6.0/7.0/vfac2lpt;
			u1 += u2LPT;
			
			double uref_2LPT[3];
			
			if( !dm_only )
			{
				if( do_CVM )
				{
					grid_hierarchy data_forIO(u1);
					for( int icoord = 0; icoord < 3; ++icoord )
					{
						the_poisson_solver->gradient(icoord, u2LPT, data_forIO );
						uref_2LPT[icoord] = kickfac * compute_finest_mean(data_forIO);
					}
					data_forIO.deallocate();
				}
				u2LPT.deallocate();	
			}
			
			double uref[3];
			
			grid_hierarchy data_forIO(u1);
			for( int icoord = 0; icoord < 3; ++icoord )
			{
				if(bdefd)
				{
					data_forIO.zero();
					*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
					poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order, 
						       data_forIO.levelmin()==data_forIO.levelmax(), decic_DM );
					*data_forIO.get_grid(data_forIO.levelmax()) /= (1<<f.levelmax());
					the_poisson_solver->gradient_add(icoord, u1, data_forIO );
				}
				else 
					the_poisson_solver->gradient(icoord, u1, data_forIO );
				
				data_forIO *= cosmo.vfact;
				
				//... velocity kick to keep refined region centered?
#if 0
				if( do_CVM )
				{
					/*uref[icoord] = kickfac * compute_finest_mean(data_forIO);
					//std::cerr << "uref_old[" << icoord << "] = " << uref[icoord] << std::endl;
					uref[icoord] = kickfac * compute_finest_mean(data_forIO) - 6./7. * uref_2LPT[icoord];
					uref_2LPT[icoord] = 6./7. * kickfac_2LPT/kickfac * uref_2LPT[icoord];
					uref[icoord] += uref_2LPT[icoord];
					//std::cerr << "uref_new[" << icoord << "] = " << uref[icoord] << std::endl; */
					uref[icoord] = compute_finest_mean(data_forIO);

					data_forIO -= uref[icoord];
				}
#endif
                if(do_CVM)
                {
                    if( !do_CVMfrompos )
                        uref[icoord] = compute_finest_mean(data_forIO);
                    else
                        uref[icoord] = -velocity_restframe[icoord];
                    
                    data_forIO -= uref[icoord];
                }
				
					
				double sigv = compute_finest_sigma( data_forIO );
				std::cerr << " - velocity component " << icoord << " : sigma = " << sigv << std::endl;
				
				coarsen_density( rh_Poisson, data_forIO, false );
                LOGUSER("Writing CDM velocities");
				the_output_plugin->write_dm_velocity(icoord, data_forIO);					
				
				if( do_baryons && !the_transfer_function_plugin->tf_has_velocities() && !bsph)
				{	
					LOGUSER("Writing baryon velocities");
					the_output_plugin->write_gas_velocity(icoord, data_forIO);				
				}
			}
			data_forIO.deallocate();
			if( !dm_only )
				u1.deallocate();
			
			
			if( do_baryons && (the_transfer_function_plugin->tf_has_velocities() || bsph) )
			{
				std::cout << "=============================================================\n";
				std::cout << "   COMPUTING BARYON VELOCITIES\n";
				std::cout << "-------------------------------------------------------------\n";
				LOGUSER("Computing baryon displacements...");
				
				GenerateDensityHierarchy(	cf, the_transfer_function_plugin, vbaryon , rh_TF, rand, f, false, bbshift );
				coarsen_density(rh_Poisson, f, bspectral_sampling);
				f.add_refinement_mask( rh_Poisson.get_coord_shift() );
                normalize_density(f);
				
				u1 = f;	u1.zero();
				
				if(bdefd)
					f2LPT=f;
				
				//... compute 1LPT term
				err = the_poisson_solver->solve(f, u1);

				LOGINFO("Writing baryon potential");
				the_output_plugin->write_gas_potential(u1);
				
				//... compute 2LPT term
				u2LPT = f; u2LPT.zero();
				
				if( !kspace2LPT )
					compute_2LPT_source(u1, f2LPT, grad_order );
				else
					compute_2LPT_source_FFT(cf, u1, f2LPT);
				
				
				err = the_poisson_solver->solve(f2LPT, u2LPT);
				
				//... if doing the hybrid step, we need a combined source term
				if( bdefd )
				{
					f2LPT*=6.0/7.0/vfac2lpt;
					f+=f2LPT;
					
					f2LPT.deallocate();
				}
				
				//... add the 2LPT contribution
				u2LPT *= 6.0/7.0/vfac2lpt;
				u1 += u2LPT;
				u2LPT.deallocate();
				
				//grid_hierarchy data_forIO(u1);
				data_forIO = u1;
				for( int icoord = 0; icoord < 3; ++icoord )
				{
					if(bdefd)
					{
						data_forIO.zero();
						*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
						poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order, 
							       data_forIO.levelmin()==data_forIO.levelmax(), decic_baryons );					
						*data_forIO.get_grid(data_forIO.levelmax()) /= (1<<f.levelmax());
						the_poisson_solver->gradient_add(icoord, u1, data_forIO );
					}
					else 
						the_poisson_solver->gradient(icoord, u1, data_forIO );
					
					data_forIO *= cosmo.vfact;
					
					//... velocity kick to keep refined region centered?
					if( do_CVM )
						data_forIO -= uref[icoord];
					
					double sigv = compute_finest_sigma( data_forIO );
					std::cerr << " - velocity component " << icoord << " : sigma = " << sigv << std::endl;
					
					coarsen_density( rh_Poisson, data_forIO, false );
					LOGUSER("Writing baryon velocities");
					the_output_plugin->write_gas_velocity(icoord, data_forIO);				
				}
				data_forIO.deallocate();
				u1.deallocate();
			}
			
			
			std::cout << "=============================================================\n";
			std::cout << "   COMPUTING DARK MATTER DISPLACEMENTS\n";
			std::cout << "-------------------------------------------------------------\n";
			LOGUSER("Computing dark matter displacements...");
			
			//... if baryons are enabled, the displacements have to be recomputed
			//... otherwise we can compute them directly from the velocities
			if( !dm_only )
			{
				// my_tf_type is cdm if do_baryons==true, total otherwise
				my_tf_type = cdm;
				if( !do_baryons || !the_transfer_function_plugin->tf_is_distinct() )
					my_tf_type = total;
				
				GenerateDensityHierarchy(	cf, the_transfer_function_plugin, my_tf_type , rh_TF, rand, f, false, false );
				coarsen_density(rh_Poisson, f, bspectral_sampling);
				f.add_refinement_mask( rh_Poisson.get_coord_shift() );
                normalize_density(f);
				
				LOGUSER("Writing CDM data");
				the_output_plugin->write_dm_density(f);
				the_output_plugin->write_dm_mass(f);
				u1 = f;	u1.zero();
				
				if(bdefd)
					f2LPT=f;
				
				//... compute 1LPT term
				err = the_poisson_solver->solve(f, u1);
				
				//... compute 2LPT term
				u2LPT = f; u2LPT.zero();
				
				if( !kspace2LPT )
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
				/*f-=f2LPT;
				the_output_plugin->write_dm_density(f);
				the_output_plugin->write_dm_mass(f);
				f+=f2LPT;*/
				
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
					poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order, 
						       data_forIO.levelmin()==data_forIO.levelmax(), decic_DM );
					*data_forIO.get_grid(data_forIO.levelmax()) /= 1<<f.levelmax();
					the_poisson_solver->gradient_add(icoord, u1, data_forIO );
				}
				else 
					the_poisson_solver->gradient(icoord, u1, data_forIO );
				
				coarsen_density( rh_Poisson, data_forIO, false );
                LOGUSER("Writing CDM displacements");
				the_output_plugin->write_dm_position(icoord, data_forIO );	
			}
			
			data_forIO.deallocate();
			u1.deallocate();
			

			if( do_baryons && !bsph )
			{	
				std::cout << "=============================================================\n";
				std::cout << "   COMPUTING BARYON DENSITY\n";
				std::cout << "-------------------------------------------------------------\n";
				LOGUSER("Computing baryon density...");
				
				GenerateDensityHierarchy(	cf, the_transfer_function_plugin, baryon , rh_TF, rand, f, true, false );
				coarsen_density(rh_Poisson, f, bspectral_sampling);
				f.add_refinement_mask( rh_Poisson.get_coord_shift() );
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
					
					if( !kspace2LPT )
						compute_2LPT_source(u1, f2LPT, grad_order );
					else
						compute_2LPT_source_FFT(cf, u1, f2LPT);
					
					err = the_poisson_solver->solve(f2LPT, u2LPT);
					u2LPT *= 3.0/7.0;
					u1 += u2LPT;
					u2LPT.deallocate();
					
					compute_LLA_density( u1, f, grad_order );
					normalize_density(f);
					
					LOGUSER("Writing baryon density");
					the_output_plugin->write_gas_density(f);
				}
			}
			else if( do_baryons && bsph )
			{
				std::cout << "=============================================================\n";
				std::cout << "   COMPUTING BARYON DISPLACEMENTS\n";
				std::cout << "-------------------------------------------------------------\n";
				LOGUSER("Computing baryon displacements...");
				
				GenerateDensityHierarchy(	cf, the_transfer_function_plugin, baryon , rh_TF, rand, f, false, bbshift );
				coarsen_density(rh_Poisson, f, bspectral_sampling);
				f.add_refinement_mask( rh_Poisson.get_coord_shift() );
                normalize_density(f);
				
				LOGUSER("Writing baryon density");
				the_output_plugin->write_gas_density(f);
				u1 = f;	u1.zero();
				
				if(bdefd)
					f2LPT=f;
				
				//... compute 1LPT term
				err = the_poisson_solver->solve(f, u1);
				
				//... compute 2LPT term
				u2LPT = f; u2LPT.zero();
				
				if( !kspace2LPT )
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
				
				data_forIO = u1;
				
				for( int icoord = 0; icoord < 3; ++icoord )
				{
					//... displacement
					if(bdefd)
					{
						data_forIO.zero();
						*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
						poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order, 
							       data_forIO.levelmin()==data_forIO.levelmax(), decic_baryons );
						*data_forIO.get_grid(data_forIO.levelmax()) /= 1<<f.levelmax();
						the_poisson_solver->gradient_add(icoord, u1, data_forIO );
					}
					else 
						the_poisson_solver->gradient(icoord, u1, data_forIO );
					
					coarsen_density( rh_Poisson, data_forIO, false );
					LOGUSER("Writing baryon displacements");
					the_output_plugin->write_gas_position(icoord, data_forIO );	
				}
			}

		}
		
		//------------------------------------------------------------------------------
		//... finish output
		//------------------------------------------------------------------------------
		
		the_output_plugin->finalize();
		delete the_output_plugin;
		
	}catch(std::runtime_error& excp){
		LOGERR("Fatal error occured. Code will exit:");
		LOGERR("Exception: %s",excp.what());
		std::cerr << " - " << excp.what() << std::endl;
		std::cerr << " - A fatal error occured. We need to exit...\n";
		bfatal = true;
	}

	std::cout << "=============================================================\n";
	
	

	if( !bfatal )
	{	
		std::cout << " - Wrote output file \'" << outfname << "\'\n     using plugin \'" << outformat << "\'...\n";
		LOGUSER("Wrote output file \'%s\'.",outfname.c_str());
	}
	
	//------------------------------------------------------------------------------
	//... clean up
	//------------------------------------------------------------------------------
	delete the_transfer_function_plugin;
	delete the_poisson_solver;

#if defined(FFTW3) and not defined(SINGLETHREAD_FFTW)
	#ifdef SINGLE_PRECISION
	fftwf_cleanup_threads();
	#else
	fftw_cleanup_threads();
	#endif
#endif
	
	
	//------------------------------------------------------------------------------
	//... we are done !
	//------------------------------------------------------------------------------
	std::cout << " - Done!" << std::endl << std::endl;
	
	ltime=time(NULL);
	
	LOGUSER("Run finished succesfully on %s",asctime( localtime(&ltime) ));
	
	cf.log_dump();
	
	
	return 0;
}
