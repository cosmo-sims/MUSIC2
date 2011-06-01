/*
 
 densities.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#include "densities.hh"
#include "convolution_kernel.hh"

//... THIS IS FOR TESTING, TODO: TAKE OUT FOR RELEASE VERSION
//... uncomment this to have a single peak in the centre and otherwise zeros
//#define SINGLE_PEAK
//#define SINGLE_OCT_PEAK
//#define OFF_OCT_PEAK
//#define DEGRADE_RAND



//TODO: this should be a larger number by default, just to maintain consistency with old default
#define DEF_RAN_CUBE_SIZE	32



/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/

void GenerateDensityUnigrid( config_file& cf, transfer_function *ptf, tf_type type, 
							refinement_hierarchy& refh, rand_gen& rand, grid_hierarchy& delta, bool smooth, bool shift )
{
	unsigned    levelmin,levelmax,levelminPoisson;
	real_t		boxlength;
	
	
	levelminPoisson	= cf.getValue<unsigned>("setup","levelmin");
	levelmin	= cf.getValueSafe<unsigned>("setup","levelmin_TF",levelminPoisson);
	levelmax	= cf.getValue<unsigned>("setup","levelmax");
	boxlength   = cf.getValue<real_t>( "setup", "boxlength" );
	
	bool kspace = cf.getValueSafe<unsigned>("setup","kspace_TF",false);
	
	unsigned	nbase	= 1<<levelmin;
	
	std::cerr << " - Running unigrid version\n";
	LOGUSER("Running unigrid density convolution...");
	
	//... select the transfer function to be used
	convolution::kernel_creator *the_kernel_creator;

	if( kspace )
	{
		std::cout << " - Using k-space transfer function kernel.\n";
		LOGUSER("Using k-space transfer function kernel.");
		
		#ifdef SINGLE_PRECISION	
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_k_float" ];
		#else
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_k_double" ];
		#endif
	}
	else
	{
		std::cout << " - Using real-space transfer function kernel.\n";
		LOGUSER("Using real-space transfer function kernel.");
		
		#ifdef SINGLE_PRECISION	
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_real_float" ];
		#else
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_real_double" ];
		#endif
	}	
		
	
	//... initialize convolution kernel
	convolution::kernel *the_tf_kernel = the_kernel_creator->create( cf, ptf, refh, type );

	//...
	std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmax << " ..." << std::endl;
	LOGUSER("Performing noise convolution on level %3d",levelmax);
	
	//... create convolution mesh
	DensityGrid<real_t> *top = new DensityGrid<real_t>( nbase, nbase, nbase );
	
	//... fill with random numbers
	rand.load( *top, levelmin );
	
#if defined(SINGLE_PEAK)
	top->zero();
	(*top)(top->size(0)/2, top->size(1)/2, top->size(2)/2) = 1.0;
#elif defined(SINGLE_OCT_PEAK)
	{
		top->zero();
		unsigned i0=top->size(0)/2, i1=top->size(1)/2, i2=top->size(2)/2;

		double weight = 1.0;
		(*top)(i0,i1,i2) = weight/8.0;			
		(*top)(i0+1,i1,i2) = weight/8.0;			
		(*top)(i0,i1+1,i2) = weight/8.0;			
		(*top)(i0,i1,i2+1) = weight/8.0;			
		(*top)(i0+1,i1+1,i2) = weight/8.0;			
		(*top)(i0+1,i1,i2+1) = weight/8.0;			
		(*top)(i0,i1+1,i2+1) = weight/8.0;			
		(*top)(i0+1,i1+1,i2+1) = weight/8.0;
	}
#endif
	
#if defined(OFF_OCT_PEAK)
	{
		top->zero();
		unsigned i0=top->size(0)/8, i1=top->size(1)/2, i2=top->size(2)/2;
		
		(*top)(i0,i1,i2) = 1.0/8.0;			
		(*top)(i0+1,i1,i2) = 1.0/8.0;			
		(*top)(i0,i1+1,i2) = 1.0/8.0;			
		(*top)(i0,i1,i2+1) = 1.0/8.0;			
		(*top)(i0+1,i1+1,i2) = 1.0/8.0;			
		(*top)(i0+1,i1,i2+1) = 1.0/8.0;			
		(*top)(i0,i1+1,i2+1) = 1.0/8.0;			
		(*top)(i0+1,i1+1,i2+1) = 1.0/8.0;
	}
#endif
		
	//... load convolution kernel
	the_tf_kernel->fetch_kernel( levelmin, false );
	
	//... perform convolution
	convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ), shift );
	
	//... clean up kernel
	delete the_tf_kernel;
	
	//... create multi-grid hierarchy
	delta.create_base_hierarchy(levelmin);
	
	//... copy convolved field to multi-grid hierarchy
	top->copy( *delta.get_grid(levelmin) );
	
	//... delete convolution grid
	delete top;
}


/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/

void GenerateDensityHierarchy(	config_file& cf, transfer_function *ptf, tf_type type, 
							  refinement_hierarchy& refh, rand_gen& rand, grid_hierarchy& delta, bool smooth, bool shift )
{
	unsigned					levelmin,levelmax,levelminPoisson;
	real_t						boxlength;
	std::vector<long>			rngseeds;
	std::vector<std::string>	rngfnames;
	bool						kspaceTF;
	
	double tstart, tend;
	tstart = omp_get_wtime();
	
	
	levelminPoisson	= cf.getValue<unsigned>("setup","levelmin");
	levelmin		= cf.getValueSafe<unsigned>("setup","levelmin_TF",levelminPoisson);
	levelmax		= cf.getValue<unsigned>("setup","levelmax");
	boxlength		= cf.getValue<real_t>( "setup", "boxlength" );
	kspaceTF		= cf.getValueSafe<bool>("setup", "kspace_TF", false);
	
	
	unsigned	nbase	= (unsigned)pow(2,levelmin);
	
	convolution::kernel_creator *the_kernel_creator;
	
	if( kspaceTF )
	{
		if( levelmin!=levelmax )
		{	
			LOGERR("K-space transfer function can only be used in unigrid density mode!");
			throw std::runtime_error("k-space transfer function can only be used in unigrid density mode");
			
		}
		
		std::cout << " - Using k-space transfer function kernel.\n";
		LOGUSER("Using k-space transfer function kernel.");
		
#ifdef SINGLE_PRECISION	
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_k_float" ];
#else
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_k_double" ];
#endif
	}
	else
	{
		std::cout << " - Using real-space transfer function kernel.\n";
		LOGUSER("Using real-space transfer function kernel.");
#ifdef SINGLE_PRECISION	
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_real_float" ];
#else
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_real_double" ];
#endif
	}	
	
	
	//... create and initialize density grids with white noise	
	PaddedDensitySubGrid<real_t>* coarse(NULL), *fine(NULL);
	
	DensityGrid<real_t>* top(NULL);
		
	convolution::kernel *the_tf_kernel = the_kernel_creator->create( cf, ptf, refh, type );
		
	/***** PERFORM CONVOLUTIONS *****/
	
	if( levelmax == levelmin )
	{
		std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmax << " ..." << std::endl;
		LOGUSER("Performing noise convolution on level %3d...",levelmax);
		
		top = new DensityGrid<real_t>( nbase, nbase, nbase );
		//rand_gen.load( *top, levelmin );
		rand.load( *top, levelmin );

#if defined(SINGLE_PEAK)
		top->zero();
		(*top)(top->size(0)/2, top->size(1)/2, top->size(2)/2) = 1.0;
#elif defined(SINGLE_OCT_PEAK)
		{
			std::cerr << ">>> setting single oct peak <<<\n";
			top->zero();
			unsigned i0=top->size(0)/2, i1=top->size(1)/2, i2=top->size(2)/2;
			
			(*top)(i0,i1,i2) = 1.0/8.0;			
			(*top)(i0+1,i1,i2) = 1.0/8.0;			
			(*top)(i0,i1+1,i2) = 1.0/8.0;			
			(*top)(i0,i1,i2+1) = 1.0/8.0;			
			(*top)(i0+1,i1+1,i2) = 1.0/8.0;			
			(*top)(i0+1,i1,i2+1) = 1.0/8.0;			
			(*top)(i0,i1+1,i2+1) = 1.0/8.0;			
			(*top)(i0+1,i1+1,i2+1) = 1.0/8.0;
		}

#endif
		
#if defined(OFF_OCT_PEAK)
		{
			top->zero();
			unsigned i0=top->size(0)/8, i1=top->size(1)/2, i2=top->size(2)/2;
			
			(*top)(i0,i1,i2) = 1.0/8.0;			
			(*top)(i0+1,i1,i2) = 1.0/8.0;			
			(*top)(i0,i1+1,i2) = 1.0/8.0;			
			(*top)(i0,i1,i2+1) = 1.0/8.0;			
			(*top)(i0+1,i1+1,i2) = 1.0/8.0;			
			(*top)(i0+1,i1,i2+1) = 1.0/8.0;			
			(*top)(i0,i1+1,i2+1) = 1.0/8.0;			
			(*top)(i0+1,i1+1,i2+1) = 1.0/8.0;
		}
#endif
				
		convolution::perform<real_t>( the_tf_kernel->fetch_kernel( levelmax ), reinterpret_cast<void*>( top->get_data_ptr() ), shift );
		the_tf_kernel->deallocate();
		
		delta.create_base_hierarchy(levelmin);
		top->copy( *delta.get_grid(levelmin) );
		delete top;
	}
		
	
	for( int i=0; i< (int)levelmax-(int)levelmin; ++i )
	{
		//.......................................................................................................//
		//... GENERATE/FILL WITH RANDOM NUMBERS .................................................................//
		//.......................................................................................................//
		
		
		if( i==0 )
		{
			top = new DensityGrid<real_t>( nbase, nbase, nbase );
			rand.load(*top,levelmin);
		}
		
		fine = new PaddedDensitySubGrid<real_t>( refh.offset(levelmin+i+1,0), refh.offset(levelmin+i+1,1), refh.offset(levelmin+i+1,2), 
												refh.size(levelmin+i+1,0), 	refh.size(levelmin+i+1,1), 	refh.size(levelmin+i+1,2) );
		rand.load(*fine,levelmin+i+1);
		
		//.......................................................................................................//
		//... PERFORM CONVOLUTIONS ..............................................................................//
		//.......................................................................................................//		
		if( i==0 )
		{
			/**********************************************************************************************************\
			 *	multi-grid: top-level grid grids .....
			 \**********************************************************************************************************/ 
			std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmin+i << " ..." << std::endl;
			LOGUSER("Performing noise convolution on level %3d",levelmin+i);
			
			LOGUSER("Creating base hierarchy...");
			delta.create_base_hierarchy(levelmin);
			
#if defined(SINGLE_PEAK) || defined(SINGLE_OCT_PEAK)
			{
				top->zero();
				(*top)(top->size(0)/2, top->size(1)/2, top->size(2)/2) = 1.0/pow(2,1.5*(levelmax-levelmin));
			}

#endif	
			
#if defined(OFF_OCT_PEAK)
			{
				top->zero();
				unsigned i0=top->size(0)/8, i1=top->size(1)/2, i2=top->size(2)/2;
				
				(*top)(i0,i1,i2) = 1.0/pow(2,1.5*(levelmax-levelmin));
			}
#endif
			
			DensityGrid<real_t> top_save( *top );

			the_tf_kernel->fetch_kernel( levelmin );
			
			//... 1) compute standard convolution for levelmin
			LOGUSER("Computing density self-contribution");
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ), shift );
			top->copy( *delta.get_grid(levelmin) );
			
			
			//... 2) compute contribution to finer grids from non-refined region
			LOGUSER("Computing long-range component for finer grid.");
			*top = top_save;
			top_save.clear();
			top->zero_subgrid(refh.offset(levelmin+i+1,0), refh.offset(levelmin+i+1,1), refh.offset(levelmin+i+1,2), 
							  refh.size(levelmin+i+1,0)/2, refh.size(levelmin+i+1,1)/2, refh.size(levelmin+i+1,2)/2 );

			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ), shift );
			the_tf_kernel->deallocate();
			
			meshvar_bnd delta_longrange( *delta.get_grid(levelmin) );
			top->copy( delta_longrange );
			delete top;			
			
			//... inject these contributions to the next level
			LOGUSER("Allocating refinement patch");
			LOGUSER("   offset=(%5d,%5d,%5d)",refh.offset(levelmin+1,0), refh.offset(levelmin+1,1), refh.offset(levelmin+1,2));
			LOGUSER("   size  =(%5d,%5d,%5d)",refh.size(levelmin+1,0), refh.size(levelmin+1,1), refh.size(levelmin+1,2));
			
			delta.add_patch( refh.offset(levelmin+1,0), refh.offset(levelmin+1,1), refh.offset(levelmin+1,2), 
							refh.size(levelmin+1,0), refh.size(levelmin+1,1), refh.size(levelmin+1,2) );
			
			LOGUSER("Injecting long range component");
			//mg_straight().prolong( delta_longrange, *delta.get_grid(levelmin+1) );
			//mg_cubic_mult().prolong( delta_longrange, *delta.get_grid(levelmin+1) );
			
			mg_cubic().prolong( delta_longrange, *delta.get_grid(levelmin+1) );
		}
		else
		{
			/**********************************************************************************************************\
			 *	multi-grid: intermediate sub-grids .....
			 \**********************************************************************************************************/ 
			
			std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmin+i << " ..." << std::endl;
			LOGUSER("Performing noise convolution on level %3d",levelmin+i);
			
			//... add new refinement patch
			LOGUSER("Allocating refinement patch");
			LOGUSER("   offset=(%5d,%5d,%5d)",refh.offset(levelmin+i+1,0), refh.offset(levelmin+i+1,1), refh.offset(levelmin+i+1,2));
			LOGUSER("   size  =(%5d,%5d,%5d)",refh.size(levelmin+i+1,0), refh.size(levelmin+i+1,1), refh.size(levelmin+i+1,2));
			
			delta.add_patch( refh.offset(levelmin+i+1,0), refh.offset(levelmin+i+1,1), refh.offset(levelmin+i+1,2), 
							refh.size(levelmin+i+1,0), refh.size(levelmin+i+1,1), refh.size(levelmin+i+1,2) );
			
			
			//... copy coarse grid long-range component to fine grid
			LOGUSER("Injecting long range component");
			//mg_straight().prolong( *delta.get_grid(levelmin+i), *delta.get_grid(levelmin+i+1) );
			mg_cubic().prolong( *delta.get_grid(levelmin+i), *delta.get_grid(levelmin+i+1) );
			
			PaddedDensitySubGrid<real_t> coarse_save( *coarse );
			the_tf_kernel->fetch_kernel( levelmin+i );
					
			//... 1) the inner region
			LOGUSER("Computing density self-contribution");
			coarse->subtract_boundary_oct_mean();
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()), shift );
			coarse->copy_add_unpad( *delta.get_grid(levelmin+i) );
			
			
			//... 2) the 'BC' for the next finer grid
			LOGUSER("Computing long-range component for finer grid.");
			*coarse = coarse_save;
			coarse->subtract_boundary_oct_mean();
			coarse->zero_subgrid(refh.offset(levelmin+i+1,0), refh.offset(levelmin+i+1,1), refh.offset(levelmin+i+1,2), 
								 refh.size(levelmin+i+1,0)/2, refh.size(levelmin+i+1,1)/2, refh.size(levelmin+i+1,2)/2 );
			
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( coarse->get_data_ptr() ), shift );
			
			//... interpolate to finer grid(s)
			meshvar_bnd delta_longrange( *delta.get_grid(levelmin+i) );
			coarse->copy_unpad( delta_longrange );
			
			LOGUSER("Injecting long range component");
			//mg_straight().prolong_add( delta_longrange, *delta.get_grid(levelmin+i+1) );
			mg_cubic().prolong_add( delta_longrange, *delta.get_grid(levelmin+i+1) );

			//... 3) the coarse-grid correction
			LOGUSER("Computing coarse grid correction");
			*coarse = coarse_save;
			coarse->subtract_oct_mean();
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()), shift );
			coarse->subtract_mean();
			coarse->upload_bnd_add( *delta.get_grid(levelmin+i-1) );
			
			//... clean up
			the_tf_kernel->deallocate();
			delete coarse;
		}
		
		
		coarse = fine;
	}
	
	//... and convolution for finest grid (outside loop)
	if( levelmax > levelmin )
	{
		/**********************************************************************************************************\
		 *	multi-grid: finest sub-grid .....
		 \**********************************************************************************************************/ 
		std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmax << " ..." << std::endl;
		LOGUSER("Performing noise convolution on level %3d",levelmax);
		
#if defined(SINGLE_PEAK) || defined(SINGLE_OCT_PEAK)
		{
			coarse->zero();
			
			int 
				i0 = (1<<(levelmax-1)) - refh.offset_abs(levelmax,0) + coarse->nx_/4,
				i1 = (1<<(levelmax-1)) - refh.offset_abs(levelmax,1) + coarse->nx_/4,
				i2 = (1<<(levelmax-1)) - refh.offset_abs(levelmax,2) + coarse->nx_/4;
#if defined(SINGLE_PEAK)
				(*coarse)(i0,i1,i2) = 1.0;
#elif defined(SINGLE_OCT_PEAK)
				(*coarse)(i0,i1,i2) = 1.0/8.0;			
				(*coarse)(i0+1,i1,i2) = 1.0/8.0;			
				(*coarse)(i0,i1+1,i2) = 1.0/8.0;			
				(*coarse)(i0,i1,i2+1) = 1.0/8.0;			
				(*coarse)(i0+1,i1+1,i2) = 1.0/8.0;			
				(*coarse)(i0+1,i1,i2+1) = 1.0/8.0;			
				(*coarse)(i0,i1+1,i2+1) = 1.0/8.0;			
				(*coarse)(i0+1,i1+1,i2+1) = 1.0/8.0;			
#endif
		}

#endif
		
#if defined(OFF_OCT_PEAK)
		coarse->zero();
#endif
		
		//... 1) grid self-contribution
		LOGUSER("Computing density self-contribution");
		PaddedDensitySubGrid<real_t> coarse_save( *coarse );
		
		//... create convolution kernel
		the_tf_kernel->fetch_kernel( levelmax );
		
		//... subtract oct mean on boundary but not in interior
		coarse->subtract_boundary_oct_mean();
		
		//... perform convolution
		convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()), shift );
		
		//... copy to grid hierarchy
		coarse->copy_add_unpad( *delta.get_grid(levelmax) );
		

		//... 2) boundary correction to top grid
		LOGUSER("Computing coarse grid correction");
		*coarse = coarse_save;
		
		//... subtract oct mean
		coarse->subtract_oct_mean();
		
		//... perform convolution
		convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()), shift );

		the_tf_kernel->deallocate();
		
		coarse->subtract_mean();
		
		//... upload data to coarser grid
		coarse->upload_bnd_add( *delta.get_grid(levelmax-1) );
			
		delete coarse;
	}
	
	delete the_tf_kernel;
	
#if 0
	// NEVER, NEVER ENABLE THE FOLLOWING
	
	//... enforce mean condition
	//for( int i=levelmin; i<(int)levelmax; ++i )
	//	enforce_mean( (*delta.get_grid(i+1)), (*delta.get_grid(i)) );
	
	/*for( unsigned i=levelmax; i>levelmin; --i )
		enforce_coarse_mean( (*delta.get_grid(i)), (*delta.get_grid(i-1)) );*/
#endif

#if 0
	//... subtract the box mean.... this will otherwise add
	//... a constant curvature term to the potential
	double sum = 0.0;
	{
		int nx,ny,nz;
		
		nx = delta.get_grid(levelmin)->size(0);
		ny = delta.get_grid(levelmin)->size(1);
		nz = delta.get_grid(levelmin)->size(2);
		
		for( int ix=0; ix<nx; ++ix )
			for( int iy=0; iy<ny; ++iy )
				for( int iz=0; iz<nz; ++iz )
					sum += (*delta.get_grid(levelmin))(ix,iy,iz);
		
		sum /= (nx*ny*nz);
	}
	
	
	std::cout << " - Top grid mean density is off by " << sum << ", correcting..." << std::endl;
	
	for( unsigned i=levelmin; i<=levelmax; ++i )
	{		
		int nx,ny,nz;
		nx = delta.get_grid(i)->size(0);
		ny = delta.get_grid(i)->size(1);
		nz = delta.get_grid(i)->size(2);
		
		for( int ix=0; ix<nx; ++ix )
			for( int iy=0; iy<ny; ++iy )
				for( int iz=0; iz<nz; ++iz )
					(*delta.get_grid(i))(ix,iy,iz) -= sum;
		
	}
#endif
		
	tend = omp_get_wtime();
	if( true )//verbosity > 1 )
		std::cout << " - Density calculation took " << tend-tstart << "s with " << omp_get_max_threads() << " threads." << std::endl;
	LOGUSER("Finished computing the density field in %fs",tend-tstart);
}


/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/

void normalize_density( grid_hierarchy& delta )
{	
	//return;
	
	long double sum = 0.0;
	unsigned levelmin = delta.levelmin(), levelmax = delta.levelmax();
	
	{
		size_t nx,ny,nz;
		
		nx = delta.get_grid(levelmin)->size(0);
		ny = delta.get_grid(levelmin)->size(1);
		nz = delta.get_grid(levelmin)->size(2);
		
		#pragma omp parallel for reduction(+:sum)
		for( int ix=0; ix<(int)nx; ++ix )
			for( size_t iy=0; iy<ny; ++iy )
				for( size_t iz=0; iz<nz; ++iz )
					sum += (*delta.get_grid(levelmin))(ix,iy,iz);
		
		sum /= (double)(nx*ny*nz);
	}
	
	std::cout << " - Top grid mean density is off by " << sum << ", correcting..." << std::endl;
	LOGUSER("Grid mean density is %g. Correcting...",sum);
	
	for( unsigned i=levelmin; i<=levelmax; ++i )
	{		
		size_t nx,ny,nz;
		nx = delta.get_grid(i)->size(0);
		ny = delta.get_grid(i)->size(1);
		nz = delta.get_grid(i)->size(2);
		
		#pragma omp parallel for
		for( int ix=0; ix<(int)nx; ++ix )
			for( size_t iy=0; iy<ny; ++iy )
				for( size_t iz=0; iz<nz; ++iz )
					(*delta.get_grid(i))(ix,iy,iz) -= sum;
	}
}


void coarsen_density( const refinement_hierarchy& rh, GridHierarchy<real_t>& u )
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

