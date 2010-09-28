/*
 
 densities.cc - This file is part of MUSIC -
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

#include "constraints.hh"
#include "densities.hh"
#include "convolution_kernel.hh"

//... uncomment this to have a single peak in the centre and otherwise zeros
//#define SINGLE_PEAK
//#define SINGLE_OCT_PEAK
//#define OFF_OCT_PEAK
//#define DEGRADE_RAND



//TODO: this should be a larger number by default, just to maintain consistency with old default
#define DEF_RAN_CUBE_SIZE	32

bool is_number(const std::string& s)
{
	for (unsigned i = 0; i < s.length(); i++)
		if (!std::isdigit(s[i])&&s[i]!='-' )
			return false;
	
	return true;
}

/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/

void GenerateDensityUnigrid( config_file& cf, transfer_function *ptf, tf_type type, 
							refinement_hierarchy& refh, grid_hierarchy& delta, bool kspace, bool bdeconvolve, bool smooth )
{
	unsigned    levelmin,levelmax,levelminPoisson;
	real_t		boxlength;
	
	
	levelminPoisson	= cf.getValue<unsigned>("setup","levelmin");
	levelmin	= cf.getValueSafe<unsigned>("setup","levelmin_TF",levelminPoisson);
	levelmax	= cf.getValue<unsigned>("setup","levelmax");
	boxlength   = cf.getValue<real_t>( "setup", "boxlength" );
	
	unsigned	nbase	= 1<<levelmin;
	
	std::cerr << " - Running unigrid version\n";
	
	//... select the transfer function to be used

	convolution::kernel_creator *the_kernel_creator;

	if( kspace )
	{
		std::cout << " - Using k-space transfer function kernel.\n";
		
		#ifdef SINGLE_PRECISION	
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_k_float" ];
		#else
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_k_double" ];
		#endif
	}
	else
	{
		std::cout << " - Using real-space transfer function kernel.\n";
		
		#ifdef SINGLE_PRECISION	
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_real_float" ];
		#else
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_real_double" ];
		#endif
	}	
		
	
	//... initialize random number generator
	random_number_generator<random_numbers<real_t>,real_t> rand_gen( cf, refh );
	
	//... initialize convolution kernel
	convolution::kernel *the_tf_kernel = the_kernel_creator->create( cf, ptf, refh, type );

	//...
	std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmax << " ..." << std::endl;
	
	//... create convolution mesh
	DensityGrid<real_t> *top = new DensityGrid<real_t>( nbase, nbase, nbase );
	
	//... fill with random numbers
	rand_gen.load( *top, levelmin );
	
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
	convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ) );
	
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
							  refinement_hierarchy& refh, grid_hierarchy& delta, bool bdeconvolve, bool smooth )
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
	
	
	random_number_generator<random_numbers<real_t>,real_t> rand_gen( cf, refh );
	
	unsigned	nbase	= (unsigned)pow(2,levelmin);
	
	convolution::kernel_creator *the_kernel_creator;
	
	if( kspaceTF )
	{
		if( levelmin!=levelmax )
			throw std::runtime_error("k-space transfer function can only be used in unigrid density mode");
		
		std::cout << " - Using k-space transfer function kernel.\n";
		
#ifdef SINGLE_PRECISION	
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_k_float" ];
#else
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_k_double" ];
#endif
	}
	else
	{
		std::cout << " - Using real-space transfer function kernel.\n";
		
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
		
		top = new DensityGrid<real_t>( nbase, nbase, nbase );
		rand_gen.load( *top, levelmin );

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
				
		convolution::perform<real_t>( the_tf_kernel->fetch_kernel( levelmax ), reinterpret_cast<void*>( top->get_data_ptr() ) );
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
			rand_gen.load(*top,levelmin);
		}
		
		fine = new PaddedDensitySubGrid<real_t>( refh.offset(levelmin+i+1,0), refh.offset(levelmin+i+1,1), refh.offset(levelmin+i+1,2), 
												refh.size(levelmin+i+1,0), 	refh.size(levelmin+i+1,1), 	refh.size(levelmin+i+1,2) );
		rand_gen.load(*fine,levelmin+i+1);
		
		//.......................................................................................................//
		//... PERFORM CONVOLUTIONS ..............................................................................//
		//.......................................................................................................//		
		if( i==0 )
		{
			/**********************************************************************************************************\
			 *	multi-grid: top-level grid grids .....
			 \**********************************************************************************************************/ 
			std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmin+i << " ..." << std::endl;
			
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
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ) );
			top->copy( *delta.get_grid(levelmin) );
			
			
			//... 2) compute contribution to finer grids from non-refined region
			*top = top_save;
			top_save.clear();
			top->zero_subgrid(refh.offset(levelmin+i+1,0), refh.offset(levelmin+i+1,1), refh.offset(levelmin+i+1,2), 
							  refh.size(levelmin+i+1,0)/2, refh.size(levelmin+i+1,1)/2, refh.size(levelmin+i+1,2)/2 );

			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ) );
			the_tf_kernel->deallocate();
			
			meshvar_bnd delta_longrange( *delta.get_grid(levelmin) );
			top->copy( delta_longrange );
			delete top;			
			
			//... restrict these contributions to the next level
			
			for( int j=1; j<=(int)levelmax-(int)levelmin; ++j )
			{
				delta.add_patch( refh.offset(levelmin+j,0), refh.offset(levelmin+j,1), refh.offset(levelmin+j,2), 
									refh.size(levelmin+j,0), refh.size(levelmin+j,1), refh.size(levelmin+j,2) );
			
				mg_cubic_mult().prolong( delta_longrange, *delta.get_grid(levelmin+j),
										refh.offset_abs(levelmin,0), refh.offset_abs(levelmin,1), refh.offset_abs(levelmin,2),
										refh.offset_abs(levelmin+j,0), refh.offset_abs(levelmin+j,1), refh.offset_abs(levelmin+j,2), j);
			}		
		}
		else
		{
			/**********************************************************************************************************\
			 *	multi-grid: intermediate sub-grids .....
			 \**********************************************************************************************************/ 
			
			std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmin+i << " ..." << std::endl;
			
			PaddedDensitySubGrid<real_t> coarse_save( *coarse );
			the_tf_kernel->fetch_kernel( levelmin+i );
					
			//... 1) the inner region
			coarse->subtract_boundary_oct_mean();
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
			coarse->copy_add_unpad( *delta.get_grid(levelmin+i) );
			
			
			//... 2) the 'BC' for the next finer grid
			*coarse = coarse_save;
			coarse->subtract_boundary_oct_mean();
			coarse->zero_subgrid(refh.offset(levelmin+i+1,0), refh.offset(levelmin+i+1,1), refh.offset(levelmin+i+1,2), 
								 refh.size(levelmin+i+1,0)/2, refh.size(levelmin+i+1,1)/2, refh.size(levelmin+i+1,2)/2 );
			
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( coarse->get_data_ptr() ) );
			
			//... interpolate to finer grid(s)
			meshvar_bnd delta_longrange( *delta.get_grid(levelmin+i) );
			coarse->copy_unpad( delta_longrange );
			
			for( int j=1; j<=(int)levelmax-(int)levelmin-i; ++j )
					mg_cubic_mult().prolong_add( delta_longrange, *delta.get_grid(levelmin+i+j),
										refh.offset_abs(levelmin+i,0), refh.offset_abs(levelmin+i,1), refh.offset_abs(levelmin+i,2),
										refh.offset_abs(levelmin+i+j,0), refh.offset_abs(levelmin+i+j,1), refh.offset_abs(levelmin+i+j,2), j);
				

			//... 3) the coarse-grid correction
			*coarse = coarse_save;
			coarse->subtract_oct_mean();
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
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
		PaddedDensitySubGrid<real_t> coarse_save( *coarse );
		
		//... create convolution kernel
		the_tf_kernel->fetch_kernel( levelmax );
		
		//... subtract oct mean on boundary but not in interior
		coarse->subtract_boundary_oct_mean();
		
		//... perform convolution
		convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
		
		//... copy to grid hierarchy
		coarse->copy_add_unpad( *delta.get_grid(levelmax) );
		

		//... 2) boundary correction to top grid
		*coarse = coarse_save;
		
		//... subtract oct mean
		coarse->subtract_oct_mean();
		
		
		//... perform convolution
		convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );

		the_tf_kernel->deallocate();
		
		//coarse->subtract_mean();
		
		//... upload data to coarser grid
		coarse->upload_bnd_add( *delta.get_grid(levelmax-1) );
		
		//delete the_tf_kernel;		
		delete coarse;
	}
	
	delete the_tf_kernel;
	
#if 0
	// NEVER, NEVER ENABLE THE FOLLOWING
	
	//... enforce mean condition
	//for( int i=levelmin; i<(int)levelmax; ++i )
	//	enforce_mean( (*delta.get_grid(i+1)), (*delta.get_grid(i)) );
	
	for( unsigned i=levelmax; i>levelmin; --i )
		enforce_coarse_mean( (*delta.get_grid(i)), (*delta.get_grid(i-1)) );
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
	
}


/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/

void normalize_density( grid_hierarchy& delta )
{	
	//return;
	
	double sum = 0.0;
	unsigned levelmin = delta.levelmin(), levelmax = delta.levelmax();
	
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
}

