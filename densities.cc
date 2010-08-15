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

// TODO: use optimized convolution routine when in unigrid mode 
void GenerateDensityUnigrid( config_file& cf, transfer_function *ptf, tf_type type, 
							refinement_hierarchy& refh, grid_hierarchy& delta, bool kspace, bool bdeconvolve, bool smooth )
{
	unsigned    levelmin,levelmax,levelminPoisson;
	real_t		boxlength;
	unsigned	ran_cube_size;
	
	std::vector<long> rngseeds;
	std::vector<std::string> rngfnames;
	
	
	levelminPoisson	= cf.getValue<unsigned>("setup","levelmin");
	levelmin	= cf.getValueSafe<unsigned>("setup","levelmin_TF",levelminPoisson);
	levelmax	= cf.getValue<unsigned>("setup","levelmax");
	boxlength   = cf.getValue<real_t>( "setup", "boxlength" );
	
	ran_cube_size = cf.getValueSafe<unsigned>("random","cubesize",DEF_RAN_CUBE_SIZE);
	
	std::cerr << " - Running unigrid version\n";
	
	//... parse random number options
	for( int i=0; i<=100; ++i )
	{
		char seedstr[128];
		std::string tempstr;
		sprintf(seedstr,"seed[%d]",i);
		if( cf.containsKey( "random", seedstr ) )
			tempstr = cf.getValue<std::string>( "random", seedstr );
		else
			tempstr = std::string("-2");
		
		if( is_number( tempstr ) )
		{	
			long ltemp;
			cf.convert( tempstr, ltemp );
			rngfnames.push_back( "" );
			rngseeds.push_back( ltemp );
		}else{
			rngfnames.push_back( tempstr );
			rngseeds.push_back(-1);
			std::cout << " - Random numbers for level " << std::setw(3) << i << " will be read from file.\n";
		}
		
	}
	
	//... parse grid setup parameters
	unsigned	nbase	= (unsigned)pow(2,levelmin);
	float lxref[3];
	std::string temp	= cf.getValue<std::string>( "setup", "ref_extent" );
	sscanf( temp.c_str(), "%g,%g,%g", &lxref[0],&lxref[1],&lxref[2] );
	
	int shift[3];
	shift[0] = cf.getValue<int>("setup","shift_x");
	shift[1] = cf.getValue<int>("setup","shift_y");
	shift[2] = cf.getValue<int>("setup","shift_z");
	
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
		
	convolution::parameters conv_param;
	conv_param.ptf = ptf;
	conv_param.pcf = &cf;
	
		
	//.. determine for which levels random seeds/random number files are given
	int lmaxread = -1, lmingiven = -1;
	for( unsigned ilevel = 0; ilevel < rngseeds.size(); ++ilevel )
	{	
		if( rngfnames[ilevel].size() > 0 )
			lmaxread = ilevel;
		if( rngseeds[ilevel] > 0 && lmingiven == -1 )
			lmingiven = ilevel;
	}
	
	if( (unsigned)lmingiven!=levelmin || levelmin!=levelminPoisson )
	{
		std::cerr <<" - Internal error: GenerateDensityUnigrid was called for a non-trivial\n"
				  <<"       problem set-up. This should not happen, GenerateDensityHierarchy\n"
				  <<"       should be called instead\n";
		throw std::runtime_error("Internal error");
	}
	
	std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmax << " ..." << std::endl;
	
	DensityGrid<real_t> *top = new DensityGrid<real_t>( nbase, nbase, nbase );
	
	double x0[3] = { refh.offset(levelmin,0), refh.offset(levelmin,1), refh.offset(levelmin,2) };
	double lx[3] = { refh.size(levelmin,0), refh.size(levelmin,1), refh.size(levelmin,2) };
	x0[0] /= pow(2,levelmin); x0[1] /= pow(2,levelmin); x0[2] /= pow(2,levelmin);
	lx[0] /= pow(2,levelmin); lx[1] /= pow(2,levelmin); lx[2] /= pow(2,levelmin);
	
	random_numbers<real_t> *rc = new random_numbers<real_t>( nbase, ran_cube_size, rngseeds[levelmin], true );//, x0, lx );


	if( shift[0]!=0||shift[1]!=0||shift[2]!=0 )
		std::cout << " - WARNING: will ignore non-zero shift in unigrid mode!\n";
	
	top->fill_rand( rc, 1.0, 0, 0, 0, true );
	//top->fill_rand( rc, x0[0], x0[1], x0[2], true );
	
	//rc->fill_all(*top);
	delete rc;
	
#if defined(SINGLE_PEAK)
	top->zero();
	(*top)(top->size(0)/2, top->size(1)/2, top->size(2)/2) = 1.0;
#elif defined(SINGLE_OCT_PEAK)
	{
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
	

	
	
	conv_param.lx = boxlength;
	conv_param.ly = boxlength;
	conv_param.lz = boxlength;
	conv_param.nx = top->nx_;
	conv_param.ny = top->ny_;
	conv_param.nz = top->nz_;
	conv_param.coarse_fact = 0;
	conv_param.deconvolve = bdeconvolve;
	conv_param.is_finest = true;
	conv_param.smooth = smooth;
	
	convolution::kernel *the_tf_kernel = the_kernel_creator->create( conv_param );
	convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ) );
	delete the_tf_kernel;
	
	delta.create_base_hierarchy(levelmin);
	top->copy( *delta.get_grid(levelmin) );
	delete top;

	for( int i=levelmax; i>0; --i )
		mg_straight().restrict( (*delta.get_grid(i)), (*delta.get_grid(i-1)) );
}

void GenerateDensityHierarchy(	config_file& cf, transfer_function *ptf, tf_type type, 
							  refinement_hierarchy& refh, grid_hierarchy& delta, bool bdeconvolve, bool smooth )
{
	unsigned					levelmin,levelmax,levelminPoisson;
	real_t						boxlength;
	std::vector<long>			rngseeds;
	std::vector<std::string>	rngfnames;
	bool force_shift(false),	kspaceTF;
	unsigned					ran_cube_size;
	
	constraint_set				contraints(cf);
	
	double tstart, tend;
	tstart = omp_get_wtime();
	
	
	levelminPoisson	= cf.getValue<unsigned>("setup","levelmin");
	levelmin		= cf.getValueSafe<unsigned>("setup","levelmin_TF",levelminPoisson);
	levelmax		= cf.getValue<unsigned>("setup","levelmax");
	boxlength		= cf.getValue<real_t>( "setup", "boxlength" );
	force_shift		= cf.getValueSafe<bool>("setup", "force_shift", force_shift );
	kspaceTF		= cf.getValueSafe<bool>("setup", "kspace_TF", false);
	ran_cube_size	= cf.getValueSafe<unsigned>("random","cubesize",DEF_RAN_CUBE_SIZE);
	
	
	//... parse random number options
	for( int i=0; i<=100; ++i )
	{
		char seedstr[128];
		std::string tempstr;
		sprintf(seedstr,"seed[%d]",i);
		if( cf.containsKey( "random", seedstr ) )
			tempstr = cf.getValue<std::string>( "random", seedstr );
		else
			tempstr = std::string("-2");

		if( is_number( tempstr ) )
		{	
			long ltemp;
			cf.convert( tempstr, ltemp );
			rngfnames.push_back( "" );
			rngseeds.push_back( ltemp );
		}else{
			rngfnames.push_back( tempstr );
			rngseeds.push_back(-1);
			std::cout << " - Random numbers for level " << std::setw(3) << i << " will be read from file.\n";
		}
			
	}
	
	//... parse grid setup parameters
	unsigned	nbase	= (unsigned)pow(2,levelmin);
	float lxref[3];
	std::string temp	= cf.getValue<std::string>( "setup", "ref_extent" );
	sscanf( temp.c_str(), "%g,%g,%g", &lxref[0],&lxref[1],&lxref[2] );
	
	
	int shift[3];
	shift[0] = cf.getValue<int>("setup","shift_x");
	shift[1] = cf.getValue<int>("setup","shift_y");
	shift[2] = cf.getValue<int>("setup","shift_z");
	
#ifdef SINGLE_PRECISION	
	convolution::kernel_creator *the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_real_float" ];
#else
	convolution::kernel_creator *the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_real_double" ];
#endif
	
	convolution::parameters conv_param;
	conv_param.ptf = ptf;
	conv_param.pcf = &cf;
	
	
	//... compute absolute grid offsets
	std::vector<int> offtotx(levelmax+1,0),offtoty(levelmax+1,0),offtotz(levelmax+1,0);
	for( unsigned ilevel = levelmin+1; ilevel<=levelmax; ++ilevel )
	{
		//... build a partial sum to get absolute offsets
		offtotx[ilevel] = 2*(offtotx[ilevel-1]+refh.offset(ilevel,0));
		offtoty[ilevel] = 2*(offtoty[ilevel-1]+refh.offset(ilevel,1));
		offtotz[ilevel] = 2*(offtotz[ilevel-1]+refh.offset(ilevel,2));
	}
	for( unsigned ilevel = levelmin+1; ilevel<=levelmax; ++ilevel )
	{
		//... the arrays are doubled in size for the isolated BCs
		offtotx[ilevel] -= refh.size(ilevel,0)/2;
		offtoty[ilevel] -= refh.size(ilevel,1)/2;
		offtotz[ilevel] -= refh.size(ilevel,2)/2;
	}
	
	//.. determine for which levels random seeds/random number files are given
	int lmaxread = -1, lmingiven = -1;
	for( unsigned ilevel = 0; ilevel < rngseeds.size(); ++ilevel )
	{	
		if( rngfnames[ilevel].size() > 0 )
			lmaxread = ilevel;
		if( rngseeds[ilevel] > 0 && lmingiven == -1 )
			lmingiven = ilevel;
	}
	
	/***** RESORT TO UNIGRID VERSION FOR SIMPLE GRID SETUPS *****/
	
	// TODO: need to make sure unigrid gets called whenever possible
	// FIXME: temporarily disabled
	if( lmingiven == (int)levelmin && levelmin == levelmax && levelmin==levelminPoisson )
	{	
		GenerateDensityUnigrid(cf,ptf,type,refh,delta,kspaceTF,bdeconvolve,smooth);
		return;
	}
	
	/***** GENERATE WHITE NOISE FIELDS FOR MULTI-LEVEL GRID *****/
	
	//... if random numbers are to be read from file, do this now
	std::vector< random_numbers<real_t>* > randc;
	randc.assign(std::max(lmaxread,std::max(lmingiven,(int)levelmax))+1,(random_numbers<real_t>*)NULL);
	
	if( lmaxread >= (int)levelmin )
	{
		randc[lmaxread] = new random_numbers<real_t>( (unsigned)pow(2,lmaxread), rngfnames[lmaxread] );
		for( int ilevel = lmaxread-1; ilevel >= (int)levelmin; --ilevel )
			randc[ilevel] = new random_numbers<real_t>( *randc[ilevel+1] );
	}
	
	
	//... if random numbers are not given for lower levels, obtain them by averaging
	//if( lmingiven >= (int)levelmin )
	if( lmingiven > (int)levelmin )
	{
		randc[lmingiven] = new random_numbers<real_t>( (unsigned)pow(2,lmingiven), ran_cube_size, rngseeds[lmingiven], true );//, x0, lx );
		
		for( int ilevel = lmingiven-1; ilevel >= (int)levelmin; --ilevel ){
			if( rngseeds[ilevel-levelmin] > 0 )
				std::cerr << " - Warning: random seed for level " << ilevel << " will be ignored.\n"
				<< "            consistency requires that it is obtained by restriction from level " << lmingiven << std::endl;
			
			randc[ilevel] = new random_numbers<real_t>( *randc[ilevel+1] );
		}
	}
	
	//... if random seeds are given for levels coarser than levelmin, use them as constraints
	if( lmingiven < (int)levelmin && !(lmaxread>=(int)levelmin) )
	{
		randc[lmingiven] = new random_numbers<real_t>( (unsigned)pow(2,lmingiven), ran_cube_size, rngseeds[lmingiven], true );//, x0, lx );
		
		for( int ilevel = lmingiven+1; ilevel <= (int)levelmin; ++ilevel )
		{
			long seed = rngseeds[ilevel];
			if( seed <= 0 )
				seed = rngseeds[lmingiven]+ilevel;
			
			randc[ilevel] = new random_numbers<real_t>( *randc[ilevel-1], ran_cube_size, seed, true );
			
			delete randc[ilevel-1];
			randc[ilevel-1]=NULL;
		}
	}	
	
	//... create and initialize density grids with white noise	
	PaddedDensitySubGrid<real_t>* coarse(NULL), *fine(NULL);
	DensityGrid<real_t>* top(NULL);
		
	//... clean up ...//
	for( unsigned i=0; i<randc.size(); ++i )
	{
		if( i<levelmin || i>levelmax )
			if( randc[i] != NULL )
			{	
				delete randc[i];
				randc[i]=NULL;
			}
	}
		
		
	/***** PERFORM CONVOLUTIONS *****/
	
	if( levelmax == levelmin )
	{
		std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmax << " ..." << std::endl;
		
		top = new DensityGrid<real_t>( nbase, nbase, nbase );
	
		random_numbers<real_t> *rc;
		
		if( levelminPoisson == levelmin && !force_shift)
		{
			if( randc[levelmin] == NULL )
				rc = new random_numbers<real_t>( nbase, ran_cube_size, rngseeds[levelmin], true );
			else
				rc = randc[levelmin];
			
			if( shift[0]!=0||shift[1]!=0||shift[2]!=0 )
				std::cout << " - WARNING: will ignore non-zero shift in unigrid mode!\n";

			top->fill_rand( rc, 1.0, 0, 0, 0, true );
		}
		else
		{	
			int lfac = (int)pow(2,levelmin-levelminPoisson);
			int x0[3] = { -shift[0]*lfac, -shift[1]*lfac, -shift[2]*lfac };
			int lx[3] = { refh.size(levelmin,0), refh.size(levelmin,1), refh.size(levelmin,2) };
			
			if( randc[levelmin] == NULL )
				rc = randc[levelmin] = new random_numbers<real_t>( nbase, ran_cube_size, rngseeds[levelmin], x0, lx );
			//
			//if( randc[levelmin] == NULL )
			//	rc = new random_numbers<real_t>( nbase, ran_cube_size, rngseeds[levelmin], true );
			else
				rc = randc[levelmin];
			
			top->fill_rand( rc, 1.0, x0[0], x0[1], x0[2], true );
		}

		delete rc;
		
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
		
		conv_param.lx = boxlength;
		conv_param.ly = boxlength;
		conv_param.lz = boxlength;
		conv_param.nx = top->nx_;
		conv_param.ny = top->ny_;
		conv_param.nz = top->nz_;
		conv_param.coarse_fact = 0;
		conv_param.deconvolve = bdeconvolve;
		conv_param.is_finest = true;
		conv_param.smooth = smooth;
		
		convolution::kernel *the_tf_kernel = the_kernel_creator->create( conv_param );
		convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ) );
		delete the_tf_kernel;
		
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
			int lfac = (int)pow(2,levelmin-levelminPoisson);
			
			top = new DensityGrid<real_t>( nbase, nbase, nbase );
			
			random_numbers<real_t> *rc;
			/*int x0[3] = {	refh.offset_abs(levelmin,0)-lfac*shift[0], 
				refh.offset_abs(levelmin,1)-lfac*shift[1], 
				refh.offset_abs(levelmin,2)-lfac*shift[2] };
			
			int lx[3] = {	refh.size(levelmin,0), 
				refh.size(levelmin,1), 
				refh.size(levelmin,2) };*/
			
			if( randc[levelmin] == NULL )
			{	
				std::cout << " - Creating new random numbers for level " << levelmin << std::endl;
				rc = randc[levelmin] = new random_numbers<real_t>( nbase, ran_cube_size, rngseeds[levelmin], true );//x0, lx );
			}
			else
				rc = randc[levelmin];
			
			top->fill_rand( rc, 1.0, -lfac*shift[0], -lfac*shift[1], -lfac*shift[2], true );
			delete rc;
			randc[levelmin] = NULL;
			
		}
		
		fine = new PaddedDensitySubGrid<real_t>(	refh.offset(levelmin+i+1,0), 
												refh.offset(levelmin+i+1,1), 
												refh.offset(levelmin+i+1,2), 
												refh.size(levelmin+i+1,0), 
												refh.size(levelmin+i+1,1), 
												refh.size(levelmin+i+1,2) );
		
		{
			random_numbers<real_t> *rc;
			int x0[3],lx[3];
			int lfac = (int)pow(2,levelmin+i+1-levelminPoisson);
			
			lx[0] = 2*refh.size(levelmin+i+1,0); 
			lx[1] = 2*refh.size(levelmin+i+1,1); 
			lx[2] = 2*refh.size(levelmin+i+1,2); 
			x0[0] = refh.offset_abs(levelmin+i+1,0)-lfac*shift[0]-lx[0]/4; 
			x0[1] = refh.offset_abs(levelmin+i+1,1)-lfac*shift[1]-lx[1]/4; 
			x0[2] = refh.offset_abs(levelmin+i+1,2)-lfac*shift[2]-lx[2]/4; 
			
			if( randc[levelmin+i+1] == NULL )
			{	
				std::cout << " - Creating new random numbers for level " << levelmin+1+i << std::endl;
				rc = randc[levelmin+i+1] = new random_numbers<real_t>((unsigned)pow(2,levelmin+i+1), ran_cube_size, rngseeds[levelmin+i+1], x0, lx);
			}
			else
				rc = randc[levelmin+i+1];
			
			fine->fill_rand( rc, 1.0, x0[0], x0[1], x0[2], false );
			
				
			if( true )//i> 0 )//i+levelmin+1 > (unsigned)lmingiven )
			{	
				if(i==0)
					fine->constrain( *top );
				else
					fine->constrain( *coarse );
			}				
			
			delete rc;
			rc = NULL;
		}
		
		
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

			conv_param.lx = boxlength;
			conv_param.ly = boxlength;
			conv_param.lz = boxlength;
			conv_param.nx = top->nx_;
			conv_param.ny = top->ny_;
			conv_param.nz = top->nz_;
			conv_param.coarse_fact = levelmax-levelmin;
			conv_param.deconvolve = bdeconvolve;
			conv_param.is_finest = false;
			conv_param.smooth = smooth;
			convolution::kernel *the_tf_kernel = the_kernel_creator->create( conv_param );
			
			
			//... 1) compute standard convolution for levelmin
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ) );
			top->copy( *delta.get_grid(levelmin) );
			
			
			//... 2) compute contribution to finer grids from non-refined region
			*top = top_save;
			top_save.clear();
			top->zero_subgrid(refh.offset(levelmin+i+1,0), refh.offset(levelmin+i+1,1), refh.offset(levelmin+i+1,2), 
							  refh.size(levelmin+i+1,0)/2, refh.size(levelmin+i+1,1)/2, refh.size(levelmin+i+1,2)/2 );

			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ) );
			
			delete the_tf_kernel;

			
			meshvar_bnd delta_longrange( *delta.get_grid(levelmin) );
			top->copy( delta_longrange );
			delete top;			
			
			//... restrict these contributions to the next level
			
			for( int j=1; j<=(int)levelmax-(int)levelmin; ++j )
			{
				int R = j;
				delta.add_patch( refh.offset(levelmin+j,0), refh.offset(levelmin+j,1), refh.offset(levelmin+j,2), 
									refh.size(levelmin+j,0), refh.size(levelmin+j,1), refh.size(levelmin+j,2) );
			
				mg_cubic_mult().prolong( delta_longrange, *delta.get_grid(levelmin+j),
										refh.offset_abs(levelmin,0), refh.offset_abs(levelmin,1), refh.offset_abs(levelmin,2),
										refh.offset_abs(levelmin+j,0), refh.offset_abs(levelmin+j,1), refh.offset_abs(levelmin+j,2), R);
										
																		  
			}
			
			
			
			//delta.add_patch( refh.offset(levelmin+1,0), refh.offset(levelmin+1,1), refh.offset(levelmin+1,2), 
			//				refh.size(levelmin+1,0), refh.size(levelmin+1,1), refh.size(levelmin+1,2) );
			
			//mg_cubic().prolong( delta_longrange, *delta.get_grid(levelmin+1) );
			
					
		}
		else
		{
			/**********************************************************************************************************\
			 *	multi-grid: intermediate sub-grids .....
			 \**********************************************************************************************************/ 
			
			std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmin+i << " ..." << std::endl;
			
			//delta.add_patch( refh.offset(levelmin+i+1,0), refh.offset(levelmin+i+1,1), refh.offset(levelmin+i+1,2), 
			//				refh.size(levelmin+i+1,0), refh.size(levelmin+i+1,1), refh.size(levelmin+i+1,2) );
			
			//mg_cubic().prolong( *delta.get_grid(levelmin+i), *delta.get_grid(levelmin+i+1) );
			
			real_t dx,lx,ly,lz;
			dx = boxlength/pow(2.0,levelmin+i);
			
			lx = dx * coarse->nx_; 
			ly = dx * coarse->ny_; 
			lz = dx * coarse->nz_; 
			
			//.. set convolution parameters
			//TODO: this needs to be changed, forgetting to set a parameter will not be warned!
			conv_param.lx = lx;
			conv_param.ly = ly;
			conv_param.lz = lz;
			conv_param.nx = coarse->nx_;
			conv_param.ny = coarse->ny_;
			conv_param.nz = coarse->nz_;
			conv_param.coarse_fact = levelmax-levelmin-i;
			conv_param.deconvolve = bdeconvolve;
			conv_param.is_finest = false;
			conv_param.smooth = smooth;
			
			convolution::kernel *the_tf_kernel = the_kernel_creator->create( conv_param );
			
			PaddedDensitySubGrid<real_t> coarse_save( *coarse );
					
			//... 1) the inner region
			coarse->subtract_boundary_oct_mean();
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
			coarse->copy_add_unpad( *delta.get_grid(levelmin+i) );
			
			
			//... 2) the 'BC' for the next finer grid
			*coarse = coarse_save;
			coarse->subtract_boundary_oct_mean();
			coarse->zero_subgrid(refh.offset(levelmin+i+1,0), refh.offset(levelmin+i+1,1), refh.offset(levelmin+i+1,2), 
								 refh.size(levelmin+i+1,0)/2, refh.size(levelmin+i+1,1)/2, refh.size(levelmin+i+1,2)/2 );
			
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
			
			meshvar_bnd delta_longrange( *delta.get_grid(levelmin+i) );
			coarse->copy_unpad( delta_longrange );
			
			for( int j=1; j<=(int)levelmax-(int)levelmin-i; ++j )
			{
				int R = j;
				//delta.add_patch( refh.offset(levelmin+i+j,0), refh.offset(levelmin+i+j,1), refh.offset(levelmin+i+j,2), 
				//				refh.size(levelmin+i+j,0), refh.size(levelmin+i+j,1), refh.size(levelmin+i+j,2) );
				
				mg_cubic_mult().prolong_add( delta_longrange, *delta.get_grid(levelmin+i+j),
										refh.offset_abs(levelmin+i,0), refh.offset_abs(levelmin+i,1), refh.offset_abs(levelmin+i,2),
										refh.offset_abs(levelmin+i+j,0), refh.offset_abs(levelmin+i+j,1), refh.offset_abs(levelmin+i+j,2), R);
				
				
			}
			
			//mg_cubic().prolong_add( delta_longrange, *delta.get_grid(levelmin+i+1) );
			

			//... 3) the coarse-grid correction
			*coarse = coarse_save;
			coarse->subtract_oct_mean();
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
			coarse->upload_bnd_add( *delta.get_grid(levelmin+i-1) );
			
			delete the_tf_kernel;
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

		real_t dx,lx,ly,lz;
		dx = boxlength/pow(2.0,levelmax);
		
		lx = dx * coarse->nx_; 
		ly = dx * coarse->ny_; 
		lz = dx * coarse->nz_; 
		
		//... set convolution parameters
		conv_param.lx = lx;
		conv_param.ly = ly;
		conv_param.lz = lz;
		conv_param.nx = coarse->nx_;
		conv_param.ny = coarse->ny_;
		conv_param.nz = coarse->nz_;
		conv_param.coarse_fact = 0;
		conv_param.deconvolve = bdeconvolve;
		conv_param.is_finest = true;	
		conv_param.smooth = smooth;
		
#if defined(SINGLE_PEAK) || defined(SINGLE_OCT_PEAK)
		{
			coarse->zero();
			
			int 
				i0 = pow(2,levelmax)/2 - refh.offset_abs(levelmax,0) + coarse->nx_/4,
				i1 = pow(2,levelmax)/2 - refh.offset_abs(levelmax,1) + coarse->nx_/4,
				i2 = pow(2,levelmax)/2 - refh.offset_abs(levelmax,2) + coarse->nx_/4;
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
		
		//... 1) grid self-contributio
		
		PaddedDensitySubGrid<real_t> coarse_save( *coarse );
		
		//... create convolution kernel
		convolution::kernel *the_tf_kernel = the_kernel_creator->create( conv_param );
		
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
		
		//coarse->subtract_mean();
		
		//... upload data to coarser grid
		coarse->upload_bnd_add( *delta.get_grid(levelmax-1) );
		
		delete the_tf_kernel;		
		delete coarse;
	}
	
	//... clean up ...//
	for( unsigned i=0; i<randc.size(); ++i )
	{
		if( i<levelmin || i>levelmax )
			if( randc[i] != NULL )
			{	
				delete randc[i];
				randc[i]=NULL;
			}
	}
	
#if 0
	// NEVER, NEVER ENABLE THE FOLLOWING
	
	//... enforce mean condition
	//for( int i=levelmin; i<(int)levelmax; ++i )
	//	enforce_mean( (*delta.get_grid(i+1)), (*delta.get_grid(i)) );
	
	//	for( unsigned i=levelmax; i>levelmin; --i )
		//	enforce_coarse_mean( (*delta.get_grid(i)), (*delta.get_grid(i-1)) );
#endif

#if 1
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
	
	
	//std::cout << " - Top grid mean density is off by " << sum << ", correcting..." << std::endl;
	
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
	
	//... fill coarser levels with data from finer ones...
	for( int i=levelmax; i>0; --i )
		mg_straight().restrict( (*delta.get_grid(i)), (*delta.get_grid(i-1)) );
	
	
	tend = omp_get_wtime();
	if( true )//verbosity > 1 )
		std::cout << " - Density calculation took " << tend-tstart << "s with " << omp_get_max_threads() << " threads." << std::endl;
	
}


void normalize_density( grid_hierarchy& delta )
{
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

