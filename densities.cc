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

#include "densities.hh"
#include "convolution_kernel.hh"

bool is_number(const std::string& s)
{
	for (unsigned i = 0; i < s.length(); i++)
		if (!std::isdigit(s[i])&&s[i]!='-' )
			return false;
	
	return true;
}

// TODO: use optimized convolution routine when in unigrid mode 
void GenerateDensityUnigrid( config_file& cf, transfer_function *ptf, tf_type type, refinement_hierarchy& refh, grid_hierarchy& delta, bool kspace )
{
	unsigned    levelmin,levelmax,levelminPoisson;
	real_t		boxlength;
	std::vector<long> rngseeds;
	std::vector<std::string> rngfnames;
	
	
	levelminPoisson	= cf.getValue<unsigned>("setup","levelmin");
	levelmin	= cf.getValueSafe<unsigned>("setup","levelminTF",levelminPoisson);
	levelmax	= cf.getValue<unsigned>("setup","levelmax");
	boxlength   = cf.getValue<real_t>( "setup", "boxlength" );
	
	std::cerr << " RUNNING UNIGRID VERSION\n";
	
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
		#ifdef SINGLE_PRECISION	
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_k_float" ];
		#else
		the_kernel_creator = convolution::get_kernel_map()[ "tf_kernel_k_double" ];
		#endif
	}
	else
	{
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
	
	random_numbers<real_t> *rc = new random_numbers<real_t>( nbase, 32, rngseeds[levelmin], true );//, x0, lx );


	if( shift[0]!=0||shift[1]!=0||shift[2]!=0 )
		std::cout << " - WARNING: will ignore non-zero shift in unigrid mode!\n";
	
	//top->fill_rand( rc, 1.0, 0, 0, 0, true );
	rc->fill_all(*top);
	delete rc;
	
	conv_param.lx = boxlength;
	conv_param.ly = boxlength;
	conv_param.lz = boxlength;
	conv_param.nx = top->nx_;
	conv_param.ny = top->ny_;
	conv_param.nz = top->nz_;
	conv_param.normalize = false;
	
	convolution::kernel *the_tf_kernel = the_kernel_creator->create( conv_param );
	convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ) );
	delete the_tf_kernel;
	
	delta.create_base_hierarchy(levelmin);
	top->copy( *delta.get_grid(levelmin) );
	delete top;

	
	for( int i=levelmax; i>0; --i )
		mg_straight().restrict( (*delta.get_grid(i)), (*delta.get_grid(i-1)) );
	
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
}

void GenerateDensityHierarchy(	config_file& cf, transfer_function *ptf, tf_type type, refinement_hierarchy& refh, grid_hierarchy& delta )
{
	unsigned    levelmin,levelmax,levelminPoisson;
	real_t		boxlength;
	std::vector<long> rngseeds;
	std::vector<std::string> rngfnames;
	bool force_shift(false);
	
	
	levelminPoisson	= cf.getValue<unsigned>("setup","levelmin");
	levelmin	= cf.getValueSafe<unsigned>("setup","levelmin_TF",levelminPoisson);
	levelmax	= cf.getValue<unsigned>("setup","levelmax");
	boxlength   = cf.getValue<real_t>( "setup", "boxlength" );
	force_shift	= cf.getValueSafe<bool>("setup", "force_shift", force_shift );
	
	// TODO: need to make sure unigrid gets called whenever possible
	/*if( levelmin == levelmax && levelmin==levelminPoisson )
	{	
		GenerateDensityUnigrid(cf,ptf,type,refh,delta);
		return;
	}*/
	
	
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
	//double      lextmin = std::min(lxref[0],std::min(lxref[1],lxref[2]));
	// alpha = 0.1*..., cutoff = 0.25 /// maybe 0.08, 0.3
	real_t		alpha	= 0.05*(real_t)nbase/boxlength;//0.45*lextmin*(double)nbase;///boxlength;
	real_t		cutoff	= 0.25;//lextmin*0.9;
	
	
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
	if( lmingiven >= (int)levelmin )
	{
		randc[lmingiven] = new random_numbers<real_t>( (unsigned)pow(2,lmingiven), 32, rngseeds[lmingiven], true );//, x0, lx );
		
		for( int ilevel = lmingiven-1; ilevel >= (int)levelmin; --ilevel ){
			if( rngseeds[ilevel-levelmin] > 0 )
				std::cerr << " - Warning: random seed for level " << ilevel << " will be ignored.\n"
				<< "            consistency requires that it is obtained by restriction from level " << lmingiven << std::endl;
			
			randc[ilevel] = new random_numbers<real_t>( *randc[ilevel+1] );
		}
	}
	
	//... if random seeds are given for levels coarser than levelmin, use them as constraints
	if( lmingiven < (int)levelmin )
	{
		throw std::runtime_error("You provided a seed for a level below levelmin, this is not supported yet.");
		
		randc[lmingiven] = new random_numbers<real_t>( (unsigned)pow(2,lmingiven), 32, rngseeds[lmingiven], true );//, x0, lx );
		
		for( int ilevel = lmingiven+1; ilevel <= (int)levelmin; ++ilevel )
		{
			long seed = rngseeds[ilevel];
			if( seed <= 0 )
				seed = rngseeds[lmingiven+ilevel];
			
			randc[ilevel] = new random_numbers<real_t>( *randc[ilevel-1] );
		}
		
		delete randc[lmingiven];
		randc[lmingiven] = NULL;
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
		
		
	//... perform convolutions ...//
	if( levelmax == levelmin )
	{
		std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmax << " ..." << std::endl;
		
		top = new DensityGrid<real_t>( nbase, nbase, nbase );
	
		random_numbers<real_t> *rc;
		
		if( levelminPoisson == levelmin && !force_shift)
		{
			if( randc[levelmin] == NULL )
				rc = new random_numbers<real_t>( nbase, 32, rngseeds[levelmin], true );
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
				rc = randc[levelmin] = new random_numbers<real_t>( nbase, 32, rngseeds[levelmin], x0, lx );
			//
			//if( randc[levelmin] == NULL )
			//	rc = new random_numbers<real_t>( nbase, 32, rngseeds[levelmin], true );
			else
				rc = randc[levelmin];
			
			top->fill_rand( rc, 1.0, x0[0], x0[1], x0[2], true );
		}

		delete rc;
		
		conv_param.lx = boxlength;
		conv_param.ly = boxlength;
		conv_param.lz = boxlength;
		conv_param.nx = top->nx_;
		conv_param.ny = top->ny_;
		conv_param.nz = top->nz_;
		conv_param.normalize = false;
		
		convolution::kernel *the_tf_kernel = the_kernel_creator->create( conv_param );
		convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ) );
		delete the_tf_kernel;
		
		delta.create_base_hierarchy(levelmin);
		top->copy( *delta.get_grid(levelmin) );
		delete top;
	}
	//}		
		
	
	for( int i=0; i< (int)levelmax-(int)levelmin; ++i )
	{
		//.......................................................................................................//
		//... GENERATE/FILL WITH RANDOM NUMBERS .................................................................//
		//.......................................................................................................//
		
		
		if( i==0 )
		{
			top = new DensityGrid<real_t>( nbase, nbase, nbase );
			
			random_numbers<real_t> *rc;
			int x0[3] = {	refh.offset_abs(levelmin,0)-shift[0], 
				refh.offset_abs(levelmin,1)-shift[1], 
				refh.offset_abs(levelmin,2)-shift[2] };
			
			int lx[3] = {	refh.size(levelmin,0), 
				refh.size(levelmin,1), 
				refh.size(levelmin,2) };
			
			if( randc[levelmin] == NULL )
				rc = randc[levelmin] = new random_numbers<real_t>( nbase, 32, rngseeds[levelmin], x0, lx );
			else
				rc = randc[levelmin];
			
			top->fill_rand( rc, 1.0, -shift[0], -shift[1], -shift[2], true );
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
			int lfac = (int)pow(2.0,i+1);
			x0[0] = refh.offset_abs(levelmin+i+1,0)-lfac*shift[0]; //((real_t)(offtotx[levelmin+i+1]+lfac*shift[0]))/pow(2.0,levelmin+i+1);
			x0[1] = refh.offset_abs(levelmin+i+1,1)-lfac*shift[1]; //((real_t)(offtoty[levelmin+i+1]+lfac*shift[1]))/pow(2.0,levelmin+i+1);
			x0[2] = refh.offset_abs(levelmin+i+1,2)-lfac*shift[2]; //((real_t)(offtotz[levelmin+i+1]+lfac*shift[2]))/pow(2.0,levelmin+i+1);
			lx[0] = refh.size(levelmin+i+1,0); // /pow(2.0,levelmin+i+1);
			lx[1] = refh.size(levelmin+i+1,1); // /pow(2.0,levelmin+i+1);
			lx[2] = refh.size(levelmin+i+1,2); // /pow(2.0,levelmin+i+1);
			
			//x0[0] -= 0.5*lx[0];  lx[0] *= 2.0;
			//x0[1] -= 0.5*lx[1];  lx[1] *= 2.0;
			//x0[2] -= 0.5*lx[2];  lx[2] *= 2.0;
			
			if( randc[levelmin+i+1] == NULL )
				rc = randc[levelmin+i+1] = new random_numbers<real_t>((unsigned)pow(2,levelmin+i+1), 32, rngseeds[levelmin+i+1], x0, lx);
			else
				rc = randc[levelmin+i+1];
			
			
			{
				/*int llfac = (int)pow(2.0,i+1);
				fine->fill_rand( rc, 1.0,	offtotx[levelmin+i+1]-llfac*shift[0], 
								offtoty[levelmin+i+1]-llfac*shift[1], 
								offtotz[levelmin+i+1]-llfac*shift[2] );*/
				
				fine->fill_rand( rc, 1.0, x0[0]-fine->nx_/4, x0[1]-fine->ny_/4, x0[2]-fine->nz_/4 );
				
				if( i+levelmin+1 > (unsigned)lmingiven )
				{	
					if(i==0)
						fine->constrain( *top );
					else
						fine->constrain( *coarse );
				}				
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
			
			DensityGrid<real_t> top_save( *top );
			
			
			conv_param.lx = boxlength;
			conv_param.ly = boxlength;
			conv_param.lz = boxlength;
			conv_param.nx = top->nx_;
			conv_param.ny = top->ny_;
			conv_param.nz = top->nz_;
			conv_param.normalize = true;
			
			convolution::kernel *the_tf_kernel = the_kernel_creator->create( conv_param );
			
			
			//... 1) compute standard convolution for levelmin
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ) );
			//convolution::perform_filtered<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ) );
			top->copy( *delta.get_grid(levelmin) );
			
			//... 2) compute contribution to finer grids from non-refined region
			*top = top_save;
			top_save.clear();
			top->zero_subgrid(refh.offset(levelmin+i+1,0), refh.offset(levelmin+i+1,1), refh.offset(levelmin+i+1,2), 
							  refh.size(levelmin+i+1,0)/2, refh.size(levelmin+i+1,1)/2, refh.size(levelmin+i+1,2)/2 );

			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ) );
			//convolution::perform_filtered<real_t>( the_tf_kernel, reinterpret_cast<void*>( top->get_data_ptr() ) );
			
			
			
			meshvar_bnd delta_longrange( *delta.get_grid(levelmin) );
			top->copy( delta_longrange );
			delete the_tf_kernel;
			delete top;
			
			//... restrict these contributions to the next level
			delta.add_patch( refh.offset(levelmin+1,0), refh.offset(levelmin+1,1), refh.offset(levelmin+1,2), 
							refh.size(levelmin+1,0), refh.size(levelmin+1,1), refh.size(levelmin+1,2) );
			
			//mg_linear().prolong( delta_longrange, *delta.get_grid(levelmin+1) );
			mg_cubic().prolong( delta_longrange, *delta.get_grid(levelmin+1) );
			//mg_lin().prolong( delta_longrange, *delta.get_grid(levelmin+1) );
			
			
		}
		else
		{
			/**********************************************************************************************************\
			 *	multi-grid: intermediate sub-grids .....
			 \**********************************************************************************************************/ 
			
			std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmin+i << " ..." << std::endl;
			
			delta.add_patch( refh.offset(levelmin+i+1,0), refh.offset(levelmin+i+1,1), refh.offset(levelmin+i+1,2), 
							refh.size(levelmin+i+1,0), refh.size(levelmin+i+1,1), refh.size(levelmin+i+1,2) );
			
			//mg_linear().prolong( *delta.get_grid(levelmin+i), *delta.get_grid(levelmin+i+1) );
			mg_cubic().prolong( *delta.get_grid(levelmin+i), *delta.get_grid(levelmin+i+1) );
			
			real_t dx,lx,ly,lz;
			dx = boxlength/pow(2.0,levelmin+i);
			
			lx = dx * coarse->nx_; 
			ly = dx * coarse->ny_; 
			lz = dx * coarse->nz_; 
			
			real_t lmin = std::min( lx, std::min(ly,lz) );
			
			
			conv_param.lx = lx;
			conv_param.ly = ly;
			conv_param.lz = lz;
			conv_param.nx = coarse->nx_;
			conv_param.ny = coarse->ny_;
			conv_param.nz = coarse->nz_;
			conv_param.normalize = false;
			
			convolution::kernel *the_tf_kernel = the_kernel_creator->create( conv_param );
			
			PaddedDensitySubGrid<real_t> coarse_save( *coarse );
					
			//... 2) the inner region
			
			coarse->zero_boundary();
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
			//convolution::perform_filtered<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
			
			
			coarse->copy_add_unpad( *delta.get_grid(levelmin+i) );
			
			
			//... 3) the 'BC' for the next finer grid
			*coarse = coarse_save;

			coarse->zero_subgrid(refh.offset(levelmin+i+1,0), refh.offset(levelmin+i+1,1), refh.offset(levelmin+i+1,2), 
								 refh.size(levelmin+i+1,0)/2, refh.size(levelmin+i+1,1)/2, refh.size(levelmin+i+1,2)/2 );
			coarse->zero_boundary();
			
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
			//convolution::perform_filtered<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
			
			
			
			meshvar_bnd delta_longrange( *delta.get_grid(levelmin+i) );
			coarse->copy_unpad( delta_longrange );
			//mg_linear().prolong_add( delta_longrange, *delta.get_grid(levelmin+i+1) );
			mg_cubic().prolong_add( delta_longrange, *delta.get_grid(levelmin+i+1) );
			
			
#if 1
			//... FFT (isolated) boundary noise contribution
			//convolution::truncate( the_tf_kernel, cutoff*lmin, 1e-8*alpha );
			convolution::truncate_sharp( the_tf_kernel, cutoff*lmin );
			*coarse = coarse_save;
			//coarse_save.clear();
			coarse->zero_subgrid(0,0,0,coarse->nx_/2,coarse->ny_/2,coarse->nz_/2);
			coarse->subtract_oct_mean();
			
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
			//convolution::perform_filtered<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
			
			coarse->copy_add_unpad( *delta.get_grid(levelmin+i) );
#endif
			
#if 0
			*coarse = coarse_save;
			coarse->zero_boundary();
			coarse->subtract_oct_mean();
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
			//coarse->zero_subgrid(0,0,0,coarse->nx_/2,coarse->ny_/2,coarse->nz_/2);
			
			coarse->upload_bnd_add( *delta.get_grid(levelmin+i-1) );
#endif
			
			/**coarse = coarse_save;
			coarse->zero_subgrid(0,0,0,coarse->nx_/2,coarse->ny_/2,coarse->nz_/2);
			coarse->set_to_oct_mean();
			
			convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
			coarse->copy_subtract_unpad( *delta.get_grid(levelmin+i) );*/
			
			
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
		real_t lmin = std::min( lx, std::min(ly,lz) );
		
		
		conv_param.lx = lx;
		conv_param.ly = ly;
		conv_param.lz = lz;
		conv_param.nx = coarse->nx_;
		conv_param.ny = coarse->ny_;
		conv_param.nz = coarse->nz_;
		conv_param.normalize = false;
		
		//... with LR/SR splitting and full subgrid		
		PaddedDensitySubGrid<real_t> coarse_save( *coarse );
		
		convolution::kernel *the_tf_kernel = the_kernel_creator->create( conv_param );
		//... 2) the inner region
		
		coarse->zero_boundary();
		convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
		coarse->copy_add_unpad( *delta.get_grid(levelmax) );
		
		
		//... compute convolution for isolated grid
		//... 1) the padded boundary region
	
	
		
		
#if 1
		//... boundary correction
		convolution::truncate( the_tf_kernel, cutoff*lmin, alpha );
		//convolution::truncate_sharp( the_tf_kernel, cutoff*lmin );
		
		
		*coarse = coarse_save;
		coarse->zero_subgrid(0,0,0,coarse->nx_/2,coarse->ny_/2,coarse->nz_/2);
		coarse->subtract_oct_mean();
		
		convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
		coarse->copy_add_unpad( *delta.get_grid(levelmax) );

#endif
#if 0
		// coarse correction
		*coarse = coarse_save;
		coarse->zero_boundary();
		coarse->subtract_oct_mean();
		convolution::perform<real_t>( the_tf_kernel, reinterpret_cast<void*> (coarse->get_data_ptr()) );
		
		coarse->upload_bnd_add( *delta.get_grid(levelmax-1) );
		
		//	meshvar_bnd delta_longrange( *delta.get_grid(levelmax) );
		//coarse->copy_unpad( delta_longrange );
		//mg_straight().restrict_add( delta_longrange, *delta.get_grid(levelmax-1) );
#endif	
		
		
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
	
	//... enforce mean condition
	//for( int i=levelmin; i<(int)levelmax; ++i )
		//	enforce_mean( (*delta.get_grid(i+1)), (*delta.get_grid(i)) );
	 
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

