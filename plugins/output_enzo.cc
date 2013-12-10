/*

output_enzo.cc - This file is part of MUSIC -
a code to generate multi-scale initial conditions
for cosmological simulations

Copyright (C) 2010  Oliver Hahn

*/

#ifdef HAVE_HDF5

#include <sys/types.h>
#include <sys/stat.h>

#include "output.hh"

#include "HDF_IO.hh"

#define MAX_SLAB_SIZE   268435456  // = 256 MBytes


class enzo_output_plugin : public output_plugin
{
protected:

  struct patch_header{
    int component_rank;
    size_t component_size;
    std::vector<int> dimensions;
    int rank;
    std::vector<int> top_grid_dims;
    std::vector<int> top_grid_end;
    std::vector<int> top_grid_start;
  };
  
  struct sim_header{
    std::vector<int> dimensions;
    std::vector<int> offset;
    float a_start;
    float dx;
    float h0;
    float omega_b;
    float omega_m;
    float omega_v;
    float vfact;
  };
  
  
  sim_header the_sim_header;
  
  void write_sim_header( std::string fname, const sim_header& h )
  {
    HDFWriteGroupAttribute( fname, "/", "Dimensions", h.dimensions );
    HDFWriteGroupAttribute( fname, "/", "Offset", h.offset );
    HDFWriteGroupAttribute( fname, "/", "a_start", h.a_start );
    HDFWriteGroupAttribute( fname, "/", "dx", h.dx );
    HDFWriteGroupAttribute( fname, "/", "h0", h.h0 );
    HDFWriteGroupAttribute( fname, "/", "omega_b", h.omega_b );
    HDFWriteGroupAttribute( fname, "/", "omega_m", h.omega_m );
    HDFWriteGroupAttribute( fname, "/", "omega_v", h.omega_v );
    HDFWriteGroupAttribute( fname, "/", "vfact", h.vfact );
  }
  
  void write_patch_header( std::string fname, std::string dsetname, const patch_header& h )
  {
    HDFWriteDatasetAttribute( fname, dsetname, "Component_Rank", h.component_rank );
    HDFWriteDatasetAttribute( fname, dsetname, "Component_Size", h.component_size );
    HDFWriteDatasetAttribute( fname, dsetname, "Dimensions", h.dimensions );
    HDFWriteDatasetAttribute( fname, dsetname, "Rank", h.rank );
    HDFWriteDatasetAttribute( fname, dsetname, "TopGridDims", h.top_grid_dims );
    HDFWriteDatasetAttribute( fname, dsetname, "TopGridEnd", h.top_grid_end );
    HDFWriteDatasetAttribute( fname, dsetname, "TopGridStart", h.top_grid_start );
  }
  
  void dump_mask( const grid_hierarchy& gh )
  {
    char enzoname[256], filename[256];
    std::string fieldname("RefinementMask");
    
    for(unsigned ilevel=levelmin_; ilevel<=levelmax_; ++ilevel )
      {
	std::vector<int> ng, ng_fortran;
	ng.push_back( gh.get_grid(ilevel)->size(0) );
	ng.push_back( gh.get_grid(ilevel)->size(1) );
	ng.push_back( gh.get_grid(ilevel)->size(2) );
	
	ng_fortran.push_back( gh.get_grid(ilevel)->size(2) );
	ng_fortran.push_back( gh.get_grid(ilevel)->size(1) );
	ng_fortran.push_back( gh.get_grid(ilevel)->size(0) );
	
	//... need to copy data because we need to get rid of the ghost zones
	//... write in slabs if data is more than MAX_SLAB_SIZE (default 128 MB)
	
	//... full 3D block size
	size_t all_data_size = (size_t)ng[0] * (size_t)ng[1] * (size_t)ng[2];
	
	//... write in slabs of MAX_SLAB_SIZE unless all_data_size is anyway smaller
	size_t max_slab_size = std::min((size_t)MAX_SLAB_SIZE/sizeof(double), all_data_size );
	
	//... but one slab hast to be at least the size of one slice
	max_slab_size = std::max(((size_t)ng[0] * (size_t)ng[1]), max_slab_size );
	
	//... number of slices in one slab
	size_t slices_in_slab = (size_t)((double)max_slab_size / ((size_t)ng[0] * (size_t)ng[1]));
	
	size_t nsz[3] = { ng[2], ng[1], ng[0] };
	
	if( levelmin_ != levelmax_ )
	  sprintf( enzoname, "%s.%d", fieldname.c_str(), ilevel-levelmin_ );
	else
	  sprintf( enzoname, "%s", fieldname.c_str() );
	
	sprintf( filename, "%s/%s", fname_.c_str(), enzoname );
	
	HDFCreateFile( filename );
	write_sim_header( filename, the_sim_header );
	
	//... create full array in file
	HDFHyperslabWriter3Ds<int> *slab_writer = new HDFHyperslabWriter3Ds<int>( filename, enzoname, nsz );
	
	//... create buffer
	int *data_buf = new int[ slices_in_slab * (size_t)ng[0] * (size_t)ng[1] ];

	//... write slice by slice
	size_t slices_written = 0;
	while( slices_written < (size_t)ng[2] )
	  {
	    slices_in_slab = std::min( (size_t)ng[2]-slices_written, slices_in_slab );
	    
            #pragma omp parallel for
	    for( int k=0; k<(int)slices_in_slab; ++k )
	      for( int j=0; j<ng[1]; ++j )
		for( int i=0; i<ng[0]; ++i )
		  {
		    int mask_val = -1;
		    
		    if( gh.is_in_mask(ilevel,i,j,k+slices_written) )
		      {
			if( gh.is_refined(ilevel,i,j,k+slices_written) )
			  mask_val = 1;
			else
			  mask_val = 0;
		      }
		    data_buf[ (size_t)(k*ng[1]+j)*(size_t)ng[0]+(size_t)i ] = mask_val;
		  }
	    
	    size_t count[3], offset[3];
	    
	    count[0] = slices_in_slab;
	    count[1] = ng[1];
	    count[2] = ng[0];
	    
	    offset[0] = slices_written;;
	    offset[1] = 0;
	    offset[2] = 0;
	    
	    slab_writer->write_slab( data_buf, count, offset );
	    slices_written += slices_in_slab;
	    
	  }
	

	delete[] data_buf;
	delete slab_writer;

	//... header data for the patch
	patch_header ph;

	ph.component_rank       = 1;
	ph.component_size       = (size_t)ng[0]*(size_t)ng[1]*(size_t)ng[2];
	ph.dimensions           = ng;
	ph.rank                 = 3;

	ph.top_grid_dims.assign(3, 1<<levelmin_);

	//... offset_abs is in units of the current level cell size
	
	double rfac = 1.0/(1<<(ilevel-levelmin_));
	
	ph.top_grid_start.push_back( (int)(gh.offset_abs(ilevel, 0)*rfac) );
	ph.top_grid_start.push_back( (int)(gh.offset_abs(ilevel, 1)*rfac) );
	ph.top_grid_start.push_back( (int)(gh.offset_abs(ilevel, 2)*rfac) );
	
	ph.top_grid_end.push_back( ph.top_grid_start[0] + (int)(ng[0]*rfac) );
	ph.top_grid_end.push_back( ph.top_grid_start[1] + (int)(ng[1]*rfac) );
	ph.top_grid_end.push_back( ph.top_grid_start[2] + (int)(ng[2]*rfac) );
	
	write_patch_header( filename, enzoname, ph );
      }
  }

  void dump_grid_data( std::string fieldname, const grid_hierarchy& gh, double factor = 1.0, double add = 0.0 )
  {
    char enzoname[256], filename[256];
    
    for(unsigned ilevel=levelmin_; ilevel<=levelmax_; ++ilevel )
      {
	std::vector<int> ng, ng_fortran;
	ng.push_back( gh.get_grid(ilevel)->size(0) );
	ng.push_back( gh.get_grid(ilevel)->size(1) );
	ng.push_back( gh.get_grid(ilevel)->size(2) );
	
	ng_fortran.push_back( gh.get_grid(ilevel)->size(2) );
	ng_fortran.push_back( gh.get_grid(ilevel)->size(1) );
	ng_fortran.push_back( gh.get_grid(ilevel)->size(0) );
	
	
	//... need to copy data because we need to get rid of the ghost zones
	//... write in slabs if data is more than MAX_SLAB_SIZE (default 128 MB)
	
	//... full 3D block size
	size_t all_data_size = (size_t)ng[0] * (size_t)ng[1] * (size_t)ng[2];
	
	//... write in slabs of MAX_SLAB_SIZE unless all_data_size is anyway smaller
	size_t max_slab_size = std::min((size_t)MAX_SLAB_SIZE/sizeof(double), all_data_size );
	
	//... but one slab hast to be at least the size of one slice
	max_slab_size = std::max(((size_t)ng[0] * (size_t)ng[1]), max_slab_size );
	
	//... number of slices in one slab
	size_t slices_in_slab = (size_t)((double)max_slab_size / ((size_t)ng[0] * (size_t)ng[1]));
	
	size_t nsz[3] = { ng[2], ng[1], ng[0] };
	
	if( levelmin_ != levelmax_ )
	  sprintf( enzoname, "%s.%d", fieldname.c_str(), ilevel-levelmin_ );
	else
	  sprintf( enzoname, "%s", fieldname.c_str() );
	
	sprintf( filename, "%s/%s", fname_.c_str(), enzoname );
	
	HDFCreateFile( filename );
	write_sim_header( filename, the_sim_header );
	
#ifdef SINGLE_PRECISION
	//... create full array in file
	HDFHyperslabWriter3Ds<float> *slab_writer = new HDFHyperslabWriter3Ds<float>( filename, enzoname, nsz );
	
	//... create buffer
	float *data_buf = new float[ slices_in_slab * (size_t)ng[0] * (size_t)ng[1] ];
#else
	//... create full array in file
	HDFHyperslabWriter3Ds<double> *slab_writer = new HDFHyperslabWriter3Ds<double>( filename, enzoname, nsz );
	
	//... create buffer
	double *data_buf = new double[ slices_in_slab * (size_t)ng[0] * (size_t)ng[1] ];
#endif
	
	//... write slice by slice
	size_t slices_written = 0;
	while( slices_written < (size_t)ng[2] )
	  {
	    slices_in_slab = std::min( (size_t)ng[2]-slices_written, slices_in_slab );
	    
            #pragma omp parallel for
	    for( int k=0; k<(int)slices_in_slab; ++k )
	      for( int j=0; j<ng[1]; ++j )
		for( int i=0; i<ng[0]; ++i )
		  data_buf[ (size_t)(k*ng[1]+j)*(size_t)ng[0]+(size_t)i ] =
		    (add+(*gh.get_grid(ilevel))(i,j,k+slices_written))*factor;
	    
	    size_t count[3], offset[3];
	    
	    count[0] = slices_in_slab;
	    count[1] = ng[1];
	    count[2] = ng[0];
	    
	    offset[0] = slices_written;;
	    offset[1] = 0;
	    offset[2] = 0;
	    
	    slab_writer->write_slab( data_buf, count, offset );
	    slices_written += slices_in_slab;
	    
	  }
	
	//... free buffer
	delete[] data_buf;
	
	//... finalize writing and close dataset
	delete slab_writer;
	
	
	//... header data for the patch
	patch_header ph;
	
	ph.component_rank       = 1;
	ph.component_size       = (size_t)ng[0]*(size_t)ng[1]*(size_t)ng[2];
	ph.dimensions           = ng;
	ph.rank                         = 3;
	
	ph.top_grid_dims.assign(3, 1<<levelmin_);
	
	//... offset_abs is in units of the current level cell size
	
	double rfac = 1.0/(1<<(ilevel-levelmin_));
	
	ph.top_grid_start.push_back( (int)(gh.offset_abs(ilevel, 0)*rfac) );
	ph.top_grid_start.push_back( (int)(gh.offset_abs(ilevel, 1)*rfac) );
	ph.top_grid_start.push_back( (int)(gh.offset_abs(ilevel, 2)*rfac) );
	
	ph.top_grid_end.push_back( ph.top_grid_start[0] + (int)(ng[0]*rfac) );
	ph.top_grid_end.push_back( ph.top_grid_start[1] + (int)(ng[1]*rfac) );
	ph.top_grid_end.push_back( ph.top_grid_start[2] + (int)(ng[2]*rfac) );
	
	write_patch_header( filename, enzoname, ph );
      }
  }
  
public:
  
  enzo_output_plugin( config_file& cf )
  : output_plugin( cf )
  {
    if( mkdir( fname_.c_str(), 0777 ) )
      {
	perror( fname_.c_str() );
	throw std::runtime_error("Error in enzo_output_plugin!");
      }
    
    bool bhave_hydro = cf_.getValue<bool>("setup","baryons");
    bool align_top                  = cf.getValueSafe<bool>( "setup", "align_top", false );
    
    if( !align_top )
      LOGWARN("Old ENZO versions may require \'align_top=true\'!");
    
    the_sim_header.dimensions.push_back( 1<<levelmin_ );
    the_sim_header.dimensions.push_back( 1<<levelmin_ );
    the_sim_header.dimensions.push_back( 1<<levelmin_ );
    
    the_sim_header.offset.push_back( 0 );
    the_sim_header.offset.push_back( 0 );
    the_sim_header.offset.push_back( 0 );
    
    the_sim_header.a_start          = 1.0/(1.0+cf.getValue<double>("setup","zstart"));
    the_sim_header.dx                       = cf.getValue<double>("setup","boxlength")/the_sim_header.dimensions[0]/(cf.getValue<double>("cosmology","H0")*0.01); // not sure?!?
    the_sim_header.h0                       = cf.getValue<double>("cosmology","H0")*0.01;
    
    if( bhave_hydro )
      the_sim_header.omega_b          = cf.getValue<double>("cosmology","Omega_b");
    else
      the_sim_header.omega_b          = 0.0;
    
    the_sim_header.omega_m          = cf.getValue<double>("cosmology","Omega_m");
    the_sim_header.omega_v          = cf.getValue<double>("cosmology","Omega_L");
    the_sim_header.vfact            = cf.getValue<double>("cosmology","vfact")*the_sim_header.h0;   //.. need to multiply by h, ENZO wants this factor for non h-1 units
    
  }
  
  ~enzo_output_plugin()
  { }
  
  void write_dm_mass( const grid_hierarchy& gh )
  {       /* do nothing, not needed */    }
  
  
  void write_dm_density( const grid_hierarchy& gh )
  {       /* write the parameter file data */
    
    bool bhave_hydro = cf_.getValue<bool>("setup","baryons");
    double refine_region_fraction  = cf_.getValueSafe<double>( "output", "enzo_refine_region_fraction", 0.8 );
    char filename[256];
    unsigned nbase = (unsigned)pow(2,levelmin_);

    // write out the refinement masks
    dump_mask( gh );

    // write out a parameter file
    
    sprintf( filename, "%s/parameter_file.txt", fname_.c_str() );
    
    std::ofstream ofs( filename, std::ios::trunc );
    
    ofs
      << "# Relevant Section of Enzo Paramter File (NOT COMPLETE!) \n"
      << "ProblemType                              = 30      // cosmology simulation\n"
      << "TopGridRank                              = 3\n"
      << "TopGridDimensions                        = " << nbase << " " << nbase << " " << nbase << "\n"
      << "SelfGravity                              = 1       // gravity on\n"
      << "TopGridGravityBoundary                   = 0       // Periodic BC for gravity\n"
      << "LeftFaceBoundaryCondition                = 3 3 3   // same for fluid\n"
      << "RightFaceBoundaryCondition               = 3 3 3\n"
      << "RefineBy                                 = 2\n"
      << "\n"
      << "#\n";
    
    if( bhave_hydro )
      ofs
	<< "CosmologySimulationOmegaBaryonNow        = " << the_sim_header.omega_b << "\n"
	<< "CosmologySimulationOmegaCDMNow           = " << the_sim_header.omega_m-the_sim_header.omega_b << "\n";
    else
      ofs
	<< "CosmologySimulationOmegaBaryonNow        = " << 0.0 << "\n"
	<< "CosmologySimulationOmegaCDMNow           = " << the_sim_header.omega_m << "\n";
    
    if( bhave_hydro )
      ofs
	<< "CosmologySimulationDensityName           = GridDensity\n"
	<< "CosmologySimulationVelocity1Name         = GridVelocities_x\n"
	<< "CosmologySimulationVelocity2Name         = GridVelocities_y\n"
	<< "CosmologySimulationVelocity3Name         = GridVelocities_z\n";
    
    ofs
      << "CosmologySimulationCalculatePositions    = 1\n"
      << "CosmologySimulationParticleVelocity1Name = ParticleVelocities_x\n"
      << "CosmologySimulationParticleVelocity2Name = ParticleVelocities_y\n"
      << "CosmologySimulationParticleVelocity3Name = ParticleVelocities_z\n"
      << "CosmologySimulationParticleDisplacement1Name = ParticleDisplacements_x\n"
      << "CosmologySimulationParticleDisplacement2Name = ParticleDisplacements_y\n"
      << "CosmologySimulationParticleDisplacement3Name = ParticleDisplacements_z\n"
      << "\n"
      << "#\n"
      << "#  define cosmology parameters\n"
      << "#\n"
      << "ComovingCoordinates                      = 1       // Expansion ON\n"
      << "CosmologyOmegaMatterNow                  = " << the_sim_header.omega_m << "\n"
      << "CosmologyOmegaLambdaNow                  = " << the_sim_header.omega_v << "\n"
      << "CosmologyHubbleConstantNow               = " << the_sim_header.h0 << "     // in 100 km/s/Mpc\n"
      << "CosmologyComovingBoxSize                 = " << cf_.getValue<double>("setup","boxlength") << "    // in Mpc/h\n"
      << "CosmologyMaxExpansionRate                = 0.015   // maximum allowed delta(a)/a\n"
      << "CosmologyInitialRedshift                 = " << cf_.getValue<double>("setup","zstart") << "      //\n"
      << "CosmologyFinalRedshift                   = 0       //\n"
      << "GravitationalConstant                    = 1       // this must be true for cosmology\n"
      << "#\n"
      << "#\n"
      << "ParallelRootGridIO                       = 1\n"
      << "ParallelParticleIO                       = 1\n"
      << "PartitionNestedGrids                     = 1\n"
      << "CosmologySimulationNumberOfInitialGrids  = " << 1+levelmax_-levelmin_ << "\n";
    
    
    int num_prec = 10;
    
    if( levelmax_ > 15 )
      num_prec = 17;
    
    //... only for additionally refined grids
    for( unsigned ilevel = 0; ilevel< levelmax_-levelmin_; ++ilevel )
      {
	double h = 1.0/(1<<(levelmin_+1+ilevel));
	
	ofs
	  
	  << "CosmologySimulationGridDimension[" << 1+ilevel << "]      = "
	  << std::setw(16) << gh.size( levelmin_+ilevel+1, 0 ) << " "
	  << std::setw(16) << gh.size( levelmin_+ilevel+1, 1 ) << " "
	  << std::setw(16) << gh.size( levelmin_+ilevel+1, 2 ) << "\n"
	  
	  << "CosmologySimulationGridLeftEdge[" << 1+ilevel << "]       = "
	  << std::setw(num_prec+6) << std::setprecision(num_prec) << h*gh.offset_abs(levelmin_+ilevel+1, 0) << " "
	  << std::setw(num_prec+6) << std::setprecision(num_prec) << h*gh.offset_abs(levelmin_+ilevel+1, 1) << " "
	  << std::setw(num_prec+6) << std::setprecision(num_prec) << h*gh.offset_abs(levelmin_+ilevel+1, 2) << "\n"
	  
	  << "CosmologySimulationGridRightEdge[" << 1+ilevel << "]      = "
	  << std::setw(num_prec+6) << std::setprecision(num_prec) << h*(gh.offset_abs(levelmin_+ilevel+1, 0)+gh.size( levelmin_+ilevel+1, 0 )) << " "
	  << std::setw(num_prec+6) << std::setprecision(num_prec) << h*(gh.offset_abs(levelmin_+ilevel+1, 1)+gh.size( levelmin_+ilevel+1, 1 )) << " "
	  << std::setw(num_prec+6) << std::setprecision(num_prec) << h*(gh.offset_abs(levelmin_+ilevel+1, 2)+gh.size( levelmin_+ilevel+1, 2 )) << "\n"
	  
	  << "CosmologySimulationGridLevel[" << 1+ilevel << "]          = " << 1+ilevel << "\n";
      }
    
    if( levelmin_ != levelmax_ )
      {
	double h = 1.0/(1<<levelmax_);
	
	
	double cen[3],le[3],re[3];
	for (int i=0;i<3;i++)
	  {
	    cen[i] = gh.offset_abs(levelmax_, i)+gh.size( levelmax_, i )/2;
	    le[i]  = cen[i]-refine_region_fraction*gh.size( levelmax_,i)/2;
	    re[i]  = le[i] +refine_region_fraction*gh.size( levelmax_, i);
	  }
	
	
	ofs
	  << "#\n"
	  << "# region allowed for further refinement\n"
	  << "#\n"
	  //                << "RefineRegionAutoAdjust                   = 1\n"
	  << "RefineRegionLeftEdge                     = "
	  << std::setw(num_prec+6) << std::setprecision(num_prec) << h*le[0] << " "
	  << std::setw(num_prec+6) << std::setprecision(num_prec) << h*le[1] << " "
	  << std::setw(num_prec+6) << std::setprecision(num_prec) << h*le[2] << "\n"
	  << "RefineRegionRightEdge                     = "
	  << std::setw(num_prec+6) << std::setprecision(num_prec) << h*re[0] << " "
	  << std::setw(num_prec+6) << std::setprecision(num_prec) << h*re[1] << " "
	  << std::setw(num_prec+6) << std::setprecision(num_prec) << h*re[2]<< "\n";
      }
    
    
    // determine density maximum and minimum location
    real_t rhomax = -1e30, rhomin = 1e30;
    double loc_rhomax[3] = {0.0,0.0,0.0}, loc_rhomin[3] = {0.0,0.0,0.0};
    int lvl_rhomax = 0, lvl_rhomin = 0;
    real_t rhomax_lm = -1e30, rhomin_lm = 1e30;
    double loc_rhomax_lm[3] = {0.0,0.0,0.0}, loc_rhomin_lm[3] = {0.0,0.0,0.0};
    
    
    for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
      for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
	for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
	  for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
	    if( ! gh.is_refined(ilevel,i,j,k) )
	      {
		real_t rho = (*gh.get_grid(ilevel))(i,j,k);
		
		if( rho > rhomax )
		  {
		    rhomax = rho;
		    lvl_rhomax = ilevel;
		    gh.cell_pos(ilevel, i, j, k, loc_rhomax);
		  }
		
		if( rho < rhomin )
		  {
		    rhomin = rho;
		    lvl_rhomin = ilevel;
		    gh.cell_pos(ilevel, i, j, k, loc_rhomin);
		  }
		
		if( ilevel == (int)gh.levelmax() )
		  {
		    if( rho > rhomax_lm )
		      {
			rhomax_lm = rho;
			gh.cell_pos(ilevel, i, j, k, loc_rhomax_lm);
		      }
		    
		    if( rho < rhomin_lm )
		      {
			rhomin_lm = rho;
			gh.cell_pos(ilevel, i, j, k, loc_rhomin_lm);
		      }
		  }
	      }
    
    double h = 1.0/(1<<levelmin_);
    double shift[3];
    shift[0] = -(double)cf_.getValue<int>( "setup", "shift_x" )*h;
    shift[1] = -(double)cf_.getValue<int>( "setup", "shift_y" )*h;
    shift[2] = -(double)cf_.getValue<int>( "setup", "shift_z" )*h;
    
    if( gh.levelmin() != gh.levelmax() )
      {
	LOGINFO("Global density extrema: ");
	LOGINFO("  minimum: delta=%f at (%f,%f,%f) (level=%d)",rhomin,loc_rhomin[0],loc_rhomin[1],loc_rhomin[2],lvl_rhomin);
	LOGINFO("       shifted back at (%f,%f,%f)",loc_rhomin[0]+shift[0],loc_rhomin[1]+shift[1],loc_rhomin[2]+shift[2]);
	LOGINFO("  maximum: delta=%f at (%f,%f,%f) (level=%d)",rhomax,loc_rhomax[0],loc_rhomax[1],loc_rhomax[2],lvl_rhomax);
	LOGINFO("       shifted back at (%f,%f,%f)",loc_rhomax[0]+shift[0],loc_rhomax[1]+shift[1],loc_rhomax[2]+shift[2]);
	
	LOGINFO("Density extrema on finest level: ");
	LOGINFO("  minimum: delta=%f at (%f,%f,%f)",rhomin_lm,loc_rhomin_lm[0],loc_rhomin_lm[1],loc_rhomin_lm[2]);
	LOGINFO("       shifted back at (%f,%f,%f)",loc_rhomin_lm[0]+shift[0],loc_rhomin_lm[1]+shift[1],loc_rhomin_lm[2]+shift[2]);
	LOGINFO("  maximum: delta=%f at (%f,%f,%f)",rhomax_lm,loc_rhomax_lm[0],loc_rhomax_lm[1],loc_rhomax_lm[2]);
	LOGINFO("       shifted back at (%f,%f,%f)",loc_rhomax_lm[0]+shift[0],loc_rhomax_lm[1]+shift[1],loc_rhomax_lm[2]+shift[2]);
	
      }else{
      LOGINFO("Global density extrema: ");
      LOGINFO("  minimum: delta=%f at (%f,%f,%f)",rhomin,loc_rhomin[0],loc_rhomin[1],loc_rhomin[2]);
      LOGINFO("       shifted back at (%f,%f,%f)",loc_rhomin[0]+shift[0],loc_rhomin[1]+shift[1],loc_rhomin[2]+shift[2]);
      LOGINFO("  maximum: delta=%f at (%f,%f,%f)",rhomax,loc_rhomax[0],loc_rhomax[1],loc_rhomax[2]);
      LOGINFO("       shifted back at (%f,%f,%f)",loc_rhomax[0]+shift[0],loc_rhomax[1]+shift[1],loc_rhomax[2]+shift[2]);
      
    }
    
  }
  
  
  void write_dm_velocity( int coord, const grid_hierarchy& gh )
  {
    char enzoname[256];
    sprintf( enzoname, "ParticleVelocities_%c", (char)('x'+coord) );
    
    double vunit = 1.0/(1.225e2*sqrt(the_sim_header.omega_m/the_sim_header.a_start));
    
    dump_grid_data( enzoname, gh, vunit );
  }
  
  
  void write_dm_position( int coord, const grid_hierarchy& gh )
  {
    char enzoname[256];
    sprintf( enzoname, "ParticleDisplacements_%c", (char)('x'+coord) );
    
    dump_grid_data( enzoname, gh );
  }
  
  void write_dm_potential( const grid_hierarchy& gh )
  { }
  
  void write_gas_potential( const grid_hierarchy& gh )
  { }
  
  
  void write_gas_velocity( int coord, const grid_hierarchy& gh )
  {
    double vunit = 1.0/(1.225e2*sqrt(the_sim_header.omega_m/the_sim_header.a_start));
    
    char enzoname[256];
    sprintf( enzoname, "GridVelocities_%c", (char)('x'+coord) );
    dump_grid_data( enzoname, gh, vunit );
  }
  
  
  void write_gas_position( int coord, const grid_hierarchy& gh )
  {
    /* do nothing, not needed */
  }
  
  
  void write_gas_density( const grid_hierarchy& gh )
  {
    
    char enzoname[256];
    sprintf( enzoname, "GridDensity" );
    dump_grid_data( enzoname, gh, the_sim_header.omega_b/the_sim_header.omega_m, 1.0 );
  }
  
  
  void finalize( void )
  {       }
  
  
};

namespace{
  output_plugin_creator_concrete<enzo_output_plugin> creator("enzo");
}

#endif


