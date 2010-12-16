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

class enzo_output_plugin : public output_plugin
{
protected:
	
	struct patch_header{
		int component_rank;
		int component_size;
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
	
	
	void dump_grid_data( std::string fieldname, const grid_hierarchy& gh, double factor = 1.0, double add = 0.0 )
	{
		char enzoname[256], filename[256];
		
		for(unsigned ilevel=levelmin_; ilevel<=levelmax_; ++ilevel )
		{
			std::vector<int> ng;
			ng.push_back( gh.get_grid(ilevel)->size(0) );
			ng.push_back( gh.get_grid(ilevel)->size(1) );
			ng.push_back( gh.get_grid(ilevel)->size(2) );
			
			//... need to copy data because we need to get rid of the ghost zones
			std::vector<double> data;
			data.reserve( ng[0]*ng[1]*ng[2] );
			
			for( int k=0; k<ng[2]; ++k )
				for( int j=0; j<ng[1]; ++j )
					for( int i=0; i<ng[0]; ++i )
						data.push_back( (add+(*gh.get_grid(ilevel))(i,j,k))*factor );
			
			sprintf( enzoname, "%s.%d", fieldname.c_str(), ilevel-levelmin_ );
			sprintf( filename, "%s/%s", fname_.c_str(), enzoname );
			
			HDFCreateFile( filename );
			write_sim_header( filename, the_sim_header );
			HDFWriteDataset3Ds( filename, enzoname, reinterpret_cast<unsigned*>(&ng[0]), data );
			
			//... header data for the patch
			patch_header ph;
			
			ph.component_rank	= 1;
			ph.component_size	= ng[0]*ng[1]*ng[2];
			ph.dimensions		= ng;
			ph.rank				= 3;
			
			ph.top_grid_dims.assign(3,(int)pow(2,levelmin_));
			
			//... offset_abs is in units of the current level cell size
			
			double rfac = 1.0/pow(2,ilevel-levelmin_);
			
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
		bool align_top			= cf.getValueSafe<bool>( "setup", "align_top", true );
		
		if( !align_top )
			throw std::runtime_error("ENZO output plug-in requires that \'align_top=true\'!");
		
		the_sim_header.dimensions.push_back( (int)pow(2,levelmin_) );
		the_sim_header.dimensions.push_back( (int)pow(2,levelmin_) );
		the_sim_header.dimensions.push_back( (int)pow(2,levelmin_) );
		
		the_sim_header.offset.push_back( 0 );
		the_sim_header.offset.push_back( 0 );
		the_sim_header.offset.push_back( 0 );
		
		the_sim_header.a_start		= 1.0/(1.0+cf.getValue<double>("setup","zstart"));
		the_sim_header.dx			= cf.getValue<double>("setup","boxlength")/the_sim_header.dimensions[0]/(cf.getValue<double>("cosmology","H0")*0.01); // not sure?!?
		the_sim_header.h0			= cf.getValue<double>("cosmology","H0")*0.01;
		
		if( bhave_hydro )
			the_sim_header.omega_b		= cf.getValue<double>("cosmology","Omega_b");
		else
			the_sim_header.omega_b		= 0.0;
		
		the_sim_header.omega_m		= cf.getValue<double>("cosmology","Omega_m");
		the_sim_header.omega_v		= cf.getValue<double>("cosmology","Omega_L");
		the_sim_header.vfact		= cf.getValue<double>("cosmology","vfact")*the_sim_header.h0;   //.. need to multiply by h, ENZO wants this factor for non h-1 units
		
	}
	
	~enzo_output_plugin()
	{ }
	
	void write_dm_mass( const grid_hierarchy& gh )
	{	/* do nothing, not needed */	}
	
	
	void write_dm_density( const grid_hierarchy& gh )
	{	/* write the parameter file data */	
	
		bool bhave_hydro = cf_.getValue<bool>("setup","baryons");
		char filename[256];
		unsigned nbase = (unsigned)pow(2,levelmin_);
		
		sprintf( filename, "%s/parameter_file.txt", fname_.c_str() );
		
		std::ofstream ofs( filename, std::ios::trunc );
		
		ofs
			<< "#\n"
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
			<< "CosmologySimulationNumberOfInitialGrids  = " << 1+levelmax_-levelmin_ << "\n";
		
		
		//... only for additionally refined grids
		for( unsigned ilevel = 0; ilevel< levelmax_-levelmin_; ++ilevel )
		{
			double h = 1.0/pow(2,levelmin_+1+ilevel);
			
			ofs
			
			<< "CosmologySimulationGridDimension[" << 1+ilevel << "]      = " 
				<< std::setw(16) << gh.size( levelmin_+ilevel+1, 0 ) << " "
				<< std::setw(16) << gh.size( levelmin_+ilevel+1, 1 ) << " "
				<< std::setw(16) << gh.size( levelmin_+ilevel+1, 2 ) << "\n"
			
			<< "CosmologySimulationGridLeftEdge[" << 1+ilevel << "]       = "
				<< std::setw(16) << std::setprecision(10) << h*gh.offset_abs(levelmin_+ilevel+1, 0) << " "
				<< std::setw(16) << std::setprecision(10) << h*gh.offset_abs(levelmin_+ilevel+1, 1) << " "
				<< std::setw(16) << std::setprecision(10) << h*gh.offset_abs(levelmin_+ilevel+1, 2) << "\n"
			
			<< "CosmologySimulationGridRightEdge[" << 1+ilevel << "]      = "
				<< std::setw(16) << std::setprecision(10) << h*(gh.offset_abs(levelmin_+ilevel+1, 0)+gh.size( levelmin_+ilevel+1, 0 )) << " "
				<< std::setw(16) << std::setprecision(10) << h*(gh.offset_abs(levelmin_+ilevel+1, 1)+gh.size( levelmin_+ilevel+1, 1 )) << " "
				<< std::setw(16) << std::setprecision(10) << h*(gh.offset_abs(levelmin_+ilevel+1, 2)+gh.size( levelmin_+ilevel+1, 2 )) << "\n"
			
			<< "CosmologySimulationGridLevel[" << 1+ilevel << "]          = " << 1+ilevel << "\n";
		}
		
		
		
		
		
	
	}
	
	
	void write_dm_velocity( int coord, const grid_hierarchy& gh )
	{
		char enzoname[256];
		sprintf( enzoname, "ParticleVelocities_%c", (char)('x'+coord) );
		
		double vunit = 1.0/(1.225e2*sqrt(the_sim_header.omega_m/the_sim_header.a_start));// /the_sim_header.h0;
		
		dump_grid_data( enzoname, gh, vunit );
	}
	
	
	void write_dm_position( int coord, const grid_hierarchy& gh )
	{
		char fieldname[256];
		sprintf( fieldname, "ParticleDisplacements_%c", (char)('x'+coord) );
		
		char enzoname[256], filename[256];
		
		for(unsigned ilevel=levelmin_; ilevel<=levelmax_; ++ilevel )
		{
			std::vector<int> ng;
			ng.push_back( gh.get_grid(ilevel)->size(0) );
			ng.push_back( gh.get_grid(ilevel)->size(1) );
			ng.push_back( gh.get_grid(ilevel)->size(2) );
			
			//... need to copy data because we need to get rid of the ghost zones
			std::vector<double> data;
			data.reserve( ng[0]*ng[1]*ng[2] );

			//double dx = 1.0/(1<<ilevel);

			for( int k=0; k<ng[2]; ++k )
				for( int j=0; j<ng[1]; ++j )
					for( int i=0; i<ng[0]; ++i )
						data.push_back( (*gh.get_grid(ilevel))(i,j,k) );
						
			//--------------------------------------------------------------------------
			
			sprintf( enzoname, "%s.%d", fieldname, ilevel-levelmin_ );
			sprintf( filename, "%s/%s", fname_.c_str(), enzoname );
			
			HDFCreateFile( filename );
			write_sim_header( filename, the_sim_header );
			HDFWriteDataset3Ds( filename, enzoname, reinterpret_cast<unsigned*>(&ng[0]), data );
			
			//--------------------------------------------------------------------------
			//... header data for the patch
			patch_header ph;
			
			ph.component_rank	= 1;
			ph.component_size	= ng[0]*ng[1]*ng[2];
			ph.dimensions		= ng;
			ph.rank				= 3;
			
			ph.top_grid_dims.assign(3,(int)pow(2,levelmin_));
			
			//... offset_abs is in units of the current level cell size
			
			double rfac = 1.0/pow(2,ilevel-levelmin_);
			
			ph.top_grid_start.push_back( (int)(gh.offset_abs(ilevel, 0)*rfac) );
			ph.top_grid_start.push_back( (int)(gh.offset_abs(ilevel, 1)*rfac) );
			ph.top_grid_start.push_back( (int)(gh.offset_abs(ilevel, 2)*rfac) );
			
			ph.top_grid_end.push_back( ph.top_grid_start[0] + (int)(ng[0]*rfac) );
			ph.top_grid_end.push_back( ph.top_grid_start[1] + (int)(ng[1]*rfac) );
			ph.top_grid_end.push_back( ph.top_grid_start[2] + (int)(ng[2]*rfac) );
			
			write_patch_header( filename, enzoname, ph );
			//--------------------------------------------------------------------------
		}
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
	{	}
		
	
};

namespace{
	output_plugin_creator_concrete<enzo_output_plugin> creator("enzo");
}

#endif

