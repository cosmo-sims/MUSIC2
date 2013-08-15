/*
 * output_arepo.cc - This file is part of MUSIC -
 * a code to generate multi-scale initial conditions 
 * for cosmological simulations
 * 
 * Copyright (C) 2010  Oliver Hahn
 * 
 * Plugin: Dylan Nelson (dnelson@cfa.harvard.edu)
 */
 
#ifndef HAVE_HDF5

#include "output.hh"

class arepo_output_plugin : public output_plugin
{
public:
	arepo_output_plugin( config_file& cf ) : output_plugin( cf ) 
	{
		std::cerr << "\n Arepo output requires HAVE_HDF5 (otherwise use gadget2 format)!\n\n";
		exit(0);
	}
	~arepo_output_plugin() {	}
};

#else

#define GAS_PARTTYPE 0
#define HIGHRES_DM_PARTTYPE 1
#define COARSE_DM_DEFAULT_PARTTYPE 2
#define STAR_PARTTYPE 4
#define NTYPES 6

#include <sstream>
#include <string>
#include <algorithm>
#include "output.hh"
#include "HDF_IO.hh"

class arepo_output_plugin : public output_plugin
{ 
protected:
	
	// header/config
	std::vector<int> nPart;
	std::vector<int> nPartTotal;
	std::vector<double> massTable;
	double time, redshift, boxSize;
	int numFiles, doublePrec;
	
	double omega0, omega_L, hubbleParam;
	
	// configuration
	double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
	double omega_b, rhoCrit;
	double posFac, velFac;
	int coarsePartType, nPartTotAllTypes;
	bool doBaryons, useLongIDs;
	
	size_t npfine, npart, npcoarse;
	std::vector<size_t> levelcounts;
	
	// parameter file hints
	int pmgrid, gridboost;
	float softening, Tini;
	
	using output_plugin::cf_;
	
	// Nx1 vector (e.g. masses,particleids)
	template< typename T >
	void writeHDF5_a( std::string fieldName, int partTypeNum, const std::vector<T> &data )
	{
    hid_t HDF_FileID, HDF_GroupID, HDF_DatasetID, HDF_DataspaceID, HDF_Type;
    hsize_t HDF_Dims;
		
		std::stringstream GrpName;
    GrpName << "PartType" << partTypeNum;
		
    HDF_FileID = H5Fopen( fname_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
    HDF_GroupID = H5Gopen( HDF_FileID, GrpName.str().c_str() );

    HDF_Type         = GetDataType<T>();
    HDF_Dims         = data.size();
    HDF_DataspaceID  = H5Screate_simple(1, &HDF_Dims, NULL);
    HDF_DatasetID    = H5Dcreate( HDF_GroupID, fieldName.c_str(), HDF_Type, HDF_DataspaceID, H5P_DEFAULT );
		
		// write and close
    H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0] );
		
    H5Dclose( HDF_DatasetID );
    H5Sclose( HDF_DataspaceID );

    H5Gclose( HDF_GroupID );
    H5Fclose( HDF_FileID );
	}
	
	// Nx3 vector (e.g. positions,velocities), where coord = index of the second dimension (writen one at a time)
	void writeHDF5_b( std::string fieldName, int coord, int partTypeNum, std::vector<float> &data, bool readFlag = false )
	{
    hid_t HDF_FileID, HDF_GroupID, HDF_DatasetID, HDF_DataspaceID, HDF_Type;
    hsize_t HDF_Dims[2], HDF_DimsMem[2];
		
		std::stringstream GrpName;
    GrpName << "PartType" << partTypeNum;

    HDF_FileID = H5Fopen( fname_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
    HDF_GroupID = H5Gopen( HDF_FileID, GrpName.str().c_str() );

    HDF_Type         = GetDataType<float>();
    HDF_Dims[0]      = data.size();
		HDF_Dims[1]      = 3;
		
		// if dataset does not yet exist, create it (on first coord call)	
		if( !(H5Lexists(HDF_GroupID, fieldName.c_str(), H5P_DEFAULT)) )
		{
		  HDF_DataspaceID = H5Screate_simple(2, HDF_Dims, NULL);
      HDF_DatasetID = H5Dcreate( HDF_GroupID, fieldName.c_str(), HDF_Type, HDF_DataspaceID, H5P_DEFAULT );
			
			H5Sclose( HDF_DataspaceID );
			H5Dclose( HDF_DatasetID );
		}
		
		// make memory space (just indicates the size/shape of data)
		HDF_DimsMem[0] = HDF_Dims[0];
		HDF_DimsMem[1] = 1;
		hid_t HDF_MemoryspaceID = H5Screate_simple(2, HDF_DimsMem, NULL);
		
		// open hyperslab
		hsize_t count[2]={1,1}, stride[2]={1,1}, offset[2]={0,0};
		
		offset[1] = coord;       // set where in the second dimension to write
		count[0]  = HDF_Dims[0]; // set size in the first dimension (num particles of this type)
		
		HDF_DatasetID   = H5Dopen(HDF_GroupID, fieldName.c_str());
		HDF_DataspaceID = H5Dget_space(HDF_DatasetID);
		
		H5Sselect_hyperslab(HDF_DataspaceID, H5S_SELECT_SET, offset, stride, count, NULL /*, HDF_Dims*/); //HDF_DimsMem
		
		// write (or read) and close
		if( readFlag )
			H5Dread( HDF_DatasetID, HDF_Type, HDF_MemoryspaceID, HDF_DataspaceID, H5P_DEFAULT, &data[0] );
		else
			H5Dwrite( HDF_DatasetID, HDF_Type, HDF_MemoryspaceID, HDF_DataspaceID, H5P_DEFAULT, &data[0] );
		
    H5Dclose( HDF_DatasetID );
    H5Gclose( HDF_GroupID );
    H5Fclose( HDF_FileID );
	}
	
	// called from finalize()
  void generateAndWriteIDs( void )
  {
		long long offset = 0;
		nPartTotAllTypes = 0;
		
		for( size_t i=0; i < nPartTotal.size(); i++ )
		{
			if( !nPartTotal[i] )
				continue;
				
			nPartTotAllTypes += nPartTotal[i];
				
		  if( !useLongIDs ) 
			{
				std::vector<int> ids = std::vector<int>(nPartTotal[i]);
				for( int j=0; j < nPartTotal[i]; j++ )
					ids[j] = offset + j;
					
				writeHDF5_a( "ParticleIDs", i, ids );
			}
			else
			{
		    std::vector<long long> ids = std::vector<long long>(nPartTotal[i]);
				for( long long j=0; j < nPartTotal[i]; j++ )
					ids[j] = offset + j;
					
				writeHDF5_a( "ParticleIDs", i, ids );
			}
			
			// make IDs of all particle types sequential (unique) = unnecessary, but consistent with gadget output format
			offset += nPartTotal[i];
		}
	}
	
	void countLeafCells( const grid_hierarchy& gh )
	{
		npfine = 0; npart = 0; npcoarse = 0;
		
		npfine = gh.count_leaf_cells(gh.levelmax(), gh.levelmax());
		npart = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());
		
		if( levelmax_ != levelmin_ ) // multimass
			npcoarse = gh.count_leaf_cells(gh.levelmin(), gh.levelmax()-1);
	}

public:
	arepo_output_plugin( config_file& cf ) : output_plugin( cf )
	{
		// ensure that everyone knows we want to do SPH, implies: bsph=1, bbshift=1, decic_baryons=1
		// -> instead of just writing gas densities (which are here ignored), the gas displacements are also written
		cf.insertValue("setup","do_SPH","yes");
		
		// init header and config parameters
		nPart      = std::vector<int>(NTYPES,0);
		nPartTotal = std::vector<int>(NTYPES,0);
		massTable  = std::vector<double>(NTYPES,0.0);
		
		coarsePartType   = cf.getValueSafe<unsigned>("output","arepo_coarsetype",COARSE_DM_DEFAULT_PARTTYPE);
		UnitLength_in_cm = cf.getValueSafe<double>("output","arepo_unitlength",3.085678e21); // 1.0 kpc
		UnitMass_in_g    = cf.getValueSafe<double>("output","arepo_unitmass",1.989e43); // 1.0e10 solar masses
		UnitVelocity_in_cm_per_s = cf.getValueSafe<double>("output","arepo_unitvel",1e5); // 1 km/sec
		
	  omega0     = cf.getValue<double>("cosmology","Omega_m");
		omega_b    = cf.getValue<double>("cosmology","Omega_b");
		omega_L    = cf.getValue<double>("cosmology","Omega_L");
		redshift   = cf.getValue<double>("setup","zstart");
		boxSize    = cf.getValue<double>("setup","boxlength");
		doBaryons  = cf.getValueSafe<bool>("setup","baryons",false);
		useLongIDs = cf.getValueSafe<bool>("output","arepo_longids",false);
		numFiles   = cf.getValueSafe<unsigned>("output","arepo_num_files",1);
		doublePrec = cf.getValueSafe<int>("output","arepo_doubleprec",0);
		
		if( numFiles != 1 )
      throw std::runtime_error("Error: arepo_num_files>1 not yet supported.");
		if( doublePrec )
			throw std::runtime_error("Error: arepo_doubleprec not yet supported.");
		
		// factors which multiply positions and velocities
		time   = 1.0/(1.0+redshift);
		posFac = 3.085678e24 / UnitLength_in_cm; // MUSIC uses Mpc internally, i.e. posFac=1e3 for kpc output
		velFac = ( 1.0f / sqrt(time) ) * boxSize; // TODO: should be normalized by posFac?
		
		// critical density
		rhoCrit = 27.7519737e-9; // in h^2 1e10 M_sol / kpc^3
		rhoCrit *= pow(UnitLength_in_cm/3.085678e21, 3.0);
		rhoCrit *= (1.989e43/UnitMass_in_g);
		
		// calculate PMGRID suggestion
		pmgrid = pow(2,levelmin_) * 2; // unigrid
		gridboost = 1;
		
		if( levelmin_ != levelmax_ )
		{
			double lxref[3];
			double pmgrid_new;
			
			std::string temp = cf.getValue<std::string>( "setup", "ref_extent" );
			std::remove_if(temp.begin(),temp.end(),isspace);
			sscanf( temp.c_str(), "%lf,%lf,%lf", &lxref[0],&lxref[1],&lxref[2] );
			
			// fraction box length of the zoom region
			lxref[0] = pow( (lxref[0]*lxref[1]*lxref[2]),0.333 );
			
			pmgrid_new = pow(2,levelmax_) * 2; // to cover entire box at highest resolution
			pmgrid_new *= lxref[0]; // only need to cover a fraction
			
			if( (gridboost=round(pmgrid_new/pmgrid)) > 1 )
				gridboost = pow(2, ceil(log(gridboost)/log(2.0))); // round to nearest, higher power of 2
		}
		
		// calculate Tini for gas
		hubbleParam = cf.getValue<double>("cosmology","H0")/100.0;
		
		double astart = 1.0/(1.0+redshift);
		double h2     = hubbleParam*hubbleParam;
		double adec   = 1.0/( 160.0*pow(omega_b*h2/0.022,2.0/5.0) );
		double Tcmb0  = 2.726;
		
		Tini = astart<adec? Tcmb0/astart : Tcmb0/astart/astart*adec;
		
		// calculate softening suggestion
		softening = (boxSize * posFac) / pow(2,levelmax_) / 40.0;
		
		// header and sanity checks
		if ( !doBaryons )
			massTable[HIGHRES_DM_PARTTYPE] = omega0 * rhoCrit * pow(boxSize*posFac,3.0)/pow(2,3*levelmax_);
		else
			massTable[HIGHRES_DM_PARTTYPE] = (omega0-omega_b) * rhoCrit * pow(boxSize*posFac,3.0)/pow(2,3*levelmax_);
		
		if ( coarsePartType == GAS_PARTTYPE || coarsePartType == HIGHRES_DM_PARTTYPE)
      throw std::runtime_error("Error: Specified illegal Arepo particle type for coarse particles.");
		if ( coarsePartType == STAR_PARTTYPE )
			LOGWARN("WARNING: Specified coarse particle type will collide with stars if USE_SFR enabled.");
		
		// create file
		HDFCreateFile(fname_);
				
		// create particle type groups
		std::stringstream GrpName;
    GrpName << "PartType" << HIGHRES_DM_PARTTYPE;

		HDFCreateGroup(fname_, GrpName.str().c_str()); // highres or unigrid DM
		
		if( doBaryons )
		{
			GrpName.str("");
			GrpName << "PartType" << GAS_PARTTYPE;
		  HDFCreateGroup(fname_, GrpName.str().c_str()); // gas
		}
		
		if( levelmax_ != levelmin_ ) // multimass
		{
			GrpName.str("");
			GrpName << "PartType" << coarsePartType;
		  HDFCreateGroup(fname_, "PartType2"); // coarse DM
		}
	}
	
	~arepo_output_plugin()
	{	}
	
	/* ------------------------------------------------------------------------------- */
	
	void write_dm_mass( const grid_hierarchy& gh )
	{
		countLeafCells(gh);
		
		// fill levelcount for header
		levelcounts = std::vector<size_t>(levelmax_-levelmin_+1);
		for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
			levelcounts[gh.levelmax()-ilevel] = gh.count_leaf_cells(ilevel, ilevel);
		
    if( levelmax_ > levelmin_ +1 ) // morethan2bnd
		{
			// DM particles will have variable masses
			size_t count = 0;
			
			std::vector<float> data(npcoarse);
			
			for( int ilevel=gh.levelmax()-1; ilevel>=(int)gh.levelmin(); --ilevel )
			{
        // baryon particles live only on finest grid, these particles here are total matter particles
				float pmass = omega0 * rhoCrit * pow(boxSize*posFac,3.0)/pow(2,3*ilevel);	
				
				for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
					for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
						for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
							if( ! gh.is_refined(ilevel,i,j,k) )
							{
								data[count++] = pmass;
							}
			}
			
			if( count != npcoarse )
				throw std::runtime_error("Internal consistency error while writing masses");
				
			writeHDF5_a( "Masses", coarsePartType, data ); // write DM
			
		}
		else
		{
			// DM particles will all have the same mass, just write to massTable
			if( levelmax_ != levelmin_ ) // multimass
			  massTable[coarsePartType] = omega0 * rhoCrit * pow(boxSize*posFac,3.0)/pow(2,3*levelmin_);
		}		
	}
	
	void write_dm_position( int coord, const grid_hierarchy& gh )
	{
		countLeafCells(gh);
		
		// update header
		nPart[HIGHRES_DM_PARTTYPE] = npfine;
		nPart[coarsePartType]      = npcoarse;
		nPartTotal[HIGHRES_DM_PARTTYPE] = npfine;
		nPartTotal[coarsePartType]      = npcoarse;
		
		// FINE: collect displacements and convert to absolute coordinates with correct units
		int ilevel = gh.levelmax();
		
		std::vector<float> data(npfine);
		size_t count = 0;
		
		for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
			for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
				for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
					if( ! gh.is_refined(ilevel,i,j,k) )
					{
						double xx[3];
						gh.cell_pos(ilevel, i, j, k, xx);
							
						xx[coord] = (xx[coord] + (*gh.get_grid(ilevel))(i,j,k)) * boxSize;
						xx[coord] = fmod( xx[coord] + boxSize,boxSize );
						
						data[count++] = (float) (xx[coord] * posFac);
					}
						
		writeHDF5_b( "Coordinates", coord, HIGHRES_DM_PARTTYPE, data );	// write fine DM
		
		if( count != npfine )
			throw std::runtime_error("Internal consistency error while writing fine DM pos");
		
		// COARSE: collect displacements and convert to absolute coordinates with correct units
		if( levelmax_ != levelmin_ ) // multimass
		{
			data = std::vector<float> (npcoarse,0.0);
			count = 0;
			
			for( int ilevel=gh.levelmax()-1; ilevel>=(int)gh.levelmin(); --ilevel )
				for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
					for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
						for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
							if( ! gh.is_refined(ilevel,i,j,k) )
							{
								double xx[3];
								gh.cell_pos(ilevel, i, j, k, xx);
								
								xx[coord] = (xx[coord] + (*gh.get_grid(ilevel))(i,j,k)) * boxSize;
								
								if ( !doBaryons ) // if so, we will handle the mod in write_gas_position
									xx[coord] = fmod( xx[coord] + boxSize,boxSize ) * posFac;
																
								data[count++] = (float) xx[coord];
							}
				
				if( count != npcoarse )
					throw std::runtime_error("Internal consistency error while writing coarse DM pos");
					
				writeHDF5_b( "Coordinates", coord, coarsePartType, data ); // write coarse DM
		}
	}
	
	void write_dm_velocity( int coord, const grid_hierarchy& gh )
	{
		countLeafCells(gh);
			
		// FINE: collect velocities and convert to correct units
		int ilevel = gh.levelmax();
		
		std::vector<float> data(npfine);
		size_t count = 0;
		
		for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
			for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
				for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
					if( ! gh.is_refined(ilevel,i,j,k) )
					{
						data[count++] = (*gh.get_grid(ilevel))(i,j,k) * velFac;
					}
						
		writeHDF5_b( "Velocities", coord, HIGHRES_DM_PARTTYPE, data ); // write fine DM
		
		if( count != npfine )
			throw std::runtime_error("Internal consistency error while writing fine DM pos");
		
		// COARSE: collect velocities and convert to correct units
		if( levelmax_ != levelmin_ ) // multimass
		{
			data = std::vector<float> (npcoarse,0.0);
			count = 0;
			
			for( int ilevel=gh.levelmax()-1; ilevel>=(int)gh.levelmin(); --ilevel )
				for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
					for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
						for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
							if( ! gh.is_refined(ilevel,i,j,k) )
							{
								data[count++] = (*gh.get_grid(ilevel))(i,j,k) * velFac;
							}
				
				if( count != npcoarse )
					throw std::runtime_error("Internal consistency error while writing coarse DM pos");
					
				writeHDF5_b( "Velocities", coord, coarsePartType, data ); // write coarse DM
		}
	
	}
	
	void write_dm_density( const grid_hierarchy& gh )
	{ /* skip */ }
	
	void write_dm_potential( const grid_hierarchy& gh )
	{ /* skip */ }
	
	/* ------------------------------------------------------------------------------- */
	
	void write_gas_velocity( int coord, const grid_hierarchy& gh )
	{	
		countLeafCells(gh);
		
		std::vector<float> gas_data(npart); // read/write gas at all levels from the gh
		size_t count = 0;
		
		for( int ilevel=levelmax_; ilevel>=(int)levelmin_; --ilevel )
			for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
				for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
					for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
						if( ! gh.is_refined(ilevel,i,j,k) )
						{
							gas_data[count++] = (*gh.get_grid(ilevel))(i,j,k) * velFac;
						}
						
		if( count != npart )
			throw std::runtime_error("Internal consistency error while writing GAS pos");
					
		// calculate modified DM velocities if: multimass and baryons present
		if( doBaryons && npcoarse )
		{
			double facb = omega_b / omega0;
			double facc = (omega0 - omega_b) / omega0;
			
			std::vector<float> dm_data(npcoarse);
			
			writeHDF5_b( "Velocities", coord, coarsePartType, dm_data, true ); // read coarse DM vels
			
			// overwrite 
			for( size_t i=0; i < npcoarse; i++ )
				dm_data[i] = facc*dm_data[i] + facb*gas_data[npfine + i];

			writeHDF5_b( "Velocities", coord, coarsePartType, dm_data ); // overwrite coarse DM vels
		} // dm_data deallocated
		
		// restrict gas_data to fine only and request write
		std::vector<float> data( gas_data.begin() + 0, gas_data.begin() + npfine );
		
		std::vector<float>().swap( gas_data ); // deallocate
		
		writeHDF5_b( "Velocities", coord, GAS_PARTTYPE, data );	 // write highres gas
	}
	
	void write_gas_position( int coord, const grid_hierarchy& gh )
	{
		countLeafCells(gh);
		
		// update header (will actually write only gas at levelmax)
		nPart[GAS_PARTTYPE] = npfine;
		nPartTotal[GAS_PARTTYPE] = npfine;
		
		std::vector<double> gas_data(npart); // read/write gas at all levels from the gh
		size_t count = 0;
		
		double h = 1.0/(1ul<<gh.levelmax());
		
		for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
			for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
				for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
					for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
						if( ! gh.is_refined(ilevel,i,j,k) )
						{
							double xx[3];
							gh.cell_pos(ilevel, i, j, k, xx);
							
							// shift particle positions (this has to be done as the same shift
							// is used when computing the convolution kernel for SPH baryons)
							xx[coord] += 0.5*h;
														
							xx[coord] = (xx[coord] + (*gh.get_grid(ilevel))(i,j,k)) * boxSize;
											
							gas_data[count++] = xx[coord];
						}
					
		if( count != npart )
			throw std::runtime_error("Internal consistency error while writing coarse DM pos");
					
		// calculate modified DM coordinates if: multimass and baryons present
		if( doBaryons && npcoarse )
		{
			double facb = omega_b / omega0;
			double facc = (omega0 - omega_b) / omega0;
			
			std::vector<float> dm_data(npcoarse);
			
			writeHDF5_b( "Coordinates", coord, coarsePartType, dm_data, true ); // read coarse DM vels
			
			// overwrite 
			for( size_t i=0; i < npcoarse; i++ ) {
				dm_data[i] = facc*dm_data[i] + facb*gas_data[npfine + i];
				dm_data[i] = fmod( dm_data[i] + boxSize, boxSize ) * posFac;
			}

			writeHDF5_b( "Coordinates", coord, coarsePartType, dm_data ); // overwrite coarse DM vels
		}
		
		// restrict gas_data to fine only and request write
		//std::vector<float> data( gas_data.begin() + 0, gas_data.begin() + npfine );
		
		std::vector<float> data(npfine);
		
		for( size_t i = 0; i < npfine; i++ )
			data[i] = (float) ( fmod( gas_data[i] + boxSize, boxSize ) * posFac );
		
		std::vector<double>().swap( gas_data ); // deallocate
		
		writeHDF5_b( "Coordinates", coord, GAS_PARTTYPE, data ); // write highres gas

	}
	
	void write_gas_density( const grid_hierarchy& gh )
	{
		// if only saving highres gas, then all gas cells have the same initial mass
		// do not write out densities as we write out displacements
		if( doBaryons )
			massTable[GAS_PARTTYPE] = omega_b * rhoCrit * pow(boxSize*posFac,3.0)/pow(2,3*levelmax_);
	}
	
	void write_gas_potential( const grid_hierarchy& gh )
	{ /* skip */ }
	
	void finalize( void )
	{		
		// generate and add contiguous IDs for each particle type we have written
		generateAndWriteIDs();
		
		// write final header (some of these fields are required, others are extra info)
		HDFCreateGroup(fname_, "Header");
		
		std::vector<unsigned int> nPartTotalHW(nPartTotal.size());
		for( size_t i=0; i < nPartTotalHW.size(); i++ )
			nPartTotalHW[i] = (unsigned)( (size_t)nPartTotal[i] >> 32 );
			
		HDFWriteGroupAttribute(fname_, "Header", "NumPart_ThisFile",       nPart );
		HDFWriteGroupAttribute(fname_, "Header", "NumPart_Total",          nPartTotal );
		HDFWriteGroupAttribute(fname_, "Header", "NumPart_Total_HighWord", nPartTotalHW );
		HDFWriteGroupAttribute(fname_, "Header", "MassTable",              massTable );
		HDFWriteGroupAttribute(fname_, "Header", "BoxSize",                boxSize );
		HDFWriteGroupAttribute(fname_, "Header", "NumFilesPerSnapshot",    numFiles );
		HDFWriteGroupAttribute(fname_, "Header", "Time",                   time );
		HDFWriteGroupAttribute(fname_, "Header", "Redshift",               redshift );
		HDFWriteGroupAttribute(fname_, "Header", "Omega0",                 omega0 );
		HDFWriteGroupAttribute(fname_, "Header", "OmegaLambda",            omega_L );
		HDFWriteGroupAttribute(fname_, "Header", "OmegaBaryon",            omega_b );
		HDFWriteGroupAttribute(fname_, "Header", "HubbleParam",            hubbleParam );
		HDFWriteGroupAttribute(fname_, "Header", "Flag_Sfr",               0 );
		HDFWriteGroupAttribute(fname_, "Header", "Flag_Cooling",           0 );
		HDFWriteGroupAttribute(fname_, "Header", "Flag_StellarAge",        0 );
		HDFWriteGroupAttribute(fname_, "Header", "Flag_Metals",            0 );
		HDFWriteGroupAttribute(fname_, "Header", "Flag_Feedback",          0 );
		HDFWriteGroupAttribute(fname_, "Header", "Flag_DoublePrecision",   doublePrec );
		HDFWriteGroupAttribute(fname_, "Header", "Music_levelmin",         levelmin_ );
		HDFWriteGroupAttribute(fname_, "Header", "Music_levelmax",         levelmax_ );
		HDFWriteGroupAttribute(fname_, "Header", "Music_levelcounts",      levelcounts );
		HDFWriteGroupAttribute(fname_, "Header", "haveBaryons",            (int)doBaryons );
		HDFWriteGroupAttribute(fname_, "Header", "longIDs",                (int)useLongIDs );
		HDFWriteGroupAttribute(fname_, "Header", "suggested_pmgrid",       pmgrid );
		HDFWriteGroupAttribute(fname_, "Header", "suggested_gridboost",    gridboost );
		HDFWriteGroupAttribute(fname_, "Header", "suggested_highressoft",  softening );
		HDFWriteGroupAttribute(fname_, "Header", "suggested_gas_Tinit",    Tini );
		
		// output particle counts
		std::cout << " - Arepo : wrote " << nPartTotAllTypes << " particles to file..." << std::endl;
		for( size_t i=0; i < nPartTotal.size(); i++ )
			std::cout << "    type [" << i << "] : " << std::setw(12) << nPartTotal[i] << std::endl;
			
		// give config/parameter file hints		
		if( useLongIDs )
			std::cout << " - Arepo: Wrote 64bit IDs, enable LONGIDS." << std::endl;
		if( NTYPES > 6 )
			std::cout << " - Arepo: Using [" << NTYPES << "] particle types, set NTYPES to match." << std::endl;
		if( doBaryons )
			std::cout << " - Arepo: Wrote gas, set REFINEMENT_HIGH_RES_GAS and GENERATE_GAS_IN_ICS with "
			          << "SPLIT_PARTICLE_TYPE=" << pow(2,coarsePartType) << "." << std::endl;
		if( levelmax_ != levelmin_ )
			std::cout << " - Arepo: Have zoom type ICs, set PLACEHIGHRESREGION=" << pow(2,HIGHRES_DM_PARTTYPE)
                << " (suggest PMGRID=" << pmgrid << " with GRIDBOOST=" << gridboost << ")." << std::endl;
		if( levelmax_ > levelmin_ + 1 )
			std::cout << " - Arepo: More than one coarse DM mass using same type, set INDIVIDUAL_GRAVITY_SOFTENING=" 
			          << pow(2,coarsePartType) << " (+" << pow(2,STAR_PARTTYPE) << " if including stars)." << std::endl;
		if( doBaryons )
			std::cout << " - Arepo: Set initial gas temperature to " << std::fixed << std::setprecision(3) << Tini << " K." << std::endl;
		std::cout << " - Arepo: Suggest grav softening = " << std::setprecision(3) << softening << " for high res DM." << std::endl;
			
	}
	
};

namespace{
	output_plugin_creator_concrete< arepo_output_plugin > creator("arepo");
}

#endif // HAVE_HDF5
