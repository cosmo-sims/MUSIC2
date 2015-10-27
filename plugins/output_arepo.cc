/*
 * output_arepo.cc - This file is part of MUSIC -
 * a code to generate multi-scale initial conditions 
 * for cosmological simulations
 * 
 * Copyright (C) 2010  Oliver Hahn
 * 
 * Plugin: Dylan Nelson (dnelson@cfa.harvard.edu)
 */
 
#ifdef HAVE_HDF5

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
  std::vector< std::vector<unsigned int> > nPart;
  std::vector<long long> nPartTotal;
  std::vector<double> massTable;
  double time, redshift, boxSize;
  unsigned int numFiles, coarsePartType;
  
  double omega0, omega_L, hubbleParam;
  
  // configuration
  double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
  double omega_b, rhoCrit;
  double posFac, velFac;
  long long nPartTotAllTypes;
  bool doBaryons, useLongIDs, doublePrec;
  
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
    hsize_t HDF_Dims, offset = 0;
    
    std::stringstream GrpName;
    GrpName << "PartType" << partTypeNum;
    
    for( unsigned i=0; i < numFiles; i++ )
    {
      std::string filename = fname_;
      HDF_Dims = data.size();
      
      // modify local filename and write size
      if( numFiles > 1 )
      {
        std::stringstream s;
        s << "." << i << ".hdf5";
        filename.replace(filename.find(".hdf5"), 5, s.str());
        
        HDF_Dims = ceil( data.size() / numFiles );
        if( i == numFiles-1 )
          HDF_Dims = data.size() - offset;
      }
    
      HDF_FileID = H5Fopen( filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
      HDF_GroupID = H5Gopen( HDF_FileID, GrpName.str().c_str() );

      HDF_Type         = GetDataType<T>();
      HDF_DataspaceID  = H5Screate_simple(1, &HDF_Dims, NULL);
      HDF_DatasetID    = H5Dcreate(HDF_GroupID, fieldName.c_str(), HDF_Type, HDF_DataspaceID, H5P_DEFAULT);
      
      // write and close
      H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[offset] );
      
      H5Dclose( HDF_DatasetID );
      H5Sclose( HDF_DataspaceID );

      H5Gclose( HDF_GroupID );
      H5Fclose( HDF_FileID );
      
      offset += HDF_Dims;
    }
  }
  
  // Nx3 vector (e.g. pos,vel), where coord = index of the second dimension (written one at a time)
  template<typename T>
  void writeHDF5_b( std::string fieldName, int coord, int partTypeNum, std::vector<T> &data, bool readFlag = false )
  {
    hid_t HDF_FileID, HDF_GroupID, HDF_DatasetID, HDF_DataspaceID, HDF_Type;
    hsize_t HDF_Dims[2], HDF_DimsMem[2], w_offset = 0;
    
    std::stringstream GrpName;
    GrpName << "PartType" << partTypeNum;

    for( unsigned i=0; i < numFiles; i++ )
    {
      std::string filename = fname_;
      HDF_Dims[0] = data.size();
      
      // modify local filename and write size
      if( numFiles > 1 )
      {
        std::stringstream s;
        s << "." << i << ".hdf5";
        filename.replace(filename.find(".hdf5"), 5, s.str());
        
        HDF_Dims[0] = ceil( data.size() / numFiles );
        if( i == numFiles-1 )
          HDF_Dims[0] = data.size() - w_offset;
      }
    
      HDF_FileID  = H5Fopen( filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
      HDF_GroupID = H5Gopen( HDF_FileID, GrpName.str().c_str() );

      HDF_Type    = GetDataType<T>();
      HDF_Dims[1] = 3;
      
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
      
      H5Sselect_hyperslab(HDF_DataspaceID, H5S_SELECT_SET, offset, stride, count, NULL);
      
      // write (or read) and close
      if( readFlag )
        H5Dread( HDF_DatasetID, HDF_Type, HDF_MemoryspaceID, HDF_DataspaceID, H5P_DEFAULT, 
                 &data[w_offset] );
      else
        H5Dwrite( HDF_DatasetID, HDF_Type, HDF_MemoryspaceID, HDF_DataspaceID, H5P_DEFAULT, 
                  &data[w_offset] );
      
      H5Dclose( HDF_DatasetID );
      H5Gclose( HDF_GroupID );
      H5Fclose( HDF_FileID );
      
      w_offset += HDF_Dims[0];
    }
  }
  
  // called from finalize()
  void generateAndWriteIDs( void )
  {
    long long offset = 1; // don't use ID==0
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

  template< typename T >
  void __write_dm_mass( const grid_hierarchy& gh )
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
      
      std::vector<T> data(npcoarse);
      
      for( int ilevel=gh.levelmax()-1; ilevel>=(int)gh.levelmin(); --ilevel )
      {
        // baryon particles live only on finest grid, these particles here are total matter particles
        T pmass = omega0 * rhoCrit * pow(boxSize*posFac,3.0)/pow(2,3*ilevel); 
        
        for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
          for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
            for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
              if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
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
  
  template< typename T >
  void __write_dm_position( int coord, const grid_hierarchy& gh )
  {
    countLeafCells(gh);
    
    // update header
    hsize_t offset_fine = 0, offset_coarse = 0;
    
    for( unsigned i=0; i < numFiles; i++ ) {
      hsize_t dims_fine   = ceil( npfine / numFiles );
      hsize_t dims_coarse = ceil( npcoarse / numFiles );      
      
      if( i == numFiles-1 ) {
        dims_fine   = npfine - offset_fine;
        dims_coarse = npcoarse - offset_coarse;
      }
        
      nPart[i][HIGHRES_DM_PARTTYPE] = dims_fine;
      nPart[i][coarsePartType]      = dims_coarse;
      
      offset_fine   += dims_fine;
      offset_coarse += dims_coarse;
    }
    
    nPartTotal[HIGHRES_DM_PARTTYPE] = npfine;
    nPartTotal[coarsePartType]      = npcoarse;
    
    // FINE: collect displacements and convert to absolute coordinates with correct units
    int ilevel = gh.levelmax();
    
    std::vector<T> data(npfine);
    size_t count = 0;
    
    for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
      for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
        for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
          if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
          {
            double xx[3];
            gh.cell_pos(ilevel, i, j, k, xx);
              
            xx[coord] = (xx[coord] + (*gh.get_grid(ilevel))(i,j,k)) * boxSize;
            xx[coord] = fmod( xx[coord] + boxSize,boxSize );
            
            data[count++] = (T) (xx[coord] * posFac);
          }
            
    writeHDF5_b( "Coordinates", coord, HIGHRES_DM_PARTTYPE, data ); // write fine DM
    
    if( count != npfine )
      throw std::runtime_error("Internal consistency error while writing fine DM pos");
    
    // COARSE: collect displacements and convert to absolute coordinates with correct units
    if( levelmax_ != levelmin_ ) // multimass
    {
      data = std::vector<T> (npcoarse,0.0);
      count = 0;
      
      for( int ilevel=gh.levelmax()-1; ilevel>=(int)gh.levelmin(); --ilevel )
        for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
          for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
            for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
              if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
              {
                double xx[3];
                gh.cell_pos(ilevel, i, j, k, xx);
                
                xx[coord] = (xx[coord] + (*gh.get_grid(ilevel))(i,j,k)) * boxSize;
                
                if ( !doBaryons ) // if so, we will handle the mod in write_gas_position
                  xx[coord] = fmod( xx[coord] + boxSize,boxSize ) * posFac;
                                
                data[count++] = (T) xx[coord];
              }
        
        if( count != npcoarse )
          throw std::runtime_error("Internal consistency error while writing coarse DM pos");
          
        writeHDF5_b( "Coordinates", coord, coarsePartType, data ); // write coarse DM
    }
  }
  
  template< typename T >
  void __write_dm_velocity( int coord, const grid_hierarchy& gh )
  {
    countLeafCells(gh);
      
    // FINE: collect velocities and convert to correct units
    int ilevel = gh.levelmax();
    
    std::vector<T> data(npfine);
    size_t count = 0;
    
    for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
      for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
        for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
          if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
          {
            data[count++] = (T) (*gh.get_grid(ilevel))(i,j,k) * velFac;
          }
            
    writeHDF5_b( "Velocities", coord, HIGHRES_DM_PARTTYPE, data ); // write fine DM
    
    if( count != npfine )
      throw std::runtime_error("Internal consistency error while writing fine DM pos");
    
    // COARSE: collect velocities and convert to correct units
    if( levelmax_ != levelmin_ ) // multimass
    {
      data = std::vector<T> (npcoarse,0.0);
      count = 0;
      
      for( int ilevel=gh.levelmax()-1; ilevel>=(int)gh.levelmin(); --ilevel )
        for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
          for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
            for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
              if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
              {
                data[count++] = (T) (*gh.get_grid(ilevel))(i,j,k) * velFac;
              }
        
        if( count != npcoarse )
          throw std::runtime_error("Internal consistency error while writing coarse DM pos");
          
        writeHDF5_b( "Velocities", coord, coarsePartType, data ); // write coarse DM
    }
  
  }

  template< typename T >
  void __write_gas_velocity( int coord, const grid_hierarchy& gh )
  { 
    countLeafCells(gh);
    
    std::vector<T> gas_data(npart); // read/write gas at all levels from the gh
    size_t count = 0;
    
    for( int ilevel=levelmax_; ilevel>=(int)levelmin_; --ilevel )
      for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
        for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
          for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
            if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
            {
              gas_data[count++] = (T) (*gh.get_grid(ilevel))(i,j,k) * velFac;
            }
            
    if( count != npart )
      throw std::runtime_error("Internal consistency error while writing GAS pos");
          
    // calculate modified DM velocities if: multimass and baryons present
    if( doBaryons && npcoarse )
    {
      double facb = omega_b / omega0;
      double facc = (omega0 - omega_b) / omega0;
      
      std::vector<T> dm_data(npcoarse);
      
      writeHDF5_b( "Velocities", coord, coarsePartType, dm_data, true ); // read coarse DM vels
      
      // overwrite 
      for( size_t i=0; i < npcoarse; i++ )
        dm_data[i] = facc*dm_data[i] + facb*gas_data[npfine + i];

      writeHDF5_b( "Velocities", coord, coarsePartType, dm_data ); // overwrite coarse DM vels
    } // dm_data deallocated
    
    // restrict gas_data to fine only and request write
    std::vector<T> data( gas_data.begin() + 0, gas_data.begin() + npfine );
    
    std::vector<T>().swap( gas_data ); // deallocate
    
    writeHDF5_b( "Velocities", coord, GAS_PARTTYPE, data );  // write highres gas
  }
  
  template< typename T >
  void __write_gas_position( int coord, const grid_hierarchy& gh )
  {
    countLeafCells(gh);
    
    // update header (will actually write only gas at levelmax)
    hsize_t offset = 0;
    
    for( unsigned i=0; i < numFiles; i++ ) {
      hsize_t dims = ceil( npfine / numFiles );        
      if( i == numFiles-1 )
        dims = npfine - offset;
        
      nPart[i][GAS_PARTTYPE] = dims;
      offset += dims;
    }
    
    nPartTotal[GAS_PARTTYPE] = npfine;
    
    std::vector<double> gas_data(npart); // read/write gas at all levels from the gh
    size_t count = 0;
    
    double h = 1.0/(1ul<<gh.levelmax());
    
    for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
      for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
        for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
          for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
            if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
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
      
      std::vector<T> dm_data(npcoarse);
      
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
    
    std::vector<T> data(npfine);
    
    for( size_t i = 0; i < npfine; i++ )
      data[i] = (T) ( fmod( gas_data[i] + boxSize, boxSize ) * posFac );
    
    std::vector<double>().swap( gas_data ); // deallocate
    
    writeHDF5_b( "Coordinates", coord, GAS_PARTTYPE, data ); // write highres gas

  }

public:
  arepo_output_plugin( config_file& cf ) : output_plugin( cf )
  {
    // ensure that everyone knows we want to do SPH, implies: bsph=1, bbshift=1, decic_baryons=1
    // -> instead of just writing gas densities (which are here ignored), the gas displacements are also written
    cf.insertValue("setup","do_SPH","yes");
    
    // init header and config parameters
    nPartTotal = std::vector<long long>(NTYPES,0);
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
    doublePrec = cf.getValueSafe<bool>("output","arepo_doubleprec",0);
    
    for( unsigned i=0; i < numFiles; i++ )
      nPart.push_back( std::vector<unsigned int>(NTYPES,0) );
    
    // factors which multiply positions and velocities
    time   = 1.0/(1.0+redshift);
    posFac = 3.085678e24 / UnitLength_in_cm; // MUSIC uses Mpc internally, i.e. posFac=1e3 for kpc output
    velFac = ( 1.0f / sqrt(time) ) * boxSize;
    
    // critical density
    rhoCrit = 27.7519737e-9; // in h^2 1e10 M_sol / kpc^3
    rhoCrit *= pow(UnitLength_in_cm/3.085678e21, 3.0);
    rhoCrit *= (1.989e43/UnitMass_in_g);
    
    // calculate PMGRID suggestion
    pmgrid = pow(2,levelmin_) * 2; // unigrid
    gridboost = 1;
    
    if( levelmin_ != levelmax_ )
    {
      double lxref[3], x0ref[3], x1ref[3];
      double pmgrid_new;
      
      the_region_generator->get_AABB(x0ref,x1ref,levelmax_); // generalized beyond box
      for (int i=0; i < 3; i++)
        lxref[i] = x1ref[i] - x0ref[i];
      
      // fraction box length of the zoom region
      lxref[0] = pow( (lxref[0]*lxref[1]*lxref[2]),0.333 );
      
      pmgrid_new = pow(2,levelmax_) * 2; // to cover entire box at highest resolution
      pmgrid_new *= lxref[0]; // only need to cover a fraction
      
      if( (gridboost=round(pmgrid_new/pmgrid)) > 1 )
        gridboost = pow(2, ceil(log(gridboost)/log(2.0))); // round to nearest, higher power of 2
			if( gridboost == 0 )
				gridboost = 1;
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
    
    // create file(s)
    for( unsigned i=0; i < numFiles; i++ )
    {
      std::string filename = fname_;
      if( numFiles > 1 )
      {
        size_t pos = filename.find(".hdf5");
        if( pos != filename.length()-5 )
          throw std::runtime_error("Error: Unexpected output filename (doesn't end in .hdf5).");
          
        std::stringstream s;
        s << "." << i << ".hdf5";
        filename.replace(pos, 5, s.str());
      }
      
      HDFCreateFile(filename);
          
      // create particle type groups
      std::stringstream GrpName;
      GrpName << "PartType" << HIGHRES_DM_PARTTYPE;

      HDFCreateGroup(filename, GrpName.str().c_str()); // highres or unigrid DM
      
      if( doBaryons )
      {
        GrpName.str("");
        GrpName << "PartType" << GAS_PARTTYPE;
        HDFCreateGroup(filename, GrpName.str().c_str()); // gas
      }
      
      if( levelmax_ != levelmin_ ) // multimass
      {
        GrpName.str("");
        GrpName << "PartType" << coarsePartType;
        HDFCreateGroup(filename, GrpName.str().c_str()); // coarse DM
      }
    }
  }
  
  ~arepo_output_plugin()
  { }
  
  /* ------------------------------------------------------------------------------- */
  void write_dm_mass( const grid_hierarchy& gh )
  {
    if(!doublePrec)
      __write_dm_mass<float>(gh);
    else
      __write_dm_mass<double>(gh);
  }
  
  void write_dm_position( int coord, const grid_hierarchy& gh )
  {
    if(!doublePrec)
      __write_dm_position<float>(coord, gh);
    else
      __write_dm_position<double>(coord, gh);
  }
  
  void write_dm_velocity( int coord, const grid_hierarchy& gh )
  {
    if(!doublePrec)
      __write_dm_velocity<float>(coord, gh);
    else
      __write_dm_velocity<double>(coord, gh);
  }
  
  void write_dm_density( const grid_hierarchy& gh )
  { /* skip */ }
  
  void write_dm_potential( const grid_hierarchy& gh )
  { /* skip */ }
  
  /* ------------------------------------------------------------------------------- */
  void write_gas_velocity( int coord, const grid_hierarchy& gh )
  {
    if(!doublePrec)
      __write_gas_velocity<float>(coord, gh);
    else
      __write_gas_velocity<double>(coord, gh);
  }
  
  void write_gas_position( int coord, const grid_hierarchy& gh )
  {
    if(!doublePrec)
      __write_gas_position<float>(coord, gh);
    else
      __write_gas_position<double>(coord, gh);
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
    
    std::vector<unsigned int> nPartTotalLW(nPartTotal.size());
    std::vector<unsigned int> nPartTotalHW(nPartTotal.size());
    for( size_t i=0; i < nPartTotalHW.size(); i++ ) {
      nPartTotalHW[i] = (unsigned)( (size_t)nPartTotal[i] >> 32 );
      nPartTotalLW[i] = (unsigned int)( (size_t)nPartTotal[i] );
    }
    
    // output particle counts
    std::cout << " - Arepo : wrote " << nPartTotAllTypes << " particles..." << std::endl;
    for( size_t i=0; i < nPartTotal.size(); i++ )
      std::cout << "    type [" << i << "] : " << std::setw(12) << nPartTotal[i] << std::endl;
    std::cout << std::endl;

    // write final header (some of these fields are required, others are extra info)
    for( unsigned i=0; i < numFiles; i++ )
    {
      std::string filename = fname_;
      if( numFiles > 1 )
      {
        std::stringstream s;
        s << "." << i << ".hdf5";
        filename.replace(filename.find(".hdf5"), 5, s.str());
        
        std::cout << "    " << filename;
        for( size_t j=0; j < nPart[i].size(); j++ )
          std::cout << " " << std::setw(10) << nPart[i][j];
        std::cout << std::endl;
      }
      
      HDFCreateGroup(filename, "Header");
      
      HDFWriteGroupAttribute(filename, "Header", "NumPart_ThisFile",       nPart[i] );
      HDFWriteGroupAttribute(filename, "Header", "NumPart_Total",          nPartTotalLW );
      HDFWriteGroupAttribute(filename, "Header", "NumPart_Total_HighWord", nPartTotalHW );
      HDFWriteGroupAttribute(filename, "Header", "MassTable",              massTable );
      HDFWriteGroupAttribute(filename, "Header", "BoxSize",                boxSize );
      HDFWriteGroupAttribute(filename, "Header", "NumFilesPerSnapshot",    numFiles );
      HDFWriteGroupAttribute(filename, "Header", "Time",                   time );
      HDFWriteGroupAttribute(filename, "Header", "Redshift",               redshift );
      HDFWriteGroupAttribute(filename, "Header", "Omega0",                 omega0 );
      HDFWriteGroupAttribute(filename, "Header", "OmegaLambda",            omega_L );
      HDFWriteGroupAttribute(filename, "Header", "OmegaBaryon",            omega_b );
      HDFWriteGroupAttribute(filename, "Header", "HubbleParam",            hubbleParam );
      HDFWriteGroupAttribute(filename, "Header", "Flag_Sfr",               0 );
      HDFWriteGroupAttribute(filename, "Header", "Flag_Cooling",           0 );
      HDFWriteGroupAttribute(filename, "Header", "Flag_StellarAge",        0 );
      HDFWriteGroupAttribute(filename, "Header", "Flag_Metals",            0 );
      HDFWriteGroupAttribute(filename, "Header", "Flag_Feedback",          0 );
      HDFWriteGroupAttribute(filename, "Header", "Flag_DoublePrecision",   (int)doublePrec );
      HDFWriteGroupAttribute(filename, "Header", "Music_levelmin",         levelmin_ );
      HDFWriteGroupAttribute(filename, "Header", "Music_levelmax",         levelmax_ );
      HDFWriteGroupAttribute(filename, "Header", "Music_levelcounts",      levelcounts );
      HDFWriteGroupAttribute(filename, "Header", "haveBaryons",            (int)doBaryons );
      HDFWriteGroupAttribute(filename, "Header", "longIDs",                (int)useLongIDs );
      HDFWriteGroupAttribute(filename, "Header", "suggested_pmgrid",       pmgrid );
      HDFWriteGroupAttribute(filename, "Header", "suggested_gridboost",    gridboost );
      HDFWriteGroupAttribute(filename, "Header", "suggested_highressoft",  softening );
      HDFWriteGroupAttribute(filename, "Header", "suggested_gas_Tinit",    Tini );
    }
      
    // give config/parameter file hints   
    if( useLongIDs )
      std::cout << " - Arepo: Wrote 64bit IDs, enable LONGIDS." << std::endl;
		if( doublePrec )
      std::cout << " - Arepo: Double precision ICs, set INPUT_IN_DOUBLEPRECISION." << std::endl;
    if( NTYPES != 6 )
      std::cout << " - Arepo: Using [" << NTYPES << "] particle types, set NTYPES to match." << std::endl;
    if( doBaryons )
      std::cout << " - Arepo: Wrote high-res gas (only), set REFINEMENT_HIGH_RES_GAS and GENERATE_GAS_IN_ICS with "
                << "SPLIT_PARTICLE_TYPE=" << pow(2,coarsePartType) << "." << std::endl;
    if( levelmax_ != levelmin_ )
      std::cout << " - Arepo: Have zoom type ICs, set PLACEHIGHRESREGION=" << pow(2,HIGHRES_DM_PARTTYPE)
                << " (suggest PMGRID=" << pmgrid << " with GRIDBOOST=" << gridboost << ")." << std::endl;
    else
      std::cout << " - Arepo: Have unigrid type ICs (suggest PMGRID=" << pmgrid << ")." << std::endl;
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
