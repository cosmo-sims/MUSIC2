/*
 
 output_gadget2.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#include <fstream>
#include <map>
#include "log.hh"
#include "region_generator.hh"
#include "output.hh"
#include "mg_interp.hh"
#include "mesh.hh"

const int empty_fill_bytes = 60;

template< typename T_store=float >
class gadget2_output_plugin : public output_plugin
{
public:
  bool do_baryons_;
  double omegab_;
  double gamma_;
  bool shift_halfcell_;
    
protected:
  
  std::ofstream ofs_;
  bool blongids_;
  bool bhave_particlenumbers_;
    
  std::map<std::string,double> units_length_;
  std::map<std::string,double> units_mass_;
  std::map<std::string,double> units_vel_;

  double unit_length_chosen_;
  double unit_mass_chosen_;
  double unit_vel_chosen_;
	

  typedef struct io_header
  {
    int npart[6];                        
    double mass[6];                      
    double time;                         
    double redshift;                     
    int flag_sfr;                        
    int flag_feedback;                   
    unsigned int npartTotal[6];          
    int flag_cooling;                    
    int num_files;                       
    double BoxSize;                      
    double Omega0;                       
    double OmegaLambda;                  
    double HubbleParam;                  
    int flag_stellarage;                 
    int flag_metals;                     
    unsigned int npartTotalHighWord[6];  
    int  flag_entropy_instead_u;         
    char fill[empty_fill_bytes];                       
  }header;                       
  
  
  header header_;
  
  std::string fname;
  
  enum iofields {
    id_dm_mass, id_dm_vel, id_dm_pos, id_gas_vel, id_gas_rho, id_gas_temp, id_gas_pos
  };
  
  size_t np_per_type_[6];
  
  size_t block_buf_size_;
  size_t npartmax_;
  unsigned nfiles_;
  
  unsigned bndparticletype_;
  bool bmorethan2bnd_;
  bool kpcunits_;
  bool msolunits_;
  double YHe_;
  bool spread_coarse_acrosstypes_;
  
  refinement_mask refmask;
  
  void distribute_particles( unsigned nfiles, std::vector< std::vector<unsigned> >& np_per_file, std::vector<unsigned>& np_tot_per_file )
  {
    np_per_file.assign( nfiles, std::vector<unsigned>( 6, 0 ) );
    np_tot_per_file.assign( nfiles, 0 );
    
    size_t n2dist[6];
    size_t ntotal = 0;
    for( int i=0; i<6; ++i )
      {
	ntotal += np_per_type_[i];
	n2dist[i] = np_per_type_[i];
      }
    
    size_t nnominal = (size_t)((double)ntotal/(double)nfiles);
    size_t nlast = ntotal - nnominal * (nfiles-1);
    
    for( unsigned i=0; i<nfiles; ++i )
      {
	size_t nthisfile = 0;
        
	size_t nmax = (i==nfiles-1)?nlast:nnominal;
        
	for( int itype=0; itype<6; ++itype )
	  {
	    if( n2dist[itype]==0 ) continue;
	    np_per_file[i][itype] = std::min( n2dist[itype], nmax-nthisfile );
	    n2dist[itype] -= np_per_file[i][itype];
	    nthisfile += np_per_file[i][itype];
            
	    if( nthisfile >= nmax ) break;
	  }
	
	np_tot_per_file[i] = nthisfile;
      }
    
    for( int i=0; i<6; ++i )
      assert(n2dist[i]==0);
  }
	

  std::ifstream& open_and_check( std::string ffname, size_t npart, size_t offset=0 )
  {
    std::ifstream ifs( ffname.c_str(), std::ios::binary );
    size_t blk;
    ifs.read( (char*)&blk, sizeof(size_t) );
    if( blk != npart*(size_t)sizeof(T_store) )
      {	
	LOGERR("Internal consistency error in gadget2 output plug-in");
	LOGERR("Expected %ld bytes in temp file but found %ld",npart*(size_t)sizeof(T_store),blk);
	throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
      }
    ifs.seekg( offset, std::ios::cur );
    
    return ifs;
  }
  
  class pistream : public std::ifstream
  {
  public:
    pistream (std::string fname, size_t npart, size_t offset=0 )
      : std::ifstream( fname.c_str(), std::ios::binary )
    {
      size_t blk;
      
      if( !this->good() )
	{	
	  LOGERR("Could not open buffer file in gadget2 output plug-in");
	  throw std::runtime_error("Could not open buffer file in gadget2 output plug-in");
	}
      
      this->read( (char*)&blk, sizeof(size_t) );
      
      if( blk != npart*sizeof(T_store) )
	{	
	  LOGERR("Internal consistency error in gadget2 output plug-in");
	  LOGERR("Expected %ld bytes in temp file but found %ld",npart*sizeof(T_store),blk);
	  throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
	}
      
      this->seekg( offset+sizeof(size_t), std::ios::beg );
    }
    
    pistream ()
    { }
    
    void open(std::string fname, size_t npart, size_t offset=0 )
    {
      std::ifstream::open( fname.c_str(), std::ios::binary );
      size_t blk;
      
      if( !this->good() )
	{	
	  LOGERR("Could not open buffer file \'%s\' in gadget2 output plug-in",fname.c_str());
	  throw std::runtime_error("Could not open buffer file in gadget2 output plug-in");
	}
      
      this->read( (char*)&blk, sizeof(size_t) );
      
      if( blk != npart*sizeof(T_store) )
	{	
	  LOGERR("Internal consistency error in gadget2 output plug-in");
	  LOGERR("Expected %ld bytes in temp file but found %ld",npart*sizeof(T_store),blk);
	  throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
	}
      
      this->seekg( offset+sizeof(size_t), std::ios::beg );
    }
  };
  
  class postream : public std::fstream
  {
  public:
    postream (std::string fname, size_t npart, size_t offset=0 )
      : std::fstream( fname.c_str(), std::ios::binary|std::ios::in|std::ios::out )
    {
      size_t blk;
      
      if( !this->good() )
	{	
	  LOGERR("Could not open buffer file in gadget2 output plug-in");
	  throw std::runtime_error("Could not open buffer file in gadget2 output plug-in");
	}
      
      this->read( (char*)&blk, sizeof(size_t) );
      
      if( blk != npart*sizeof(T_store) )
	{	
	  LOGERR("Internal consistency error in gadget2 output plug-in");
	  LOGERR("Expected %ld bytes in temp file but found %ld",npart*sizeof(T_store),blk);
	  throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
	}
      
      this->seekg( offset, std::ios::cur );
      this->seekp( offset+sizeof(size_t), std::ios::beg );
    }
    
    postream ()
    { }
    
    void open(std::string fname, size_t npart, size_t offset=0 )
    {
      if( is_open() )
	this->close();
      
      std::fstream::open( fname.c_str(), std::ios::binary|std::ios::in|std::ios::out );
      size_t blk;
      
      if( !this->good() )
	{	
	  LOGERR("Could not open buffer file \'%s\' in gadget2 output plug-in",fname.c_str());
	  throw std::runtime_error("Could not open buffer file in gadget2 output plug-in");
	}
      
      this->read( (char*)&blk, sizeof(size_t) );
      
      if( blk != npart*sizeof(T_store) )
	{	
	  LOGERR("Internal consistency error in gadget2 output plug-in");
	  LOGERR("Expected %ld bytes in temp file but found %ld",npart*sizeof(T_store),blk);
	  throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
	}
      
      this->seekg( offset, std::ios::cur );
      this->seekp( offset+sizeof(size_t), std::ios::beg );
    }
  };
  
  void combine_components_for_coarse( void )
  {        
    const size_t
      nptot = np_per_type_[1]+np_per_type_[2]+np_per_type_[3]+np_per_type_[4]+np_per_type_[5],
      npfine = np_per_type_[1],
      npcoarse = nptot-npfine;
    
    std::vector<T_store> tmp1, tmp2;
    
    tmp1.assign(block_buf_size_,0.0);
    tmp2.assign(block_buf_size_,0.0);
    
    double facb = omegab_/header_.Omega0, facc = (header_.Omega0-omegab_)/header_.Omega0;
    
    
    for( int icomp=0; icomp < 3; ++icomp )
      {
	char fc[256], fb[256];
	postream iffs1, iffs2;
        
	/*** positions ***/
	
	sprintf( fc, "___ic_temp_%05d.bin", 100*id_dm_pos+icomp );
	sprintf( fb, "___ic_temp_%05d.bin", 100*id_gas_pos+icomp );
        
	iffs1.open( fc, nptot, npfine*sizeof(T_store) );
	iffs2.open( fb, nptot, npfine*sizeof(T_store) );
        
	size_t npleft = npcoarse;
	size_t n2read = std::min((size_t)block_buf_size_,npleft);
	while( n2read > 0ul )
	  {
	    std::streampos sp = iffs1.tellg();
	    iffs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
	    iffs2.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
            
	    for( size_t i=0; i<n2read; ++i )
	      {
		tmp1[i] = facc*tmp1[i] + facb*tmp2[i];
	      }
	    
	    iffs1.seekp( sp );
	    iffs1.write( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
            
	    npleft -= n2read;
	    n2read = std::min( (size_t)block_buf_size_,npleft );
	  }
	
	iffs1.close();
	iffs2.close();
        
	/*** velocities ***/
        
	sprintf( fc, "___ic_temp_%05d.bin", 100*id_dm_vel+icomp );
	sprintf( fb, "___ic_temp_%05d.bin", 100*id_gas_vel+icomp );
        
	iffs1.open( fc, nptot, npfine*sizeof(T_store) );
	iffs2.open( fb, nptot, npfine*sizeof(T_store) );
        
	npleft = npcoarse;
	n2read = std::min( (size_t)block_buf_size_,npleft);
	
	while( n2read > 0ul )
	  {
	    std::streampos sp = iffs1.tellg();
	    iffs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
	    iffs2.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
            
	    for( size_t i=0; i<n2read; ++i )
	      {
		tmp1[i] = facc*tmp1[i] + facb*tmp2[i];
	      }
	    
	    iffs1.seekp( sp );
	    iffs1.write( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
            
	    npleft -= n2read;
	    n2read = std::min( (size_t)block_buf_size_,npleft );
	  }
	
	iffs1.close();
	iffs2.close();
      }
    
  }
  
  void assemble_gadget_file( void )
  {
    
    if( do_baryons_ )
      combine_components_for_coarse();
    
    
    
    //............................................................................
    //... copy from the temporary files, interleave the data and save ............
		
    char fnx[256],fny[256],fnz[256],fnvx[256],fnvy[256],fnvz[256],fnm[256];
    char fnbx[256], fnby[256], fnbz[256], fnbvx[256], fnbvy[256], fnbvz[256];
    
    sprintf( fnx,  "___ic_temp_%05d.bin", 100*id_dm_pos+0 );
    sprintf( fny,  "___ic_temp_%05d.bin", 100*id_dm_pos+1 );
    sprintf( fnz,  "___ic_temp_%05d.bin", 100*id_dm_pos+2 );
    sprintf( fnvx, "___ic_temp_%05d.bin", 100*id_dm_vel+0 );
    sprintf( fnvy, "___ic_temp_%05d.bin", 100*id_dm_vel+1 );
    sprintf( fnvz, "___ic_temp_%05d.bin", 100*id_dm_vel+2 );
    sprintf( fnm,  "___ic_temp_%05d.bin", 100*id_dm_mass  );
    
    sprintf( fnbx,  "___ic_temp_%05d.bin", 100*id_gas_pos+0 );
    sprintf( fnby,  "___ic_temp_%05d.bin", 100*id_gas_pos+1 );
    sprintf( fnbz,  "___ic_temp_%05d.bin", 100*id_gas_pos+2 );
    sprintf( fnbvx, "___ic_temp_%05d.bin", 100*id_gas_vel+0 );
    sprintf( fnbvy, "___ic_temp_%05d.bin", 100*id_gas_vel+1 );
    sprintf( fnbvz, "___ic_temp_%05d.bin", 100*id_gas_vel+2 );
    
    
    pistream iffs1, iffs2, iffs3;
    
    const size_t 
      nptot = np_per_type_[0]+np_per_type_[1]+np_per_type_[2]+np_per_type_[3]+np_per_type_[4]+np_per_type_[5],
      //npgas = np_fine_gas_,
      npcdm = nptot-np_per_type_[0];
    
    size_t
      wrote_coarse = 0,
      wrote_gas  = 0,
      wrote_dm   = 0;
    
    size_t
      npleft = nptot, 
      n2read = std::min((size_t)block_buf_size_,npleft);
    
    std::cout << " - Gadget2 : writing " << nptot << " particles to file...\n";
    for( int i=0; i<6; ++i )
      if( np_per_type_[i] > 0 )
	LOGINFO("      type   %d : %12llu [m=%g]", i, np_per_type_[i], header_.mass[i] );
    
    bool bbaryons = np_per_type_[0] > 0;
    
    std::vector<T_store> adata3;
    adata3.reserve( 3*block_buf_size_ );
    T_store *tmp1, *tmp2, *tmp3;
    
    tmp1 = new T_store[block_buf_size_];
    tmp2 = new T_store[block_buf_size_];
    tmp3 = new T_store[block_buf_size_];
    
    //... for multi-file output
    //int fileno = 0;
    //size_t npart_left = nptot;
    
    //std::vector<unsigned> nfdm_per_file, nfgas_per_file, nc_per_file;
        
    std::vector< std::vector<unsigned> > np_per_file;
    std::vector<unsigned> np_tot_per_file;
    
    distribute_particles( nfiles_, np_per_file, np_tot_per_file );
    
    if( nfiles_ > 1 )
      {
	LOGINFO("Gadget2 : distributing particles to %d files", nfiles_ );
	//<< "                 " << std::setw(12) << "type 0" << "," << std::setw(12) << "type 1" << "," << std::setw(12) << "type " << bndparticletype_ << std::endl;
	for( unsigned i=0; i<nfiles_; ++i )
	  LOGINFO("      file %i : %12llu", i, np_tot_per_file[i], header_.mass[i] );
      }
    
    
    size_t curr_block_buf_size = block_buf_size_;
    
    size_t idcount = 0;
    bool bneed_long_ids = blongids_;
    if( nptot >= 1ul<<32 && !bneed_long_ids )
      {
	bneed_long_ids = true;
	LOGWARN("Need long particle IDs, will write 64bit, make sure to enable in Gadget!");
      }   
    
    for( unsigned ifile=0; ifile<nfiles_; ++ifile )
      {
	
	if( nfiles_ > 1 )
	  {
	    char ffname[256];
	    sprintf(ffname,"%s.%d",fname_.c_str(), ifile);
	    ofs_.open(ffname, std::ios::binary|std::ios::trunc );
	  }else{
	  ofs_.open(fname_.c_str(), std::ios::binary|std::ios::trunc );
	}
	
        
	size_t np_this_file = np_tot_per_file[ifile];
	
	int blksize = sizeof(header);
	
	//... write the header .......................................................
	
	header this_header( header_ );
	for( int i=0; i<6; ++i ){
	  this_header.npart[i] = np_per_file[ifile][i];
	  this_header.npartTotal[i] = (unsigned)np_per_type_[i];
	  this_header.npartTotalHighWord[i] = (unsigned)(np_per_type_[i]>>32);
	}
            
	ofs_.write( (char *)&blksize, sizeof(int) );
	ofs_.write( (char *)&this_header, sizeof(header) );
	ofs_.write( (char *)&blksize, sizeof(int) );
	
	
	//... particle positions ..................................................
	blksize = 3ul*np_this_file*sizeof(T_store);
	ofs_.write( (char *)&blksize, sizeof(int) );
	
	if( bbaryons && np_per_file[ifile][0] > 0ul )
	  {
	    
	    iffs1.open( fnbx, npcdm, wrote_gas*sizeof(T_store) );
	    iffs2.open( fnby, npcdm, wrote_gas*sizeof(T_store) );
	    iffs3.open( fnbz, npcdm, wrote_gas*sizeof(T_store) );
	    
	    npleft = np_per_file[ifile][0];
	    n2read = std::min(curr_block_buf_size,npleft);
	    while( n2read > 0ul )
	      {
		iffs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
		iffs2.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
		iffs3.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
		
		for( size_t i=0; i<n2read; ++i )
		  {
		    adata3.push_back( fmod(tmp1[i]+header_.BoxSize,header_.BoxSize) );
		    adata3.push_back( fmod(tmp2[i]+header_.BoxSize,header_.BoxSize) );
		    adata3.push_back( fmod(tmp3[i]+header_.BoxSize,header_.BoxSize) );
		  }
		ofs_.write( reinterpret_cast<char*>(&adata3[0]), 3*n2read*sizeof(T_store) );
		
		adata3.clear();
		npleft -= n2read;
		n2read = std::min( curr_block_buf_size,npleft );
	      }
	    iffs1.close();
	    iffs2.close();
	    iffs3.close();
	    
            
	  }
	
	npleft = np_this_file - np_per_file[ifile][0];
	n2read = std::min(curr_block_buf_size,npleft);
	
	iffs1.open( fnx, npcdm, wrote_dm*sizeof(T_store) );
	iffs2.open( fny, npcdm, wrote_dm*sizeof(T_store) );
	iffs3.open( fnz, npcdm, wrote_dm*sizeof(T_store) );
	
	while( n2read > 0ul )
	  {
	    iffs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
	    iffs2.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
	    iffs3.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
	    
	    for( size_t i=0; i<n2read; ++i )
	      {
		adata3.push_back( fmod(tmp1[i]+header_.BoxSize,header_.BoxSize) );
		adata3.push_back( fmod(tmp2[i]+header_.BoxSize,header_.BoxSize) );
		adata3.push_back( fmod(tmp3[i]+header_.BoxSize,header_.BoxSize) );
	      }
	    ofs_.write( reinterpret_cast<char*>(&adata3[0]), 3*n2read*sizeof(T_store) );
	    
	    adata3.clear();
	    npleft -= n2read;
	    n2read = std::min( curr_block_buf_size,npleft );
	  }
	ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
	
	iffs1.close();
	iffs2.close();
	iffs3.close();
	
	
	//... particle velocities ..................................................
	blksize = 3ul*np_this_file*sizeof(T_store);
	ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
	
	
	if( bbaryons && np_per_file[ifile][0] > 0ul )
	  {
	    iffs1.open( fnbvx, npcdm, wrote_gas*sizeof(T_store) );
	    iffs2.open( fnbvy, npcdm, wrote_gas*sizeof(T_store) );
	    iffs3.open( fnbvz, npcdm, wrote_gas*sizeof(T_store) );
	    
	    npleft = np_per_file[ifile][0];
	    n2read = std::min(curr_block_buf_size,npleft);
	    while( n2read > 0ul )
	      {
		iffs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
		iffs2.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
		iffs3.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
		
		for( size_t i=0; i<n2read; ++i )
		  {
		    adata3.push_back( tmp1[i] );
		    adata3.push_back( tmp2[i] );
		    adata3.push_back( tmp3[i] );
		  }
		
		ofs_.write( reinterpret_cast<char*>(&adata3[0]), 3*n2read*sizeof(T_store) );
		
		adata3.clear();
		npleft -= n2read;
		n2read = std::min( curr_block_buf_size,npleft );
	      }
	    
	    iffs1.close();
	    iffs2.close();
	    iffs3.close();
	    
	    
	  }
	
	iffs1.open( fnvx, npcdm, wrote_dm*sizeof(T_store) );
	iffs2.open( fnvy, npcdm, wrote_dm*sizeof(T_store) );
	iffs3.open( fnvz, npcdm, wrote_dm*sizeof(T_store) );
	
	npleft = np_this_file - np_per_file[ifile][0];
	n2read = std::min(curr_block_buf_size,npleft);
	while( n2read > 0ul )
	  {
	    iffs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
	    iffs2.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
	    iffs3.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
	    
	    for( size_t i=0; i<n2read; ++i )
	      {
		adata3.push_back( tmp1[i] );
		adata3.push_back( tmp2[i] );
		adata3.push_back( tmp3[i] );
	      }
	    
	    ofs_.write( reinterpret_cast<char*>(&adata3[0]), 3*n2read*sizeof(T_store) );
	    
	    adata3.clear();
	    npleft -= n2read;
	    n2read = std::min( curr_block_buf_size,npleft );
	  }
	ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
	
	iffs1.close();
	iffs2.close();
	iffs3.close();
	
	//... particle IDs ..........................................................
	std::vector<unsigned> short_ids;
	std::vector<size_t> long_ids;
        
	if( bneed_long_ids )
	  long_ids.assign(curr_block_buf_size,0);
	else
	  short_ids.assign(curr_block_buf_size,0);
	
	npleft	= np_this_file;
	n2read	= std::min(curr_block_buf_size,npleft);
	blksize = sizeof(unsigned)*np_this_file;
        
	if( bneed_long_ids )
	  blksize = sizeof(size_t)*np_this_file;
	
	
	//... generate contiguous IDs and store in file ..
	ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
	while( n2read > 0ul )
	  {
	    if( bneed_long_ids )
	      {
                   for( size_t i=0; i<n2read; ++i )
		     long_ids[i] = idcount++;
                   ofs_.write( reinterpret_cast<char*>(&long_ids[0]), n2read*sizeof(size_t) );
                }else{
	      for( size_t i=0; i<n2read; ++i )
		short_ids[i] = idcount++;
	      ofs_.write( reinterpret_cast<char*>(&short_ids[0]), n2read*sizeof(unsigned) );
	    }
	    npleft -= n2read;
	    n2read = std::min( curr_block_buf_size,npleft );
	  }
	ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
	
	std::vector<unsigned>().swap( short_ids );
	std::vector<size_t>().swap( long_ids );
	
	
	//... particle masses .......................................................
	if( bmorethan2bnd_ )//bmultimass_ && bmorethan2bnd_ && nc_per_file[ifile] > 0ul)
	  {
	    unsigned npcoarse = np_per_file[ifile][bndparticletype_];// nc_per_file[ifile];//header_.npart[5];
	    iffs1.open( fnm, np_per_type_[bndparticletype_], wrote_coarse*sizeof(T_store) );
	    
	    npleft  = npcoarse;
	    n2read  = std::min(curr_block_buf_size,npleft);
	    blksize = npcoarse*sizeof(T_store);
	    
	    ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
	    while( n2read > 0ul )
	      {
		iffs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
		ofs_.write( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
		
		npleft -= n2read;
		n2read = std::min( curr_block_buf_size,npleft );
		
	      }
	    ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
	    
	    iffs1.close();
	    
	  }
	
	
	//... initial internal energy for gas particles
	if( bbaryons && np_per_file[ifile][0] > 0ul )
	  {
	    
	    std::vector<T_store> eint(curr_block_buf_size,0.0);
	    
	    const double astart = 1./(1.+header_.redshift);
	    const double npol  = (fabs(1.0-gamma_)>1e-7)? 1.0/(gamma_-1.) : 1.0;
	    const double unitv = 1e5;
	    const double h2    = header_.HubbleParam*header_.HubbleParam;//*0.0001;
	    const double adec  = 1.0/(160.*pow(omegab_*h2/0.022,2.0/5.0));
	    const double Tcmb0 = 2.726;
	    const double Tini  = astart<adec? Tcmb0/astart : Tcmb0/astart/astart*adec;
	    const double mu    = (Tini>1.e4) ? 4.0/(8.-5.*YHe_) : 4.0/(1.+3.*(1.-YHe_));
	    const double ceint = 1.3806e-16/1.6726e-24 * Tini * npol / mu / unitv / unitv;
	    
	    npleft	= np_per_file[ifile][0];
	    n2read	= std::min(curr_block_buf_size,npleft);
	    blksize = sizeof(T_store)*np_per_file[ifile][0]; //*npgas
	    
	    ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
	    while( n2read > 0ul )
	      {
		for( size_t i=0; i<n2read; ++i )
		  eint[i] = ceint;
		ofs_.write( reinterpret_cast<char*>(&eint[0]), n2read*sizeof(T_store) );
		npleft -= n2read;
		n2read = std::min( curr_block_buf_size,npleft );
	      }
	    ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
	    
	    static bool bdisplayed = false;
	    if( !bdisplayed )
	      {
		LOGINFO("Gadget2 : set initial gas temperature to %.2f K/mu",Tini/mu);
		bdisplayed = true;
	      }
	  }
	
	
	ofs_.flush();
	ofs_.close();
	
	wrote_gas       += np_per_file[ifile][0];
	wrote_dm        += np_this_file-np_per_file[ifile][0];
	wrote_coarse    += np_per_file[ifile][5];
	
        
      }
    
    delete[] tmp1;
    delete[] tmp2;
    delete[] tmp3;
    
    remove( fnbx );
    remove( fnby );
    remove( fnbz );
    remove( fnx );
    remove( fny );
    remove( fnz );
    remove( fnbvx );
    remove( fnbvy );
    remove( fnbvz );
    remove( fnvx );
    remove( fnvy );
    remove( fnvz );
    remove( fnm );
    
  }
  
  void determine_particle_numbers( const grid_hierarchy& gh )
  {
    if( !bhave_particlenumbers_ )
      {
	bhave_particlenumbers_ = true;

	double rhoc = 27.7519737; // in h^2 1e10 M_sol / Mpc^3
    
	/*if( kpcunits_ )
	  rhoc *= 1e-9; // in h^2 1e10 M_sol / kpc^3
	
	if( msolunits_ )
	rhoc *= 1e10; // in h^2 M_sol / kpc^3*/

	rhoc /= unit_mass_chosen_ / (unit_length_chosen_*unit_length_chosen_*unit_length_chosen_);
	
	// only type 1 are baryons
	if( !do_baryons_ )
	  header_.mass[1] = header_.Omega0 * rhoc * pow(header_.BoxSize,3.)/pow(2,3*levelmax_);
	else{
	  header_.mass[0] = (omegab_) * rhoc * pow(header_.BoxSize,3.)/pow(2,3*levelmax_);
	  header_.mass[1] = (header_.Omega0-omegab_) * rhoc * pow(header_.BoxSize,3.)/pow(2,3*levelmax_);
	}

	//...
	for( int i=0; i<6; ++i )
	  np_per_type_[i] = 0;
	
	// determine how many particles per type exist, determine their mass
	for( int ilevel=(int)gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
	  {
	    int itype = std::min<int>((int)gh.levelmax()-ilevel+1,5);
	    np_per_type_[itype] += gh.count_leaf_cells(ilevel,ilevel);
	    if( itype > 1 )
	      header_.mass[itype] = header_.Omega0 * rhoc * pow(header_.BoxSize,3.)/pow(2,3*ilevel);
	  }
	
	// if coarse particles should not be spread across types, assign them all to type bndparticletype
	if( !spread_coarse_acrosstypes_ )
	  {
	    if( gh.levelmax() > gh.levelmin()+1 )
	      bmorethan2bnd_ = true;
	    else 
	      bmorethan2bnd_ = false;
	    
	    for( unsigned itype=2; itype<6; ++itype )
	      {
		if( itype == bndparticletype_ ) continue;
		np_per_type_[bndparticletype_] += np_per_type_[itype];
		np_per_type_[itype] = 0;
		header_.mass[itype] = 0.;
	      }
	  }
	
	if( do_baryons_ )
	  np_per_type_[0] = np_per_type_[1];
      }
  }
  
public:
  gadget2_output_plugin( config_file& cf )
  : output_plugin( cf )
  {

    units_mass_.insert( std::pair<std::string,double>( "1e10Msol", 1.0 ) );         // 1e10 M_o/h (default)
    units_mass_.insert( std::pair<std::string,double>( "Msol", 1.0e-10 ) );         // 1 M_o/h
    units_mass_.insert( std::pair<std::string,double>( "Mearth", 3.002e-16 ) );     // 1 M_earth/h

    units_length_.insert( std::pair<std::string,double>( "Mpc", 1.0 ) );            // 1 Mpc/h (default)
    units_length_.insert( std::pair<std::string,double>( "kpc", 1.0e-3 ) );         // 1 kpc/h
    units_length_.insert( std::pair<std::string,double>( "pc", 1.0e-6 ) );          // 1 pc/h

    units_vel_.insert( std::pair<std::string,double>( "km/s", 1.0 ) );              // 1 km/s (default)
    units_vel_.insert( std::pair<std::string,double>( "m/s", 1.0e-3 ) );            // 1 m/s
    units_vel_.insert( std::pair<std::string,double>( "cm/s", 1.0e-5 ) );           // 1 cm/s
      
    block_buf_size_ = cf_.getValueSafe<unsigned>("output","gadget_blksize",1048576);
    
    //... ensure that everyone knows we want to do SPH
    cf.insertValue("setup","do_SPH","yes");
    
    //bbndparticles_  = !cf_.getValueSafe<bool>("output","gadget_nobndpart",false);
    npartmax_ = 1<<30;
    
    nfiles_ = cf.getValueSafe<unsigned>("output","gadget_num_files",1);
    
    blongids_ = cf.getValueSafe<bool>("output","gadget_longids",false);
    
    shift_halfcell_ = cf.getValueSafe<bool>("output","gadget_cell_centered",false);
    
    //if( nfiles_ < (int)ceil((double)npart/(double)npartmax_) )
    //	LOGWARN("Should use more files.");
    
    if (nfiles_ > 1 ) 
      {
	for( unsigned ifile=0; ifile<nfiles_; ++ifile )
	  {
	    char ffname[256];
	    sprintf(ffname,"%s.%d",fname_.c_str(), ifile);
	    ofs_.open(ffname, std::ios::binary|std::ios::trunc );
	    if(!ofs_.good())
	      {	
		LOGERR("gadget-2 output plug-in could not open output file \'%s\' for writing!",ffname);
		throw std::runtime_error(std::string("gadget-2 output plug-in could not open output file \'")+std::string(ffname)+"\' for writing!\n");
	      }
	    ofs_.close();	
	  }
      }else{
      ofs_.open(fname_.c_str(), std::ios::binary|std::ios::trunc );
      if(!ofs_.good())
	{	
	  LOGERR("gadget-2 output plug-in could not open output file \'%s\' for writing!",fname_.c_str());
	  throw std::runtime_error(std::string("gadget-2 output plug-in could not open output file \'")+fname_+"\' for writing!\n");
	}
      ofs_.close();
    }

    bhave_particlenumbers_ = false;
    
    bmorethan2bnd_ = false;
    if( levelmax_ > levelmin_ +4)
      bmorethan2bnd_ = true;
    
    for( int i=0; i<6; ++i )
      {
	header_.npart[i] = 0;
	header_.npartTotal[i] = 0;
	header_.npartTotalHighWord[i] = 0;
	header_.mass[i] = 0.0;
      }
    
    YHe_ = cf.getValueSafe<double>("cosmology","YHe",0.248);
    gamma_ = cf.getValueSafe<double>("cosmology","gamma",5.0/3.0);
    
    do_baryons_ = cf.getValueSafe<bool>("setup","baryons",false);
    omegab_ = cf.getValueSafe<double>("cosmology","Omega_b",0.045);

    //... new way
    std::string lunitstr = cf.getValueSafe<std::string>("output","gadget_lunit","Mpc");
    std::string munitstr = cf.getValueSafe<std::string>("output","gadget_munit","1e10Msol");
    std::string vunitstr = cf.getValueSafe<std::string>("output","gadget_vunit","km/s");


    std::map<std::string,double>::iterator mapit;

    if( (mapit=units_length_.find( lunitstr ))!=units_length_.end() )
      unit_length_chosen_ = (*mapit).second;
    else{
      LOGERR("Gadget: length unit \'%s\' unknown in gadget_lunit",lunitstr.c_str() );
      throw std::runtime_error("Unknown length unit specified for Gadget output plugin");
    }

    if( (mapit=units_mass_.find( munitstr ))!=units_mass_.end() )
      unit_mass_chosen_ = (*mapit).second;
    else{
      LOGERR("Gadget: mass unit \'%s\' unknown in gadget_munit",munitstr.c_str() );
      throw std::runtime_error("Unknown mass unit specified for Gadget output plugin");
    }

    if( (mapit=units_vel_.find( vunitstr ))!=units_vel_.end() )
      unit_vel_chosen_ = (*mapit).second;
    else{
      LOGERR("Gadget: velocity unit \'%s\' unknown in gadget_vunit",munitstr.c_str() );
      throw std::runtime_error("Unknown velocity unit specified for Gadget output plugin");
    }

    //... maintain compatibility with old way of setting units
    if( cf.containsKey("output","gadget_usekpc") )
      {
	kpcunits_ = cf.getValueSafe<bool>("output","gadget_usekpc",false);
	if( kpcunits_ )
	  unit_length_chosen_ = 1e-3;
	LOGWARN("Deprecated option \'gadget_usekpc\' may override unit selection. Use \'gadget_lunit\' instead.");
      }
    if( cf.containsKey("output","gadget_usemsol") )
      {
	msolunits_ = cf.getValueSafe<bool>("output","gadget_usemsol",false);
	if( msolunits_ )
	  unit_mass_chosen_ = 1e-10;
	LOGWARN("Deprecated option \'gadget_usemsol\' may override unit selection. Use \'gadget_munit\' instead.");
      }

    //... coarse particle properties...

    spread_coarse_acrosstypes_ = cf.getValueSafe<bool>("output","gadget_spreadcoarse",false);
    bndparticletype_ = 5;

    if( !spread_coarse_acrosstypes_ )
      {
	bndparticletype_ = cf.getValueSafe<unsigned>("output","gadget_coarsetype",5);
    
	if( bndparticletype_ == 0 || //bndparticletype_ == 1 || bndparticletype_ == 4 ||
	    bndparticletype_ > 5 )
	  {
	    LOGERR("Coarse particles cannot be of Gadget particle type %d in output plugin.", bndparticletype_);
	    throw std::runtime_error("Specified illegal Gadget particle type for coarse particles");
	  }
      }
    else
      {
	if( cf.getValueSafe<unsigned>("output","gadget_coarsetype",5) != 5 )
	  LOGWARN("Gadget: Option \'gadget_spreadcoarse\' forces \'gadget_coarsetype=5\'! Will override.");
      }
    
    //... set time ......................................................
    header_.redshift = cf.getValue<double>("setup","zstart");
    header_.time = 1.0/(1.0+header_.redshift);
    
    //... SF flags
    header_.flag_sfr = 0;
    header_.flag_feedback = 0;
    header_.flag_cooling = 0;
    
    //... 
    header_.num_files = nfiles_;//1;
    header_.BoxSize = cf.getValue<double>("setup","boxlength");
    header_.Omega0 = cf.getValue<double>("cosmology","Omega_m");
    header_.OmegaLambda = cf.getValue<double>("cosmology","Omega_L");
    header_.HubbleParam = cf.getValue<double>("cosmology","H0")/100.0;
    
    header_.flag_stellarage = 0;
    header_.flag_metals = 0;
    
    header_.flag_entropy_instead_u = 0;
    
    //if( kpcunits_ )
    //  header_.BoxSize *= 1000.0;
    header_.BoxSize /= unit_length_chosen_;
    
    for( int i=0; i<empty_fill_bytes; ++i )
      header_.fill[i] = 0;
  }
  
  void write_dm_mass( const grid_hierarchy& gh )
  { 
    determine_particle_numbers( gh );

    double rhoc = 27.7519737; // in h^2 1e10 M_sol / Mpc^3

    // adjust units
    rhoc /= unit_mass_chosen_ / (unit_length_chosen_*unit_length_chosen_*unit_length_chosen_);
    
    /*if( kpcunits_ )
      rhoc *= 1e-9; // in h^2 1e10 M_sol / kpc^3
    
    if( msolunits_ )
      rhoc *= 1e10; // in h^2 M_sol / kpc^3
    */
    
    // if there are more than one kind of coarse particle assigned to the same type,
    // we have to explicitly store their masses
    if( bmorethan2bnd_ )
      {
	header_.mass[bndparticletype_] = 0.;
        
	size_t npcoarse = np_per_type_[bndparticletype_];
	size_t nwritten = 0;
	
	std::vector<T_store> temp_dat;
	temp_dat.reserve(block_buf_size_);
	
	char temp_fname[256];
	sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_mass );
	std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
	
	size_t blksize = sizeof(T_store)*npcoarse;
	
	ofs_temp.write( (char *)&blksize, sizeof(size_t) );
          
	int levelmaxcoarse = gh.levelmax()-4;
	if( !spread_coarse_acrosstypes_ )
	  levelmaxcoarse = gh.levelmax()-1;
	
	for( int ilevel=levelmaxcoarse; ilevel>=(int)gh.levelmin(); --ilevel )
	  {
	    // baryon particles live only on finest grid
	    // these particles here are total matter particles
	    double pmass = header_.Omega0 * rhoc * pow(header_.BoxSize,3.)/pow(2,3*ilevel);	
	    
	    for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
	      for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
		for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
		  if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
		    {
		      if( temp_dat.size() <  block_buf_size_ )
			temp_dat.push_back( pmass );	
		      else
			{
			  ofs_temp.write( (char*)&temp_dat[0], sizeof(T_store)*block_buf_size_ );	
			  nwritten += block_buf_size_;
			  temp_dat.clear();
			  temp_dat.push_back( pmass );	
			}
		    }
	  }
	
	if( temp_dat.size() > 0 )
	  {	
	    ofs_temp.write( (char*)&temp_dat[0], sizeof(T_store)*temp_dat.size() );		
	    nwritten+=temp_dat.size();
	  }
	
	if( nwritten != npcoarse ){
	  LOGERR("nwritten = %llu != npcoarse = %llu\n",nwritten,npcoarse);
	  throw std::runtime_error("Internal consistency error while writing temporary file for masses");
	}

	ofs_temp.write( (char *)&blksize, sizeof(size_t) );
	
	if( ofs_temp.bad() )
	  throw std::runtime_error("I/O error while writing temporary file for masses");
	
      }
  }
	
	
  void write_dm_position( int coord, const grid_hierarchy& gh )
  {
    //... count number of leaf cells ...//
    determine_particle_numbers( gh );

    size_t npart = 0;
    for( int i=1; i<6; ++i )
      npart += np_per_type_[i];
    
    //... determine if we need to shift the coordinates back
    double *shift = NULL;
    
    if( shift_halfcell_ )
      {
	double h = 1.0/(1<<(levelmin_+1));
	shift = new double[3];
	shift[0] = shift[1] = shift[2] = -h;
      }
      
    size_t nwritten = 0;
    //... collect displacements and convert to absolute coordinates with correct
    //... units
    std::vector<T_store> temp_data;
    temp_data.reserve( block_buf_size_ );
    
    char temp_fname[256];
    sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_pos+coord );
    std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
    
    size_t blksize = sizeof(T_store)*npart;
    ofs_temp.write( (char *)&blksize, sizeof(size_t) );
    
    double xfac = header_.BoxSize;
    
    for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
      for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
	for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
	  for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
	    if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
	      {
		double xx[3];
		gh.cell_pos(ilevel, i, j, k, xx);
		if( shift != NULL )
		  xx[coord] += shift[coord];
		
		xx[coord] = (xx[coord]+(*gh.get_grid(ilevel))(i,j,k))*xfac;
		
		if( temp_data.size() < block_buf_size_ )
		  temp_data.push_back( xx[coord] );
		else
		  {
		    ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*block_buf_size_ );
		    nwritten += block_buf_size_;
		    temp_data.clear();
		    temp_data.push_back( xx[coord] );
		  }
	      }
    
    if( temp_data.size() > 0 )
      {	
	ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*temp_data.size() );
	nwritten += temp_data.size();
      }
    
    if( nwritten != npart )
      throw std::runtime_error("Internal consistency error while writing temporary file for positions");
    
    //... dump to temporary file
    ofs_temp.write( (char *)&blksize, sizeof(size_t) );
    
    if( ofs_temp.bad() )
      throw std::runtime_error("I/O error while writing temporary file for positions");
    
    ofs_temp.close();
    
    if( shift != NULL )
      delete[] shift;
    
  }
  
  void write_dm_velocity( int coord, const grid_hierarchy& gh )
  {
    //... count number of leaf cells ...//
    determine_particle_numbers( gh );

    size_t npart = 0;
    for( int i=1; i<6; ++i )
      npart += np_per_type_[i];
      
    //... collect displacements and convert to absolute coordinates with correct
    //... units
    std::vector<T_store> temp_data;
    temp_data.reserve( block_buf_size_ );
    
    float isqrta = 1.0f/sqrt(header_.time);
    float vfac = isqrta*header_.BoxSize;
    
    //if( kpcunits_ )
    //  vfac /= 1000.0;
    vfac *= unit_length_chosen_ / unit_vel_chosen_;
    
    size_t nwritten = 0;
    
    char temp_fname[256];
    sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_vel+coord );
    std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
    
    size_t blksize = sizeof(T_store)*npart;
    ofs_temp.write( (char *)&blksize, sizeof(size_t) );
    
    for( int ilevel=levelmax_; ilevel>=(int)levelmin_; --ilevel )
      for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
	for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
	  for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
	    if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
	      {
		if( temp_data.size() < block_buf_size_ )
		  temp_data.push_back( (*gh.get_grid(ilevel))(i,j,k) * vfac );
		else 
		  {
		    ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*block_buf_size_ );
		    nwritten += block_buf_size_;
		    temp_data.clear();
		    temp_data.push_back( (*gh.get_grid(ilevel))(i,j,k) * vfac );
		  }
		
	      }
    if( temp_data.size() > 0 )
      {	
	ofs_temp.write( (char*)&temp_data[0], temp_data.size()*sizeof(T_store) );
	nwritten += temp_data.size();
      }
    
    if( nwritten != npart )
      throw std::runtime_error("Internal consistency error while writing temporary file for velocities");
      
    ofs_temp.write( (char *)&blksize, sizeof(int) );
    
    if( ofs_temp.bad() )
      throw std::runtime_error("I/O error while writing temporary file for velocities");
    
    ofs_temp.close();
  }
  
  void write_dm_density( const grid_hierarchy& gh )
  {
    //... we don't care about DM density for Gadget  
  }
  
  void write_dm_potential( const grid_hierarchy& gh )
  {
    //... we don't care about DM potential for Gadget  
  }
	
  void write_gas_potential( const grid_hierarchy& gh )
  {
    //... we don't care about gas potential for Gadget  
  }
  
  //... write data for gas -- don't do this
  void write_gas_velocity( int coord, const grid_hierarchy& gh )
  {	
    determine_particle_numbers( gh );
    size_t npart = 0;
    for( int i=1; i<6; ++i )
      npart += np_per_type_[i];
    
    //... collect velocities and convert to absolute coordinates with correct
    //... units
    std::vector<T_store> temp_data;
    temp_data.reserve( block_buf_size_ );
    
    float isqrta = 1.0f/sqrt(header_.time);
    float vfac = isqrta*header_.BoxSize;
    
    //if( kpcunits_ )
    //  vfac /= 1000.0;
    vfac *= unit_length_chosen_ / unit_vel_chosen_;
    
    //size_t npart = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());;;
    size_t nwritten = 0;
    
    char temp_fname[256];
    sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_gas_vel+coord );
    std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
    
    size_t blksize = sizeof(T_store)*npart;
    ofs_temp.write( (char *)&blksize, sizeof(size_t) );
    
    
    for( int ilevel=levelmax_; ilevel>=(int)levelmin_; --ilevel )
      for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
	for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
	  for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
	    if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
	      {
		if( temp_data.size() < block_buf_size_ )
		  temp_data.push_back( (*gh.get_grid(ilevel))(i,j,k) * vfac );
		else 
		  {
		    ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*block_buf_size_ );
		    nwritten += block_buf_size_;
		    temp_data.clear();
		    temp_data.push_back( (*gh.get_grid(ilevel))(i,j,k) * vfac );
		  }
		
	      }
    
    
    if( temp_data.size() > 0 )
      {	
	ofs_temp.write( (char*)&temp_data[0], temp_data.size()*sizeof(T_store) );
	nwritten += temp_data.size();
      }
    
    if( nwritten != npart )
      throw std::runtime_error("Internal consistency error while writing temporary file for gas velocities");
    
    ofs_temp.write( (char *)&blksize, sizeof(int) );
    
    if( ofs_temp.bad() )
      throw std::runtime_error("I/O error while writing temporary file for gas velocities");
    
    ofs_temp.close();
  }
  
  
  //... write only for fine level
  void write_gas_position( int coord, const grid_hierarchy& gh )
  {	
    //... count number of leaf cells ...//
    determine_particle_numbers( gh );

    size_t npart = 0;
    for( int i=1; i<6; ++i )
      npart += np_per_type_[i];
    
    //... determine if we need to shift the coordinates back
    double *shift = NULL;
    
    if( shift_halfcell_ )
      {
	double h = 1.0/(1<<(levelmin_+1));
	shift = new double[3];
	shift[0] = shift[1] = shift[2] = -h;
      }
    
    size_t nwritten = 0;
    
    //...
    //... collect displacements and convert to absolute coordinates with correct
    //... units
    std::vector<T_store> temp_data;
    temp_data.reserve( block_buf_size_ );
    
    char temp_fname[256];
    sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_gas_pos+coord );
    std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
    
    size_t blksize = sizeof(T_store)*npart;
    ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
    double xfac = header_.BoxSize;
    
    double h = 1.0/(1ul<<gh.levelmax());
    
    for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
      {	
   	for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
	  for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
	    for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
	      //if( ! gh.is_refined(ilevel,i,j,k) )
	      if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
		{
		  double xx[3];
		  gh.cell_pos(ilevel, i, j, k, xx);
		  if( shift != NULL )
		    xx[coord] += shift[coord];
		  
		  //... shift particle positions (this has to be done as the same shift
		  //... is used when computing the convolution kernel for SPH baryons)
		  xx[coord] += 0.5*h;
                  
		  xx[coord] = (xx[coord]+(*gh.get_grid(ilevel))(i,j,k))*xfac;
		  
		  if( temp_data.size() < block_buf_size_ )
		    temp_data.push_back( xx[coord] );
		  else
		    {
		      ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*block_buf_size_ );
		      nwritten += block_buf_size_;
		      temp_data.clear();
		      temp_data.push_back( xx[coord] );
		    }
		}
      }
    
    if( temp_data.size() > 0 )
      {	
	ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*temp_data.size() );
	nwritten += temp_data.size();
      }
    
    if( nwritten != npart )
      throw std::runtime_error("Internal consistency error while writing temporary file for gas positions");
    
    //... dump to temporary file
    ofs_temp.write( (char *)&blksize, sizeof(size_t) );
    
    if( ofs_temp.bad() )
      throw std::runtime_error("I/O error while writing temporary file for gas positions");
    
    ofs_temp.close();
    
    if( shift != NULL )
      delete[] shift;
  }
  
  void write_gas_density( const grid_hierarchy& gh )
  {	
    //do nothing as we write out positions
  }
  
  void finalize( void )
  {	
    this->assemble_gadget_file();
  }
};

namespace{
	output_plugin_creator_concrete< gadget2_output_plugin<float> > creator1("gadget2");
#ifndef SINGLE_PRECISION
	output_plugin_creator_concrete< gadget2_output_plugin<double> > creator2("gadget2_double");
#endif
}
