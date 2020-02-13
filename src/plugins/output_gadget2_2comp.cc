/*
 
 output_gadget2.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#include <fstream>
#include "log.hh"
#include "output.hh"
#include "mg_interp.hh"
#include "mesh.hh"

template< typename T_store=float >
class gadget2_2comp_output_plugin : public output_plugin
{
protected:
	
	std::ofstream ofs_;
	bool bmultimass_;
	
	
	typedef struct io_header
	{
		unsigned int npart[6];                        
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
                int flag_doubleprecision;
                int flag_ic_info; 
		char fill[52];                       
	}header;                       
	
	
	header header_;
	
	std::string fname;
	
	bool do_glass_;
	std::string fname_glass_baryon_, fname_glass_cdm_;
	
	enum iofields {
		id_dm_mass, id_dm_vel, id_dm_pos, id_gas_vel, id_gas_rho, id_gas_temp, id_gas_pos
	};
	
	size_t np_fine_gas_, np_fine_dm_, np_coarse_dm_;
	
	size_t block_buf_size_;
	unsigned long long npartmax_;
	unsigned nfiles_;
	
	//bool bbndparticles_;
	bool bmorethan2bnd_;
	bool kpcunits_;
	double YHe_;
	
	void distribute_particles( unsigned nfiles, size_t nfine_dm, size_t nfine_gas, size_t ncoarse, 
                              std::vector<unsigned>& nfdm_pf, std::vector<unsigned>& nfgas_pf, std::vector<unsigned>& nc_pf )
    {
        nfdm_pf.assign( nfiles, 0 );
        nfgas_pf.assign( nfiles, 0 );
        nc_pf.assign( nfiles, 0 );
        
        size_t ntotal = nfine_dm + nfine_gas + ncoarse;
        size_t nnominal = (size_t)((double)ntotal/(double)nfiles);
        
        size_t nf_dm_assigned = 0, nf_gas_assigned = 0, nc_assigned = 0;
        
        for( unsigned i=0; i<nfiles; ++i )
        {
            if( nfine_gas > 0 )
            {
                nfdm_pf[i] = std::min( nnominal/2ul, nfine_dm-nf_dm_assigned );
                nf_dm_assigned += nfdm_pf[i];
                nfgas_pf[i] = std::min( nnominal/2ul, nfine_gas-nf_gas_assigned );
                nf_gas_assigned += nfgas_pf[i];
                
            }else{
                nfdm_pf[i] = std::min( nnominal, nfine_dm-nf_dm_assigned );
                nf_dm_assigned += nfdm_pf[i];
            }
            
            // once all fine particles are assigned, start with the coarse
            if( nf_dm_assigned+nf_gas_assigned == nfine_dm+nfine_gas )
            {
                nc_pf[i] = std::min( nnominal-(size_t)(nfdm_pf[i]+nfgas_pf[i]), ncoarse-nc_assigned );
                nc_assigned += nc_pf[i];
            }
            
        }
        
        // make sure all particles are assigned
        nfdm_pf[ nfiles-1 ]     += nfine_dm-nf_dm_assigned;
        nfgas_pf[ nfiles-1 ]    += nfine_gas-nf_gas_assigned;
        nc_pf[ nfiles-1 ]       += ncoarse-nc_assigned;
        
    }
	
	std::ifstream& open_and_check( std::string ffname, size_t npart )
	{
		std::ifstream ifs( ffname.c_str(), std::ios::binary );
		unsigned long long blk, expected;
		ifs.read( (char*)&blk, sizeof(unsigned long long) );
                expected = ((unsigned long long) npart*(unsigned long long)sizeof(T_store));
		if( blk != expected )
		{	
			LOGERR("Internal consistency error in gadget2 output plug-in, open_and_check");
			LOGERR("Expected %d particles (%lld bytes) in temp file %s but found %lld",npart, expected ,ffname.c_str(), blk);
			//throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
		}
		
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
		{
			
		}
		
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
		{
			
		}
		
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
	
	void assemble_gadget_file( void )
	{
		
		
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
		
		/*const size_t 
			nptot = header_.npart[1]+header_.npart[2]+header_.npart[5],
			npgas = header_.npart[2],
			npcdm = nptot-npgas;*/
		
		const size_t 
			nptot = np_fine_gas_+np_fine_dm_+np_coarse_dm_,
			//npgas = np_fine_gas_,
			npcdm = np_fine_dm_+np_coarse_dm_;
			
		size_t
			wrote_coarse = 0,
			wrote_gas  = 0,
			wrote_dm   = 0;
		
		size_t
			npleft = nptot, 
			n2read = std::min(block_buf_size_,npleft);
		
		if( header_.npart[5] > 0 )
			LOGERR("Multi-resolution setup not supported for 2comp hack");
		
		std::cout << " - Writing " << nptot << " particles to Gadget file...\n"
				  << "      type 1 : " << header_.npart[1] << "\n"
				  << "      type 2 : " << header_.npart[2] << "\n"
				  << "      type 5 : " << header_.npart[5] << "\n";
		
		bool bbaryons = np_fine_gas_ > 0;
				
		std::vector<T_store> adata3;
		adata3.reserve( 3*block_buf_size_ );
		T_store *tmp1, *tmp2, *tmp3;
		
		tmp1 = new T_store[block_buf_size_];
		tmp2 = new T_store[block_buf_size_];
		tmp3 = new T_store[block_buf_size_];
		
		std::vector<unsigned> nfdm_per_file, nfgas_per_file, nc_per_file;
        distribute_particles( nfiles_, np_fine_dm_, np_fine_gas_, np_coarse_dm_,
                             nfdm_per_file, nfgas_per_file, nc_per_file );
		
		
		if( nfiles_ > 1 )
		{
			std::cout << " - Gadget2 : distributing particles to " << nfiles_ << " files\n"
			<< "                 " << std::setw(12) << "type 1" << "," << std::setw(12) << "type 2" << "," << std::setw(12) << "type 5" << std::endl;
			for( unsigned i=0; i<nfiles_; ++i )
			{
				std::cout << "      file " << std::setw(3) << i << " : " 
				<< std::setw(12) << nfdm_per_file[i] << "," 
				<< std::setw(12) << nfgas_per_file[i] << "," 
				<< std::setw(12) << nc_per_file[i] << std::endl;
			}			
		}
		
		
		size_t curr_block_buf_size = block_buf_size_;
        
        size_t idcount = 0;
        bool bneed_long_ids = false;
        if( nptot >= 1ul<<32 )
        {
            bneed_long_ids = true;
            LOGWARN("Need long particle IDs, make sure to enable in Gadget!");
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
			
            
			size_t np_this_file = nfgas_per_file[ifile] + nfdm_per_file[ifile] + nc_per_file[ifile];
			
			int blksize = sizeof(header);
			
			//... write the header .......................................................
			
			header this_header( header_ );
            this_header.npart[1] = nfdm_per_file[ifile];
            this_header.npart[2] = nfgas_per_file[ifile];
            this_header.npart[5] = nc_per_file[ifile];
			
			
			ofs_.write( (char *)&blksize, sizeof(int) );
			ofs_.write( (char *)&this_header, sizeof(header) );
			ofs_.write( (char *)&blksize, sizeof(int) );
			
			
			//... particle positions ..................................................
			blksize = 3ul*np_this_file*sizeof(T_store);
			ofs_.write( (char *)&blksize, sizeof(int) );
			
			npleft = nfdm_per_file[ifile];//+nc_per_file[ifile];
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
			
			iffs1.close();
			iffs2.close();
			iffs3.close();
			
			if( bbaryons && nfgas_per_file[ifile] > 0ul )
			{
				
				iffs1.open( fnbx, npcdm, wrote_gas*sizeof(T_store) );
				iffs2.open( fnby, npcdm, wrote_gas*sizeof(T_store) );
				iffs3.open( fnbz, npcdm, wrote_gas*sizeof(T_store) );
				
				npleft = nfgas_per_file[ifile];
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
			
			ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
			
			
			
			
			//... particle velocities ..................................................
			blksize = 3ul*np_this_file*sizeof(T_store);
			ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
			
			iffs1.open( fnvx, npcdm, wrote_dm*sizeof(T_store) );
			iffs2.open( fnvy, npcdm, wrote_dm*sizeof(T_store) );
			iffs3.open( fnvz, npcdm, wrote_dm*sizeof(T_store) );
			
			npleft = nfdm_per_file[ifile];//+nc_per_file[ifile];
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
			
			if( bbaryons && nfgas_per_file[ifile] > 0ul )
			{
				iffs1.open( fnbvx, npcdm, wrote_gas*sizeof(T_store) );
				iffs2.open( fnbvy, npcdm, wrote_gas*sizeof(T_store) );
				iffs3.open( fnbvz, npcdm, wrote_gas*sizeof(T_store) );
				
				npleft = nfgas_per_file[ifile];
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
			
			ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
			
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
			// multi-mass not supported here
			
			
			ofs_.flush();
			ofs_.close();
			
            wrote_gas       += nfgas_per_file[ifile];
            wrote_dm        += nfdm_per_file[ifile] + nc_per_file[ifile];
            wrote_coarse    += nc_per_file[ifile];
			
            
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
	
	void get_cic_displacement( size_t icoord, const float* ppos, size_t np, float l, const grid_hierarchy& gh, T_store* valp )
	{
		size_t N = gh.size(gh.levelmax(), 0);
		
		float facconv   = 1.f / l * (float)N/(float)(1ul<<levelmax_);
		float suboffset = (float)(gh.offset_abs(levelmax_, icoord))/((float)(1ul<<levelmax_));
		
		for( size_t ip=0; ip < np; ++ip )
		{
			float u,v,w;
			
			u = ppos[3*ip+0] / l * (float)N;
			v = ppos[3*ip+1] / l * (float)N;
			w = ppos[3*ip+2] / l * (float)N;
			
			int i,j,k;
			
			i = (((int)u)+N)%N;
			j = (((int)v)+N)%N;
			k = (((int)w)+N)%N;
			
			u -= (float)i;
			v -= (float)j;
			w -= (float)k;
			
			int i1,j1,k1;
			i1 = (i+1+N)%N;
			j1 = (j+1+N)%N;
			k1 = (k+1+N)%N;
			
			double f1,f2,f3,f4,f5,f6,f7,f8;
			
			f1 = (1.f - u) * (1.f - v) * (1.f - w);
			f2 = (1.f - u) * (1.f - v) * (w);
			f3 = (1.f - u) * (v) * (1.f - w);
			f4 = (1.f - u) * (v) * (w);
			f5 = (u) * (1.f - v) * (1.f - w);
			f6 = (u) * (1.f - v) * (w); 
			f7 = (u) * (v) * (1.f - w);
			f8 = (u) * (v) * (w);
			
			float disp = 0.0f;
			
			//disp += add_term;
			disp += suboffset;
			disp += ppos[3*ip+icoord] * facconv;
			
			disp += f1*(*gh.get_grid(levelmax_))(i,j,k);
			disp += f2*(*gh.get_grid(levelmax_))(i,j,k1);
			disp += f3*(*gh.get_grid(levelmax_))(i,j1,k);
			disp += f4*(*gh.get_grid(levelmax_))(i,j1,k1);
			disp += f5*(*gh.get_grid(levelmax_))(i1,j,k);
			disp += f6*(*gh.get_grid(levelmax_))(i1,j,k1);
			disp += f7*(*gh.get_grid(levelmax_))(i1,j1,k);
			disp += f8*(*gh.get_grid(levelmax_))(i1,j1,k1);
			
			
			
			disp = fmodf( (1.0f+disp)*header_.BoxSize, header_.BoxSize );
			
			valp[ip] = disp;
			
		}
	}
	
	void get_cic_velocity( const float* ppos, size_t np, float l, const grid_hierarchy& gh, T_store* valp )
	{
		float isqrta = 1.0f/sqrt(header_.time);
		float vfac = isqrta*header_.BoxSize;
		
		if( kpcunits_ )
			vfac /= 1000.0;
		
		size_t N = gh.size(gh.levelmax(), 0);
		//float facconv   = 1.f / l * (float)N/(float)(1ul<<levelmax_);
		
		for( size_t ip=0; ip < np; ++ip )
		{
			float u,v,w;
			
			u = ppos[3*ip+0] / l * (float)N;
			v = ppos[3*ip+1] / l * (float)N;
			w = ppos[3*ip+2] / l * (float)N;
			
			int i,j,k;
			
			i = (((int)u)+N)%N;
			j = (((int)v)+N)%N;
			k = (((int)w)+N)%N;
			
			u -= (float)i;
			v -= (float)j;
			w -= (float)k;
			
			int i1,j1,k1;
			i1 = (i+1+N)%N;
			j1 = (j+1+N)%N;
			k1 = (k+1+N)%N;
			
			double f1,f2,f3,f4,f5,f6,f7,f8;
			
			f1 = (1.f - u) * (1.f - v) * (1.f - w);
			f2 = (1.f - u) * (1.f - v) * (w);
			f3 = (1.f - u) * (v) * (1.f - w);
			f4 = (1.f - u) * (v) * (w);
			f5 = (u) * (1.f - v) * (1.f - w);
			f6 = (u) * (1.f - v) * (w); 
			f7 = (u) * (v) * (1.f - w);
			f8 = (u) * (v) * (w);
			
			float vel = 0.0f;
			
			vel += f1*(*gh.get_grid(levelmax_))(i,j,k);
			vel += f2*(*gh.get_grid(levelmax_))(i,j,k1);
			vel += f3*(*gh.get_grid(levelmax_))(i,j1,k);
			vel += f4*(*gh.get_grid(levelmax_))(i,j1,k1);
			vel += f5*(*gh.get_grid(levelmax_))(i1,j,k);
			vel += f6*(*gh.get_grid(levelmax_))(i1,j,k1);
			vel += f7*(*gh.get_grid(levelmax_))(i1,j1,k);
			vel += f8*(*gh.get_grid(levelmax_))(i1,j1,k1);
			
			vel *= vfac;			
			
			valp[ip] = vel;
			
		}
	}
	
	
public:
	
	bool do_baryons_;
	double omegab_, omegac_, omegam_;
	double gamma_;
	
	gadget2_2comp_output_plugin( config_file& cf )
	: output_plugin( cf )//, ofs_( fname_.c_str(), std::ios::binary|std::ios::trunc )	
	{
		block_buf_size_ = cf_.getValueSafe<unsigned>("output","gadget_blksize",2*1048576);
		
		//... ensure that everyone knows we want to do SPH
		cf.insertValue("setup","do_SPH","yes");
		
		//bbndparticles_  = !cf_.getValueSafe<bool>("output","gadget_nobndpart",false);
		npartmax_ = 1<<30;
		
		nfiles_ = cf.getValueSafe<unsigned>("output","gadget_num_files",1);
		
		
		
		/****************************************/
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
		
		bmorethan2bnd_ = false;
		if( levelmax_ > levelmin_ +1)
			bmorethan2bnd_ = true;

		bmultimass_ = true;
		if( levelmax_ == levelmin_ )
			bmultimass_ = false;
			
		
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
		
		//... write displacements in kpc/h rather than Mpc/h?
		kpcunits_ = cf.getValueSafe<bool>("output","gadget_usekpc",false);
		
		do_glass_ = cf.getValueSafe<bool>("output","glass", false);
		if( do_glass_ )
		{
			LOGINFO("Will use provided glass rather than Cartesian mesh for particle placement.");
			
			fname_glass_cdm_ = cf.getValue<std::string>("output","glass_file_cdm");
			
			if( do_baryons_ )
				fname_glass_baryon_ = fname_glass_cdm_;//cf.getValue<std::string>("output","glass_file_baryon");
		}
		
		
		//... set time ......................................................
		header_.redshift = cf.getValue<double>("setup","zstart");
		header_.time = 1.0/(1.0+header_.redshift);
		
		//... SF flags
		header_.flag_sfr = 0;
		header_.flag_feedback = 0;
		header_.flag_cooling = 0;
		
		//... 
		header_.num_files = nfiles_;
		header_.BoxSize = cf.getValue<double>("setup","boxlength");
		header_.Omega0 = cf.getValue<double>("cosmology","Omega_m");
                omegam_ = header_.Omega0;
                omegac_ = omegam_ - omegab_;
        
		header_.OmegaLambda = cf.getValue<double>("cosmology","Omega_L");
		header_.HubbleParam = cf.getValue<double>("cosmology","H0");
		
		header_.flag_stellarage = 0;
		header_.flag_metals = 0;
		
	        	
		header_.flag_entropy_instead_u = 0;
#ifdef SINGLE_PRECISION
                header_.flag_doubleprecision = 0;
#else
                header_.flag_doubleprecision = 1; 
#endif
                std::cout << "header_.flag_doubleprecision " << header_.flag_doubleprecision << "\n";
                header_.flag_ic_info = 0;
		
		if( kpcunits_ )
			header_.BoxSize *= 1000.0;
	}
	
	
	void write_dm_mass( const grid_hierarchy& gh )
	{
		double rhoc = 27.7519737; // in h^2 1e10 M_sol / Mpc^3
		
		if( kpcunits_ )
			rhoc *= 10.0; // in h^2 M_sol / kpc^3
		
		if(! do_glass_ )
		{
			if( !do_baryons_ )
				header_.mass[1] = header_.Omega0 * rhoc * pow(header_.BoxSize,3.)/pow(2,3*levelmax_);
			else
				header_.mass[1] = (header_.Omega0-omegab_) * rhoc * pow(header_.BoxSize,3.)/pow(2,3*levelmax_);			
		}
				
		if( bmorethan2bnd_ )
		{
			unsigned long long npcoarse = gh.count_leaf_cells(gh.levelmin(), gh.levelmax()-1);
			unsigned long long nwritten = 0;
			
			std::vector<T_store> temp_dat;
			temp_dat.reserve(block_buf_size_);
			
			char temp_fname[256];
			sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_mass );
			std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
			
			unsigned long long blksize = sizeof(T_store)*npcoarse;
			
			ofs_temp.write( (char *)&blksize, sizeof(unsigned long long) );
			
			for( int ilevel=gh.levelmax()-1; ilevel>=(int)gh.levelmin(); --ilevel )
			{
				double pmass = 0.0;
				
				if( !do_baryons_ )
					pmass = header_.Omega0 * rhoc * pow(header_.BoxSize,3.)/pow(2,3*ilevel);		
				else
					pmass = (header_.Omega0-omegab_) * rhoc * pow(header_.BoxSize,3.)/pow(2,3*ilevel);
					
				for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
					for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
						for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
							if( ! gh.is_refined(ilevel,i,j,k) )
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
			
			if( nwritten != npcoarse )
				throw std::runtime_error("Internal consistency error while writing temporary file for masses");
			
			ofs_temp.write( (char *)&blksize, sizeof(unsigned long long) );
			
			if( ofs_temp.bad() )
				throw std::runtime_error("I/O error while writing temporary file for masses");
			
		}
		else if( gh.levelmax() != gh.levelmin() )
		{
			header_.mass[5] = header_.Omega0 * rhoc * pow(header_.BoxSize,3.)/pow(2,3*levelmin_);
		}
	}
	
	
	void write_dm_position( int coord, const grid_hierarchy& gh )
	{
		//... count number of leaf cells ...//
		unsigned long long npcoarse = 0, npfine = 0;
		
		npfine   = gh.count_leaf_cells(gh.levelmax(), gh.levelmax());
		if( bmultimass_ )
			npcoarse = gh.count_leaf_cells(gh.levelmin(), gh.levelmax()-1);
		
		
		np_fine_dm_   = npfine;
        np_fine_gas_  = do_baryons_? npfine : 0ul;
        np_coarse_dm_ = npcoarse;
		
		//... determine if we need to shift the coordinates back
		double *shift = NULL;
		
		if( cf_.getValueSafe<bool>("output","shift_back",false ) )
		{
			if( coord == 0 )
				std::cout << " - gadget2 output plug-in will shift particle positions back...\n";
			
			double h = 1.0/(1<<levelmin_);
			shift = new double[3];
			shift[0] = -(double)cf_.getValue<int>( "setup", "shift_x" )*h;
			shift[1] = -(double)cf_.getValue<int>( "setup", "shift_y" )*h;
			shift[2] = -(double)cf_.getValue<int>( "setup", "shift_z" )*h;
		}
		
		size_t npart = npfine+npcoarse;
		size_t nwritten = 0;
		size_t blksize;
		
		//...
		header_.npart[1] = npfine;
		header_.npart[5] = npcoarse;
		header_.npartTotal[1] = (unsigned)npfine;
		header_.npartTotal[5] = (unsigned)npcoarse;
		header_.npartTotalHighWord[1] = (unsigned)(npfine>>32);
		header_.npartTotalHighWord[5] = (unsigned)(npfine>>32);
		
		//header_.num_files = (int)ceil((double)npart/(double)npartmax_);
		
		//... collect displacements and convert to absolute coordinates with correct
		//... units
		std::vector<T_store> temp_data;
		temp_data.reserve( block_buf_size_ );
		
		double xfac = header_.BoxSize;
		
		char temp_fname[256];
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_pos+coord );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
        
        //... if baryons are present, then stagger the two fields
        if( do_baryons_ && !do_glass_ )
        {
			
			
            double h = 1. / (1<<gh.levelmax());
            
            if( shift == NULL )
            {
                shift = new double[3];
                shift[0] = 0.0;
                shift[1] = 0.0;
                shift[2] = 0.0;
                
            }
            shift[0] -= 0.5 * h * omegab_/omegam_;
            shift[1] -= 0.5 * h * omegab_/omegam_;
            shift[2] -= 0.5 * h * omegab_/omegam_;
        }
		
		if(! do_glass_ )
		{
			
			blksize = sizeof(T_store)*npart;
			ofs_temp.write( (char *)&blksize, sizeof(unsigned long long) );
			
			for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
				for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
					for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
						for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
							if( ! gh.is_refined(ilevel,i,j,k) )
							{
								double xx[3];
								gh.cell_pos(ilevel, i, j, k, xx);
								if( shift != NULL )
									xx[coord] += shift[coord];
								
								xx[coord] = fmod( (xx[coord]+(*gh.get_grid(ilevel))(i,j,k))*xfac + header_.BoxSize, header_.BoxSize );
								
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
		else
		{
			
			std::ifstream ofg( fname_glass_cdm_.c_str(), std::ios::binary );
			
			if( !ofg.good() )
				LOGERR("could not open glass input file \'%s\'",fname_glass_cdm_.c_str());
			
			io_header glasshead;
			unsigned blksz;
			
			ofg.read( reinterpret_cast<char*>(&blksz), sizeof(unsigned) );
			assert( blksz == sizeof(io_header) );
			
			ofg.read( reinterpret_cast<char*>(&glasshead), sizeof( io_header ) );
			
			//size_t nreq = gh.size(gh.levelmax(), 0)*gh.size(gh.levelmax(), 1)*gh.size(gh.levelmax(), 2);
			/*if( nreq != (size_t)glasshead.npart[1] )
			{
				LOGERR("glass file contains %d particles, but should contain %ld",glasshead.npart[1],nreq);
				throw std::runtime_error("glass file does not contain the right amount of particles");
			}*/
			
			ofg.read( reinterpret_cast<char*>(&blksz), sizeof(unsigned) );
			ofg.read( reinterpret_cast<char*>(&blksz), sizeof(unsigned) );
			//assert( blksz == glasshead.npart[1]*sizeof(float)*3 );
			
			
			float lglass = glasshead.BoxSize;
			
			
			blksize = sizeof(T_store)*glasshead.npart[1];
			//ofs_temp.write( (char *)&blksize, sizeof(int) );
			ofs_temp.write( (char *)&blksize, sizeof(unsigned long long) );
						
			header_.npart[1] = glasshead.npart[1];
			header_.npartTotal[1] = glasshead.npartTotal[1];
			header_.npartTotalHighWord[1] = 0;

			double rhoc = 27.7519737;
			if( kpcunits_ )
				rhoc *= 10.0; // in h^2 M_sol / kpc^3
			
            if( do_baryons_ )
                header_.mass[1] = omegac_ * rhoc * pow(header_.BoxSize,3.)/(glasshead.npart[1]);
			else
                header_.mass[1] = omegam_ * rhoc * pow(header_.BoxSize,3.)/(glasshead.npart[1]);
            
			// read glass, do interpolation and write
			size_t npartdone=0;
			size_t npinter = glasshead.npart[1];
			
			blksize = sizeof(T_store)*npinter;
			ofs_temp.write( (char *)&blksize, sizeof(unsigned long long) );
			
			float *pos_tmp = new float[3*block_buf_size_];
			temp_data.assign( block_buf_size_, 0.0 );
			
			while( npartdone < npinter )
			{
				size_t npart2read = std::min(npinter-npartdone,block_buf_size_);
				
				ofg.read( reinterpret_cast<char*>(&pos_tmp[0]), npart2read*sizeof(float)*3 );
				get_cic_displacement( coord, pos_tmp, npart2read, lglass,  gh, &temp_data[0] );
				ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*npart2read );

                        //        std::cout << "npart2read " << npart2read << "\n";
                        //        std::cout << "pos temp " << temp_data[0] << " " << temp_data[1] << " " << temp_data[2] << "\n";
                        //        std::cout << "pos temp " << temp_data[npart2read-1] << " " << temp_data[npart2read-2] << " " << temp_data[npart2read-3] <<"\n";
				
				npartdone += npart2read;
				nwritten += npart2read;

			}

			delete[] pos_tmp;
			temp_data.clear();
								
			// do all lower levels with standard cartesian grid
			for( int ilevel=gh.levelmax()-1; ilevel>=(int)gh.levelmin(); --ilevel )
				for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
					for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
						for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
							if( ! gh.is_refined(ilevel,i,j,k) )
							{
								double xx[3];
								gh.cell_pos(ilevel, i, j, k, xx);
								if( shift != NULL )
									xx[coord] += shift[coord];
								
								xx[coord] = fmod( (xx[coord]+(*gh.get_grid(ilevel))(i,j,k))*xfac + header_.BoxSize, header_.BoxSize );
								
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
			throw std::runtime_error("Internal consistency error while writing temporary file for positions");

		//... dump to temporary file
		ofs_temp.write( (char *)&blksize, sizeof(unsigned long long) );
		
		if( ofs_temp.bad() )
			throw std::runtime_error("I/O error while writing temporary file for positions");
		
		ofs_temp.close();
				
		if( shift != NULL )
			delete[] shift;
		
	}
	
	void write_dm_velocity( int coord, const grid_hierarchy& gh )
	{
		//... count number of leaf cells ...//
		size_t npcoarse = 0, npfine = 0;
		
		npfine   = gh.count_leaf_cells(gh.levelmax(), gh.levelmax());
		if( bmultimass_ )
			npcoarse = gh.count_leaf_cells(gh.levelmin(), gh.levelmax()-1);
		
		header_.npart[1] = npfine;
		header_.npart[5] = npcoarse;
		header_.npartTotal[1] = npfine;
		header_.npartTotal[5] = npcoarse;
		header_.npartTotalHighWord[1] = 0;
		header_.npartTotalHighWord[5] = 0;
		
		//... collect displacements and convert to absolute coordinates with correct
		//... units
		std::vector<T_store> temp_data;
		temp_data.reserve( block_buf_size_ );
		
		float isqrta = 1.0f/sqrt(header_.time);
		float vfac = isqrta*header_.BoxSize;
		
		if( kpcunits_ )
			vfac /= 1000.0;
		
		unsigned npart = npfine+npcoarse;
		unsigned nwritten = 0;
		unsigned long long blksize;
		
		char temp_fname[256];
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_vel+coord );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
		
		
		if( !do_glass_ )
		{
			
			blksize = sizeof(T_store)*npart;
			ofs_temp.write( (char *)&blksize, sizeof(unsigned long long) );
			
			for( int ilevel=levelmax_; ilevel>=(int)levelmin_; --ilevel )
				for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
					for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
						for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
							if( ! gh.is_refined(ilevel,i,j,k) )
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
		}
		else
		{
			
			std::ifstream ofg( fname_glass_cdm_.c_str(), std::ios::binary );
			
			if( !ofg.good() )
				LOGERR("could not open glass input file \'%s\'",fname_glass_cdm_.c_str());
			
			io_header glasshead;
			unsigned blksz;
			
			ofg.read( reinterpret_cast<char*>(&blksz), sizeof(unsigned) );
			assert( blksz == sizeof(io_header) );
			
			ofg.read( reinterpret_cast<char*>(&glasshead), sizeof( io_header ) );
			
			ofg.read( reinterpret_cast<char*>(&blksz), sizeof(unsigned) );
			ofg.read( reinterpret_cast<char*>(&blksz), sizeof(unsigned) );
			//assert( blksz == glasshead.npart[1]*sizeof(float)*3 );
			
			header_.npart[1] = glasshead.npart[1];
			header_.npartTotal[1] = glasshead.npartTotal[1];
			header_.npartTotalHighWord[1] = glasshead.npartTotalHighWord[1];
			float lglass = glasshead.BoxSize;
			
			// read glass, do interpolation and write
			size_t npartdone=0;
			size_t npinter = glasshead.npart[1];
			
			blksize = sizeof(T_store)*npinter;
			ofs_temp.write( (char *)&blksize, sizeof(unsigned long long) );
			
			float *pos_tmp = new float[3*block_buf_size_];
			temp_data.assign( block_buf_size_, 0.0 );
			
			while( npartdone < npinter )
			{
				size_t npart2read = std::min(npinter-npartdone,block_buf_size_);
				
				ofg.read( reinterpret_cast<char*>(&pos_tmp[0]), npart2read*sizeof(float)*3 );
				get_cic_velocity( pos_tmp, npart2read, lglass,  gh, &temp_data[0] );
				ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*npart2read );
				
				npartdone += npart2read;
				nwritten += npart2read;
				
			}
			
			delete[] pos_tmp;
			temp_data.clear();
			
			for( int ilevel=levelmax_-1; ilevel>=(int)levelmin_; --ilevel )
				for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
					for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
						for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
							if( ! gh.is_refined(ilevel,i,j,k) )
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
	{ }
	
	void write_gas_potential( const grid_hierarchy& gh )
	{ }
	
	
	
	//... write data for gas -- don't do this
	void write_gas_velocity( int coord, const grid_hierarchy& gh )
	{	
		//... count number of leaf cells ...//
		size_t npcoarse = 0, npfine = 0;
		
		npfine   = gh.count_leaf_cells(gh.levelmax(), gh.levelmax());
		
		header_.npart[2] = npfine;
		header_.npartTotal[2] = (unsigned)npfine;
		header_.npartTotalHighWord[2] = (unsigned)(npfine>>32);
		
		//... collect displacements and convert to absolute coordinates with correct
		//... units
		std::vector<T_store> temp_data;
		temp_data.reserve( block_buf_size_ );
		
		float isqrta = 1.0f/sqrt(header_.time);
		float vfac = isqrta*header_.BoxSize;
		
		if( kpcunits_ )
			vfac /= 1000.0;
		
		unsigned npart = npfine+npcoarse;
		unsigned nwritten = 0;
		
		char temp_fname[256];
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_gas_vel+coord );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );

		unsigned long long blksize;
		
		
		if(!do_glass_)
		{
			blksize = sizeof(T_store)*npart;
			ofs_temp.write( (char *)&blksize, sizeof(unsigned long long) );
			
			
			const unsigned ilevel = gh.levelmax();
			const unsigned 
			nx = gh.get_grid(ilevel)->size(0),
			ny = gh.get_grid(ilevel)->size(1),
			nz = gh.get_grid(ilevel)->size(2);
			
			for( unsigned i=0; i<nx; ++i )
				for( unsigned j=0; j<ny; ++j )
					for( unsigned k=0; k<nz; ++k )
					{	
						double v = (*gh.get_grid(ilevel))(i,j,k);
						
						if( temp_data.size() < block_buf_size_ )
							temp_data.push_back( v * vfac );
						else 
						{
							ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*block_buf_size_ );
							nwritten += block_buf_size_;
							temp_data.clear();
							temp_data.push_back( v * vfac );
						}
						
					}
		}else{
			
			std::ifstream ofg( fname_glass_baryon_.c_str(), std::ios::binary );
			
			if( !ofg.good() )
				LOGERR("could not open glass input file \'%s\'",fname_glass_cdm_.c_str());
			
			io_header glasshead;
			unsigned blksz;
			
			ofg.read( reinterpret_cast<char*>(&blksz), sizeof(unsigned) );
			assert( blksz == sizeof(io_header) );
			
			ofg.read( reinterpret_cast<char*>(&glasshead), sizeof( io_header ) );
			
			//size_t nreq = gh.size(gh.levelmax(), 0)*gh.size(gh.levelmax(), 1)*gh.size(gh.levelmax(), 2);
			/*if( nreq != (size_t)glasshead.npart[1] )
			{
				LOGERR("glass file contains %d particles, but should contain %ld",glasshead.npart[1],nreq);
				throw std::runtime_error("glass file does not contain the right amount of particles");
			}*/
			
			ofg.read( reinterpret_cast<char*>(&blksz), sizeof(unsigned) );
			ofg.read( reinterpret_cast<char*>(&blksz), sizeof(unsigned) );
			//assert( blksz == (glasshead.npart[1]+glasshead.npart[2])*sizeof(float)*3 );
			//ofg.seekg( sizeof(float)*3*glasshead.npart[1], std::ios_base::cur );
			
			// do the highest level with the glass
			float lglass = glasshead.BoxSize;
			
			header_.npart[2] = glasshead.npart[2];
			header_.npartTotal[2] = glasshead.npartTotal[2];
			header_.npartTotalHighWord[2] = glasshead.npartTotalHighWord[2];
			
			// read glass, do interpolation and write
			size_t npartdone=0;
			size_t npinter = (size_t)glasshead.npart[2];
			
			
			blksize = sizeof(T_store)*npinter;
			ofs_temp.write( (char *)&blksize, sizeof(unsigned long long) );
			
			float *pos_tmp = new float[3*block_buf_size_];
			temp_data.assign( block_buf_size_, 0.0 );
			
			while( npartdone < npinter )
			{
				size_t npart2read = std::min(npinter-npartdone,block_buf_size_);
				
				ofg.read( reinterpret_cast<char*>(&pos_tmp[0]), npart2read*sizeof(float)*3 );
				get_cic_velocity( pos_tmp, npart2read, lglass,  gh, &temp_data[0] );
				ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*npart2read );
				
				npartdone += npart2read;
				nwritten += npart2read;
			}
			
			delete[] pos_tmp;
			temp_data.clear();
			
			for( int ilevel=levelmax_-1; ilevel>=(int)levelmin_; --ilevel )
				for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
					for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
						for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
							if( ! gh.is_refined(ilevel,i,j,k) )
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
		unsigned long long npfine = 0;
		
		npfine   = gh.count_leaf_cells(gh.levelmax(), gh.levelmax());
		
		//... determine if we need to shift the coordinates back
		double *shift = NULL;
		
		if( cf_.getValueSafe<bool>("output","shift_back",false ) )
		{
			if( coord == 0 )
				std::cout << " - gadget2 output plug-in will shift particle positions back...\n";
			
			double h = 1.0/(1<<levelmin_);
			shift = new double[3];
			shift[0] = -(double)cf_.getValue<int>( "setup", "shift_x" )*h;
			shift[1] = -(double)cf_.getValue<int>( "setup", "shift_y" )*h;
			shift[2] = -(double)cf_.getValue<int>( "setup", "shift_z" )*h;
		}
		
		unsigned long long npart = npfine;
		unsigned long long nwritten = 0;
		
		//...
		header_.npart[2] = npfine;
		header_.npartTotal[2] = (unsigned)npfine;
		header_.npartTotalHighWord[2] = (unsigned)(npfine>>32);
		
		//header_.num_files = (int)ceil((double)npart/(double)npartmax_);
		
		//... collect displacements and convert to absolute coordinates with correct
		//... units
		std::vector<T_store> temp_data;
		temp_data.reserve( block_buf_size_ );
		
		
		char temp_fname[256];
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_gas_pos+coord );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
		
		unsigned long long blksize;
		
		double xfac = header_.BoxSize;
        
        //... shift particle positions (this has to be done as the same shift
        //... is used when computing the convolution kernel for SPH baryons)
        if( do_baryons_ )
        {
            double h = 1. / (1<<gh.levelmax());
            
            if( shift == NULL )
            {
                shift = new double[3];
                shift[0] = 0.0;
                shift[1] = 0.0;
                shift[2] = 0.0;
                
            }
            shift[0] += 0.5 * h * omegac_/omegam_;
            shift[1] += 0.5 * h * omegac_/omegam_;
            shift[2] += 0.5 * h * omegac_/omegam_;
        }
		
		//... only do finest grid
		if( !do_glass_ )
		{
			
			blksize = sizeof(T_store)*npart;
			ofs_temp.write( (char *)&blksize, sizeof(unsigned long long) );
			
			const unsigned ilevel = gh.levelmax();
			const unsigned 
				nx = gh.get_grid(ilevel)->size(0),
				ny = gh.get_grid(ilevel)->size(1),
				nz = gh.get_grid(ilevel)->size(2);
			
			for( unsigned i=0; i<nx; ++i )
				for( unsigned j=0; j<ny; ++j )
					for( unsigned k=0; k<nz; ++k )
					{	
						double xx[3];
						gh.cell_pos(ilevel, i, j, k, xx);
						if( shift != NULL )
							xx[coord] += shift[coord];
						
						double v = (*gh.get_grid(ilevel))(i,j,k);
						
						xx[coord] = fmod( (xx[coord]+v)*xfac + header_.BoxSize, header_.BoxSize );
						
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
		}else{
			
			std::ifstream ofg( fname_glass_baryon_.c_str(), std::ios::binary );
			
			if( !ofg.good() )
				LOGERR("could not open glass input file \'%s\'",fname_glass_cdm_.c_str());
			
			io_header glasshead;
			unsigned blksz;
			
			ofg.read( reinterpret_cast<char*>(&blksz), sizeof(unsigned) );
			assert( blksz == sizeof(io_header) );
			
			ofg.read( reinterpret_cast<char*>(&glasshead), sizeof( io_header ) );
			
			//size_t nreq = gh.size(gh.levelmax(), 0)*gh.size(gh.levelmax(), 1)*gh.size(gh.levelmax(), 2);
			/*if( nreq != (size_t)glasshead.npart[1] )
			{
				LOGERR("glass file contains %d particles, but should contain %ld",glasshead.npart[1],nreq);
				throw std::runtime_error("glass file does not contain the right amount of particles");
			}*/
			
			ofg.read( reinterpret_cast<char*>(&blksz), sizeof(unsigned) );
			ofg.read( reinterpret_cast<char*>(&blksz), sizeof(unsigned) );
			//assert( blksz == (glasshead.npart[1]+glasshead.npart[2])*sizeof(float)*3 );
			//ofg.seekg( sizeof(float)*3*glasshead.npart[1], std::ios_base::cur );
			
			float lglass = glasshead.BoxSize;
			
			
			blksize = sizeof(T_store)*glasshead.npart[1];
			ofs_temp.write( (char *)&blksize, sizeof(unsigned long long) );
			
			header_.npart[2] = glasshead.npart[2];
			header_.npartTotal[2] = glasshead.npartTotal[2];
			header_.npartTotalHighWord[2] = glasshead.npartTotalHighWord[2];
			
			double rhoc = 27.7519737;
			if( kpcunits_ )
				rhoc *= 10.0; // in h^2 M_sol / kpc^3
			
			header_.mass[2] = omegab_ * rhoc * pow(header_.BoxSize,3.)/(glasshead.npart[2]);
			
			// read glass, do interpolation and write
			size_t npartdone=0;
			size_t npinter = glasshead.npart[2];
			
			blksize = sizeof(T_store)*npinter;
			ofs_temp.write( (char *)&blksize, sizeof(unsigned long long) );
			
			float *pos_tmp = new float[3*block_buf_size_];
			temp_data.assign( block_buf_size_, 0.0 );
			
			while( npartdone < npinter )
			{
				size_t npart2read = std::min(npinter-npartdone,block_buf_size_);
				
				ofg.read( reinterpret_cast<char*>(&pos_tmp[0]), npart2read*sizeof(float)*3 );
				get_cic_displacement( coord, pos_tmp, npart2read, lglass,  gh, &temp_data[0] );
				ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*npart2read );
				
				npartdone += npart2read;
				nwritten += npart2read;
				
			}
			
			delete[] pos_tmp;
			temp_data.clear();		
		}
		
		if( temp_data.size() > 0 )
		{	
			ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*temp_data.size() );
			nwritten += temp_data.size();
		}
		
		if( nwritten != npart )
			throw std::runtime_error("Internal consistency error while writing temporary file for gas positions");
		
		//... dump to temporary file
		ofs_temp.write( (char *)&blksize, sizeof(unsigned long long) );
		
		if( ofs_temp.bad() )
			throw std::runtime_error("I/O error while writing temporary file for gas positions");
		
		ofs_temp.close();
		
		if( shift != NULL )
			delete[] shift;
	}
	
	void write_gas_density( const grid_hierarchy& gh )
	{	
		double rhoc = 27.7519737; // h^2 1e10 M_sol / Mpc^3
		
		if( kpcunits_ )
			rhoc *= 10.0; // in h^2 M_sol / kpc^3
		
		if( do_baryons_ && !do_glass_ )
			header_.mass[2] = omegab_ * rhoc * pow(header_.BoxSize,3.)/pow(2,3*levelmax_);
	}
	
	void finalize( void )
	{	
		this->assemble_gadget_file();
	}
};



namespace{
	output_plugin_creator_concrete< gadget2_2comp_output_plugin<float> > creator1("gadget2_2c");
#ifndef SINGLE_PRECISION
	output_plugin_creator_concrete< gadget2_2comp_output_plugin<double> > creator2("gadget2_2c_double");
#endif
}

