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
class gadget2_output_plugin : public output_plugin
{
protected:
	
	std::ofstream ofs_;
	bool bmultimass_;
	
	
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
		char fill[60];                       
	}header;                       
	
	
	header header_;
	
	std::string fname;
	
	enum iofields {
		id_dm_mass, id_dm_vel, id_dm_pos, id_gas_vel, id_gas_rho, id_gas_temp, id_gas_pos
	};
	
	unsigned block_buf_size_;
	unsigned long long npartmax_;
	unsigned nfiles_;
	
	//bool bbndparticles_;
	bool bmorethan2bnd_;
	bool kpcunits_;
	double YHe_;
	
	std::ifstream& open_and_check( std::string ffname, size_t npart )
	{
		std::ifstream ifs( ffname.c_str(), std::ios::binary );
		long long blk;
		ifs.read( (char*)&blk, sizeof(long long) );
		if( blk != npart*(long long)sizeof(T_store) )
		{	
			LOGERR("Internal consistency error in gadget2 output plug-in");
			LOGERR("Expected %d bytes in temp file but found %d",npart*(unsigned)sizeof(T_store),blk);
			throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
		}
		
		return ifs;
	}
	
	class pistream : public std::ifstream
	{
	public:
		pistream (std::string fname, size_t npart )
		: std::ifstream( fname.c_str(), std::ios::binary )
		{
			long long blk;
			
			if( !this->good() )
			{	
				LOGERR("Could not open buffer file in gadget2 output plug-in");
				throw std::runtime_error("Could not open buffer file in gadget2 output plug-in");
			}
			
			this->read( (char*)&blk, sizeof(long long) );
			
			if( blk != (long long)(npart*sizeof(T_store)) )
			{	
				LOGERR("Internal consistency error in gadget2 output plug-in");
				LOGERR("Expected %d bytes in temp file but found %d",npart*(unsigned)sizeof(T_store),blk);
				throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
			}
		}
		
		pistream ()
		{
			
		}
		
		void open(std::string fname, size_t npart )
		{
			std::ifstream::open( fname.c_str(), std::ios::binary );
			long long blk;
			
			if( !this->good() )
			{	
				LOGERR("Could not open buffer file \'%s\' in gadget2 output plug-in",fname.c_str());
				throw std::runtime_error("Could not open buffer file in gadget2 output plug-in");
			}
			
			this->read( (char*)&blk, sizeof(long long) );
			
			if( blk != (long long)(npart*sizeof(T_store)) )
			{	
				LOGERR("Internal consistency error in gadget2 output plug-in");
				LOGERR("Expected %d bytes in temp file but found %d",npart*(unsigned)sizeof(T_store),blk);
				throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
			}
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
		
		const unsigned 
			nptot = header_.npart[0]+header_.npart[1]+header_.npart[5],
			npgas = header_.npart[0],
			npcdm = nptot-npgas;
			
		unsigned
			npleft = nptot, 
			n2read = std::min((unsigned)block_buf_size_,npleft);
		
		std::cout << " - Writing " << nptot << " particles to Gadget file...\n"
				  << "      type 0 : " << header_.npart[0] << "\n"
				  << "      type 1 : " << header_.npart[1] << "\n"
				  << "      type 5 : " << header_.npart[5] << "\n";
		
		bool bbaryons = header_.npart[0] > 0;
				
		std::vector<T_store> adata3;
		adata3.reserve( 3*block_buf_size_ );
		T_store *tmp1, *tmp2, *tmp3;
		
		tmp1 = new T_store[block_buf_size_];
		tmp2 = new T_store[block_buf_size_];
		tmp3 = new T_store[block_buf_size_];
		
		//... for multi-file output
		//int fileno = 0;
		//size_t npart_left = nptot;
		
		while( true )
		{
			int blksize = sizeof(header);
			
			//... write the header .......................................................
			
			header this_header( header_ );
			
			ofs_.write( (char *)&blksize, sizeof(int) );
			ofs_.write( (char *)&this_header, sizeof(header) );
			ofs_.write( (char *)&blksize, sizeof(int) );
			
			
			//... particle positions ..................................................
			blksize = 3*nptot*sizeof(T_store);
			ofs_.write( (char *)&blksize, sizeof(int) );
			
			if( bbaryons )
			{
				
				
				iffs1.open( fnbx, npgas );
				iffs2.open( fnby, npgas );
				iffs3.open( fnbz, npgas );
				
				npleft = npgas;
				n2read = std::min(block_buf_size_,npleft);
				while( n2read > 0 )
				{
					iffs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
					iffs2.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
					iffs3.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
					
					for( unsigned i=0; i<n2read; ++i )
					{
						adata3.push_back( tmp1[i] );
						adata3.push_back( tmp2[i] );
						adata3.push_back( tmp3[i] );
					}
					ofs_.write( reinterpret_cast<char*>(&adata3[0]), 3*n2read*sizeof(T_store) );
					
					adata3.clear();
					npleft -= n2read;
					n2read = std::min( block_buf_size_,npleft );
				}
				iffs1.close();
				iffs2.close();
				iffs3.close();
				remove( fnbx );
				remove( fnby );
				remove( fnbz );
			}
			
			npleft = npcdm;
			n2read = std::min(block_buf_size_,npleft);
			
			iffs1.open( fnx, npcdm );
			iffs2.open( fny, npcdm );
			iffs3.open( fnz, npcdm );
			
			while( n2read > 0 )
			{
				iffs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
				iffs2.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
				iffs3.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
				
				for( unsigned i=0; i<n2read; ++i )
				{
					adata3.push_back( tmp1[i] );
					adata3.push_back( tmp2[i] );
					adata3.push_back( tmp3[i] );
				}
				ofs_.write( reinterpret_cast<char*>(&adata3[0]), 3*n2read*sizeof(T_store) );
				
				adata3.clear();
				npleft -= n2read;
				n2read = std::min( block_buf_size_,npleft );
			}
			ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
			
			iffs1.close();
			iffs2.close();
			iffs3.close();
			remove( fnx );
			remove( fny );
			remove( fnz );
			
			//... particle velocities ..................................................
			blksize = 3*nptot*sizeof(T_store);
			ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
			
			
			if( bbaryons )
			{
				iffs1.open( fnbvx, npgas );
				iffs2.open( fnbvy, npgas );
				iffs3.open( fnbvz, npgas );
				
				npleft = npgas;
				n2read = std::min(block_buf_size_,npleft);
				while( n2read > 0 )
				{
					iffs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
					iffs2.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
					iffs3.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
					
					for( unsigned i=0; i<n2read; ++i )
					{
						adata3.push_back( tmp1[i] );
						adata3.push_back( tmp2[i] );
						adata3.push_back( tmp3[i] );
					}
					
					ofs_.write( reinterpret_cast<char*>(&adata3[0]), 3*n2read*sizeof(T_store) );
					
					adata3.clear();
					npleft -= n2read;
					n2read = std::min( block_buf_size_,npleft );
				}
				
				iffs1.close();
				iffs2.close();
				iffs3.close();
				remove( fnbvx );
				remove( fnbvy );
				remove( fnbvz );
				
			}
			
			iffs1.open( fnvx, npcdm );
			iffs2.open( fnvy, npcdm );
			iffs3.open( fnvz, npcdm );
			
			npleft = npcdm;
			n2read = std::min(block_buf_size_,npleft);
			while( n2read > 0 )
			{
				iffs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
				iffs2.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
				iffs3.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
				
				for( unsigned i=0; i<n2read; ++i )
				{
					adata3.push_back( tmp1[i] );
					adata3.push_back( tmp2[i] );
					adata3.push_back( tmp3[i] );
				}
				
				ofs_.write( reinterpret_cast<char*>(&adata3[0]), 3*n2read*sizeof(T_store) );
				
				adata3.clear();
				npleft -= n2read;
				n2read = std::min( block_buf_size_,npleft );
			}
			ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
			
			iffs1.close();
			iffs2.close();
			iffs3.close();
			remove( fnvx );
			remove( fnvy );
			remove( fnvz );
			
			
			delete[] tmp2;
			delete[] tmp3;
			
			
			//... particle IDs ..........................................................
			std::vector<unsigned> ids(block_buf_size_,0);
			
			unsigned idcount = 0;
			npleft	= nptot;
			n2read	= std::min(block_buf_size_,npleft);
			blksize = sizeof(unsigned)*nptot;
			
			//... generate contiguous IDs and store in file .......................//
			ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
			while( n2read > 0 )
			{
				for( unsigned i=0; i<n2read; ++i )
					ids[i] = idcount++;
				ofs_.write( reinterpret_cast<char*>(&ids[0]), n2read*sizeof(unsigned) );
				ids.clear();
				npleft -= n2read;
				n2read = std::min( block_buf_size_,npleft );
			}
			ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
			
			std::vector<unsigned>().swap( ids );
			
			
			//... particle masses .......................................................
			if( bmultimass_ && bmorethan2bnd_ )
			{
				unsigned npcoarse = header_.npart[5];
				
				iffs1.open( fnm, npcoarse );
				
				npleft  = npcoarse;
				n2read  = std::min(block_buf_size_,npleft);
				blksize = npcoarse*sizeof(T_store);
				
				ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
				while( n2read > 0 )
				{
					iffs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
					ofs_.write( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
					
					npleft -= n2read;
					n2read = std::min( block_buf_size_,npleft );

				}
				ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
				
				iffs1.close();
				remove( fnm );
			}
			
			//... initial internal energy for gas particles
			if( bbaryons )
			{
				
				std::vector<T_store> eint(block_buf_size_,0.0);
				
				const double astart = 1./(1.+header_.redshift);
				const double npol  = (fabs(1.0-gamma_)>1e-7)? 1.0/(gamma_-1.) : 1.0;
				const double unitv = 1e5;
				const double h2    = header_.HubbleParam*header_.HubbleParam*0.0001;
				const double adec  = 1.0/(160.*pow(omegab_*h2/0.022,2.0/5.0));
				const double Tcmb0 = 2.726;
				const double Tini  = astart<adec? Tcmb0/astart : Tcmb0/astart/astart*adec;
				const double mu    = (Tini>1.e4) ? 4.0/(8.-5.*YHe_) : 4.0/(1.+3.*(1.-YHe_));
				const double ceint = 1.3806e-16/1.6726e-24 * Tini * npol / mu / unitv / unitv;
				
				npleft	= npgas;
				n2read	= std::min(block_buf_size_,npleft);
				blksize = sizeof(T_store)*npgas;
				
				ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
				while( n2read > 0 )
				{
					for( unsigned i=0; i<n2read; ++i )
						eint[i] = ceint;
					ofs_.write( reinterpret_cast<char*>(&eint[0]), n2read*sizeof(T_store) );
					ids.clear();
					npleft -= n2read;
					n2read = std::min( block_buf_size_,npleft );
				}
				ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
				
				LOGINFO("Gadget2 : set initial gas temperature to %.2f K/mu",Tini/mu);
			}
			
			
			delete[] tmp1;
			ofs_.flush();
			
			break;
		}
		
	}
	
	
public:
	
	bool do_baryons_;
	double omegab_;
	double gamma_;
	
	gadget2_output_plugin( config_file& cf )
	: output_plugin( cf ), ofs_( fname_.c_str(), std::ios::binary|std::ios::trunc )	
	{
		block_buf_size_ = cf_.getValueSafe<unsigned>("output","gadget_blksize",1048576);
		
		//... ensure that everyone knows we want to do SPH
		cf.insertValue("setup","do_SPH","yes");
		
		//bbndparticles_  = !cf_.getValueSafe<bool>("output","gadget_nobndpart",false);
		npartmax_ = 1<<30;
		
		if(!ofs_.good())
		{	
			LOGERR("gadget-2 output plug-in could not open output file \'%s\' for writing!",fname_.c_str());
			throw std::runtime_error(std::string("gadget-2 output plug-in could not open output file \'")+fname_+"\' for writing!\n");
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
		
		//... set time ......................................................
		header_.redshift = cf.getValue<double>("setup","zstart");
		header_.time = 1.0/(1.0+header_.redshift);
		
		//... SF flags
		header_.flag_sfr = 0;
		header_.flag_feedback = 0;
		header_.flag_cooling = 0;
		
		//... 
		header_.num_files = 1;
		header_.BoxSize = cf.getValue<double>("setup","boxlength");
		header_.Omega0 = cf.getValue<double>("cosmology","Omega_m");
		header_.OmegaLambda = cf.getValue<double>("cosmology","Omega_L");
		header_.HubbleParam = cf.getValue<double>("cosmology","H0");
		
		header_.flag_stellarage = 0;
		header_.flag_metals = 0;
		
		
		header_.flag_entropy_instead_u = 0;
		
		if( kpcunits_ )
			header_.BoxSize *= 1000.0;
	}
	
	
	void write_dm_mass( const grid_hierarchy& gh )
	{
		double rhoc = 27.7519737; // in h^2 1e10 M_sol / Mpc^3
		
		if( kpcunits_ )
			rhoc *= 10.0; // in h^2 M_sol / kpc^3
		
		
		if( !do_baryons_ )
			header_.mass[1] = header_.Omega0 * rhoc * pow(header_.BoxSize,3.)/pow(2,3*levelmax_);
		else
			header_.mass[1] = (header_.Omega0-omegab_) * rhoc * pow(header_.BoxSize,3.)/pow(2,3*levelmax_);
				
		if( bmorethan2bnd_ )
		{
			unsigned long long npcoarse = gh.count_leaf_cells(gh.levelmin(), gh.levelmax()-1);
			unsigned long long nwritten = 0;
			
			std::vector<T_store> temp_dat;
			temp_dat.reserve(block_buf_size_);
			
			char temp_fname[256];
			sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_mass );
			std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
			
			long long blksize = sizeof(T_store)*npcoarse;
			
			ofs_temp.write( (char *)&blksize, sizeof(long long) );
			
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
			
			ofs_temp.write( (char *)&blksize, sizeof(long long) );
			
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
		
		unsigned long long npart = npfine+npcoarse;
		unsigned long long nwritten = 0;
		
		//...
		header_.npart[1] = npfine;
		header_.npart[5] = npcoarse;
		header_.npartTotal[1] = (unsigned)npfine;
		header_.npartTotal[5] = (unsigned)npcoarse;
		header_.npartTotalHighWord[1] = (unsigned)(npfine>>32);
		header_.npartTotalHighWord[5] = (unsigned)(npfine>>32);
		
		header_.num_files = (int)ceil((double)npart/(double)npartmax_);
		
		//... collect displacements and convert to absolute coordinates with correct
		//... units
		std::vector<T_store> temp_data;
		temp_data.reserve( block_buf_size_ );
		
		
		char temp_fname[256];
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_pos+coord );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
		
		long long blksize = sizeof(T_store)*npart;
		ofs_temp.write( (char *)&blksize, sizeof(long long) );
		
		double xfac = header_.BoxSize;
		
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
		
		if( temp_data.size() > 0 )
		{	
			ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*temp_data.size() );
			nwritten += temp_data.size();
		}
		
		if( nwritten != npart )
			throw std::runtime_error("Internal consistency error while writing temporary file for positions");
		
		//... dump to temporary file
		ofs_temp.write( (char *)&blksize, sizeof(long long) );
		
		if( ofs_temp.bad() )
			throw std::runtime_error("I/O error while writing temporary file for positions");
		
		ofs_temp.close();
				
		if( shift != NULL )
			delete[] shift;
		
	}
	
	void write_dm_velocity( int coord, const grid_hierarchy& gh )
	{
		//... count number of leaf cells ...//
		unsigned npcoarse = 0, npfine = 0;
		
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
		
		char temp_fname[256];
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_vel+coord );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
		
		long long blksize = sizeof(T_store)*npart;
		ofs_temp.write( (char *)&blksize, sizeof(long long) );
		
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
		unsigned npcoarse = 0, npfine = 0;
		
		npfine   = gh.count_leaf_cells(gh.levelmax(), gh.levelmax());
		
		header_.npart[0] = npfine;
		header_.npartTotal[0] = npfine;
		header_.npartTotalHighWord[0] = 0;
		
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

		long long blksize = sizeof(T_store)*npart;
		ofs_temp.write( (char *)&blksize, sizeof(long long) );
		
		{
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
		header_.npart[0] = npfine;
		header_.npartTotal[0] = (unsigned)npfine;
		header_.npartTotalHighWord[0] = (unsigned)(npfine>>32);
		
		header_.num_files = (int)ceil((double)npart/(double)npartmax_);
		
		//... collect displacements and convert to absolute coordinates with correct
		//... units
		std::vector<T_store> temp_data;
		temp_data.reserve( block_buf_size_ );
		
		
		char temp_fname[256];
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_gas_pos+coord );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
		
		long long blksize = sizeof(T_store)*npart;
		ofs_temp.write( (char *)&blksize, sizeof(long long) );
		
		double xfac = header_.BoxSize;
		
		//... only do finest grid
		{
			const unsigned ilevel = gh.levelmax();
			const double h = 1.0/(1<<ilevel);
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
						
						//... shift particle positions (this has to be done as the same shift
						//... is used when computing the convolution kernel for SPH baryons)
						xx[coord] += 0.5*h;

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
		}
		
		if( temp_data.size() > 0 )
		{	
			ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*temp_data.size() );
			nwritten += temp_data.size();
		}
		
		if( nwritten != npart )
			throw std::runtime_error("Internal consistency error while writing temporary file for gas positions");
		
		//... dump to temporary file
		ofs_temp.write( (char *)&blksize, sizeof(long long) );
		
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
		
		if( do_baryons_ )
			header_.mass[0] = omegab_ * rhoc * pow(header_.BoxSize,3.)/pow(2,3*levelmax_);

		//do nothing as we write out positions
		
		//std::cout << " - WARNING: Gadget-2 output plug-in does not support baryons yet!\n"
		//		  << "            Baryon data is not written to file!" << std::endl;
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

