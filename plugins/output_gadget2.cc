/*
 
 output_gadget2.cc - This file is part of MUSIC -
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

#include <fstream>
#include "output.hh"

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
		id_dm_mass, id_dm_vel, id_dm_pos, id_gas_vel, id_gas_rho, id_gas_temp
	};
	
	unsigned block_buf_size_;
	//bool bbndparticles_;
	bool bmorethan2bnd_;
	
	
	void assemble_scalar( unsigned nptot, std::string ifname )
	{
		T_store *tmp;
		
		std::ifstream 
		ifs( ifname.c_str(), std::ios::binary );
		
		int
		npleft = nptot, 
		n2read = std::min((int)block_buf_size_,npleft);
		
		tmp = new T_store[block_buf_size_];
		
		std::vector<T_store> adata;
		adata.reserve( block_buf_size_ );
		
		unsigned blk;
		ifs.read( (char *)&blk, sizeof(int) );
		if( blk != nptot*sizeof(T_store) ){
			delete[] tmp;
			throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
		}
		
		npleft  = nptot;
		n2read  = std::min((int)block_buf_size_,npleft);
		int blksize = nptot*sizeof(T_store);
		
		ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
		while( n2read > 0 )
		{
			ifs.read( reinterpret_cast<char*>(&tmp[0]), n2read*sizeof(T_store) );
			ofs_.write( reinterpret_cast<char*>(&tmp[0]), n2read*sizeof(T_store) );
			npleft -= n2read;
			n2read = std::min( (int)block_buf_size_,npleft );
		}
		ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
		
		remove( fname.c_str() );
		
		
	}
	
	void assemble_vector( unsigned nptot, std::string fname1, std::string fname2,std::string fname3 )
	{
		T_store *tmp1, *tmp2, *tmp3;

		std::ifstream 
			ifs1( fname1.c_str(), std::ios::binary ),	
			ifs2( fname2.c_str(), std::ios::binary ),	
			ifs3( fname3.c_str(), std::ios::binary );
		
		int
			npleft = nptot, 
			n2read = std::min((int)block_buf_size_,npleft);
		
		
		tmp1 = new T_store[block_buf_size_];
		tmp2 = new T_store[block_buf_size_];
		tmp3 = new T_store[block_buf_size_];
		
		std::vector<T_store> adata3;
		adata3.reserve( 3*block_buf_size_ );
		
		unsigned blk;
		ifs1.read( (char *)&blk, sizeof(int) );
		if( blk != nptot*sizeof(T_store) ){
			delete[] tmp1;
			delete[] tmp2;
			delete[] tmp3;
			throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
		}
			
		
		ifs2.read( (char *)&blk, sizeof(int) );
		if( blk != nptot*sizeof(T_store) )
		{	
			delete[] tmp1;
			delete[] tmp2;
			delete[] tmp3;
			throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
		}
		
		ifs3.read( (char *)&blk, sizeof(int) );
		if( blk != nptot*sizeof(T_store) )
		{	
			delete[] tmp1;
			delete[] tmp2;
			delete[] tmp3;
			throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
		}
		
		
		//int blksize = 3*nptot*sizeof(T_store);
		//ofs_.write( (char *)&blksize, sizeof(int) );
		
		while( n2read > 0 )
		{
			//ifs1.read( (char*)&tmp1[0], n2read*sizeof(T_store) );
			ifs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
			ifs2.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
			ifs3.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
			
			for( int i=0; i<n2read; ++i )
			{
				adata3.push_back( tmp1[i] );
				adata3.push_back( tmp2[i] );
				adata3.push_back( tmp3[i] );
			}
			ofs_.write( reinterpret_cast<char*>(&adata3[0]), 3*n2read*sizeof(T_store) );
			
			adata3.clear();
			npleft -= n2read;
			n2read = std::min( (int)block_buf_size_,npleft );
		}
		//ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
		
		remove( fname1.c_str() );
		remove( fname2.c_str() );
		remove( fname3.c_str() );
		
		
		delete[] tmp1;
		delete[] tmp2;
		delete[] tmp3;
		
	}
	
	
	void assemble_gadget_file( void )
	{
		int blksize = sizeof(header);
		
		//... write the header .......................................................
		ofs_.write( (char *)&blksize, sizeof(int) );
		ofs_.write( (char *)&header_, sizeof(header) );
		ofs_.write( (char *)&blksize, sizeof(int) );
		
		//............................................................................
		//... copy from the temporary files, interleave the data and save ............
		
		char fnx[256],fny[256],fnz[256],fnvx[256],fnvy[256],fnvz[256],fnm[256];
		
		sprintf( fnx,  "___ic_temp_%05d.bin", 100*id_dm_pos+0 );
		sprintf( fny,  "___ic_temp_%05d.bin", 100*id_dm_pos+1 );
		sprintf( fnz,  "___ic_temp_%05d.bin", 100*id_dm_pos+2 );
		sprintf( fnvx, "___ic_temp_%05d.bin", 100*id_dm_vel+0 );
		sprintf( fnvy, "___ic_temp_%05d.bin", 100*id_dm_vel+1 );
		sprintf( fnvz, "___ic_temp_%05d.bin", 100*id_dm_vel+2 );
		sprintf( fnm,  "___ic_temp_%05d.bin", 100*id_dm_mass  );
		
		std::ifstream 
		ifs1( fnx, std::ios::binary ),	
		ifs2( fny, std::ios::binary ),	
		ifs3( fnz, std::ios::binary );
		
		
		const int 
		nptot = header_.npart[1]+header_.npart[5];
		int
		npleft = nptot, 
		n2read = std::min((int)block_buf_size_,npleft);
		
		std::cout << " - Writing " << nptot << " particles to Gadget file...\n"
				  << "      type 1 : " << header_.npart[1] << "\n"
				  << "      type 5 : " << header_.npart[5] << "\n";
		
		
		//... particle coordinates ..................................................
		unsigned blk;
		ifs1.read( (char *)&blk, sizeof(int) );
		if( blk != nptot*sizeof(T_store) )
			throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
		
		ifs2.read( (char *)&blk, sizeof(int) );
		if( blk != nptot*sizeof(T_store) )
			throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
		
		ifs3.read( (char *)&blk, sizeof(int) );
		if( blk != nptot*sizeof(T_store) )
			throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
		
		std::vector<T_store> adata3;
		adata3.reserve( 3*block_buf_size_ );
		T_store *tmp1, *tmp2, *tmp3;
		
		tmp1 = new T_store[block_buf_size_];
		tmp2 = new T_store[block_buf_size_];
		tmp3 = new T_store[block_buf_size_];
		
		blksize = 3*nptot*sizeof(T_store);
		ofs_.write( (char *)&blksize, sizeof(int) );
		
		while( n2read > 0 )
		{
			//ifs1.read( (char*)&tmp1[0], n2read*sizeof(T_store) );
			ifs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
			ifs2.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
			ifs3.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
			
			for( int i=0; i<n2read; ++i )
			{
				adata3.push_back( tmp1[i] );
				adata3.push_back( tmp2[i] );
				adata3.push_back( tmp3[i] );
			}
			ofs_.write( reinterpret_cast<char*>(&adata3[0]), 3*n2read*sizeof(T_store) );
			
			adata3.clear();
			npleft -= n2read;
			n2read = std::min( (int)block_buf_size_,npleft );
		}
		ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
		
		remove( fnx );
		remove( fny );
		remove( fnz );
		
		//... particle velocities ..................................................
		ifs1.close(); ifs1.open( fnvx, std::ios::binary );
		ifs2.close(); ifs2.open( fnvy, std::ios::binary );
		ifs3.close(); ifs3.open( fnvz, std::ios::binary );
		
		ifs1.read( (char *)&blk, sizeof(int) );
		if( blk != nptot*sizeof(T_store)){
			delete[] tmp1;
			delete[] tmp2;
			delete[] tmp3;
			throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
		}
		
		ifs2.read( (char *)&blk, sizeof(int) );
		if( blk != nptot*sizeof(T_store)){
			delete[] tmp1;
			delete[] tmp2;
			delete[] tmp3;
			throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
		}
		
		ifs3.read( (char *)&blk, sizeof(int) );
		if( blk != nptot*sizeof(T_store)){
			delete[] tmp1;
			delete[] tmp2;
			delete[] tmp3;
			throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
		}
		
		npleft = nptot;
		n2read = std::min((int)block_buf_size_,npleft);
		
		
		ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
		
		while( n2read > 0 )
		{
			ifs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
			ifs2.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
			ifs3.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
			
			for( int i=0; i<n2read; ++i )
			{
				adata3.push_back( tmp1[i] );
				adata3.push_back( tmp2[i] );
				adata3.push_back( tmp3[i] );
			}
			
			ofs_.write( reinterpret_cast<char*>(&adata3[0]), 3*n2read*sizeof(T_store) );
			
			adata3.clear();
			npleft -= n2read;
			n2read = std::min( (int)block_buf_size_,npleft );
		}
		ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
		remove( fnvx );
		remove( fnvy );
		remove( fnvz );
		
		
		delete[] tmp2;
		delete[] tmp3;
		
		
		//... particle IDs ..........................................................
		std::vector<unsigned> ids(block_buf_size_,0);
		
		unsigned idcount = 0;
		npleft	= nptot;
		n2read	= std::min((int)block_buf_size_,npleft);
		blksize = sizeof(unsigned)*nptot;
		
		//... generate contiguous IDs and store in file .......................//
		ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
		while( n2read > 0 )
		{
			for( int i=0; i<n2read; ++i )
				ids[i] = idcount++;
			ofs_.write( reinterpret_cast<char*>(&ids[0]), n2read*sizeof(unsigned) );
			ids.clear();
			npleft -= n2read;
			n2read = std::min( (int)block_buf_size_,npleft );
		}
		ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
		
		
		//... particle masses .......................................................
		if( bmultimass_ && bmorethan2bnd_ )
		{
			unsigned npcoarse = header_.npart[5];
			
			ifs1.close();
			ifs1.open( fnm, std::ios::binary );
			
			if( !ifs1.good() )
				throw std::runtime_error("Could not open buffer file in gadget2 output plug-in");
			
			ifs1.read( (char *)&blk, sizeof(int) );
			if( blk != npcoarse*sizeof(T_store)){
				delete[] tmp1;
				throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
			}
			
			
			npleft  = npcoarse;
			n2read  = std::min((int)block_buf_size_,npleft);
			blksize = npcoarse*sizeof(T_store);
			
			ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
			while( n2read > 0 )
			{
				ifs1.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
				ofs_.write( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
				npleft -= n2read;
				n2read = std::min( (int)block_buf_size_,npleft );
			}
			ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
			
			remove( fnm );
		}
		
		delete[] tmp1;
		ofs_.flush();
	}
	
	
public:
	
	gadget2_output_plugin( config_file& cf )//std::string afname, Cosmology cosm, Parameters param, unsigned block_buf_size = 100000 )
	: output_plugin( cf ), ofs_( fname_.c_str(), std::ios::binary|std::ios::trunc )	
	{
		block_buf_size_ = cf_.getValueSafe<unsigned>("output","gadget_blksize",1048576);
		//block_buf_size_ = cf_.getValueSafe<unsigned>("output","gadget_blksize",100000);
		//bbndparticles_  = !cf_.getValueSafe<bool>("output","gadget_nobndpart",false);
		
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
	}
	
	
	void write_dm_mass( const grid_hierarchy& gh )
	{
		double rhoc = 27.7519737; // h^2 1e10 M_sol / Mpc^3
		
		header_.mass[1] = header_.Omega0 * rhoc * pow(header_.BoxSize,3.)/pow(2,3*levelmax_);
		
		if( bmorethan2bnd_ )
		{
			
			std::vector<T_store> temp_dat;
			temp_dat.clear();
						
			for( int ilevel=gh.levelmax()-1; ilevel>=(int)gh.levelmin(); --ilevel )
			{
				double pmass = header_.Omega0 * rhoc * pow(header_.BoxSize,3.)/pow(2,3*ilevel);			
				
				for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
					for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
						for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
							if( ! gh.is_refined(ilevel,i,j,k) )
								temp_dat.push_back( pmass );
			}
			
			char temp_fname[256];
			sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_mass );
			std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
			
			int blksize = sizeof(T_store)*temp_dat.size();
			ofs_temp.write( (char *)&blksize, sizeof(int) );
			ofs_temp.write( (char*)&temp_dat[0], blksize );		
			ofs_temp.write( (char *)&blksize, sizeof(int) );
			
			
		}
		else 
		{
			header_.mass[5] = header_.Omega0 * rhoc * pow(header_.BoxSize,3.)/pow(2,3*levelmin_);
		}
	}
	
	
	void write_dm_position( int coord, const grid_hierarchy& gh )
	{
		//... count number of leaf cells ...//
		unsigned npcoarse = 0, npfine = 0;
		
		npfine   = gh.count_leaf_cells(gh.levelmax(), gh.levelmax());
		if( bmultimass_ )
			npcoarse = gh.count_leaf_cells(gh.levelmin(), gh.levelmax()-1);
		
		
		//... determine if we need to shift the coordinates back
		double *shift = NULL;
		
		if( cf_.getValueSafe<bool>("output","shift_back",false ) )
		{
			if( coord == 0 )
				std::cout << " - gadget2 output plug-in will shift particle positions back...\n";
			
			double h = 1.0/pow(2,levelmin_);
			shift = new double[3];
			shift[0] = -(double)cf_.getValue<int>( "setup", "shift_x" )*h;
			shift[1] = -(double)cf_.getValue<int>( "setup", "shift_y" )*h;
			shift[2] = -(double)cf_.getValue<int>( "setup", "shift_z" )*h;
		}
		
		/*if( !cf_.getValueSafe<bool>("output","stagger_particles",false ) )
		{
			double h = 1.0/pow(2,levelmax_);
			if( shift==NULL )
			{
				shift = new double[3];
				shift[0] = -0.5*h;
				shift[1] = -0.5*h;
				shift[2] = -0.5*h;	
			}else{
				shift[0] -= 0.5*h;
				shift[1] -= 0.5*h;
				shift[2] -= 0.5*h;
			}
			
		}*/
		
		//...
		header_.npart[1] = npfine;
		header_.npart[5] = npcoarse;
		header_.npartTotal[1] = npfine;
		header_.npartTotal[5] = npcoarse;
		header_.npartTotalHighWord[1] = 0;
		header_.npartTotalHighWord[5] = 0;
		
		
		//... collect displacements and convert to absolute coordinates with correct
		//... units
		std::vector<T_store> temp_data;
		temp_data.reserve( npfine );
		
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
							temp_data.push_back( xx[coord] );
						}
		
		
		//... dump to temporary file
		char temp_fname[256];
		
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_pos+coord );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
		
		int blksize = sizeof(T_store)*temp_data.size();
		ofs_temp.write( (char *)&blksize, sizeof(int) );
		ofs_temp.write( (char*)&temp_data[0], blksize );
		ofs_temp.write( (char *)&blksize, sizeof(int) );
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
		temp_data.reserve( npfine+npcoarse );
		
		float isqrta = 1.0f/sqrt(header_.time);
		//float vfac = cosm_.astart/100.0*isqrta*param_.boxlength;
		float vfac = isqrta*header_.BoxSize;
		
		for( int ilevel=levelmax_; ilevel>=(int)levelmin_; --ilevel )
			for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
				for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
					for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
						if( ! gh.is_refined(ilevel,i,j,k) )
							temp_data.push_back( (*gh.get_grid(ilevel))(i,j,k) * vfac );
							
		//... dump to temporary file
		char temp_fname[256];
		
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_vel+coord );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
		
		int blksize = sizeof(T_store)*temp_data.size();
		ofs_temp.write( (char *)&blksize, sizeof(int) );
		ofs_temp.write( (char*)&temp_data[0], blksize );
		ofs_temp.write( (char *)&blksize, sizeof(int) );
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
		std::cout << " - WARNING: Gadget-2 output plug-in does not support baryons yet!\n"
				  << "            Baryon data is not written to file!" << std::endl;
	}
	
	void write_gas_position( int coord, const grid_hierarchy& gh )
	{	
		std::cout << " - WARNING: Gadget-2 output plug-in does not support baryons yet!\n"
				  << "            Baryon data is not written to file!" << std::endl;
	}
	
	void write_gas_density( const grid_hierarchy& gh )
	{	
		std::cout << " - WARNING: Gadget-2 output plug-in does not support baryons yet!\n"
				  << "            Baryon data is not written to file!" << std::endl;
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

