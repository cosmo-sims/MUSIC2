//
//  output_tipsy.cc
//  MUSIC
//
//  Created by Oliver Hahn on 4/6/11.
//  Copyright 2011 KIPAC/SLAC. All rights reserved.
//


#include <stdio.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <fstream>

#include "output.hh"


template< typename T_store=float >
class tipsy_output_plugin : public output_plugin
{
protected:
	
	std::ofstream ofs_;
	
	typedef T_store Real;
	
	struct gas_particle {
		Real mass;
		Real pos[3];
		Real vel[3];
		Real rho;
		Real temp;
		Real hsmooth;
		Real metals ;
		Real phi ;
	};
	
	struct gas_particle *gas_particles;
	
	struct dark_particle {
		Real mass;
		Real pos[3];
		Real vel[3];
		Real eps;
		Real phi ;
	};
	
	struct dark_particle *dark_particles;
	
	struct star_particle {
		Real mass;
		Real pos[3];
		Real vel[3];
		Real metals ;
		Real tform ;
		Real eps;
		Real phi ;
	};
	
	struct star_particle *star_particles;
	
	struct dump {
		double time ;
		int nbodies ;
		int ndim ;
		int nsph ;
		int ndark ;
		int nstar ;
	};
	
	enum iofields {
		id_dm_mass, id_dm_vel, id_dm_pos, id_gas_vel, id_gas_rho, id_gas_temp, id_gas_pos
	};
	
	dump header_;
	FILE *fp_;
	unsigned block_buf_size_;
	size_t npartmax_;
	bool bmorethan2bnd_;
	double epsfac_;
	double boxsize_;
	double astart_;
	double omegam_;
    double H0_;
	double BoxSize_;

	class pistream : public std::ifstream
	{
	public:
		pistream (std::string fname, size_t npart )
		: std::ifstream( fname.c_str(), std::ios::binary )
		{
			size_t blk;
			
			if( !this->good() )
			{	
				LOGERR("Could not open buffer file in gadget2 output plug-in");
				throw std::runtime_error("Could not open buffer file in gadget2 output plug-in");
			}
			
			this->read( (char*)&blk, sizeof(size_t) );
			
			if( blk != (size_t)(npart*sizeof(T_store)) )
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
			size_t blk;
			
			if( !this->good() )
			{	
				LOGERR("Could not open buffer file \'%s\' in gadget2 output plug-in",fname.c_str());
				throw std::runtime_error("Could not open buffer file in gadget2 output plug-in");
			}
			
			this->read( (char*)&blk, sizeof(size_t) );
			
			if( blk != (size_t)(npart*sizeof(T_store)) )
			{	
				LOGERR("Internal consistency error in gadget2 output plug-in");
				LOGERR("Expected %d bytes in temp file but found %d",npart*(unsigned)sizeof(T_store),blk);
				throw std::runtime_error("Internal consistency error in gadget2 output plug-in");
			}
		}
	};
	
	int xdr_dump( XDR *xdrs, T_store *fp )
	{ return 0; }
	
	int convert_header_XDR( XDR *pxdrs, struct dump* ph )
	{
		int pad = 0;
		
		if (!xdr_double(pxdrs,&ph->time)) return 0;
		if (!xdr_int(pxdrs,&ph->nbodies)) return 0;
		if (!xdr_int(pxdrs,&ph->ndim)) return 0;
		if (!xdr_int(pxdrs,&ph->nsph)) return 0;
		if (!xdr_int(pxdrs,&ph->ndark)) return 0;
		if (!xdr_int(pxdrs,&ph->nstar)) return 0;
		if (!xdr_int(pxdrs,&pad)) return 0;
		return 1;
		
	}
	
	inline T_store mass2eps( T_store& m )
	{
		return pow(m/omegam_,0.333333333333)*epsfac_;
	}
	
	void assemble_tipsy_file( void )
	{
		
		fp_ = fopen( fname_.c_str(), "w+" );
					
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
		
		
		pistream ifs_x, ifs_y, ifs_z, ifs_vx, ifs_vy, ifs_vz, ifs_m;
		pistream ifs_bx, ifs_by, ifs_bz, ifs_bvx, ifs_bvy, ifs_bvz;
		
		
		const unsigned 
			nptot = header_.nbodies,
		//npgas = header_.nsph ,
			npcdm = header_.ndark ;
		
		unsigned
			npleft = nptot, 
			n2read = std::min((unsigned)block_buf_size_,npleft);
			
		std::cout << " - Writing " << nptot << " particles to tipsy file...\n";
		
		//bool bbaryons = npgas > 0;
		
		std::vector<T_store> adata3;
		adata3.reserve( 3*block_buf_size_ );
		T_store *tmp1, *tmp2, *tmp3, *tmp4, *tmp5, *tmp6, *tmp7;
		
		tmp1 = new T_store[block_buf_size_];
		tmp2 = new T_store[block_buf_size_];
		tmp3 = new T_store[block_buf_size_];
		tmp4 = new T_store[block_buf_size_];
		tmp5 = new T_store[block_buf_size_];
		tmp6 = new T_store[block_buf_size_];
		tmp7 = new T_store[block_buf_size_];
		
		ifs_x.open( fnx, npcdm );
		ifs_y.open( fny, npcdm );
		ifs_z.open( fnz, npcdm );
		ifs_vx.open( fnvx, npcdm );
		ifs_vy.open( fnvy, npcdm );
		ifs_vz.open( fnvz, npcdm );
		ifs_m.open( fnm, npcdm );
		
		T_store zero = (T_store)0.0;
		
		while( true )
		{
			
			//... write the header .......................................................
			XDR xdrs;
			xdrstdio_create(&xdrs, fp_, XDR_ENCODE);
			convert_header_XDR( &xdrs, &header_ );
			
			std::vector<T_store> dump_store ( 9*block_buf_size_, (T_store)0.0 );
			
			npleft = npcdm;
			n2read = std::min(block_buf_size_,npleft);
			while( n2read > 0 )
			{
				ifs_x.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
				ifs_y.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
				ifs_z.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
				ifs_vx.read( reinterpret_cast<char*>(&tmp4[0]), n2read*sizeof(T_store) );
				ifs_vy.read( reinterpret_cast<char*>(&tmp5[0]), n2read*sizeof(T_store) );
				ifs_vz.read( reinterpret_cast<char*>(&tmp6[0]), n2read*sizeof(T_store) );
				ifs_m.read( reinterpret_cast<char*>(&tmp7[0]), n2read*sizeof(T_store) );
				
				for( size_t i=0; i<n2read; ++i )
				{
					xdr_dump(&xdrs,&tmp7[i]); // mass
					xdr_dump(&xdrs,&tmp1[i]); // x
					xdr_dump(&xdrs,&tmp2[i]); // y
					xdr_dump(&xdrs,&tmp3[i]); // z
					xdr_dump(&xdrs,&tmp4[i]); // vx
					xdr_dump(&xdrs,&tmp5[i]); // vy
					xdr_dump(&xdrs,&tmp6[i]); // vz
					
					T_store eps = mass2eps( tmp7[i] );
					
					xdr_dump(&xdrs, &eps ); // epsilon
					xdr_dump(&xdrs, &zero ); //potential
					
					/*dump_store[9*i+0] = tmp7[i];
					dump_store[9*i+1] = tmp1[i];
					dump_store[9*i+2] = tmp2[i];
					dump_store[9*i+3] = tmp3[i];
					dump_store[9*i+4] = tmp4[i];
					dump_store[9*i+5] = tmp5[i];
					dump_store[9*i+6] = tmp6[i];
					dump_store[9*i+7] = mass2eps( tmp7[i] );
					dump_store[9*i+8] = zero;*/
				}
				
				
				
				npleft -= n2read;
				n2read = std::min( block_buf_size_,npleft );
			}
			
			break;
		}
		
		
		fclose( fp_ );
	}
	
	
	
public:
	
	tipsy_output_plugin( config_file& cf )
	: output_plugin( cf ), ofs_( fname_.c_str(), std::ios::binary|std::ios::trunc )	
	{
		block_buf_size_ = cf_.getValueSafe<unsigned>("output","tipsy_blksize",1048576);
		
		//... ensure that everyone knows we want to do SPH
		cf.insertValue("setup","do_SPH","yes");
		
		//bbndparticles_  = !cf_.getValueSafe<bool>("output","gadget_nobndpart",false);
		npartmax_ = 1<<30;
		
		if(!ofs_.good())
		{	
			LOGERR("tipsy output plug-in could not open output file \'%s\' for writing!",fname_.c_str());
			throw std::runtime_error(std::string("tipsy output plug-in could not open output file \'")+fname_+"\' for writing!\n");
		}
		ofs_.close();
		
		double zstart = cf.getValue<double>("setup","zstart");
		astart_ = 1.0/(1.0+zstart);
		omegam_  = cf.getValue<double>("cosmology","Omega_m");
		boxsize_ = cf.getValue<double>("setup","boxlength");
		epsfac_ = cf.getValueSafe<double>("output","tipsy_eps",0.05);
        H0_ = cf.getValue<double>("cosmology","H0");

	}
	
	void write_dm_mass( const grid_hierarchy& gh )
	{
		//.. store header data
		header_.nbodies = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());
		header_.nsph	= 0;
		header_.ndark	= header_.nbodies;
		header_.nstar	= 0;
		header_.ndim	= 3;
		header_.time	= astart_;
		
		//... write data
		size_t nptot = header_.nbodies;
		
		std::vector<T_store> temp_dat;
		temp_dat.reserve(block_buf_size_);
		
		char temp_fname[256];
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_mass );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
		
		
		size_t blksize = sizeof(T_store)*nptot;
		ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
		size_t nwritten = 0;
		for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
		{
		  double pmass = omegam_/(1ul<<(3*ilevel));
			
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
		
		if( nwritten != nptot )
			throw std::runtime_error("Internal consistency error while writing temporary file for masses");
		
		ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
		if( ofs_temp.bad() )
			throw std::runtime_error("I/O error while writing temporary file for masses");
		
	}
	
	
	void write_dm_position( int coord, const grid_hierarchy& gh )
	{
		size_t nptot = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());
		
		std::vector<T_store> temp_data;
		temp_data.reserve( block_buf_size_ );
		
		
		char temp_fname[256];
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_pos+coord );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
		
		size_t blksize = sizeof(T_store)*nptot;
		ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
		size_t nwritten = 0;
		for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
			for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
				for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
					for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
						if( ! gh.is_refined(ilevel,i,j,k) )
						{
							double xx[3];
							gh.cell_pos(ilevel, i, j, k, xx);
							
							//xx[coord] = fmod( (xx[coord]+(*gh.get_grid(ilevel))(i,j,k)) + 1.0, 1.0 ) - 0.5;
							xx[coord] = (xx[coord]+(*gh.get_grid(ilevel))(i,j,k)) - 0.5;
							
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
		
		if( nwritten != nptot )
			throw std::runtime_error("Internal consistency error while writing temporary file for positions");
		
		//... dump to temporary file
		ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
		if( ofs_temp.bad() )
			throw std::runtime_error("I/O error while writing temporary file for positions");
		
		ofs_temp.close();
	}
	
	void write_dm_velocity( int coord, const grid_hierarchy& gh )
	{
		size_t nptot = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());
		
		std::vector<T_store> temp_data;
		temp_data.reserve( block_buf_size_ );
		
		double vfac = 2.894405/(100.0 * astart_); 

		char temp_fname[256];
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_vel+coord );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
		
		size_t blksize = sizeof(T_store)*nptot;
		ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
		size_t nwritten = 0;
		for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
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
			ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*temp_data.size() );
			nwritten += temp_data.size();
		}
		
		if( nwritten != nptot )
			throw std::runtime_error("Internal consistency error while writing temporary file for positions");
		
		//... dump to temporary file
		ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
		if( ofs_temp.bad() )
			throw std::runtime_error("I/O error while writing temporary file for positions");
		
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
	{ }
	
	
	//... write only for fine level
	void write_gas_position( int coord, const grid_hierarchy& gh )
	{ }
	
	void write_gas_density( const grid_hierarchy& gh )
	{ }
	
	void finalize( void )
	{	
		this->assemble_tipsy_file();
	}
};
	
template<>
int tipsy_output_plugin<float>::xdr_dump( XDR *xdrs, float*p )
{
	return xdr_float(xdrs,p);
}
	
template<>
int tipsy_output_plugin<double>::xdr_dump( XDR *xdrs, double*p )
{
	return xdr_double(xdrs,p);
}
	

namespace{
	output_plugin_creator_concrete< tipsy_output_plugin<float> > creator1("tipsy");
#ifndef SINGLE_PRECISION
	output_plugin_creator_concrete< tipsy_output_plugin<double> > creator2("tipsy_double");
#endif
}
