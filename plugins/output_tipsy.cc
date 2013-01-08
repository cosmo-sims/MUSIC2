//
//  output_tipsy.cc
//  MUSIC
//
//  Created by Oliver Hahn on 4/6/11.
//  Credits go to Kyle Stewart and Shea Garrison-Kimmel.
//  Copyright 2011 KIPAC/SLAC. All rights reserved.
//

#include <unistd.h>
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
	Real hsmooth; // force softening
	Real metals ;
	Real phi ;
    };
    
    struct gas_particle *gas_particles;
    
    struct dark_particle {
	Real mass;
	Real pos[3];
	Real vel[3];
	Real eps;     // force softening
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
	int pad ; // add the pad to make it 8-byte aligned
    };
    
    enum iofields {
	id_dm_mass, id_dm_vel, id_dm_pos, id_gas_vel, id_gas_rho, id_gas_temp, id_gas_pos, id_gas_mass
    };
    
    dump header_;
    FILE *fp_;
    unsigned block_buf_size_;
    size_t npartmax_;
    bool bmorethan2bnd_;
    bool bmultimass_;
    double epsfac_;
    double boxsize_;
    double astart_;
    double omegam_;
    double omegab_;
    double H0_;
    double YHe_;
    double gamma_;
    double BoxSize_;
    
    bool with_baryons_;
    size_t np_fine_gas_, np_fine_dm_, np_coarse_dm_;
    
    bool native_;
    
    class pistream : public std::ifstream
    {
    public:
	pistream (std::string fname, size_t npart )
	    : std::ifstream( fname.c_str(), std::ios::binary )
	{
	    size_t blk;
	    
	    if( !this->good() )
	    {	
		LOGERR("Could not open buffer file in TIPSY output plug-in");
		throw std::runtime_error("Could not open buffer file in TIPSY output plug-in");
	    }
	    
	    this->read( (char*)&blk, sizeof(size_t) );
	    
	    if( blk != (size_t)(npart*sizeof(T_store)) )
	    {	
		LOGERR("Internal consistency error in TIPSY output plug-in");
		LOGERR("Expected %d bytes in temp file but found %d",npart*(unsigned)sizeof(T_store),blk);
		throw std::runtime_error("Internal consistency error in TIPSY output plug-in");
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
		LOGERR("Could not open buffer file \'%s\' in TIPSY output plug-in",fname.c_str());
		throw std::runtime_error("Could not open buffer file in TIPSY output plug-in");
	    }
	    
	    this->read( (char*)&blk, sizeof(size_t) );
	    
	    if( blk != (size_t)(npart*sizeof(T_store)) )
	    {	
		LOGERR("Internal consistency error in TIPSY output plug-in");
		LOGERR("Expected %d bytes in temp file but found %d",npart*(unsigned)sizeof(T_store),blk);
		throw std::runtime_error("Internal consistency error in TIPSY output plug-in");
	    }
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
		LOGERR("Could not open buffer file in TIPSY output plug-in");
		throw std::runtime_error("Could not open buffer file in TIPSY output plug-in");
	    }
	    
            this->read( (char*)&blk, sizeof(size_t) );
	    
	    if( blk != npart*sizeof(T_store) )
	    {	
		LOGERR("Internal consistency error in TIPSY output plug-in");
		LOGERR("Expected %ld bytes in temp file but found %ld",npart*sizeof(T_store),blk);
		throw std::runtime_error("Internal consistency error in TIPSY output plug-in");
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
		LOGERR("Could not open buffer file \'%s\' in TIPSY output plug-in",fname.c_str());
		throw std::runtime_error("Could not open buffer file in TIPSY output plug-in");
	    }
	    
            this->read( (char*)&blk, sizeof(size_t) );
	    
	    if( blk != npart*sizeof(T_store) )
	    {	
		LOGERR("Internal consistency error in TIPSY output plug-in");
		LOGERR("Expected %ld bytes in temp file but found %ld",npart*sizeof(T_store),blk);
		throw std::runtime_error("Internal consistency error in TIPSY output plug-in");
	    }
	    
            this->seekg( offset, std::ios::cur );
            this->seekp( offset+sizeof(size_t), std::ios::beg );
	}
    };
    
    int xdr_dump( XDR *xdrs, T_store *fp )
    { return 0; }
    
    int convert_header_XDR( XDR *pxdrs, struct dump* ph )
    {
	if (!xdr_double(pxdrs,&ph->time)) return 0;
	if (!xdr_int(pxdrs,&ph->nbodies)) return 0;
	if (!xdr_int(pxdrs,&ph->ndim)) return 0;
	if (!xdr_int(pxdrs,&ph->nsph)) return 0;
	if (!xdr_int(pxdrs,&ph->ndark)) return 0;
	if (!xdr_int(pxdrs,&ph->nstar)) return 0;
	if (!xdr_int(pxdrs,&ph->pad)) return 0;
	return 1;
	
    }
    
    inline T_store mass2eps( T_store& m )
    {
	return pow(m/omegam_,0.333333333333)*epsfac_;
    }
    
    void combine_components_for_coarse( void )
    {        
        const size_t 
	    nptot = np_fine_dm_+np_coarse_dm_,
	    npfine = np_fine_dm_,
	    npcoarse = np_coarse_dm_;
	
        std::vector<T_store> tmp1, tmp2;
	
        tmp1.assign(block_buf_size_,0.0);
        tmp2.assign(block_buf_size_,0.0);
	
        double facb = omegab_/omegam_, facc = (omegam_-omegab_)/omegam_;
	
	
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
    
    void assemble_tipsy_file( void )
    {
	
	if( with_baryons_ && bmultimass_ )
	    combine_components_for_coarse();
	
	
	fp_ = fopen( fname_.c_str(), "w+" );
	
	//............................................................................
	//... copy from the temporary files, interleave the data and save ............
	
	char fnx[256],fny[256],fnz[256],fnvx[256],fnvy[256],fnvz[256],fnm[256];
	char fnbx[256], fnby[256], fnbz[256], fnbvx[256], fnbvy[256], fnbvz[256], fnbm[256];
	
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
	sprintf( fnbm,  "___ic_temp_%05d.bin", 100*id_gas_mass  );
	
	
	pistream ifs_x, ifs_y, ifs_z, ifs_vx, ifs_vy, ifs_vz, ifs_m;
	pistream ifs_bx, ifs_by, ifs_bz, ifs_bvx, ifs_bvy, ifs_bvz;
	
	
	const unsigned 
	    nptot = header_.nbodies,
            npgas = header_.nsph ,
	    npcdm = header_.ndark ;
	
	unsigned
	    npleft = nptot, 
	    n2read = std::min((unsigned)block_buf_size_,npleft);
	
	//std::cout << " - Writing " << nptot << " particles to tipsy file...\n";
	
	LOGINFO("TIPSY : output plugin will write:\n       DM particles   : %d\n       SPH particles  : %d",header_.ndark, header_.nsph);
	
	
	
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
	
	
	T_store zero = (T_store)0.0;
	
	while( true )
	{
            //... write the header .......................................................
	    XDR xdrs;
	    
	    if( native_ ) {

		fwrite( &header_, sizeof(dump), 1, fp_ );

	    } else {

		xdrstdio_create(&xdrs, fp_, XDR_ENCODE);
		convert_header_XDR( &xdrs, &header_ );
		
		std::vector<T_store> dump_store ( 9*block_buf_size_, (T_store)0.0 );
	    }

            //... sph particles ..................................................
            if( with_baryons_ )
            {
		
		LOGINFO("TIPSY : writing baryon data");
		
                // compute gas temperature
		
                const double astart = astart_;
		//const double npol  = (fabs(1.0-gamma_)>1e-7)? 1.0/(gamma_-1.) : 1.0;
		//const double unitv = 1e5;
		const double h2    = H0_ * H0_ * 0.0001;
		const double adec  = 1.0/(160.*pow(omegab_*h2/0.022,2.0/5.0));
		const double Tcmb0 = 2.726;
		const double Tini  = astart<adec? Tcmb0/astart : Tcmb0/astart/astart*adec;
                const double mu    = (Tini>1.e4) ? 4.0/(8.-5.*YHe_) : 4.0/(1.+3.*(1.-YHe_));
		//const double ceint = 1.3806e-16/1.6726e-24 * Tini * npol / mu / unitv / unitv;
		
                T_store temperature = (T_store) Tini;
                LOGINFO("TIPSY : set initial gas temperature to %.2f K (mu = %.2f)",Tini,mu);
		
		
                // write
                ifs_x.open( fnbx, npcdm );
                ifs_y.open( fnby, npcdm );
                ifs_z.open( fnbz, npcdm );
                ifs_vx.open( fnbvx, npcdm );
                ifs_vy.open( fnbvy, npcdm );
                ifs_vz.open( fnbvz, npcdm );
                ifs_m.open( fnbm, npgas );
		
                npleft = npgas;
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
		    
		    if( native_ ) {

			for( size_t i=0; i<n2read; ++i )
			{                        
			    fwrite( &tmp7[i], sizeof(tmp7[i]), 1, fp_); // mass
			    fwrite( &tmp1[i], sizeof(tmp1[i]), 1, fp_); // x
			    fwrite( &tmp2[i], sizeof(tmp2[i]), 1, fp_); // y
			    fwrite( &tmp3[i], sizeof(tmp3[i]), 1, fp_); // z
			    fwrite( &tmp4[i], sizeof(tmp4[i]), 1, fp_); // vx
			    fwrite( &tmp5[i], sizeof(tmp5[i]), 1, fp_); // vy
			    fwrite( &tmp6[i], sizeof(tmp6[i]), 1, fp_); // vz
			    fwrite( &zero, sizeof(zero), 1, fp_);       // rho
			    fwrite( &temperature, sizeof(temperature), 1, fp_); // temp
			    
			    T_store eps = mass2eps( tmp7[i] );
			    
			    fwrite( &eps, sizeof(eps), 1, fp_);
			    fwrite( &zero, sizeof(zero), 1, fp_);
			    fwrite( &zero, sizeof(zero), 1, fp_);
			}

		    } else {

			for( size_t i=0; i<n2read; ++i )
			{                        
			    xdr_dump(&xdrs,&tmp7[i]); // mass
			    xdr_dump(&xdrs,&tmp1[i]); // x
			    xdr_dump(&xdrs,&tmp2[i]); // y
			    xdr_dump(&xdrs,&tmp3[i]); // z
			    xdr_dump(&xdrs,&tmp4[i]); // vx
			    xdr_dump(&xdrs,&tmp5[i]); // vy
			    xdr_dump(&xdrs,&tmp6[i]); // vz
			    xdr_dump(&xdrs, &zero );  // rho
			    xdr_dump(&xdrs, &temperature );  // temp
			    
			    T_store eps = mass2eps( tmp7[i] );
			    
			    xdr_dump(&xdrs, &eps );   // epsilon / hsmooth
			    xdr_dump(&xdrs, &zero );  // metals
			    xdr_dump(&xdrs, &zero );  //potential
			}
		    
		    }
		    
		    
                    npleft -= n2read;
                    n2read = std::min( block_buf_size_,npleft );
                }
		
                ifs_x.close();
                ifs_y.close();
                ifs_z.close();
                ifs_vx.close();
                ifs_vy.close();
                ifs_vz.close();
                ifs_m.close();
		
            }
	    
	    
	    
	    //... dark matter particles ..................................................
	    LOGINFO("TIPSY : writing DM data");
	    
            ifs_x.open( fnx, npcdm );
            ifs_y.open( fny, npcdm );
            ifs_z.open( fnz, npcdm );
            ifs_vx.open( fnvx, npcdm );
            ifs_vy.open( fnvy, npcdm );
            ifs_vz.open( fnvz, npcdm );
            ifs_m.open( fnm, npcdm );
	    
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

		
		if( native_ ) {
		    
		    for( size_t i=0; i<n2read; ++i )
		    {                        
			fwrite( &tmp7[i], sizeof(tmp7[i]), 1, fp_); // mass
			fwrite( &tmp1[i], sizeof(tmp1[i]), 1, fp_); // x
			fwrite( &tmp2[i], sizeof(tmp2[i]), 1, fp_); // y
			fwrite( &tmp3[i], sizeof(tmp3[i]), 1, fp_); // z
			fwrite( &tmp4[i], sizeof(tmp4[i]), 1, fp_); // vx
			fwrite( &tmp5[i], sizeof(tmp5[i]), 1, fp_); // vy
			fwrite( &tmp6[i], sizeof(tmp6[i]), 1, fp_); // vz

			T_store eps = mass2eps( tmp7[i] );
			
			fwrite( &eps, sizeof(eps), 1, fp_);
			fwrite( &zero, sizeof(zero), 1, fp_);
		    }
		    
		} else {
		    
		    
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
			
		    }

		}		
		
		npleft -= n2read;
		n2read = std::min( block_buf_size_,npleft );
	    }
	    
            ifs_x.close();
            ifs_y.close();
            ifs_z.close();
            ifs_vx.close();
            ifs_vy.close();
            ifs_vz.close();
            ifs_m.close();
	    
	    
	    break;
	}
	
	
	fclose( fp_ );
	
	// clean up temp files
	unlink(fnx);
	unlink(fny);
	unlink(fnz);
	unlink(fnvx);
	unlink(fnvy);
	unlink(fnvz);
	unlink(fnm);
	
	if( with_baryons_ )
	{
	    unlink(fnbx);
	    unlink(fnby);
	    unlink(fnbz);
	    unlink(fnbvx);
	    unlink(fnbvy);
	    unlink(fnbvz);
	    unlink(fnbm);	
	}
	
	
	
	LOGINFO("TIPSY : done writing.");
    }
    
    
    
public:
    
    tipsy_output_plugin( config_file& cf )
	: output_plugin( cf ), ofs_( fname_.c_str(), std::ios::binary|std::ios::trunc )	
    {
	block_buf_size_ = cf_.getValueSafe<unsigned>("output","tipsy_blksize",10485760); // default buffer size is 10 MB
	
	//... ensure that everyone knows we want to do SPH
	cf.insertValue("setup","do_SPH","yes");
        with_baryons_ = cf_.getValue<bool>("setup","baryons");
	
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
        omegab_  = cf.getValue<double>("cosmology","Omega_b");
	boxsize_ = cf.getValue<double>("setup","boxlength");
	epsfac_ = cf.getValueSafe<double>("output","tipsy_eps",0.05);
        H0_ = cf.getValue<double>("cosmology","H0");
        YHe_ = cf.getValueSafe<double>("cosmology","YHe",0.248);
        gamma_ = cf.getValueSafe<double>("cosmology","gamma",5.0/3.0);
		
	native_ = cf.getValueSafe<bool>("output","tipsy_native",false);

	
    }
    
    void write_dm_mass( const grid_hierarchy& gh )
    {
	//.. get meta data about coarse/fine particle number
	
	size_t npcoarse = 0, npfine = 0;
	
	bmultimass_ = gh.levelmax() != gh.levelmin();
	
	npfine   = gh.count_leaf_cells(gh.levelmax(), gh.levelmax());
	if( bmultimass_ )
	    npcoarse = gh.count_leaf_cells(gh.levelmin(), gh.levelmax()-1);
	
	np_fine_dm_   = npfine;
        np_fine_gas_  = with_baryons_? npfine : 0ul;
        np_coarse_dm_ = npcoarse;
	
	
	//.. store header data
	header_.nbodies = npfine + npcoarse;
	header_.ndark	= header_.nbodies;
	header_.nsph	= 0;
	
        if( with_baryons_ )
        {
            header_.nsph	= npfine;
            header_.nbodies += header_.nsph;
        }
	
	
	header_.nstar	= 0;
	header_.ndim	= 3;
	header_.time	= astart_;
	
	
	//... write data for dark matter......
	size_t nptot = header_.ndark;
	
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
        
            if( with_baryons_ && ilevel == (int)gh.levelmax() )
                pmass *= (omegam_-omegab_)/omegam_;
	    
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
	{
	    LOGERR("TIPSY output plugin wrote %ld, should have %ld", nwritten, nptot);
	    throw std::runtime_error("Internal consistency error while writing temporary file for DM masses");
	}
	ofs_temp.write( (char *)&blksize, sizeof(size_t) );
	
	if( ofs_temp.bad() )
	    throw std::runtime_error("I/O error while writing temporary file for DM masses");
	
	
        ofs_temp.close();
	
        //... write data for baryons......
        if( with_baryons_ )
        {
            nptot = header_.nsph;
	    
            temp_dat.clear();
            temp_dat.reserve(block_buf_size_);
	    
            char temp_fnameb[256];
            sprintf( temp_fnameb, "___ic_temp_%05d.bin", 100*id_gas_mass );
            ofs_temp.open( temp_fnameb, std::ios::binary|std::ios::trunc );
	    
	    
            blksize = sizeof(T_store)*nptot;
            ofs_temp.write( (char *)&blksize, sizeof(size_t) );
	    
            nwritten = 0;
            int ilevel=gh.levelmax();
	    double pmass = omegam_/(1ul<<(3*ilevel));
	    
            pmass *= omegab_/omegam_;
	    
            for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
                for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
                    for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
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
	    
            if( temp_dat.size() > 0 )
            {	
                ofs_temp.write( (char*)&temp_dat[0], sizeof(T_store)*temp_dat.size() );		
                nwritten+=temp_dat.size();
            }
	    
            if( nwritten != nptot )
                throw std::runtime_error("Internal consistency error while writing temporary file for baryon masses");
	    
            ofs_temp.write( (char *)&blksize, sizeof(size_t) );
	    
            if( ofs_temp.bad() )
                throw std::runtime_error("I/O error while writing temporary file for baryon masses");
        }
	
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
	{
	    LOGERR("TIPSY output plugin wrote %ld, should have %ld", nwritten, nptot);
	    throw std::runtime_error("Internal consistency error while writing temporary file for positions");
	}
	
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
	{
	    LOGERR("TIPSY output plugin wrote %ld, should have %ld", nwritten, nptot);
	    throw std::runtime_error("Internal consistency error while writing temporary file for DM velocities");
	}
	
	//... dump to temporary file
	ofs_temp.write( (char *)&blksize, sizeof(size_t) );
	
	if( ofs_temp.bad() )
	    throw std::runtime_error("I/O error while writing temporary file for DM velocities");
	
	ofs_temp.close();
    }
    
    void write_dm_density( const grid_hierarchy& gh )
    {
	//... we don't care about DM density for TIPSY
    }
    
    void write_dm_potential( const grid_hierarchy& gh )
    {
        //... we don't care about DM potential for TIPSY
    }
    
    void write_gas_potential( const grid_hierarchy& gh )
    {
        //... we don't care about baryon potential for TIPSY
    }
    
    //... write data for gas -- don't do this
    void write_gas_velocity( int coord, const grid_hierarchy& gh )
    {
        //size_t npgas = gh.count_leaf_cells(gh.levelmax(), gh.levelmax());
	size_t npart = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());
	
	std::vector<T_store> temp_data;
	temp_data.reserve( block_buf_size_ );
	
	double vfac = 2.894405/(100.0 * astart_); 
	
	char temp_fname[256];
	sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_gas_vel+coord );
	std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
	
	size_t blksize = sizeof(T_store)*npart;
	ofs_temp.write( (char *)&blksize, sizeof(size_t) );
	
	size_t nwritten = 0;
	
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
	    ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*temp_data.size() );
	    nwritten += temp_data.size();
	}
	
	if( nwritten != npart )
	{
	    LOGERR("TIPSY output plugin wrote %ld, should have %ld", nwritten, npart);
	    throw std::runtime_error("Internal consistency error while writing temporary file for baryon velocities");
	}
	
	//... dump to temporary file
	ofs_temp.write( (char *)&blksize, sizeof(size_t) );
	
	if( ofs_temp.bad() )
	    throw std::runtime_error("I/O error while writing temporary file for baryon velocities");
	
	ofs_temp.close();
    }
    
    
    //... write only for fine level
    void write_gas_position( int coord, const grid_hierarchy& gh )
    { 
        //size_t npgas = gh.count_leaf_cells(gh.levelmax(), gh.levelmax());
	size_t npart = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());
	
	std::vector<T_store> temp_data;
	temp_data.reserve( block_buf_size_ );
	
	
	char temp_fname[256];
	sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_gas_pos+coord );
	std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
	
	size_t blksize = sizeof(T_store)*npart;
	ofs_temp.write( (char *)&blksize, sizeof(size_t) );
	
        double h = 1.0/(1ul<<gh.levelmax());
	
	
	size_t nwritten = 0;
	for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
	{
	    for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
		for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
		    for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
			if( ! gh.is_refined(ilevel,i,j,k) )
			{
			    double xx[3];
			    gh.cell_pos(ilevel, i, j, k, xx);
			    
			    //... shift particle positions (this has to be done as the same shift
			    //... is used when computing the convolution kernel for SPH baryons)
			    xx[coord] += 0.5*h;
			    
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
	}
	
	if( temp_data.size() > 0 )
	{	
	    ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*temp_data.size() );
	    nwritten += temp_data.size();
	}
	
	if( nwritten != npart )
	{
	    LOGERR("TIPSY output plugin wrote %ld, should have %ld", nwritten, npart);
	    throw std::runtime_error("Internal consistency error while writing temporary file for baryon positions");
	}
	
	//... dump to temporary file
	ofs_temp.write( (char *)&blksize, sizeof(size_t) );
	
	if( ofs_temp.bad() )
	    throw std::runtime_error("I/O error while writing temporary file for baryon positions");
	
	ofs_temp.close();
    }
    
    void write_gas_density( const grid_hierarchy& gh )
    { 
        //... we don't care about gas density for TIPSY
    }
    
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
