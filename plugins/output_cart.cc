/*

 output_cart.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions
 for cosmological simulations

 Copyright (C) 2012  Jose Onorbe & Oliver Hahn

 */
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>

#include "output.hh"

template<typename T>
inline T bytereorder(T v )
{
	T rval;
	(reinterpret_cast<unsigned char*>(&rval))[3] = (reinterpret_cast<unsigned char*>(&v))[0];
	(reinterpret_cast<unsigned char*>(&rval))[2] = (reinterpret_cast<unsigned char*>(&v))[1];
	(reinterpret_cast<unsigned char*>(&rval))[1] = (reinterpret_cast<unsigned char*>(&v))[2];
	(reinterpret_cast<unsigned char*>(&rval))[0] = (reinterpret_cast<unsigned char*>(&v))[3];
	return rval;
}

#define NFILL 75
template< typename T_store=double >
class cart_output_plugin : public output_plugin
{
	public:
		bool do_baryons_;
		bool swap_endianness_;
		double omegab_, omegam_;
		double gamma_;
		double astart_;
		double zstart_;
		size_t npcdm_;
		int hsize_;

	protected:

		enum iofields {
			id_dm_pos, id_dm_vel, id_gas_pos, id_gas_vel, id_gas_pma //add fields here
		};

		typedef struct io_header
		{
			char head[45];
			float aexpN; // current expansion factor
			float aexp0; // initial expansion factor
			float amplt; // Amplitude of density fluctuations
			float astep; // Delta a -> time step.
			// This value is also stored in pt.dat (binary 1 float)
			// It is recalculated by art for the next steps so just a small value should work
			int istep; // step (=0 in IC)
			float partw; // mass of highest res particle.
			float TINTG; //=0 in IC
			float EKIN; //SUM 0.5 * m_i*(v_i**2) in code units
			float EKIN1; //=0 in IC
			float EKIN2; //=0 in IC
			float AU0; //=0 in IC
			float AEU0; //=0 in IC
			int NROWC; // Number of particles in 1 dim (number of particles per page = NROW**2)
			int NGRIDC; // Number of cells in 1 dim
			int nspecies; // number of dm species
			int Nseed; // random number used ( 0 for MUSIC? or set the random number used in the lowest level?)
			float Om0; //Omega_m
			float Oml0; //Omega_L
			float hubble; //hubble constant h=H/100
			float Wp5; //
			float Ocurv; //Omega_k
			float Omb0; // this parameter only appears in header in hydro runs
			float wpart[10]; // extras[0-9] particle masses from high res to low res (normalized to low res particle)
			//  Mass of smallest particle=wpart[0]*0m0*2.746e+11*(Box/NGRID)**3 -> Msun/h
			//  Mass of largest  particle=wpart[nspecies-1]*0m0*2.746e+11*(Box/NGRID)**3 -> Msun/h
			int lpart[10]; // extras[10-19] number of particles from high res to low res cumulative!!!
			//(i.e., lpart[0]=Nhigh res particles; lpart[1]=lpart[0]+N_this_level; etc) so lpart[nspecies-1]=N total
			float magic1;
			float DelDC;
			float abox;
			float Hbox;
			float magic2;
			float extras[NFILL]; //extras[25-99]
			//extras[74]=Lbox (Mpc/h)
		}header;

		typedef struct io_ptf
		{
			float astep;
		}ptf;

		header header_;
		ptf ptf_;
		std::string fname;
		size_t np_fine_gas_, np_fine_dm_, np_coarse_dm_;
		size_t block_buf_size_;
		size_t npartmax_;

		double YHe_;


		// helper class to read temp files
		class pistream : public std::ifstream
	{
		public:
			pistream (std::string fname, size_t npart )
				: std::ifstream( fname.c_str(), std::ios::binary )
			{
				size_t blk;

				if( !this->good() )
				{
					LOGERR("Could not open buffer file in CART output plug-in");
					throw std::runtime_error("Could not open buffer file in CART output plug-in");
				}

				this->read( (char*)&blk, sizeof(size_t) );

				if( blk != (size_t)(npart*sizeof(T_store)) )
				{
					LOGERR("Internal consistency error in CART output plug-in");
					LOGERR("Expected %d bytes in temp file but found %d",npart*(unsigned)sizeof(T_store),blk);
					throw std::runtime_error("Internal consistency error in CART output plug-in");
				}
			}

			pistream (){}

			void open(std::string fname, size_t npart )
			{
				std::ifstream::open( fname.c_str(), std::ios::binary );
				size_t blk;

				if( !this->good() )
				{
					LOGERR("Could not open buffer file \'%s\' in CART output plug-in",fname.c_str());
					throw std::runtime_error("Could not open buffer file in CART output plug-in");
				}

				this->read( (char*)&blk, sizeof(size_t) );

				if( blk != (size_t)(npart*sizeof(T_store)) )
				{
					LOGERR("Internal consistency error in CART output plug-in");
					LOGERR("Expected %d bytes in temp file but found %d",npart*(unsigned)sizeof(T_store),blk);
					throw std::runtime_error("Internal consistency error in CART output plug-in");
				}
			}
	};


		// non-public member functions
		void write_header_file( void ) //PMcrd.DAT
		{
			std::string headfname;
			if(do_baryons_){
				headfname = fname_ + "/music_H.mdh";
			}else{
				headfname = fname_ + "/music_D.mdh";
			}
			//std::string headfname = fname_ + "/PMcrd.DAT";
			std::ofstream ofs( headfname.c_str(), std::ios::trunc );
			//ofs.open(fname_.c_str(), std::ios::binary|std::ios::trunc );
			header this_header(header_);
			//Should be 529 in a dm only run; 533 in a baryon run
			//but not working for alignment so it must be written one by one:
			int blksize = hsize_;
			if( swap_endianness_ )
			{
				LOGINFO("CART : swap_endianness option enabled");
				blksize = bytereorder( blksize );
				this_header.aexpN = bytereorder( this_header.aexpN );
				this_header.aexp0 = bytereorder( this_header.aexp0 );
				this_header.amplt = bytereorder( this_header.amplt );
				this_header.astep = bytereorder( this_header.astep );
				this_header.istep = bytereorder( this_header.istep );
				this_header.partw = bytereorder( this_header.partw );
				this_header.TINTG = bytereorder( this_header.TINTG );
				this_header.EKIN = bytereorder( this_header.EKIN );
				this_header.EKIN1 = bytereorder( this_header.EKIN1 );
				this_header.EKIN2 = bytereorder( this_header.EKIN2 );
				this_header.AEU0 = bytereorder( this_header.AEU0 );
				this_header.AEU0 = bytereorder( this_header.AEU0 );
				this_header.NROWC = bytereorder( this_header.NROWC );
				this_header.NGRIDC = bytereorder( this_header.NGRIDC );
				this_header.nspecies = bytereorder( this_header.nspecies );
				this_header.Nseed = bytereorder( this_header.Nseed );
				this_header.Om0 = bytereorder( this_header.Om0);
				this_header.Oml0 = bytereorder( this_header.Oml0 );
				this_header.hubble = bytereorder( this_header.hubble );
				this_header.Wp5 = bytereorder( this_header.Wp5 );
				this_header.Ocurv = bytereorder( this_header.Ocurv );
				this_header.Omb0 = bytereorder( this_header.Omb0 );
				for( int i=0; i<10; ++i )
				{
					this_header.wpart[i] = bytereorder( this_header.wpart[i] );
					this_header.lpart[i] = bytereorder( this_header.lpart[i] );
				}
				for( int i=0; i<NFILL; ++i ) /* CCART fill is 75 long */
				{
					this_header.extras[i] = bytereorder( this_header.extras[i] );
				}
			}
			ofs.write( (char *)&blksize, sizeof(int) );
			//ofs.write( (char *)&this_header,sizeof(header));  //Not working because struct aligment, so:
			ofs.write( (char *)&this_header.head,sizeof(this_header.head));
			ofs.write( (char *)&this_header.aexpN,sizeof(this_header.aexpN));
			ofs.write( (char *)&this_header.aexp0,sizeof(this_header.aexp0));
			ofs.write( (char *)&this_header.amplt,sizeof(this_header.amplt));
			ofs.write( (char *)&this_header.astep,sizeof(this_header.astep));
			ofs.write( (char *)&this_header.istep,sizeof(this_header.istep));
			ofs.write( (char *)&this_header.partw,sizeof(this_header.partw));
			ofs.write( (char *)&this_header.TINTG,sizeof(this_header.TINTG));
			ofs.write( (char *)&this_header.EKIN,sizeof(this_header.EKIN));
			ofs.write( (char *)&this_header.EKIN1,sizeof(this_header.EKIN1));
			ofs.write( (char *)&this_header.EKIN2,sizeof(this_header.EKIN2));
			ofs.write( (char *)&this_header.AEU0,sizeof(this_header.AEU0));
			ofs.write( (char *)&this_header.AEU0,sizeof(this_header.AEU0));
			ofs.write( (char *)&this_header.NROWC,sizeof(this_header.NROWC));
			ofs.write( (char *)&this_header.NGRIDC,sizeof(this_header.NGRIDC));
			ofs.write( (char *)&this_header.nspecies,sizeof(this_header.nspecies));
			ofs.write( (char *)&this_header.Nseed,sizeof(this_header.Nseed));
			ofs.write( (char *)&this_header.Om0,sizeof(this_header.Om0));
			ofs.write( (char *)&this_header.Oml0,sizeof(this_header.Oml0));
			ofs.write( (char *)&this_header.hubble,sizeof(this_header.hubble));
			ofs.write( (char *)&this_header.Wp5,sizeof(this_header.Wp5));
			ofs.write( (char *)&this_header.Ocurv,sizeof(this_header.Ocurv));
			ofs.write( (char *)&this_header.Omb0,sizeof(this_header.Omb0));   /*Omb0 decides whether the header is a baryon or no-baryon header*/
			ofs.write( (char *)&this_header.wpart,sizeof(this_header.wpart));   //mass[10]
			ofs.write( (char *)&this_header.lpart,sizeof(this_header.lpart));     //num[10]
			ofs.write( (char *)&this_header.magic1,sizeof(this_header.magic1));
			ofs.write( (char *)&this_header.DelDC,sizeof(this_header.DelDC));
			ofs.write( (char *)&this_header.abox,sizeof(this_header.abox));
			ofs.write( (char *)&this_header.Hbox,sizeof(this_header.Hbox));
			ofs.write( (char *)&this_header.magic2,sizeof(this_header.magic2));
			ofs.write( (char *)&this_header.extras,sizeof(this_header.extras));
			ofs.write( (char *)&blksize, sizeof(int) );
			ofs.close();
			LOGINFO("CART : done writing header file.");
		}

		void write_pt_file( void ) //pt.dat
		{
			// 	    std::string partfname = fname_ + "/pt.dat";
			// 	    std::ofstream ofs( partfname.c_str(), std::ios::trunc );
			// 	    //ofs.open(fname_.c_str(), std::ios::binary|std::ios::trunc );
			// 	    ptf this_ptf(ptf_);
			// 	    int blksize = sizeof(ptf); //4
			// 	    if( swap_endianness_ )
			// 		{
			// 		    blksize = bytereorder( blksize );
			// 		    this_ptf = bytereorder( this_ptf );
			// 		}
			// 	    ofs.write( (char *)&blksize, sizeof(int) );
			// 	    ofs.write( (char *)&this_ptf,sizeof(ptf));
			// 	    ofs.write( (char *)&blksize, sizeof(int) );
			// 	    ofs.close();
			// 	    LOGINFO("CART : done writing pt file.");
		}


		void adjust_buf_endianness( T_store* buf )
		{
			if( swap_endianness_ )
			{
				for( size_t i=0; i<block_buf_size_; ++i )
					buf[i] = bytereorder<T_store>( buf[i] );
			}
		}

		/*
		   The direct format write the particle data in pages. Each page of particles is read into a common block,
		   which has the structure: X(Npage),Y(Npage),Z(Npage),Vx(Npage),Vy(Npage),Vz(Npage).
		   There are NO Fortran size blocks pre or after these blocks!!

		   The number of particles in each page (Npage) is Npage = Nrow**2; Npages = (N_particles -1)/NPAGE +1
		   so in last page sometimes can be tricky (zooms): N_in_last=N_particles -NPAGE*(Npages-1)
		   But keep in mind that CART expects all pages to be written in full regarding of the actual number of particles
		   you care about.

*/
		void assemble_DM_file( void ) //PMcrs0.DAT
		{
			// file name
			std::string partfname;
			if(do_baryons_){
				partfname = fname_ + "/music_H.mdxv";
			}else{
				partfname = fname_ + "/music_D.mdxv";
			}
			//std::string partfname = fname_ + "/PMcrs0.DAT";
			std::ofstream ofs( partfname.c_str(), std::ios::trunc );

			// generate all temp file names
			char fnx[256],fny[256],fnz[256],fnvx[256],fnvy[256],fnvz[256];
			sprintf( fnx,  "___ic_temp_%05d.bin", 100*id_dm_pos+0 );
			sprintf( fny,  "___ic_temp_%05d.bin", 100*id_dm_pos+1 );
			sprintf( fnz,  "___ic_temp_%05d.bin", 100*id_dm_pos+2 );
			sprintf( fnvx, "___ic_temp_%05d.bin", 100*id_dm_vel+0 );
			sprintf( fnvy, "___ic_temp_%05d.bin", 100*id_dm_vel+1 );
			sprintf( fnvz, "___ic_temp_%05d.bin", 100*id_dm_vel+2 );

			// create buffers for temporary data
			T_store *tmp1, *tmp2, *tmp3, *tmp4, *tmp5, *tmp6;

			tmp1 = new T_store[block_buf_size_];
			tmp2 = new T_store[block_buf_size_];
			tmp3 = new T_store[block_buf_size_];
			tmp4 = new T_store[block_buf_size_];
			tmp5 = new T_store[block_buf_size_];
			tmp6 = new T_store[block_buf_size_];


			// read in the data from the temporary files in slabs and write it to the output file
			size_t npleft, n2read;
			size_t npcdm = npcdm_;

			LOGINFO("writing DM data to CART format file");
			//ofs.open(fname_.c_str(), std::ios::binary|std::ios::trunc );

			pistream ifs_x, ifs_y, ifs_z, ifs_vx, ifs_vy, ifs_vz;

			ifs_x.open( fnx, npcdm );
			ifs_y.open( fny, npcdm );
			ifs_z.open( fnz, npcdm );
			ifs_vx.open( fnvx, npcdm );
			ifs_vy.open( fnvy, npcdm );
			ifs_vz.open( fnvz, npcdm );

			npleft = npcdm;
			n2read = std::min(block_buf_size_,npleft);
			while( n2read > 0 )
			{
				// To make sure last page in zooms have 0s in non-relevant values
				// NOT MANDATORY. Can be commented if makes things slow
				// but I do not like the idea of writting data in the file
				//  that could be interpreted as real.
				if(n2read<block_buf_size_)
				{
					for (uint i = 0; i < block_buf_size_; i++)
					{
						tmp1[i]=0.0;tmp2[i]=0.0;tmp3[i]=0.0;tmp4[i]=0.0;
						tmp5[i]=0.0;tmp6[i]=0.0;
					}
				}
				ifs_x.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
				ifs_y.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
				ifs_z.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
				ifs_vx.read( reinterpret_cast<char*>(&tmp4[0]), n2read*sizeof(T_store) );
				ifs_vy.read( reinterpret_cast<char*>(&tmp5[0]), n2read*sizeof(T_store) );
				ifs_vz.read( reinterpret_cast<char*>(&tmp6[0]), n2read*sizeof(T_store) );

				adjust_buf_endianness( tmp1 );
				adjust_buf_endianness( tmp2 );
				adjust_buf_endianness( tmp3 );
				adjust_buf_endianness( tmp4 );
				adjust_buf_endianness( tmp5 );
				adjust_buf_endianness( tmp6 );

				ofs.write( reinterpret_cast<char*>(&tmp1[0]), block_buf_size_*sizeof(T_store) );
				ofs.write( reinterpret_cast<char*>(&tmp2[0]), block_buf_size_*sizeof(T_store) );
				ofs.write( reinterpret_cast<char*>(&tmp3[0]), block_buf_size_*sizeof(T_store) );
				ofs.write( reinterpret_cast<char*>(&tmp4[0]), block_buf_size_*sizeof(T_store) );
				ofs.write( reinterpret_cast<char*>(&tmp5[0]), block_buf_size_*sizeof(T_store) );
				ofs.write( reinterpret_cast<char*>(&tmp6[0]), block_buf_size_*sizeof(T_store) );

				npleft -= n2read;
				n2read = std::min( block_buf_size_,npleft );
			}

			ifs_x.close();
			ifs_y.close();
			ifs_z.close();
			ifs_vx.close();
			ifs_vy.close();
			ifs_vz.close();
			ofs.close();

			// clean up temp files
			unlink(fnx);
			unlink(fny);
			unlink(fnz);
			unlink(fnvx);
			unlink(fnvy);
			unlink(fnvz);

			delete[] tmp1;
			delete[] tmp2;
			delete[] tmp3;
			delete[] tmp4;
			delete[] tmp5;
			delete[] tmp6;

			LOGINFO("CART : done writing DM file.");

		}


		/*
		   CART users currently create the baryon grid structure from the dark matter data file.
		   Therefore they have decided that the best way to implement baryons for CART in MUSIC was
		   by creating a file with the same dm format but using the baryon displacements and velocities.
		   From this file they will create the actual grid suign their tools.

		   So here we have just to re-create the dark matter file format but using the baryon data.
		   */
		void assemble_gas_file( void ) //PMcrs0IC.DAT
		{
			// file name
			//	    std::string hydrofname = fname_ + "/PMcrs0IC.DAT";
			std::string hydrofname;
			if(do_baryons_){
				hydrofname = fname_ + "/music_H.md";
			}else{
				hydrofname = fname_ + "/music_D.md";
			}
			std::ofstream ofs( hydrofname.c_str(), std::ios::trunc );

			// generate all temp file names
			char fnx[256],fny[256],fnz[256],fnvx[256],fnvy[256],fnvz[256],fnpma[256]; //add fields here
			sprintf( fnx,  "___ic_temp_%05d.bin", 100*id_gas_pos+0 );
			sprintf( fny,  "___ic_temp_%05d.bin", 100*id_gas_pos+1 );
			sprintf( fnz,  "___ic_temp_%05d.bin", 100*id_gas_pos+2 );
			sprintf( fnvx, "___ic_temp_%05d.bin", 100*id_gas_vel+0 );
			sprintf( fnvy, "___ic_temp_%05d.bin", 100*id_gas_vel+1 );
			sprintf( fnvz, "___ic_temp_%05d.bin", 100*id_gas_vel+2 );
			sprintf( fnpma,  "___ic_temp_%05d.bin", 100*id_gas_pma ); //add fields here

			// create buffers for temporary data
			T_store *tmp1, *tmp2, *tmp3, *tmp4, *tmp5, *tmp6, *tmp7; //add fields here

			tmp1 = new T_store[block_buf_size_];
			tmp2 = new T_store[block_buf_size_];
			tmp3 = new T_store[block_buf_size_];
			tmp4 = new T_store[block_buf_size_];
			tmp5 = new T_store[block_buf_size_];
			tmp6 = new T_store[block_buf_size_];
			tmp7 = new T_store[block_buf_size_]; //add fields here


			// read in the data from the temporary files in slabs and write it to the output file
			size_t npleft, n2read;
			size_t npcgas = npcdm_; // # of gas elemets should be equal to # of dm elements

			LOGINFO("writing gas data to CART format file");
			//ofs.open(fname_.c_str(), std::ios::binary|std::ios::trunc );

			pistream ifs_x, ifs_y, ifs_z, ifs_vx, ifs_vy, ifs_vz, ifs_pma;


			ifs_x.open( fnx, npcgas );
			ifs_y.open( fny, npcgas );
			ifs_z.open( fnz, npcgas );
			ifs_vx.open( fnvx, npcgas );
			ifs_vy.open( fnvy, npcgas );
			ifs_vz.open( fnvz, npcgas );
			ifs_pma.open( fnpma, npcgas );

			npleft = npcgas;
			n2read = std::min(block_buf_size_,npleft);
			while( n2read > 0 )
			{
				// To make sure last page in zooms have 0s in non-relevant values
				// NOT MANDATORY. Can be commented if makes things slow
				// but I do not like the idea of writting data in the file
				//  that could be interpreted as real.
				if(n2read<block_buf_size_)
				{
					for (uint i = 0; i < block_buf_size_; i++)
					{
						tmp1[i]=0.0;tmp2[i]=0.0;tmp3[i]=0.0;tmp4[i]=0.0;
						tmp5[i]=0.0;tmp6[i]=0.0;
					}
				}
				ifs_x.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
				ifs_y.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
				ifs_z.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
				ifs_vx.read( reinterpret_cast<char*>(&tmp4[0]), n2read*sizeof(T_store) );
				ifs_vy.read( reinterpret_cast<char*>(&tmp5[0]), n2read*sizeof(T_store) );
				ifs_vz.read( reinterpret_cast<char*>(&tmp6[0]), n2read*sizeof(T_store) );
				ifs_pma.read( reinterpret_cast<char*>(&tmp7[0]), n2read*sizeof(T_store) );

				adjust_buf_endianness( tmp1 );
				adjust_buf_endianness( tmp2 );
				adjust_buf_endianness( tmp3 );
				adjust_buf_endianness( tmp4 );
				adjust_buf_endianness( tmp5 );
				adjust_buf_endianness( tmp6 );
				adjust_buf_endianness( tmp7 );

				ofs.write( reinterpret_cast<char*>(&tmp1[0]), block_buf_size_*sizeof(T_store) );
				ofs.write( reinterpret_cast<char*>(&tmp2[0]), block_buf_size_*sizeof(T_store) );
				ofs.write( reinterpret_cast<char*>(&tmp3[0]), block_buf_size_*sizeof(T_store) );
				ofs.write( reinterpret_cast<char*>(&tmp4[0]), block_buf_size_*sizeof(T_store) );
				ofs.write( reinterpret_cast<char*>(&tmp5[0]), block_buf_size_*sizeof(T_store) );
				ofs.write( reinterpret_cast<char*>(&tmp6[0]), block_buf_size_*sizeof(T_store) );
				ofs.write( reinterpret_cast<char*>(&tmp7[0]), block_buf_size_*sizeof(T_store) );

				npleft -= n2read;
				n2read = std::min( block_buf_size_,npleft );
			}

			ifs_x.close();
			ifs_y.close();
			ifs_z.close();
			ifs_vx.close();
			ifs_vy.close();
			ifs_vz.close();
			ifs_pma.close();
			ofs.close();

			// clean up temp files
			unlink(fnx);
			unlink(fny);
			unlink(fnz);
			unlink(fnvx);
			unlink(fnvy);
			unlink(fnvz);
			unlink(fnpma);

			delete[] tmp1;
			delete[] tmp2;
			delete[] tmp3;
			delete[] tmp4;
			delete[] tmp5;
			delete[] tmp6;
			delete[] tmp7;

			LOGINFO("CART : done writing gas file.");

		}

	public:


		explicit cart_output_plugin ( config_file& cf )
			: output_plugin( cf )
		{
			mkdir( fname_.c_str(), 0777 );

			do_baryons_ = cf.getValue<bool>("setup","baryons");
			hsize_ = 533; // dm & hydro run (omega_b is included in header -- 529 for oldstyle)

			omegab_  = cf.getValue<double>("cosmology","Omega_b");
			omegam_  = cf.getValue<double>("cosmology","Omega_m");
			zstart_  = cf.getValue<double>("setup","zstart");
			astart_ = 1.0/(1.0+zstart_);

			//snl this doesn't corrently swap particle endianness and you check on the CART end anyway
			swap_endianness_ = cf.getValueSafe<bool>("output","art_swap_endian",false);

			int levelmin = cf.getValue<unsigned>("setup","levelmin");
			int levelmax = cf.getValue<unsigned>("setup","levelmax");
			block_buf_size_ = (size_t) (pow(pow(2,levelmax),2)); //Npage=nrow^2; Number of particles in each page

			YHe_ = cf.getValueSafe<double>("cosmology","YHe",0.248);
			gamma_ = cf.getValueSafe<double>("cosmology","gamma",5.0/3.0);
			// Set header
			std::string thead;
			thead=cf.getValueSafe<std::string>("output","header","ICs generated using MUSIC");
			strcpy(header_.head,thead.c_str()); // text for the header; any easy way to add also the version?
			std::string ws = " "; // Filling with blanks. Any better way?
			for (int i=thead.size(); i<45;i++)
			{
				header_.head[i]=ws[0];
			}
			header_.aexpN = astart_;
			header_.aexp0 = header_.aexpN;
			header_.amplt = 0.0; // Amplitude of density fluctuations
			header_.astep = 0.0; //cf.getValue<double>("output","astep");
			ptf_.astep=header_.astep; // to write pt file
			header_.istep = 0; // step (=0 in IC)
			header_.partw = 0.0; // mass of highest res particle. SEE BELOW

			/* primordial gas state */
			//	    const double npol  = (fabs(1.0-gamma_)>1e-7)? 1.0/(gamma_-1.) : 1.0;
			//	    const double unitv = 1e5;
			const double h2    = header_.hubble*header_.hubble*0.0001;
			const double adec  = 1.0/(160.*pow(omegab_*h2/0.022,2.0/5.0));
			const double Tcmb0 = 2.726;
			const double Tini  = astart_<adec? Tcmb0/astart_ : Tcmb0/astart_/astart_*adec;
			//	    const double mu    = (Tini>1.e4) ? 4.0/(8.-5.*YHe_) : 4.0/(1.+3.*(1.-YHe_));
			//	    const double ceint = 1.3806e-16/1.6726e-24 * Tini * npol / mu / unitv / unitv;
			header_.TINTG = Tini; //=0 in IC

			header_.EKIN = 0.0;
			header_.EKIN1 = 0;
			header_.EKIN2 = 0;
			header_.AU0 = 0;
			header_.AEU0 = 0;

			header_.NROWC = (int) pow(2,levelmax); // Number of particles in 1 dim (number of particles per page = NROW**2)
			header_.NGRIDC = (int) pow(2,levelmin); // Number of cells in 1 dim

			header_.nspecies = 0; // number of dm particle species
			for( int ilevel=levelmax; ilevel>=(int)levelmin; --ilevel )
			{
				header_.nspecies+=1;
			}

			header_.Nseed = 0; // random number used ( 0 for MUSIC? or set the random number used in the lowest level?)
			header_.Omb0 = cf.getValue<double>("cosmology","Omega_b");; // this parameter only appears in header in hydro runs
			header_.Om0 = cf.getValue<double>("cosmology","Omega_m"); //Omega_m
			header_.Oml0 = cf.getValue<double>("cosmology","Omega_L"); //Omega_L
			header_.hubble = cf.getValue<double>("cosmology","H0")/100.; //hubble constant h=H/100
			header_.Wp5 = 0.0;
			header_.Ocurv = 1.0 - header_.Oml0 - header_.Om0;

			for (int i=0;i<10;i++)
			{
				header_.wpart[i] = 0.0; // extras[0-9] part. masses from high res to low res (normalized to low res particle)
				header_.lpart[i] = 0; // extras[10-19] # particles from high res to low res cumulative!!!
			}
			for (int i=0;i<header_.nspecies;i++)
			{
				header_.wpart[i] = 1.0/pow(8.0,(header_.nspecies-i-1)); //from high res to lo res // 8 should be changed for internal variable?
				if(do_baryons_){
					header_.wpart[i] *= (header_.Om0 - header_.Omb0)/header_.Om0 ; //from high res to lo res // 8 should be changed for internal variable?
				}
			}
			header_.partw = header_.wpart[0]; // mass of highest res particle.
			for (int i=0;i<NFILL;i++)
			{
				header_.extras[i] = 0.0; //extras[20-99]
			}
			header_.extras[13] = cf.getValueSafe<double>("cosmology","Omega_b",0.0);
			header_.extras[14] = cf.getValue<double>("cosmology","sigma_8");
			header_.extras[15] = cf.getValue<double>("cosmology","nspec"); //Slope of the Power spectrum
			header_.magic1 = 0.1234f ;
			header_.DelDC = 0.0;
			header_.abox = astart_;
			header_.Hbox = 0.0;  // not used
			header_.magic2 = 0.1234f ;


			//header_.extras[NFILL-2] = Tini;
			header_.extras[NFILL-1] = cf.getValue<double>("setup","boxlength");

			LOGINFO("CART : done header info.");

		}



		void write_dm_mass( const grid_hierarchy& gh ){}

		void write_dm_position( int coord, const grid_hierarchy& gh )
		{
			size_t nptot = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());
			//... store all the meta data about the grid hierarchy in header variables
			npcdm_ = nptot;
			for (int i=0;i<header_.nspecies;i++)
			{
				header_.lpart[i] = gh.count_leaf_cells(gh.levelmax()-i, gh.levelmax()); //cumulative!!
			}

			// Now, let us write the dm particle info
			std::vector<T_store> temp_data;
			temp_data.reserve( block_buf_size_ );


			//coordinates are in the range 1 - (NGRID+1)
			// so scale factor is  scaleX = Box/NGRID -> to Mpc/h (Box in Mpc/h)
			double xfac = (double) header_.NGRIDC;

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
							if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
							{
								double xx[3];
								gh.cell_pos(ilevel, i, j, k, xx);

								xx[coord] = fmod( (xx[coord]+(*gh.get_grid(ilevel))(i,j,k)) + 1.0, 1.0 ) ;
								xx[coord] = (xx[coord]*xfac)+1.0;
								//xx[coord] = ((xx[coord]+(*gh.get_grid(ilevel))(i,j,k)));

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

			// t0_internal = 2 * aexpn^2/(100*h*sqrt(Om0))
			// r0_internal^{-1} = Ng/(boxh/h*aexpn)
			// v0^{-1} = t0/r0 = aexpn * Ng / {50*sqrt(Om0) * boxh}

			// v_code = v_kms * v0^{-1}

			// v_kms = v_music * boxh           [from gadget2 vfac ]
			// v_code = v_music * boxh * v0^{-1}

			// v_code = v_music * aexpn * Ng / {50*sqrt(Om0)}
			double vfac ;
			if(do_baryons_){
				vfac =  (header_.aexpN*header_.NGRIDC) / (50.0*sqrt(header_.Om0) );
			}else{
				vfac =  (header_.aexpN*header_.NGRIDC) / (50.0*sqrt(header_.Om0) );
			}

			//snl	    std::cout << "snl ae" <<  header_.aexpN << " ng "<< header_.NGRIDC << " h " <<  header_.hubble <<" om0 "<< sqrt(header_.Om0) << " boxh "<<header_.extras[NFILL-1] << "\n" ;
			//snl	    std::cout << "snl " << 1.0 / ( (header_.aexpN*header_.NGRIDC) / (50.0*header_.hubble*sqrt(header_.Om0) * header_.extras[NFILL-1] ) ) << "\n" ;
			//snl	    exit(1);

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
							if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
							{
								if( temp_data.size() < block_buf_size_ ){
									//snl					std::cout << "coord " << coord<< " "<< i <<" " << j << " " << k << " " << (*gh.get_grid(ilevel))(i,j,k) * header_.extras[NFILL-1] << "\n" ; //snl
									temp_data.push_back( (*gh.get_grid(ilevel))(i,j,k) * vfac );
								}else
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
				throw std::runtime_error("Internal consistency error while writing temporary file for DM velocities");

			//... dump to temporary file
			ofs_temp.write( (char *)&blksize, sizeof(size_t) );

			if( ofs_temp.bad() )
				throw std::runtime_error("I/O error while writing temporary file for DM velocities");

			ofs_temp.close();
		}

		void write_dm_density( const grid_hierarchy& gh ){} //... we don't care about DM density for art
		void write_dm_potential( const grid_hierarchy& gh ){}

		void write_gas_position( int coord, const grid_hierarchy& gh ) {}

		void write_gas_velocity( int coord, const grid_hierarchy& gh )
		{
			size_t nptot = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());

			std::vector<T_store> temp_data;
			temp_data.reserve( block_buf_size_ );

			// 	    // t0_internal = 2 * aexpn^2/(100*h*sqrt(Om0))
			// 	    // r0_internal^{-1} = Ng/(boxh*aexpn)
			// 	    // v0^{-1} = t0/r0 = aexpn * Ng / {50*h*sqrt(Om0) * boxh}
			// 	    // v_kms * v0^{-1} = v_code
			// 	    // v_kms = v_music * boxh/sqrt(aexpn)          [from gadget2 vfac ]
			// 	    // v_code = v_music * boxh/sqrt(aexpn) * v0^{-1}
			// 	    // v_code = v_music * sqrt(aexpn) * Ng / {50*h*sqrt(Om0)}

			//Jose says  	    // v_kms = v_music * boxh          [from gadget2 vfac ]  ....
			double vfac ;
			if(do_baryons_){
				vfac =  ((header_.aexpN)*header_.NGRIDC) / (50.0*header_.hubble*sqrt(header_.Om0) );
			}else{
				vfac =  ((header_.aexpN)*header_.NGRIDC) / (50.0*header_.hubble*sqrt(header_.Om0) );
			}

			char temp_fname[256];
			sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_gas_vel+coord );
			std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );

			size_t blksize = sizeof(T_store)*nptot;
			ofs_temp.write( (char *)&blksize, sizeof(size_t) );

			size_t nwritten = 0;
			for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
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
				ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*temp_data.size() );
				nwritten += temp_data.size();
			}

			if( nwritten != nptot )
				throw std::runtime_error("Internal consistency error while writing temporary file for gas velocities");

			//... dump to temporary file
			ofs_temp.write( (char *)&blksize, sizeof(size_t) );

			if( ofs_temp.bad() )
				throw std::runtime_error("I/O error while writing temporary file for gas velocities");

			ofs_temp.close();
		}

		void write_gas_density( const grid_hierarchy& gh ){
			size_t nptot = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());
			char temp_fname[256];
			std::vector<T_store> temp_data;
			temp_data.reserve( block_buf_size_ );
			size_t blksize = sizeof(T_store)*nptot;

			double xfac = (double) header_.NGRIDC;

			// write gas positions to cell centers
			for (int coord=0; coord < 3; coord++ ) {
				sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_gas_pos+coord );
				std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
				ofs_temp.write( (char *)&blksize, sizeof(size_t) );

				size_t nwritten = 0;
				for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
					for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
						for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
							for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
								if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
								{
									double xx[3];
									gh.cell_pos(ilevel, i, j, k, xx);

									// for gas positions just leave particle centered on the grid cell (no gh shift)
									xx[coord] = (xx[coord]*xfac)+1.0;

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
					temp_data.clear();
				}

				if( nwritten != nptot )
					throw std::runtime_error("Internal consistency error while writing temporary file for gas positions");

				//... dump to temporary file
				ofs_temp.write( (char *)&blksize, sizeof(size_t) );

				if( ofs_temp.bad() )
					throw std::runtime_error("I/O error while writing temporary file for gas positions");

				ofs_temp.close();
			}

			// write gas densities
			{
				double pmafac = header_.Omb0 / header_.Om0 ;
				double pma;
				sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_gas_pma);
				std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
				ofs_temp.write( (char *)&blksize, sizeof(size_t) );

				size_t nwritten = 0;
				for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
					for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
						for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
							for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
								if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
								{
									pma = ( 1 + (*gh.get_grid(ilevel))(i,j,k) ) * pmafac * pow(8.0, -1.0*(ilevel-gh.levelmin())) ;
									if( temp_data.size() < block_buf_size_ )
										temp_data.push_back( pma );
									else
									{
										ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*block_buf_size_ );
										nwritten += block_buf_size_;
										temp_data.clear();
										temp_data.push_back( pma );
									}
								}

				if( temp_data.size() > 0 )
				{
					ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*temp_data.size() );
					nwritten += temp_data.size();
				}

				if( nwritten != nptot )
					throw std::runtime_error("Internal consistency error while writing temporary file for gas densities");

				//... dump to temporary file
				ofs_temp.write( (char *)&blksize, sizeof(size_t) );

				if( ofs_temp.bad() )
					throw std::runtime_error("I/O error while writing temporary file for gas densities");

				ofs_temp.close();
			}
		}


		void write_gas_potential( const grid_hierarchy& gh ){}

		void finalize( void )
		{
			this->write_header_file();
			//	    this->write_pt_file();
			this->assemble_DM_file();
			if(do_baryons_)
			{
				this->assemble_gas_file();
			}
		}
};

namespace{ output_plugin_creator_concrete<cart_output_plugin<double> > creator("cart");}
