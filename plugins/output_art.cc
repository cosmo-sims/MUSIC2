/*
 
 output_art.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2012  Jose Onorbe & Oliver Hahn
 
 */
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>

#include "output.hh"

template< typename T_store=float >
class art_output_plugin : public output_plugin
{
public:
	bool do_baryons_;
	double omegab_, omegam_;
	double gamma_;
    double astart_;
    
protected:
    
    enum iofields {
		id_dm_mass, id_dm_vel, id_dm_pos
	};
	
	typedef struct io_header
	{
		char head[45];
		float aexpN; // current expansion factor
        float aexp0; // initial expansion factor
		float amplt; // Amplitude of density fluctuations
		float astep; // Delta a -> time step. 
				// This value is also stored in pt.dat (binary 1 float)
				// It is recalculated by art so just a small value should work
		int istep; // step (=0 in IC)
		int partw; // mass of highest res particle.
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
	        //float Omb0; // this parameter only appears in header in hydro runs
		float wpart[10]; // extras[0-9] particle masses from high res to low res (normalized to low res particle)
		//  Mass of smallest particle=wpart[0]*0m0*2.746e+11*(Box/NGRID)**3 -> Msun/h
		//  Mass of largest  particle=wpart[nspecies-1]*0m0*2.746e+11*(Box/NGRID)**3 -> Msun/h
		int lpart[10]; // extras[10-19] number of particles from high res to low res cumulative!!! 
		//(i.e., lpart[0]=Nhigh res particles; lpart[1]=lpart[0]+N_this_level; etc) so lpart[nspecies-1]=N total
	        float extras[80]; //extras[20-99] 
		     //extras[9]=iLblock ->0 in IC 
                     //extras[10]=LevMin  ->0 in IC
                     //extras[11]=LevSmall ->0 in IC
                     //extras[12]=LevLarge ->0 in IC
                     //extras[13]=Omegab  ->0 in IC; fix it?
                     //extras[14]=sig8    ->0 in IC; fix it?
                     //extras[15]=Spslope ->0 in IC; fix it? Slope of the Power spectrum
                     //extras[16]=iDEswtch ->0 in IC; DE Flag=0:LCDM 1:w 2:RP 3:SUGRA
                     //extras[17]=DEw0    ->0 in IC; w0 for DE z=0
                     //extras[18]=DEwprime ->0 in IC; DE parameter
		     //extras[59]= 0 or 1; is used as switch for random numbers generators [do not apply in music use 0?]
		     //extras[60]= lux - level of luxury  [do not apply in music use 0?]
		     //extras[79]=Lbox (Mpc/h)

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

public:


	explicit art_output_plugin ( config_file& cf )
	: output_plugin( cf )
	{
		
		if( mkdir( fname_.c_str(), 0777 ) )
                {
                        perror( fname_.c_str() );
                        throw std::runtime_error("Error in art_output_plugin!");
                }

		do_baryons_ = cf.getValueSafe<bool>("setup","baryons",false);
		omegab_  = cf.getValueSafe<double>("cosmology","Omega_b",0.045);
        omegam_  = cf.getValue<double>("cosmology","Omega_m");
        astart_  = cf.getValue<double>("cosmology","astart");

		YHe_ = cf.getValueSafe<double>("cosmology","YHe",0.248);
                gamma_ = cf.getValueSafe<double>("cosmology","gamma",5.0/3.0);
		//... set time ......................................................
        header_.aexpN = astart_;//1.0/(1.0+header_.redshift);
        header_.aexp0 = header_.aexpN;
		//etc, etc
        
        
        

	}

	void write_header_file() //PMcrd.DAT
	{
        char filename[256];
		sprintf( filename, "%s/PMcrd.DAT", fname_.c_str() );
        std::ofstream ofs_;
        ofs_.open(fname_.c_str(), std::ios::binary|std::ios::trunc );
		header this_header(header_);
		int blksize = sizeof(header); //529 in a dm only run; 533 in a baryon run
		ofs_.write( (char *)&blksize, sizeof(int) );
		ofs_.write( (char *)&this_header,sizeof(header));
		ofs_.write( (char *)&blksize, sizeof(int) );
		ofs_.close();
		
	}

	void write_pt_file() //pt.dat
	{
		char filename[256];
		sprintf( filename, "%s/pt.dat", fname_.c_str() );
        std::ofstream ofs_;
        ofs_.open(fname_.c_str(), std::ios::binary|std::ios::trunc );
		ptf this_ptf(ptf_);
		int blksize = sizeof(ptf); //4
		ofs_.write( (char *)&blksize, sizeof(int) );
		ofs_.write( (char *)&this_ptf,sizeof(ptf));
		ofs_.write( (char *)&blksize, sizeof(int) );
		ofs_.close();
		
	}

	void write_dm_pages()
	{
		//The direct format write the particle data in pages.
		// Each page of particles is read into a common block,
		// which has the structure: X(Npage),Y(Npage),Z(Npage),
		// Vx(Npage),Vy(Npage),Vz(Npage). 
		///The number of particles in each page (Npage) is Npage = Nrow**2
		// Npages = (N_particles -1)/NPAGE +1
		// so in last page sometimes can be tricky
	        // N_in_last=N_particles -NPAGE*(Npages-1)
		// There are NO Fortran size blocks pre or after these blocks!!
		//coordinates are in the range 1 - (NGRID+1)
		// so scale factor is  scaleX = Box/NGRID -> to Mpc/h (Box in Mpc/h) 
                //velocities are P = a_expansion*V_pec/(x_0H_0) where x_0 = comoving cell_size=Box/Ngrid;H_0 = Hubble at z=0
		// so scale factor is scaleV = BoxV/AEXPN/NGRID -> to km/s (BoxV is Box*100; aexpn=current expansion factor)
		// Contradiction with documentation?? one file for each type of particle
		// however Daniel sent me just one file for a zoom all particle info together. 

	}
    
    void write_dm_mass( const grid_hierarchy& gh )
	{
        //... write data for dark matter......
		size_t nptot = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());
		
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
            double pmass = omegam_/(1ul<<(3*ilevel)); // this needs to be adjusted to have the right units
            
            if( do_baryons_ )
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
			throw std::runtime_error("Internal consistency error while writing temporary file for DM masses");
		
		ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
		if( ofs_temp.bad() )
			throw std::runtime_error("I/O error while writing temporary file for DM masses");
        
        
        ofs_temp.close();
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
			throw std::runtime_error("Internal consistency error while writing temporary file for DM velocities");
		
		//... dump to temporary file
		ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
		if( ofs_temp.bad() )
			throw std::runtime_error("I/O error while writing temporary file for DM velocities");
		
		ofs_temp.close();
    }
    
    void write_dm_density( const grid_hierarchy& gh )
	{
		//... we don't care about DM density for art
	}
    
    void write_dm_potential( const grid_hierarchy& gh )
	{ }
	
	void write_gas_potential( const grid_hierarchy& gh )
	{ }
	
	void write_gas_velocity( int coord, const grid_hierarchy& gh )
	{
        
    }
	
    void write_gas_position( int coord, const grid_hierarchy& gh )
	{
        //... we don't care about gas positions in art
    }
    
    void write_gas_density( const grid_hierarchy& gh )
	{
        
    }

	void finalize( void )
	{ 	}
};

namespace{
	output_plugin_creator_concrete<art_output_plugin<float> > creator("art");
}
