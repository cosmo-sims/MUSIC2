/*
 
 output_gadget2.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#include <fstream>
#include <cstdlib>
#include <iostream>
#include <bitset>
#include <map>
#include <cstring>
#include <algorithm>

#include "log.hh"
#include "region_generator.hh"
#include "output.hh"
#include "mg_interp.hh"
#include "mesh.hh"


#define get_lagrange_coord(x,i) ( ((x)>>(20*(2-(i))))&1048575ll)
#define get_lagrange_ID(ix,iy,iz) ((((long long)(ix))<<40)+(((long long)(iy))<<20)+((long long)(iz)))

#define REMOVE_REFINEMENT_COUNTER(x) ((x)&(~(15ll<<60)))
#define ADD_REFINEMENT_COUNTERS(x,y) ((x)+((y)&(15ll<<60)))  // adds the counter of y to x, maintains rest of x
#define SET_REFINEMENT_COUNTER(x,n) ((((x)&(~(15ll<<60)))+(((long long)n)<<60)))
#define GET_REFINEMENT_COUNTER(x) (((x)>>60)&(15ll))

#define REMOVE_FREE_PARTICLE_BIT(x) ((x)&(~(1ll<<63)))
#define SET_FREE_PARTICLE_BIT(x) ((x) |= (1ll<<63))
#define GET_FREE_PARTICLE_BIT(x) ((x>>63)&1ll)

#define REMOVE_DECORATION_BITS(x) ((x)&(~(1ll<<63))&(~(15ll<<60)))


typedef float MyFloat;
typedef long long MyIDType;

const int vert[8][3] = { {0,0,0}, {0,0,1}, {0,1,0}, {0,1,1}, {1,0,0}, {1,0,1}, {1,1,0}, {1,1,1} };
const MyIDType vert_long[8][3] = { {0,0,0}, {0,0,1}, {0,1,0}, {0,1,1}, {1,0,0}, {1,0,1}, {1,1,0}, {1,1,1} };
const int conn[][4] = { {1,0,2,4}, {3,1,2,4}, {3,5,1,4}, {3,6,5,4}, {3,2,6,4}, {3,7,5,6} };


std::map<MyIDType,size_t> idmap;
typedef std::map<MyIDType,size_t>::iterator idmap_it;

const size_t num_p_realloc_blocksize = 1000000;
const size_t newnum_p_per_split = 19;

int tetgrid_levelmax = 0;


struct particle
{
    //MyIDType ID;
    MyIDType Lagrange_ID;
    int Type;
    int Level;
    
    
    bool operator<( const particle& p ) const
    {
        if( REMOVE_DECORATION_BITS(Lagrange_ID) < REMOVE_DECORATION_BITS(p.Lagrange_ID) )
            return true;
        if( REMOVE_DECORATION_BITS(Lagrange_ID) > REMOVE_DECORATION_BITS(p.Lagrange_ID) )
            return false;
        
        if( Type < p.Type )
            return true;
        
        return false;
    }
    
    particle& operator=( const particle& p )
    {
        memcpy( this, &p, sizeof(particle) );
        return *this;
    }
    
    MyIDType get_vertex( int idx, int isub=0 )
    {
        int rfac = 20-Level-isub;
        //MyIDType nfull = ((MyIDType)1)<<20;
        MyIDType lid  = REMOVE_DECORATION_BITS(Lagrange_ID);
        MyIDType lidx = get_lagrange_coord( lid, 0 );
        MyIDType lidy = get_lagrange_coord( lid, 1 );
        MyIDType lidz = get_lagrange_coord( lid, 2 );
        
        //MyIDType lidold = lid;
        
        lidx = (lidx+(vert_long[idx][0] << (rfac))) % (1ull<<20);
        lidy = (lidy+(vert_long[idx][1] << (rfac))) % (1ull<<20);
        lidz = (lidz+(vert_long[idx][2] << (rfac))) % (1ull<<20);
        
        lid = (lidx<<40) + (lidy<<20) + lidz;
        
        return lid;
    }
    
    bool can_refine( void )
    {
        MyIDType lid = REMOVE_DECORATION_BITS( Lagrange_ID );
        MyIDType xl[3] = {get_lagrange_coord(lid,0),get_lagrange_coord(lid,1),get_lagrange_coord(lid,2)};
        MyIDType xm = (1<<20) - (1<<(20-Level));
        
        if( xl[0]>0 && xl[1]>0 && xl[2] > 0 && xl[0] < xm && xl[1] < xm && xl[2] < xm )
            return true;
        
        return false;
        
    }
    
};

bool sortP_bytype( const particle& p1, const particle& p2)
{  return p1.Type < p2.Type;  }


particle *P;
size_t num_p, num_p_alloc;
int TetRefinementIDFactor = 0;



inline real_t get_cic( const grid_hierarchy & gh, int ilevel, real_t u, real_t v, real_t w )
{
    int i,j,k;
    //int levelmax_ = gh.levelmax();
    ilevel = std::min<int>( ilevel, (int)gh.levelmax() );
    
    int nx,ny,nz, ox, oy, oz;
    nx = (int)gh.size(ilevel,0);
    ny = (int)gh.size(ilevel,1);
    nz = (int)gh.size(ilevel,2);
    
    ox = (int)gh.offset_abs(ilevel,0);
    oy = (int)gh.offset_abs(ilevel,1);
    oz = (int)gh.offset_abs(ilevel,2);
    
    real_t usave(u), vsave(v), wsave(w);
    
    u = u*(1ull<<ilevel) - (float)ox;
    v = v*(1ull<<ilevel) - (float)oy;
    w = w*(1ull<<ilevel) - (float)oz;
    
    i = (int)u;
    j = (int)v;
    k = (int)w;
    
    
    if( ilevel>(int)gh.levelmin() && (i<0||i>=(int)gh.size(ilevel,0)-1 ||
       j<0||j>=(int)gh.size(ilevel,1)-1 ||
       k<0||k>=(int)gh.size(ilevel,2)-1) )
        return get_cic( gh, ilevel-1, usave,vsave,wsave);
    
    
    u -= (float)i;
    v -= (float)j;
    w -= (float)k;
    
    int i1,j1,k1;
    i1 = i+1;
    j1 = j+1;
    k1 = k+1;
    
    if( i>= nx ) i = nx-1;
    if( j>= ny ) j = ny-1;
    if( k>= nz ) k = nz-1;
    if( i1>=nx ) i1= nx-1;
    if( j1>=ny ) j1= ny-1;
    if( k1>=nz ) k1= nz-1;
    
    
    float f1,f2,f3,f4,f5,f6,f7,f8;
    
    f1 = (1.f - u) * (1.f - v) * (1.f - w);
    f2 = (1.f - u) * (1.f - v) * (w);
    f3 = (1.f - u) * (v) * (1.f - w);
    f4 = (1.f - u) * (v) * (w);
    f5 = (u) * (1.f - v) * (1.f - w);
    f6 = (u) * (1.f - v) * (w);
    f7 = (u) * (v) * (1.f - w);
    f8 = (u) * (v) * (w);
    
    real_t val = 0.0f;
    
    val += f1*(*gh.get_grid(ilevel))(i,j,k);
    val += f2*(*gh.get_grid(ilevel))(i,j,k1);
    val += f3*(*gh.get_grid(ilevel))(i,j1,k);
    val += f4*(*gh.get_grid(ilevel))(i,j1,k1);
    val += f5*(*gh.get_grid(ilevel))(i1,j,k);
    val += f6*(*gh.get_grid(ilevel))(i1,j,k1);
    val += f7*(*gh.get_grid(ilevel))(i1,j1,k);
    val += f8*(*gh.get_grid(ilevel))(i1,j1,k1);
    
    return val;
}


MyIDType compute_midpoint( const size_t *connect )
{
    
    MyIDType lid1, lid2, newlid, lcoord1[3], lcoord2[3];
    
    idmap_it it;
    
    if( (it=idmap.find( REMOVE_DECORATION_BITS(connect[0]) )) == idmap.end() )
    { lid1 = 0; std::cerr << "1111 -- UWEHFFB YAUDASBOHASDJ" << std::endl; }
    else
        lid1 = connect[0];
    
    if( (it=idmap.find( connect[1] )) == idmap.end() )
    { lid2 = 0; std::cerr << "1111 -- UWEHFFB YAUDASBOHASDJ" << std::endl; }
    else
        lid2 = connect[1];
    
    // take care of periodic boundary conditions
    for( int k=0; k<3; ++k )
    {
        lcoord1[k] = get_lagrange_coord(REMOVE_DECORATION_BITS(lid1),k);
        lcoord2[k] = get_lagrange_coord(REMOVE_DECORATION_BITS(lid2),k);
        
        
        long long dx = lcoord2[k]-lcoord1[k];
        if( dx < -(1ll<<19) ) dx += (1ll<<20);
        if( dx > (1ll<<19)  ) dx -= (1ll<<20);
        
        lcoord1[k] = (lcoord1[k]+(dx>>1)+(1ll<<20)) % (1ll<<20);
    }
    
    newlid = SET_REFINEMENT_COUNTER(get_lagrange_ID(lcoord1[0],lcoord1[1],lcoord1[2]),1);
    return newlid;
}



void split_lagrange_cube( size_t ip )
{
    int k, count;
    
    size_t newi = num_p;//, new_tracer_id = num_p;
    
    const int massc_edge[][2] = {{0,1},{0,2},{1,2},{0,4},{1,4},{2,4},{3,4}};
    const int nmassc = 7; // +1 for the already existing 0
                          //                              8    9    10     11    12    13    14    15    16   17     18   19
    const int massl_edge[][2] = {{1,3},{1,5},{2,3},{2,6},{3,5},{3,6},{3,7},{4,5},{4,6},{5,6},{5,7},{6,7}};
    const int nmassl = 12; // the massless dangling vertices
    
    /*const int new_conn[][8] =
    {
        { 0, 8, 9,10,11,12,13,14},
        { 8, 1,10,15,12,16,14,19},
        { 9,10, 2,17,13,14,18,20},
        {10,15,17, 3,14,19,20,21},
        {11,12,13,14, 4,22,23,24},
        {12,16,14,19,22, 5,24,25},
        {13,14,18,20,23,24, 6,26},
        {14,19,20,21,24,25,26, 7},
    };*/
    
    //////////////////////
    /** spawn new particles **/
    
    //float newmass = P[ip].Mass * 0.125;
    //float zeromass = newmass; // this should eventually be set to zero
    //int masslesstype = 2;
    int this_level = P[ip].Level + 1;
    
    tetgrid_levelmax = std::max(this_level,tetgrid_levelmax);
    
    count = 0;
    
    
    // insert mass carying "true" particles
    for( k=0; k<nmassc; ++k )
    {
        size_t edge_ids[2] = { P[ip].get_vertex( massc_edge[k][0] ), P[ip].get_vertex( massc_edge[k][1] ) };
        
        P[newi+count].Lagrange_ID = compute_midpoint( edge_ids );
        
        assert( REMOVE_DECORATION_BITS(P[newi+count].Lagrange_ID) == REMOVE_DECORATION_BITS(P[ip].get_vertex(k+1,1)));
        
        //P[newi+count].ID = new_tracer_id+count;
        //P[newi+count].Mass = newmass;
        P[newi+count].Type = 1;
        P[newi+count].Level = this_level;
        count++;
    }
    // the innermost particle is a freely moving particle to begin with
    SET_FREE_PARTICLE_BIT( P[newi+count-1].Lagrange_ID );
    
    // insert massless helper particles
    for( k=0; k<nmassl; ++k )
    {
        size_t edge_ids[2] = { P[ip].get_vertex( massl_edge[k][0] ), P[ip].get_vertex( massl_edge[k][1] ) };
        
        //P[newi+count].ID = new_tracer_id+count;
        P[newi+count].Lagrange_ID = compute_midpoint( edge_ids );
        
        //P[newi+count].Mass = zeromass;
        P[newi+count].Type = 2;//masslesstype;
        P[newi+count].Level = this_level;
        count++;
    }
    
    // update spawning particle
    //P[ip].Mass = newmass;
    P[ip].Level = this_level;
    
    
    num_p += count;
}

void init_base_grid( int ilevel, size_t prealloc_particles = 0 )
{
    size_t nbase = 1<<ilevel, lid;
    TetRefinementIDFactor = 20-ilevel;
    
    prealloc_particles = std::max( prealloc_particles, nbase*nbase*nbase );
    P = (particle*)malloc( sizeof(particle) * prealloc_particles );
    num_p_alloc = prealloc_particles;
    
    //double dx = 1.0 / nbase;
    
    num_p = 0;
    
    for( size_t ix=0; ix<nbase; ++ix )
        for( size_t iy=0; iy<nbase; ++iy )
            for( size_t iz=0; iz<nbase; ++iz )
            {
                size_t idx = (ix*nbase+iy)*nbase+iz;
                
                P[idx].Lagrange_ID =
                SET_REFINEMENT_COUNTER((ix<<(TetRefinementIDFactor+40))
                                       +(iy<<(TetRefinementIDFactor+20))
                                       +(iz<<TetRefinementIDFactor),1);
                
                //P[idx].ID = idx;
                P[idx].Type = 1;
                P[idx].Level = ilevel;
                
                SET_FREE_PARTICLE_BIT( P[idx].Lagrange_ID );
                
                //... insert into access map
                lid = REMOVE_FREE_PARTICLE_BIT( REMOVE_REFINEMENT_COUNTER( P[idx].Lagrange_ID ) );
                idmap[ REMOVE_DECORATION_BITS(P[idx].Lagrange_ID) ] = num_p;
                
                num_p++;
            }
    
    tetgrid_levelmax = ilevel;
}



void delete_duplicates( void )
{
    std::sort<particle*>( P, P+num_p );
    
    idmap.clear();
    idmap[ REMOVE_DECORATION_BITS(P[0].Lagrange_ID) ] = 0;
    
    size_t j=1;
    for( size_t i=1; i<num_p; ++i )
    {
        if( REMOVE_DECORATION_BITS(P[i].Lagrange_ID) != REMOVE_DECORATION_BITS(P[i-1].Lagrange_ID) )
        {
            memcpy( &P[j], &P[i], sizeof(particle) );
            idmap[ REMOVE_DECORATION_BITS(P[j].Lagrange_ID) ] = j;
            ++j;
        }
    }
    //size_t ndeleted = num_p-j;
    num_p = j;
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

const int empty_fill_bytes = 60;

template< typename T_store=float >
class gadget_tetmesh_output_plugin : public output_plugin
{
public:
  bool do_baryons_;
  double omegab_;
  double gamma_;
  bool shift_halfcell_;
    
protected:
	
	std::ofstream ofs_;
	bool bmultimass_;
    bool blongids_;
	
	
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
		id_dm_mass, id_dm_vel, id_dm_pos, id_gas_vel, id_gas_rho, id_gas_temp, id_gas_pos, id_dm_conn, id_dm_level, id_dm_lagrangeid
	};
    
    size_t np_type1_, np_type2_, np_type5_;
    
	size_t block_buf_size_;
	size_t npartmax_;
	unsigned nfiles_;
	
	//bool bbndparticles_;
	bool kpcunits_;
    bool msolunits_;
    double YHe_;
    
    refinement_mask refmask;
    
    void distribute_particles( unsigned nfiles, size_t n_t1, size_t n_t2, size_t n_t5,
                              std::vector<unsigned>& nf_t1, std::vector<unsigned>& nf_t2,
                              std::vector<unsigned>& nf_t5 )
    {
        nf_t1.assign( nfiles, 0 );
        nf_t2.assign( nfiles, 0 );
        nf_t5.assign( nfiles, 0 );
        
        
        size_t ntotal = n_t1 + n_t2 + n_t5;
        size_t nnominal = (size_t)((double)ntotal/(double)nfiles);
        
        size_t n_t1_assigned = 0, n_t2_assigned = 0, n_t5_assigned = 0;
        
        for( unsigned i=0; i<nfiles; ++i )
        {
            if( n_t1_assigned < n_t1 )
            {
                nf_t1[i] = std::min( nnominal, n_t1-n_t1_assigned );
                n_t1_assigned += nf_t1[i];
            }
            if( n_t1_assigned == n_t1 )
            {
                nf_t2[i] = std::min( nnominal, n_t2-n_t2_assigned );
                n_t2_assigned += nf_t2[i];
                
                
                if( n_t2_assigned == n_t2 )
                {
                    nf_t5[i] = std::min( nnominal, n_t5-n_t5_assigned );
                    n_t5_assigned += nf_t5[i];
                }
            }
        }
        
        // make sure all particles are assigned
        nf_t1[ 0 ]          += n_t1-n_t1_assigned;
        nf_t2[ nfiles-1 ]   += n_t2-n_t2_assigned;
        nf_t5[ nfiles-1 ]   += n_t5-n_t5_assigned;
        
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
                
		//if( do_baryons_ )
        //    combine_components_for_coarse();
		
        
        
		//............................................................................
		//... copy from the temporary files, interleave the data and save ............
		
		char fnx[256],fny[256],fnz[256],fnvx[256],fnvy[256],fnvz[256],fnm[256];
		char fnc[256], fnl[256], fnlid[256];
        
		sprintf( fnx,  "___ic_temp_%05d.bin", 100*id_dm_pos+0 );
		sprintf( fny,  "___ic_temp_%05d.bin", 100*id_dm_pos+1 );
		sprintf( fnz,  "___ic_temp_%05d.bin", 100*id_dm_pos+2 );
		sprintf( fnvx, "___ic_temp_%05d.bin", 100*id_dm_vel+0 );
		sprintf( fnvy, "___ic_temp_%05d.bin", 100*id_dm_vel+1 );
		sprintf( fnvz, "___ic_temp_%05d.bin", 100*id_dm_vel+2 );
		sprintf( fnm,  "___ic_temp_%05d.bin", 100*id_dm_mass  );

		sprintf( fnc,  "___ic_temp_%05d.bin", 100*id_dm_conn );
		sprintf( fnl,  "___ic_temp_%05d.bin", 100*id_dm_level );
		sprintf( fnlid,  "___ic_temp_%05d.bin", 100*id_dm_lagrangeid );
		
        
        
        
		pistream iffs1, iffs2, iffs3;
		//np_type1_ = np_type1;
        //np_type2_ = np_type2;
        
        
		const size_t 
            nptot = np_type1_+np_type2_+np_type5_;
		
        size_t wrote_p  = 0;
		size_t
			npleft = nptot, 
			n2read = std::min((size_t)block_buf_size_,npleft);
		
		std::cout << " - Gadget2 : writing " << nptot << " particles to file...\n"
			<< "      type   1 : " << std::setw(12) << np_type1_ << "\n"
			<< "      type   2 : " << std::setw(12) << np_type2_ << "\n"
            << "      type   5 : " << std::setw(12) << np_type5_ << "\n";
        
			
				
		std::vector<T_store> adata3;
		adata3.reserve( 3*block_buf_size_ );
		T_store *tmp1, *tmp2, *tmp3;
		
		tmp1 = new T_store[block_buf_size_];
		tmp2 = new T_store[block_buf_size_];
		tmp3 = new T_store[block_buf_size_];
		
		//... for multi-file output
		//int fileno = 0;
		//size_t npart_left = nptot;
        
        std::vector<unsigned> nftype1_per_file, nftype2_per_file, nftype5_per_file;
        distribute_particles( nfiles_, np_type1_, np_type2_, np_type5_, nftype1_per_file, nftype2_per_file, nftype5_per_file );
        
		if( nfiles_ > 1 )
		{
			std::cout << " - Gadget2 : distributing particles to " << nfiles_ << " files\n"
                << "                 " << std::setw(12) << "type 1" << "," << std::setw(12)
                << "type 2" << "," << std::setw(12) << "type 5" << std::endl;
			for( unsigned i=0; i<nfiles_; ++i )
			{
				std::cout << "      file " << std::setw(3) << i << " : " 
				<< std::setw(12) << nftype1_per_file[i] << ","
				<< std::setw(12) << nftype2_per_file[i] << ","
                << std::setw(12) << nftype5_per_file[i] << std::endl;
			}			
		}
        
        size_t curr_block_buf_size = block_buf_size_;
        
        size_t idcount = 0;
        LOGWARN("Need long particle IDs, will write 64bit, make sure to enable in Gadget!");
    	
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
			
            
			size_t np_this_file = nftype1_per_file[ifile] + nftype2_per_file[ifile] + nftype5_per_file[ifile];
			int blksize = sizeof(header);
			
			//... write the header .......................................................
			
			header this_header( header_ );
            this_header.npart[1] = nftype1_per_file[ifile];
            this_header.npart[2] = nftype2_per_file[ifile];
            this_header.npart[5] = nftype5_per_file[ifile];
            
            
			ofs_.write( (char *)&blksize, sizeof(int) );
			ofs_.write( (char *)&this_header, sizeof(header) );
			ofs_.write( (char *)&blksize, sizeof(int) );
			
			
			//... particle positions ..................................................
			blksize = 3ul*np_this_file*sizeof(T_store);
			ofs_.write( (char *)&blksize, sizeof(int) );
			
			npleft = np_this_file;
			n2read = std::min(curr_block_buf_size,npleft);
			
			iffs1.open( fnx, num_p, wrote_p*sizeof(T_store) );
			iffs2.open( fny, num_p, wrote_p*sizeof(T_store) );
			iffs3.open( fnz, num_p, wrote_p*sizeof(T_store) );
			
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
			
			iffs1.open( fnvx, num_p, wrote_p*sizeof(T_store) );
			iffs2.open( fnvy, num_p, wrote_p*sizeof(T_store) );
			iffs3.open( fnvz, num_p, wrote_p*sizeof(T_store) );
			
			npleft = nftype1_per_file[ifile]+nftype2_per_file[ifile];
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
			{
                std::vector<size_t> ids;
                ids.assign(curr_block_buf_size,0);
                npleft	= np_this_file;
                n2read	= std::min(curr_block_buf_size,npleft);
                blksize = sizeof(size_t)*np_this_file;
                
                
                //... generate contiguous IDs and store in file ..
                ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
                while( n2read > 0ul )
                {
                    for( size_t i=0; i<n2read; ++i )
                        ids[i] = idcount++;
                    ofs_.write( reinterpret_cast<char*>(&ids[0]), n2read*sizeof(size_t) );
                    npleft -= n2read;
                    n2read = std::min( curr_block_buf_size,npleft );
                }
                ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
                
                std::vector<size_t>().swap( ids );
			}
			
			//... particle masses .......................................................
			if( true )//bmultimass_ && bmorethan2bnd_ && nc_per_file[ifile] > 0ul)
			{
				iffs1.open( fnm, num_p, wrote_p*sizeof(T_store) );
				
				npleft = nftype1_per_file[ifile]+nftype2_per_file[ifile]+nftype5_per_file[ifile];
                n2read  = std::min(curr_block_buf_size,npleft);
				blksize = np_this_file*sizeof(T_store);
                
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
			
            
            // vert_ids
            {
                std::vector<long long> itemp;
                itemp.assign(curr_block_buf_size,0);
                
				iffs1.open( fnc, 8*num_p, 8*wrote_p*sizeof(size_t) );
				
				npleft = 8*(nftype1_per_file[ifile]+nftype2_per_file[ifile]+nftype5_per_file[ifile]);
                n2read  = std::min(curr_block_buf_size,npleft);
				blksize = 8*np_this_file*sizeof(long long);
                
				ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
				while( n2read > 0ul )
				{
					iffs1.read( reinterpret_cast<char*>(&itemp[0]), n2read*sizeof(long long) );
					ofs_.write( reinterpret_cast<char*>(&itemp[0]), n2read*sizeof(long long) );
					
					npleft -= n2read;
					n2read = std::min( curr_block_buf_size,npleft );
                    
				}
				ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
				
				iffs1.close();
				
			}
            
            // lagrange_id
			{
				std::vector<size_t> itemp;
                itemp.assign(curr_block_buf_size,0);
                iffs1.open( fnlid, num_p, wrote_p*sizeof(size_t) );
				
				npleft = nftype1_per_file[ifile]+nftype2_per_file[ifile]+nftype5_per_file[ifile];
                n2read  = std::min(curr_block_buf_size,npleft);
				blksize = np_this_file*sizeof(size_t);
                
				ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
				while( n2read > 0ul )
				{
					iffs1.read( reinterpret_cast<char*>(&itemp[0]), n2read*sizeof(size_t) );
					ofs_.write( reinterpret_cast<char*>(&itemp[0]), n2read*sizeof(size_t) );
					
					npleft -= n2read;
					n2read = std::min( curr_block_buf_size,npleft );
                    
				}
				ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
				
				iffs1.close();
				
			}
            
			// level
			{
				std::vector<size_t> itemp;
                itemp.assign(curr_block_buf_size,0);
                iffs1.open( fnl, num_p, wrote_p*sizeof(int) );
				
				npleft = nftype1_per_file[ifile]+nftype2_per_file[ifile]+nftype5_per_file[ifile];
                n2read  = std::min(curr_block_buf_size,npleft);
				blksize = np_this_file*sizeof(int);
                
				ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
				while( n2read > 0ul )
				{
					iffs1.read( reinterpret_cast<char*>(&itemp[0]), n2read*sizeof(int) );
					ofs_.write( reinterpret_cast<char*>(&itemp[0]), n2read*sizeof(int) );
					
					npleft -= n2read;
					n2read = std::min( curr_block_buf_size,npleft );
                    
				}
				ofs_.write( reinterpret_cast<char*>(&blksize), sizeof(int) );
				
				iffs1.close();
				
			}
			
			ofs_.flush();
			ofs_.close();
			
            wrote_p += nftype1_per_file[ifile] + nftype2_per_file[ifile] + nftype5_per_file[ifile];
		}
        
        delete[] tmp1;
		delete[] tmp2;
        delete[] tmp3;
        
        remove( fnx );
        remove( fny );
        remove( fnz );
        remove( fnvx );
        remove( fnvy );
        remove( fnvz );
        remove( fnm );
		remove( fnc );
		remove( fnl );
		remove( fnlid );
		
	}
	
	
public:
	
	
	
	gadget_tetmesh_output_plugin( config_file& cf )
	: output_plugin( cf )
	{
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
        /*
		bmorethan2bnd_ = false;
		if( levelmax_ > levelmin_ +1)
			bmorethan2bnd_ = true;

		bmultimass_ = true;
		if( levelmax_ == levelmin_ )
			bmultimass_ = false;
		*/
        //bmorethan2bnd_ = true;
        bmultimass_ = true;
         
		
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
        msolunits_ = cf.getValueSafe<bool>("output","gadget_usemsol",false);
        /*bndparticletype_ = cf.getValueSafe<unsigned>("output","gadget_coarsetype",5);
        
        if( bndparticletype_ == 0 || bndparticletype_ == 1 || bndparticletype_ == 4 ||
           bndparticletype_ > 5 )
        {
            LOGERR("Coarse particles cannot be of Gadget particle type %d in output plugin.", bndparticletype_);
            throw std::runtime_error("Specified illegal Gadget particle type for coarse particles");
        }*/
        
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
		
		if( kpcunits_ )
			header_.BoxSize *= 1000.0;
        
        
		for( int i=0; i<empty_fill_bytes; ++i )
		  header_.fill[i] = 0;
	}
    
    ~gadget_tetmesh_output_plugin()
    {
        free(P);
    }
    
	void write_dm_mass( const grid_hierarchy& gh )
	{
        // generate the refinement hierarchy using the tets
        
        int ibaselvl = gh.levelmin();
        
        init_base_grid( ibaselvl );

        
        for( int ilevel=gh.levelmin(); ilevel<(int)gh.levelmax(); ++ilevel )
        {
            int rfac = 20-ilevel;
            
            MyIDType off[3];
            
            off[0] = gh.offset_abs(ilevel, 0) * (1<<rfac);// + (1<<(rfac-1));
            off[1] = gh.offset_abs(ilevel, 1) * (1<<rfac);// + (1<<(rfac-1));
            off[2] = gh.offset_abs(ilevel, 2) * (1<<rfac);// + (1<<(rfac-1));
            
            
            size_t num_p_before_refinement_sweep = num_p;
            
            for( size_t ip=0; ip<num_p_before_refinement_sweep; ++ip )
            {
                if( P[ip].Level != ilevel || !P[ip].can_refine() )
                    continue;

                MyIDType lidsum = 0, xc[3] = {0,0,0};
                
                bool foundall = true;
                bool dorefine = true;
                
                if( P[ip].Type == 2 ) continue;
                
                for( int i=0; i<8; ++i )
                {
                    idmap_it it = idmap.find( REMOVE_DECORATION_BITS( P[ip].get_vertex(i) ) );
                    
                    if( it==idmap.end() ){
                        foundall = false;
                        LOGERR("This should not happen : Lagrange ID %llu not found!", REMOVE_DECORATION_BITS( P[ip].get_vertex(i) ));
                        throw std::runtime_error("FATAL");
                        break;
                    }
                    
                    
                    
                    size_t cid = (*it).second;
                    MyIDType lid = P[cid].Lagrange_ID;
                    
                    lid = REMOVE_DECORATION_BITS(lid);
                    
                    xc[0] = get_lagrange_coord( lid, 0 );//; - (1<<(rfac-1)); // -dx/2
                    xc[1] = get_lagrange_coord( lid, 1 );// - (1<<(rfac-1)); // -dx/2
                    xc[2] = get_lagrange_coord( lid, 2 );// - (1<<(rfac-1)); // -dx/2
                
                    if( xc[0] < off[0] || xc[1] < off[1] || xc[2] < off[2] ){
                        dorefine = false;
                        break;
                    }
                    
                    int ix, iy, iz;
                    
                    ix = (xc[0]-off[0])>>rfac;//(int)((double)(xc[0] - off[0])*dlm);
                    iy = (xc[1]-off[1])>>rfac;//(int)((double)(xc[1] - off[1])*dlm);
                    iz = (xc[2]-off[2])>>rfac;//(int)((double)(xc[2] - off[2])*dlm);
                    
                    
                    if( ix >= (int)gh.size(ilevel,0) || iy >= (int)gh.size(ilevel,1) || iz >= (int)gh.size(ilevel,2) )
                    {
                        dorefine = false;
                        break;
                    }
                    
                    
                    if( !gh.is_refined(ilevel,ix,iy,iz) )
                    {
                        dorefine = false;
                        break;
                    }
                    
                    
                    
                    lidsum += lid;
                }
                
                if( !foundall ) continue;
                if( !dorefine ) continue;
                
                if( num_p + newnum_p_per_split > num_p_alloc )
                {
                    P = (particle*) realloc( P, (num_p_alloc+=num_p_realloc_blocksize)*sizeof(particle) );
                    LOGINFO("reallocated particle buffer. new size = %llu MBytes.",  num_p_alloc * sizeof(particle)/1024/1024 );
                }
                split_lagrange_cube( ip );
            }
            
            delete_duplicates();
            LOGINFO("refined tet mesh to level %d : now have %lld particles", ilevel+1, num_p );
            
        }
        
        // make all particles that are type 1 but not at highest ref level type 5
        for( size_t ip = 0; ip < num_p; ++ip )
        {
            if( P[ip].Type == 1 && P[ip].Level < tetgrid_levelmax )
                P[ip].Type = 5;
        }
        
        
        /// now we sort all particles by type
        std::sort<particle*>( P, P+num_p, sortP_bytype );

        
        
        
        ////////////////////////////////////////////////////////////////
        
		double rhoc = 27.7519737; // in h^2 1e10 M_sol / Mpc^3
		
		if( kpcunits_ )
            rhoc *= 1e-9; // in h^2 1e10 M_sol / kpc^3
        
		//	rhoc *= 10.0; // in h^2 M_sol / kpc^3
        
        if( msolunits_ )
            rhoc *= 1e10;
        
        header_.mass[1] = 0.0;
        
        idmap.clear();
        
        size_t num_p_t1 = 0, num_p_t2 = 0, num_p_t5 = 0;
        for( size_t ip=0; ip<num_p; ++ip )
        {
            if( P[ip].Type == 1 ) num_p_t1++;
            else if( P[ip].Type == 2 ) num_p_t2++;
            else if( P[ip].Type == 5 ) num_p_t5++;
            
            
            idmap[ REMOVE_DECORATION_BITS( P[ip].Lagrange_ID ) ] = ip;
        }
        
        np_type1_ = num_p_t1;
        np_type2_ = num_p_t2;
        np_type5_ = num_p_t5;
        
        
        
        header_.npart[1] = num_p_t1;
        header_.npartTotal[1] = num_p_t1;
        
        header_.npart[2] = num_p_t2;
        header_.npartTotal[2] = num_p_t2;
        
        header_.npart[5] = num_p_t5;
        header_.npartTotal[5] = num_p_t5;
        
        LOGINFO("   active HR particles  (type 1) : %llu", num_p_t1 );
        LOGINFO("   active LR particles  (type 5) : %llu", num_p_t5 );
        LOGINFO("   passive particles    (type 2) : %llu", num_p_t2 );
        
        
        // write all particle masses
        {
            size_t nwritten = 0;//, ncount = 0;
            std::vector<T_store> temp_dat;
			temp_dat.reserve(block_buf_size_);
            
            char temp_fname[256];
			sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_mass );
			std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
            
            double mfac = header_.Omega0 * rhoc * pow(header_.BoxSize,3.);
            
            size_t blksize = sizeof(T_store)*num_p;
			ofs_temp.write( (char *)&blksize, sizeof(size_t) );
            
            for( size_t ip=0; ip<num_p; ++ip )
            {
                double pmass = mfac / (1ull << (3*P[ip].Level));
                
                //if( P[ip].Type == 2 )
                //    pmass = 0.0;
            
                if( temp_dat.size() < block_buf_size_ )
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
            
            if( nwritten != num_p )
            {
                LOGERR("Internal consistency error while writing temporary file for masses");
                throw std::runtime_error("Internal consistency error while writing temporary file for masses");
            }
            
            ofs_temp.write( (char *)&blksize, sizeof(size_t) );
            
            if( ofs_temp.bad() ){
                LOGERR("I/O error while writing temporary file for masses");
                throw std::runtime_error("I/O error while writing temporary file for masses");
            }
        }
        
        
        
        // write all particle connectivities
        {
            size_t nwritten = 0;
            std::vector<long long> temp_dat;
            
            if( block_buf_size_%8 != 0 )
                block_buf_size_ = (block_buf_size_/8+1)*8;
			temp_dat.reserve(block_buf_size_);
            
            char temp_fname[256];
			sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_conn );
			std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
            
            size_t blksize = sizeof(T_store)*num_p*8;
			ofs_temp.write( (char *)&blksize, sizeof(size_t) );
            
            for( size_t ip=0; ip<num_p; ++ip )
            for( size_t j=0; j<8; ++j )
            {
                size_t lid;
                
                if( P[ip].Type == 2 )
                    lid = -1;
                else
                    lid = (long long)idmap[ REMOVE_DECORATION_BITS(P[ip].get_vertex(j)) ];
                
                if( temp_dat.size() < block_buf_size_ )
                    temp_dat.push_back( lid );
                else
                {
                    ofs_temp.write( (char*)&temp_dat[0], sizeof(long long)*block_buf_size_ );
                    nwritten += block_buf_size_;
                    temp_dat.clear();
                
                    temp_dat.push_back( lid );
                }
            }
            
            if( temp_dat.size() > 0 )
            {
                ofs_temp.write( (char*)&temp_dat[0], sizeof(long long)*temp_dat.size() );
                nwritten+=temp_dat.size();
            }
            
            if( nwritten != 8*num_p )
            {
                LOGERR("Internal consistency error while writing temporary file for connectivities");
                throw std::runtime_error("Internal consistency error while writing temporary file for connectivities");
            }
            
            ofs_temp.write( (char *)&blksize, sizeof(size_t) );
            
            if( ofs_temp.bad() ){
                LOGERR("I/O error while writing temporary file for connectivities");
                throw std::runtime_error("I/O error while writing temporary file for connectivities");
            }
        }
        
        // write all particle levels
        {
            size_t nwritten = 0;
            std::vector<unsigned> temp_dat;
			temp_dat.reserve(block_buf_size_);
            
            char temp_fname[256];
			sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_level );
			std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
            
            size_t blksize = sizeof(T_store)*num_p;
			ofs_temp.write( (char *)&blksize, sizeof(size_t) );
            
            for( size_t ip=0; ip<num_p; ++ip )
            {
                if( temp_dat.size() < block_buf_size_ )
                    temp_dat.push_back( P[ip].Level );
                else
                {
                    ofs_temp.write( (char*)&temp_dat[0], sizeof(T_store)*block_buf_size_ );
                    nwritten += block_buf_size_;
                    temp_dat.clear();
                    temp_dat.push_back( P[ip].Level );
                }
            }
            
            if( temp_dat.size() > 0 )
            {
                ofs_temp.write( (char*)&temp_dat[0], sizeof(T_store)*temp_dat.size() );
                nwritten+=temp_dat.size();
            }
            
            if( nwritten != num_p )
            {
                LOGERR("Internal consistency error while writing temporary file for masses");
                throw std::runtime_error("Internal consistency error while writing temporary file for masses");
            }
            
            ofs_temp.write( (char *)&blksize, sizeof(size_t) );
            
            if( ofs_temp.bad() ){
                LOGERR("I/O error while writing temporary file for masses");
                throw std::runtime_error("I/O error while writing temporary file for masses");
            }
        }
        
        // write all particle lagrange_ids
        {
            size_t nwritten = 0;
            std::vector<size_t> temp_dat;
			temp_dat.reserve(block_buf_size_);
            
            char temp_fname[256];
			sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_lagrangeid );
			std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
            
            size_t blksize = sizeof(T_store)*num_p;
			ofs_temp.write( (char *)&blksize, sizeof(size_t) );
            
            for( size_t ip=0; ip<num_p; ++ip )
            {
                if( temp_dat.size() < block_buf_size_ )
                    temp_dat.push_back( P[ip].Lagrange_ID );
                else
                {
                    ofs_temp.write( (char*)&temp_dat[0], sizeof(T_store)*block_buf_size_ );
                    nwritten += block_buf_size_;
                    temp_dat.clear();
                    temp_dat.push_back( P[ip].Lagrange_ID );
                }
            }
            
            if( temp_dat.size() > 0 )
            {
                ofs_temp.write( (char*)&temp_dat[0], sizeof(T_store)*temp_dat.size() );
                nwritten+=temp_dat.size();
            }
            
            if( nwritten != num_p )
            {
                LOGERR("Internal consistency error while writing temporary file for masses");
                throw std::runtime_error("Internal consistency error while writing temporary file for masses");
            }
            
            ofs_temp.write( (char *)&blksize, sizeof(size_t) );
            
            if( ofs_temp.bad() ){
                LOGERR("I/O error while writing temporary file for masses");
                throw std::runtime_error("I/O error while writing temporary file for masses");
            }
        }
        
	}
	
	
	void write_dm_position( int coord, const grid_hierarchy& gh )
	{
        size_t nwritten = 0;
        std::vector<T_store> temp_dat;
        temp_dat.reserve(block_buf_size_);
        
        double xfac = header_.BoxSize;
        
        char temp_fname[256];
        sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_pos+coord );
        std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
        
        // write all particle masses
        {
            size_t blksize = sizeof(T_store)*num_p;
			ofs_temp.write( (char *)&blksize, sizeof(size_t) );
            
            for( size_t ip = 0; ip < num_p; ++ip )
            {
                double px[3];
                px[0] = get_lagrange_coord( REMOVE_DECORATION_BITS( P[ip].Lagrange_ID ), 0 );
                px[1] = get_lagrange_coord( REMOVE_DECORATION_BITS( P[ip].Lagrange_ID ), 1 );
                px[2] = get_lagrange_coord( REMOVE_DECORATION_BITS( P[ip].Lagrange_ID ), 2 );
                
                px[0] /= (double)(1<<20);
                px[1] /= (double)(1<<20);
                px[2] /= (double)(1<<20);
                
                double dd = get_cic( gh, P[ip].Level, px[0], px[1], px[2] );
                
                if( temp_dat.size() < block_buf_size_ )
                    temp_dat.push_back( (px[coord]+dd)*xfac );
                else
                {
                    ofs_temp.write( (char*)&temp_dat[0], sizeof(T_store)*block_buf_size_ );
                    nwritten += block_buf_size_;
                    temp_dat.clear();
                    temp_dat.push_back( (px[coord]+dd)*xfac );
                }
            }
            
            if( temp_dat.size() > 0 )
            {
                ofs_temp.write( (char*)&temp_dat[0], sizeof(T_store)*temp_dat.size() );
                nwritten+=temp_dat.size();
            }
            
            if( nwritten != num_p )
            {
                LOGERR("Internal consistency error while writing temporary file for masses");
                throw std::runtime_error("Internal consistency error while writing temporary file for masses");
            }
            
            ofs_temp.write( (char *)&blksize, sizeof(size_t) );
            
            if( ofs_temp.bad() ){
                LOGERR("I/O error while writing temporary file for masses");
                throw std::runtime_error("I/O error while writing temporary file for masses");
            }
        }
	}
	
	void write_dm_velocity( int coord, const grid_hierarchy& gh )
	{
        float isqrta = 1.0f/sqrt(header_.time);
		float vfac = isqrta*header_.BoxSize;
		
		if( kpcunits_ )
			vfac /= 1000.0;
		
        size_t nwritten = 0;
        std::vector<T_store> temp_dat;
        temp_dat.reserve(block_buf_size_);
        
        char temp_fname[256];
        sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_vel+coord );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
        
        // write all particle masses
        {
            size_t blksize = sizeof(T_store)*num_p;
			ofs_temp.write( (char *)&blksize, sizeof(size_t) );
            
            for( size_t ip = 0; ip < num_p; ++ip )
            {
                double px[3];
                px[0] = get_lagrange_coord( REMOVE_DECORATION_BITS( P[ip].Lagrange_ID ), 0 );
                px[1] = get_lagrange_coord( REMOVE_DECORATION_BITS( P[ip].Lagrange_ID ), 1 );
                px[2] = get_lagrange_coord( REMOVE_DECORATION_BITS( P[ip].Lagrange_ID ), 2 );
                
                px[0] /= (double)(1<<20);
                px[1] /= (double)(1<<20);
                px[2] /= (double)(1<<20);
                
                double v = get_cic( gh, P[ip].Level, px[0], px[1], px[2] );
                
                if( temp_dat.size() < block_buf_size_ )
                    temp_dat.push_back( v*vfac );
                else
                {
                    ofs_temp.write( (char*)&temp_dat[0], sizeof(T_store)*block_buf_size_ );
                    nwritten += block_buf_size_;
                    temp_dat.clear();
                    temp_dat.push_back( v*vfac );
                }
            }
            
            if( temp_dat.size() > 0 )
            {
                ofs_temp.write( (char*)&temp_dat[0], sizeof(T_store)*temp_dat.size() );
                nwritten+=temp_dat.size();
            }
            
            if( nwritten != num_p )
            {
                LOGERR("Internal consistency error while writing temporary file for velocities");
                throw std::runtime_error("Internal consistency error while writing temporary file for vleocities");
            }
            
            ofs_temp.write( (char *)&blksize, sizeof(size_t) );
            
            if( ofs_temp.bad() ){
                LOGERR("I/O error while writing temporary file for velocities");
                throw std::runtime_error("I/O error while writing temporary file for velocities");
            }
        }
        
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
		this->assemble_gadget_file();
	}
};



namespace{
	output_plugin_creator_concrete< gadget_tetmesh_output_plugin<float> > creator1("gadget_tetmesh");
#ifndef SINGLE_PRECISION
	output_plugin_creator_concrete< gadget_tetmesh_output_plugin<double> > creator2("gadget_tetmesh_double");
#endif
}

