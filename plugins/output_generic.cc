/*
 
 output_generic.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010-13  Oliver Hahn
 
 */


#ifdef HAVE_HDF5

#include "output.hh"
#include "HDF_IO.hh"


class generic_output_plugin : public output_plugin
{
protected:
	
	using output_plugin::cf_;
		
	template< typename Tt >
	void write2HDF5( std::string fname, std::string dname, const MeshvarBnd<Tt>& data )
	{
		int n0 = data.size(0), n1 = data.size(1), n2 = data.size(2), nb = data.m_nbnd;
		std::vector<Tt> vdata;
		vdata.reserve((unsigned)(n0+2*nb)*(n1+2*nb)*(n2+2*nb));
		for(int i=-nb; i<n0+nb; ++i )
			for(int j=-nb; j<n1+nb; ++j )
				for(int k=-nb; k<n2+nb; ++k )
					vdata.push_back( data(i,j,k) );
		
		unsigned nd[3] = { (unsigned)(n0+2*nb),(unsigned)(n1+2*nb),(unsigned)(n2+2*nb)	};
		HDFWriteDataset3D( fname, dname, nd, vdata);
	}
	
public:
	generic_output_plugin( config_file& cf )//std::string fname, Cosmology cosm, Parameters param )
	: output_plugin( cf )//fname, cosm, param )
	{

		HDFCreateFile(fname_);
		
		HDFCreateGroup(fname_, "header");

		HDFWriteDataset(fname_,"/header/grid_off_x",offx_);
		HDFWriteDataset(fname_,"/header/grid_off_y",offy_);
		HDFWriteDataset(fname_,"/header/grid_off_z",offz_);
		
		HDFWriteDataset(fname_,"/header/grid_len_x",sizex_);
		HDFWriteDataset(fname_,"/header/grid_len_y",sizey_);
		HDFWriteDataset(fname_,"/header/grid_len_z",sizez_);
		
		HDFWriteGroupAttribute(fname_, "header", "levelmin", levelmin_ );
		HDFWriteGroupAttribute(fname_, "header", "levelmax", levelmax_ );
	}
	
	~generic_output_plugin()
	{	}
	
	void write_dm_mass( const grid_hierarchy& gh )
	{	}
	
	void write_dm_velocity( int coord, const grid_hierarchy& gh )
	{
		char sstr[128];
		
		for( unsigned ilevel=0; ilevel<=levelmax_; ++ilevel )
		{
			if( coord == 0 )
				sprintf(sstr,"level_%03d_DM_vx",ilevel);
			else if( coord == 1 )
				sprintf(sstr,"level_%03d_DM_vy",ilevel);
			else if( coord == 2 )
				sprintf(sstr,"level_%03d_DM_vz",ilevel);
			
			write2HDF5( fname_, sstr, *gh.get_grid(ilevel) );
		}
	}
	
	void write_dm_position( int coord, const grid_hierarchy& gh )
	{
		char sstr[128];
		
		for( unsigned ilevel=0; ilevel<=levelmax_; ++ilevel )
		{
			if( coord == 0 )
				sprintf(sstr,"level_%03d_DM_dx",ilevel);
			else if( coord == 1 )
				sprintf(sstr,"level_%03d_DM_dy",ilevel);
			else if( coord == 2 )
				sprintf(sstr,"level_%03d_DM_dz",ilevel);
			
			write2HDF5( fname_, sstr, *gh.get_grid(ilevel) );
		}
	}
	
	void write_dm_density( const grid_hierarchy& gh )
	{
		char sstr[128];
		
		for( unsigned ilevel=0; ilevel<=levelmax_; ++ilevel )
		{
			sprintf(sstr,"level_%03d_DM_rho",ilevel);
			write2HDF5( fname_, sstr, *gh.get_grid(ilevel) );
		}



		// determine density maximum and minimum location
		real_t rhomax = -1e30, rhomin = 1e30;
		double loc_rhomax[3] = {0.0,0.0,0.0}, loc_rhomin[3] = {0.0,0.0,0.0};
		int lvl_rhomax = 0, lvl_rhomin = 0;
		real_t rhomax_lm = -1e30, rhomin_lm = 1e30;
		double loc_rhomax_lm[3] = {0.0,0.0,0.0}, loc_rhomin_lm[3] = {0.0,0.0,0.0};
        
        
		for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
		  for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
		    for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
		      for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
			if( ! gh.is_refined(ilevel,i,j,k) )
			  {
			    real_t rho = (*gh.get_grid(ilevel))(i,j,k);
                            
                            if( rho > rhomax )
			      {
                                rhomax = rho;
                                lvl_rhomax = ilevel;
                                gh.cell_pos(ilevel, i, j, k, loc_rhomax);
			      }
                            
                            if( rho < rhomin )
			      {
                                rhomin = rho;
                                lvl_rhomin = ilevel;
                                gh.cell_pos(ilevel, i, j, k, loc_rhomin);
			      }
                            
                            if( ilevel == (int)gh.levelmax() )
			      {
                                if( rho > rhomax_lm )
				  {
                                    rhomax_lm = rho;
                                    gh.cell_pos(ilevel, i, j, k, loc_rhomax_lm);
				  }
                                
                                if( rho < rhomin_lm )
				  {
                                    rhomin_lm = rho;
                                    gh.cell_pos(ilevel, i, j, k, loc_rhomin_lm);
				  }
			      }
			  }

		double h = 1.0/(1<<levelmin_);
		double shift[3];
		shift[0] = -(double)cf_.getValue<int>( "setup", "shift_x" )*h;
		shift[1] = -(double)cf_.getValue<int>( "setup", "shift_y" )*h;
		shift[2] = -(double)cf_.getValue<int>( "setup", "shift_z" )*h;
			
		if( gh.levelmin() != gh.levelmax() )
		  {
		    LOGINFO("Global density extrema: ");
		    LOGINFO("  minimum: delta=%f at (%f,%f,%f) (level=%d)",rhomin,loc_rhomin[0],loc_rhomin[1],loc_rhomin[2],lvl_rhomin);
		    LOGINFO("       shifted back at (%f,%f,%f)",loc_rhomin[0]+shift[0],loc_rhomin[1]+shift[1],loc_rhomin[2]+shift[2]);
		    LOGINFO("  maximum: delta=%f at (%f,%f,%f) (level=%d)",rhomax,loc_rhomax[0],loc_rhomax[1],loc_rhomax[2],lvl_rhomax);
		    LOGINFO("       shifted back at (%f,%f,%f)",loc_rhomax[0]+shift[0],loc_rhomax[1]+shift[1],loc_rhomax[2]+shift[2]);
		    
		    LOGINFO("Density extrema on finest level: ");
		    LOGINFO("  minimum: delta=%f at (%f,%f,%f)",rhomin_lm,loc_rhomin_lm[0],loc_rhomin_lm[1],loc_rhomin_lm[2]);
		    LOGINFO("       shifted back at (%f,%f,%f)",loc_rhomin_lm[0]+shift[0],loc_rhomin_lm[1]+shift[1],loc_rhomin_lm[2]+shift[2]);
		    LOGINFO("  maximum: delta=%f at (%f,%f,%f)",rhomax_lm,loc_rhomax_lm[0],loc_rhomax_lm[1],loc_rhomax_lm[2]);
		    LOGINFO("       shifted back at (%f,%f,%f)",loc_rhomax_lm[0]+shift[0],loc_rhomax_lm[1]+shift[1],loc_rhomax_lm[2]+shift[2]);
		    
		  }else{
		  LOGINFO("Global density extrema: ");
		  LOGINFO("  minimum: delta=%f at (%f,%f,%f)",rhomin,loc_rhomin[0],loc_rhomin[1],loc_rhomin[2]);
		  LOGINFO("       shifted back at (%f,%f,%f)",loc_rhomin[0]+shift[0],loc_rhomin[1]+shift[1],loc_rhomin[2]+shift[2]);
		  LOGINFO("  maximum: delta=%f at (%f,%f,%f)",rhomax,loc_rhomax[0],loc_rhomax[1],loc_rhomax[2]);
		  LOGINFO("       shifted back at (%f,%f,%f)",loc_rhomax[0]+shift[0],loc_rhomax[1]+shift[1],loc_rhomax[2]+shift[2]);
            
		}
	}
	
	void write_dm_potential( const grid_hierarchy& gh )
	{ 
		char sstr[128];
		
		for( unsigned ilevel=0; ilevel<=levelmax_; ++ilevel )
		{
			sprintf(sstr,"level_%03d_DM_potential",ilevel);
			write2HDF5( fname_, sstr, *gh.get_grid(ilevel) );
		}
	}
	
	void write_gas_potential( const grid_hierarchy& gh )
	{ 
		char sstr[128];
		
		for( unsigned ilevel=0; ilevel<=levelmax_; ++ilevel )
		{
			sprintf(sstr,"level_%03d_BA_potential",ilevel);
			write2HDF5( fname_, sstr, *gh.get_grid(ilevel) );
		}
	}
	
	
	
	void write_gas_velocity( int coord, const grid_hierarchy& gh )
	{	
		char sstr[128];
		
		for( unsigned ilevel=0; ilevel<=levelmax_; ++ilevel )
		{
			if( coord == 0 )
				sprintf(sstr,"level_%03d_BA_vx",ilevel);
			else if( coord == 1 )
				sprintf(sstr,"level_%03d_BA_vy",ilevel);
			else if( coord == 2 )
				sprintf(sstr,"level_%03d_BA_vz",ilevel);
			
			write2HDF5( fname_, sstr, *gh.get_grid(ilevel) );
		}
	}
	
	void write_gas_position( int coord, const grid_hierarchy& gh )
	{	}
	
	void write_gas_density( const grid_hierarchy& gh )
	{	
		char sstr[128];
		
		for( unsigned ilevel=0; ilevel<=levelmax_; ++ilevel )
		{
			sprintf(sstr,"level_%03d_BA_rho",ilevel);
			write2HDF5( fname_, sstr, *gh.get_grid(ilevel) );
		}
	}
	
	void finalize( void )
	{	}
};



namespace{
	output_plugin_creator_concrete< generic_output_plugin > creator("generic");
}


#endif

