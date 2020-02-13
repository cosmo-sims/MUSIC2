/*
 
 output_grafic2.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include "output.hh"


//! Implementation of class grafic2_output_plugin 
/*!
 This class implements a grafic-2 (cf. Bertschinger 2001) compatible
 output format. With some RAMSES extras.
*/
class grafic2_output_plugin : public output_plugin
{
protected:
	
	
	typedef struct{
		int n1, n2, n3;
		float dxini0;
		float xoff10,xoff20,xoff30;
		float astart0,omega_m0,omega_l0,h00;
		
	}header;
	
	bool bhavehydro_;
  //float metal_floor_;
    int passive_variable_index_;
    float passive_variable_value_;
	
	void write_file_header( std::ofstream& ofs, unsigned ilevel, const grid_hierarchy& gh )
	{
		header loc_head;
		
		double 
			boxlength	= cf_.getValue<double>("setup","boxlength"),
			H0			= cf_.getValue<double>("cosmology","H0"),
			zstart		= cf_.getValue<double>("setup","zstart"),
			astart		= 1.0/(1.0+zstart),
			omegam		= cf_.getValue<double>("cosmology","Omega_m"),
			omegaL		= cf_.getValue<double>("cosmology","Omega_L");
		
		loc_head.n1 = gh.get_grid(ilevel)->size(0);
		loc_head.n2 = gh.get_grid(ilevel)->size(1);
		loc_head.n3 = gh.get_grid(ilevel)->size(2);
		
		loc_head.dxini0 = boxlength / (H0*0.01) / pow(2.0,ilevel);
		
		loc_head.xoff10 = gh.offset_abs(ilevel,0) * loc_head.dxini0;
		loc_head.xoff20 = gh.offset_abs(ilevel,1) * loc_head.dxini0;
		loc_head.xoff30 = gh.offset_abs(ilevel,2) * loc_head.dxini0;
		
		loc_head.astart0 = astart;
		loc_head.omega_m0 = omegam;
		loc_head.omega_l0 = omegaL;
		loc_head.h00 = H0;
		
		
		int blksz = sizeof(header);
		ofs.write( reinterpret_cast<char*> (&blksz), sizeof(int) );
		ofs.write( reinterpret_cast<char*> (&loc_head), blksz );
		ofs.write( reinterpret_cast<char*> (&blksz), sizeof(int) );
		
	}
	
	void write_sliced_array( std::ofstream& ofs, unsigned ilevel, const grid_hierarchy& gh, float fac = 1.0f )
	{
		unsigned n1,n2,n3;
		n1 = gh.get_grid(ilevel)->size(0);
		n2 = gh.get_grid(ilevel)->size(1);
		n3 = gh.get_grid(ilevel)->size(2);
		
		std::vector<float> data(n1*n2,0.0f);
		
		for( unsigned i=0; i<n3; ++i )
		{
			
			data.clear();
			
			for( unsigned j=0; j<n2; ++j )
				for( unsigned k=0; k<n1; ++k )
					data[j*n1+k] = (*gh.get_grid(ilevel))(k,j,i) * fac;
			
			unsigned blksize = n1*n2*sizeof(float);
			
			ofs.write( reinterpret_cast<char*> (&blksize), sizeof(unsigned) );
			ofs.write( reinterpret_cast<char*> (&data[0]), blksize );
			ofs.write( reinterpret_cast<char*> (&blksize), sizeof(unsigned) );
			
		}
	}
    
    size_t restrict_mask( size_t n1, size_t n2, size_t n3, size_t o1, size_t o2, size_t o3,
                        size_t n1c, size_t n2c, size_t n3c, const float* finemask, float* coarsemask )
    {
        //unsigned n1p = n1/2, n2p = n2/2, n3p = n3/2;
        
        for( size_t i=0; i<n1c*n2c*n3c; ++i )
            coarsemask[i] = 0.0f;
        
        for( size_t i=0; i<n1; ++i )
        {
            size_t ii=i/2+o1;
            for( size_t j=0; j<n2; ++j )
            {
                size_t jj=j/2+o2;
                for( size_t k=0; k<n3; ++k )
                {
                    size_t kk=k/2+o3;
                    if( finemask[ (i*n2+j)*n3+k ] )
                        coarsemask[(ii*n2c+jj)*n3c+kk] += 1.0f;
                }
            }
        }
        
        size_t count_ref = 0;
        for( size_t i=0; i<n1c*n2c*n3c; ++i )
            if( coarsemask[i] > 0.1f )
            {
                coarsemask[i] = 1.0f;
                ++count_ref;
            }
        return count_ref;
        
        
    }
    
    void write_refinement_mask( const grid_hierarchy& gh )
    {
        
        // generate mask for highest level
        char ff[256];
       
        size_t n1,n2,n3;
            n1 = gh.get_grid(gh.levelmax())->size(0);
            n2 = gh.get_grid(gh.levelmax())->size(1);
            n3 = gh.get_grid(gh.levelmax())->size(2);
        
        std::vector<float> data(n1*n2*n3,0.0f);
        
        // do finest level
        {
            // get mask for levelmax
            for( size_t i=0; i<n1; ++i )
                for( size_t j=0; j<n2; ++j )
                    for( size_t k=0; k<n3; ++k )
                        if( gh.is_in_mask(gh.levelmax(),i,j,k) && !gh.is_refined(gh.levelmax(),i,j,k) )
                            data[(i*n2+j)*n3+k] = 1.0;
                        else
                            data[(i*n2+j)*n3+k] = 0.0;
            
            // write mask
            sprintf(ff,"%s/level_%03d/ic_refmap",fname_.c_str(), gh.levelmax() );
            std::ofstream ofs(ff,std::ios::binary|std::ios::trunc);
	    write_file_header( ofs, gh.levelmax(), gh );

	    std::ofstream ofs_metals;
	    
	    if( passive_variable_value_ > 0.0f )
	      {
            sprintf(ff,"%s/level_%03d/ic_pvar_%05d",fname_.c_str(), gh.levelmax(), passive_variable_index_ );
            ofs_metals.open(ff,std::ios::binary|std::ios::trunc);
            write_file_header( ofs_metals, gh.levelmax(), gh );
	      }

           
            
            std::vector<float> block(n1*n2,0.0f);
            for( unsigned k=0; k<n3; ++k )
            {
                for( unsigned j=0; j<n2; ++j )
                    for( unsigned i=0; i<n1; ++i )
                        block[j*n1+i] = data[(i*n2+j)*n3+k];
                
                unsigned blksize = n1*n2*sizeof(float);
                
                ofs.write( reinterpret_cast<char*> (&blksize), sizeof(unsigned) );
                ofs.write( reinterpret_cast<char*> (&block[0]), blksize );
                ofs.write( reinterpret_cast<char*> (&blksize), sizeof(unsigned) );


		if( passive_variable_value_ > 0.0f ){

		  for( unsigned j=0; j<n2; ++j )
		    for( unsigned i=0; i<n1; ++i )
		      block[j*n1+i] = data[(i*n2+j)*n3+k] * passive_variable_value_;
                
		  unsigned blksize = n1*n2*sizeof(float);
                
		  ofs_metals.write( reinterpret_cast<char*> (&blksize), sizeof(unsigned) );
		  ofs_metals.write( reinterpret_cast<char*> (&block[0]), blksize );
		  ofs_metals.write( reinterpret_cast<char*> (&blksize), sizeof(unsigned) );
		}
            }
        }
        
        // do all coarser levels
        for( unsigned ilevel=levelmax_-1; ilevel>=levelmin_; --ilevel )
        {
            size_t n1c,n2c,n3c,o1,o2,o3;
            n1c = gh.get_grid(ilevel)->size(0);
            n2c = gh.get_grid(ilevel)->size(1);
            n3c = gh.get_grid(ilevel)->size(2);
            
            n1 = gh.get_grid(ilevel+1)->size(0);
            n2 = gh.get_grid(ilevel+1)->size(1);
            n3 = gh.get_grid(ilevel+1)->size(2);
            
            o1 = gh.get_grid(ilevel+1)->offset(0);
            o2 = gh.get_grid(ilevel+1)->offset(1);
            o3 = gh.get_grid(ilevel+1)->offset(2);
            
            std::vector<float> data_coarse( n1c*n2c*n3c, 0.0f );
            
            /*if( ilevel <= levelmax_-2 )
            {
                for( size_t i=0; i<n1*n2*n3; ++i )
                    data[i] = 0.0;
                
                for( unsigned i=2; i<n1-2; ++i )
                    for( unsigned j=2; j<n2-2; ++j )
                        for( unsigned k=2; k<n3-2; ++k )
                            data[(i*n2+j)*n3+k] = 1.0;
			    }*/
            
            size_t nref;
            nref = restrict_mask( n1, n2, n3, o1, o2, o3, n1c, n2c, n3c, &data[0], &data_coarse[0] );
            
            LOGINFO("%f of cells on level %d are refined",(double)nref/(n1c*n2c*n3c),ilevel);
            
            sprintf(ff,"%s/level_%03d/ic_refmap",fname_.c_str(), ilevel );
            std::ofstream ofs(ff,std::ios::binary|std::ios::trunc);
            write_file_header( ofs, ilevel, gh );

	    std::ofstream ofs_metals;
	    if( passive_variable_value_ > 0.0f )
	      {
		sprintf(ff,"%s/level_%03d/ic_pvar_%05d",fname_.c_str(), ilevel, passive_variable_index_ );
		ofs_metals.open(ff,std::ios::binary|std::ios::trunc);
		write_file_header( ofs_metals, ilevel, gh );
	      }


            std::vector<float> block(n1c*n2c,0.0f);
            for( unsigned i=0; i<n3c; ++i )
            {
                for( unsigned j=0; j<n2c; ++j )
                    for( unsigned k=0; k<n1c; ++k )
                        block[j*n1c+k] = data_coarse[(k*n2c+j)*n3c+i];
                
                unsigned blksize = n1c*n2c*sizeof(float);
                
                ofs.write( reinterpret_cast<char*> (&blksize), sizeof(unsigned) );
                ofs.write( reinterpret_cast<char*> (&block[0]), blksize );
                ofs.write( reinterpret_cast<char*> (&blksize), sizeof(unsigned) );

		if( passive_variable_value_ > 0.0f ){

		  for( unsigned j=0; j<n2c; ++j )
                    for( unsigned k=0; k<n1c; ++k )
                        block[j*n1c+k] = data_coarse[(k*n2c+j)*n3c+i] * passive_variable_value_;
                
		  unsigned blksize = n1c*n2c*sizeof(float);
                
		  ofs_metals.write( reinterpret_cast<char*> (&blksize), sizeof(unsigned) );
		  ofs_metals.write( reinterpret_cast<char*> (&block[0]), blksize );
		  ofs_metals.write( reinterpret_cast<char*> (&blksize), sizeof(unsigned) );
		}
            }
            
            data.swap( data_coarse );
            
        }
    }
    
    void write_ramses_namelist( const grid_hierarchy& gh )
	{
		//... also write the refinement options to a dummy namelist file
		char ff[256];
		sprintf(ff,"%s/ramses.nml",fname_.c_str() );
		
		std::ofstream ofst(ff,std::ios::trunc);
      
        // -- RUN_PARAMS -- //
        ofst
          << "&RUN_PARAMS\n"
          << "cosmo=.true.\n"
          << "pic=.true.\n"
          << "poisson=.true.\n";
      
        if( bhavehydro_ )
          ofst << "hydro=.true.\n";
        else
          ofst << "hydro=.false.\n";
          
        ofst
          << "nrestart=0\n"
          << "nremap=1\n"
          << "nsubcycle=";
      
        for( unsigned ilevel=gh.levelmin(); ilevel<=gh.levelmax(); ++ilevel )
          ofst << "1,";
        ofst << "1,2\n";
      
        ofst
          << "ncontrol=1\n"
          << "verbose=.false.\n/\n\n";
		
        // -- INIT_PARAMS -- //
        ofst
          << "&INIT_PARAMS\n"
          << "filetype=\'grafic\'\n";
		for( unsigned i=gh.levelmin();i<=gh.levelmax(); ++i)
		{
			sprintf(ff,"initfile(%d)=\'%s/level_%03d\'\n",i-gh.levelmin()+1,fname_.c_str(), i );
			ofst << std::string(ff);
		}
		ofst << "/\n\n";
		
		
        unsigned naddref = 8; // initialize with settings for 10 additional levels of refinement
        unsigned nexpand = (cf_.getValue<unsigned>("setup","padding")-1)/2;
        
        // -- AMR_PARAMS -- //
        ofst << "&AMR_PARAMS\n"
            << "levelmin=" << gh.levelmin() << "\n"
            << "levelmax=" << gh.levelmax()+naddref << "\n"
            << "nexpand=";
      
        if( gh.levelmax() == gh.levelmin() )
          ofst << "1";
        else 
        {
          for( unsigned ilevel=gh.levelmin(); ilevel<gh.levelmax()-1; ++ilevel )
            ofst << nexpand << ",";
          ofst << "1,1";
          
        }
      
        ofst << "\n"
             << "ngridtot=2000000\n"
             << "nparttot=3000000\n"
             << "/\n\n";
      
        ofst << "&REFINE_PARAMS\n"
            << "m_refine=" << gh.levelmax()-gh.levelmin()+1+naddref << "*8.,\n";
      
        if( bhavehydro_ )
          ofst << "ivar_refine=" << 5+passive_variable_index_ << "\n"
               << "var_cut_refine=" << passive_variable_value_*0.01 << "\n";
        else
          ofst << "ivar_refine=0\n";
      
        ofst << "mass_cut_refine=" << 2.0/pow(2,3*gh.levelmax()) << "\n"
             << "interpol_var=1\n"
             << "interpol_type=0\n"
             << "/\n\n";
        
        
		LOGINFO("The grafic2 output plug-in wrote the grid data to a partial");
		LOGINFO("   RAMSES namelist file \'%s\'",fname_.c_str() );
    }
	
	void write_ramses_namelist_old( const grid_hierarchy& gh )
	{
		//... also write the refinement options to a dummy namelist file
		char ff[256];
		sprintf(ff,"%s/ramses.nml",fname_.c_str() );
		
		std::ofstream ofst(ff,std::ios::trunc);
		
		ofst 
		<< "&INIT_PARAMS\n"
		<< "filetype=\'grafic\'\n";
		for( unsigned i=gh.levelmin();i<=gh.levelmax(); ++i)
		{
			sprintf(ff,"initfile(%d)=\'%s/level_%03d\'\n",i-gh.levelmin()+1,fname_.c_str(), i );
			ofst << std::string(ff);
		}
		ofst << "/\n\n";
		
		
		double xc,yc,zc,l;
		
		ofst
		<< "&AMR_PARAMS\n"
		<< "levelmin=" << gh.levelmin() << "\n"
		<< "levelmax=" << gh.levelmax() << "\n"
		<< "ngridtot=2000000\n"
		<< "nparttot=3000000\n"
		<< "nexpand=1\n/\n\n";
		
		const size_t fprec = 12, fwid = 16;
		
		if( gh.levelmax() > gh.levelmin() )
		{
			l = (double)(1l<<(gh.levelmin()+1));
			xc = ((double)gh.offset_abs(gh.levelmin()+1,0)+0.5*(double)gh.size(gh.levelmin()+1,0))/l;
			yc = ((double)gh.offset_abs(gh.levelmin()+1,1)+0.5*(double)gh.size(gh.levelmin()+1,1))/l;
			zc = ((double)gh.offset_abs(gh.levelmin()+1,2)+0.5*(double)gh.size(gh.levelmin()+1,2))/l;	
		
			ofst << "&REFINE_PARAMS\n"
			<< "m_refine=  "<< std::setw(fwid) << std::setprecision(fprec) << 0.0;
			
			
			for( unsigned i=gh.levelmin()+1;i<gh.levelmax(); ++i)
				ofst << "," << std::setw(fwid) << std::setprecision(fprec) << 0.0;
			ofst << "\nx_refine=  "<< std::setw(fwid) << std::setprecision(fprec) << xc;
			for( unsigned i=gh.levelmin()+1;i<gh.levelmax(); ++i)
			{	
				l = (double)(1l<<(i+1));
				xc = ((double)gh.offset_abs(i+1,0)+0.5*(double)gh.size(i+1,0))/l;
				ofst << ","<< std::setw(fwid) << std::setprecision(fprec) << xc;
			}
			ofst << "\ny_refine=  "<< std::setw(fwid) << std::setprecision(fprec) << yc;
			for( unsigned i=gh.levelmin()+1;i<gh.levelmax(); ++i)
			{	
				l = (double)(1l<<(i+1));
				yc = ((double)gh.offset_abs(i+1,1)+0.5*(double)gh.size(i+1,1))/l;
				ofst << ","<< std::setw(fwid) << std::setprecision(fprec) << yc;
			}
			ofst << "\nz_refine=  "<< std::setw(fwid) << std::setprecision(fprec) << zc;
			for( unsigned i=gh.levelmin()+1;i<gh.levelmax(); ++i)
			{	
				l = (double)(1l<<(i+1));
				zc = ((double)gh.offset_abs(i+1,2)+0.5*(double)gh.size(i+1,2))/l;
				ofst << ","<< std::setw(fwid) << std::setprecision(fprec) << zc;
			}
			
			ofst << "\nr_refine=  ";
			for(unsigned i=gh.levelmin();i<gh.levelmax(); ++i )
			{
				size_t nmax = std::min(gh.size(i+1,0),std::min(gh.size(i+1,1),gh.size(i+1,2)));
				
				double r = (nmax-4.0)/(double)(1l<<(i+1));
				if( i==gh.levelmin() )
					ofst << std::setw(fwid) << std::setprecision(fprec) << r;
				else
					ofst << "," << std::setw(fwid) << std::setprecision(fprec) << r;
			}
			ofst << "\nexp_refine=" << std::setw(fwid) << std::setprecision(fprec) << 10.0;
			for( unsigned i=gh.levelmin()+1;i<gh.levelmax(); ++i)
				ofst << "," << std::setw(fwid) << std::setprecision(fprec) << 10.0;
			ofst << "\n/\n";
		}
		
		sprintf(ff,"%s/ramses.nml",fname_.c_str() );
		std::cout	<< " - The grafic2 output plug-in wrote the grid data to a partial\n"
		<< "   RAMSES namelist file \'" << ff << "\'\n"; 
	}
	
public:
	
	grafic2_output_plugin( config_file& cf )
	: output_plugin( cf )
	{
		// create directory structure
		remove( fname_.c_str() );
		mkdir( fname_.c_str(), 0777 );
		for(unsigned ilevel=levelmin_; ilevel<=levelmax_; ++ilevel )
		{
			char fp[256];
			sprintf(fp,"%s/level_%03d",fname_.c_str(), ilevel );
			mkdir( fp, 0777 );
		}
		
		
		bhavehydro_ = cf.getValue<bool>("setup","baryons");
      //metal_floor_ = cf.getValueSafe<float>("output","ramses_metal_floor",1e-5);
        passive_variable_index_ = cf.getValueSafe<int>("output","ramses_pvar_idx",1);
        passive_variable_value_ = cf.getValueSafe<float>("output","ramses_pvar_val",1.0f);
	}
	
	/*~grafic2_output_plugin()
	 { }*/
	
	
	void write_dm_position( int coord, const grid_hierarchy& gh  )
	{
		double 
		boxlength	= cf_.getValue<double>("setup","boxlength");
		
		for(unsigned ilevel=levelmin_; ilevel<=levelmax_; ++ilevel )
		{
			
			char ff[256];
			sprintf(ff,"%s/level_%03d/ic_posc%c",fname_.c_str(), ilevel, (char)('x'+coord) );
			
			std::ofstream ofs(ff,std::ios::binary|std::ios::trunc);
			
			write_file_header( ofs, ilevel, gh );
			write_sliced_array( ofs, ilevel, gh, boxlength );
		}
	}
	
	void write_dm_velocity( int coord, const grid_hierarchy& gh )
	{
		double 
		boxlength	= cf_.getValue<double>("setup","boxlength");
		
		for(unsigned ilevel=levelmin_; ilevel<=levelmax_; ++ilevel )
		{
			
			char ff[256];
			sprintf(ff,"%s/level_%03d/ic_velc%c",fname_.c_str(), ilevel, (char)('x'+coord) );
			
			std::ofstream ofs(ff,std::ios::binary|std::ios::trunc);
			
			write_file_header( ofs, ilevel, gh );
			write_sliced_array( ofs, ilevel, gh, boxlength );
		}
	}
	
	void write_gas_velocity( int coord, const grid_hierarchy& gh )
	{	
		double 
		boxlength	= cf_.getValue<double>("setup","boxlength");
		
		for(unsigned ilevel=levelmin_; ilevel<=levelmax_; ++ilevel )
		{
			
			char ff[256];
			sprintf(ff,"%s/level_%03d/ic_velb%c",fname_.c_str(), ilevel, (char)('x'+coord) );
			
			std::ofstream ofs(ff,std::ios::binary|std::ios::trunc);
			
			write_file_header( ofs, ilevel, gh );
			write_sliced_array( ofs, ilevel, gh, boxlength );
		}
	}
	
	void write_gas_density( const grid_hierarchy& gh )
	{	
		for(unsigned ilevel=levelmin_; ilevel<=levelmax_; ++ilevel )
		{
			
			char ff[256];
			sprintf(ff,"%s/level_%03d/ic_deltab",fname_.c_str(), ilevel );
			
			std::ofstream ofs(ff,std::ios::binary|std::ios::trunc);
			
			write_file_header( ofs, ilevel, gh );
			write_sliced_array( ofs, ilevel, gh );
		}
		
	}
	
	
	void write_dm_density( const grid_hierarchy& gh )
	{	
		if(! bhavehydro_ )
			write_gas_density(gh);
		
		if( cf_.getValueSafe<bool>("output","ramses_nml",true) )
			write_ramses_namelist(gh);
        else if( cf_.getValueSafe<bool>("output","ramses_old_nml",false) )
			write_ramses_namelist_old(gh);
      
        if( gh.levelmin() != gh.levelmax() )
          write_refinement_mask( gh );
	}
	
	void write_dm_mass( const grid_hierarchy& gh )
	{	/* do nothing, not used... */ }
	
	void write_dm_potential( const grid_hierarchy& gh )
	{	/* do nothing, not used... */ }
	
	void write_gas_potential( const grid_hierarchy& gh )
	{	/* do nothing, not used... */ }
	
	void write_gas_position( int coord, const grid_hierarchy& gh )
	{	/* do nothing, not used... */ }
	
	void finalize( void )
	{	}
	
};

namespace{
	output_plugin_creator_concrete<grafic2_output_plugin> creator("grafic2");
}

