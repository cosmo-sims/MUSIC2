/*
 
 output_grafic2.cc - This file is part of MUSIC -
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
		
		for( unsigned i=0; i<n1; ++i )
		{
			
			data.clear();
			
			for( unsigned j=0; j<n2; ++j )
				for( unsigned k=0; k<n3; ++k )
					data[j*n3+k] = (*gh.get_grid(ilevel))(k,j,i) * fac;
			
			unsigned blksize = n2*n3*sizeof(float);
			
			ofs.write( reinterpret_cast<char*> (&blksize), sizeof(unsigned) );
			ofs.write( reinterpret_cast<char*> (&data[0]), blksize );
			ofs.write( reinterpret_cast<char*> (&blksize), sizeof(unsigned) );
			
		}
	}
	
	void write_ramses_namelist( const grid_hierarchy& gh )
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
		
		
		if( gh.levelmax() > gh.levelmin() )
		{
			l = pow(2.0,gh.levelmin()+1);
			xc = ((double)gh.offset_abs(gh.levelmin()+1,0)+0.5*(double)gh.size(gh.levelmin()+1,0))/l;
			yc = ((double)gh.offset_abs(gh.levelmin()+1,1)+0.5*(double)gh.size(gh.levelmin()+1,1))/l;
			zc = ((double)gh.offset_abs(gh.levelmin()+1,2)+0.5*(double)gh.size(gh.levelmin()+1,2))/l;	
		
			ofst << "&REFINE_PARAMS\n"
				<< "m_refine=0.";
			
			
			for( unsigned i=gh.levelmin()+1;i<gh.levelmax(); ++i)
				ofst << ",0.";
			ofst << "\nx_refine="<< std::setw(12) << xc;
			for( unsigned i=gh.levelmin()+1;i<gh.levelmax(); ++i)
			{	
				l = pow(2.0,i+1);
				xc = ((double)gh.offset_abs(i+1,0)+0.5*(double)gh.size(i+1,0))/l;
				ofst << ","<< std::setw(12) << xc;
			}
			ofst << "\ny_refine="<< std::setw(12)<< yc;
			for( unsigned i=gh.levelmin()+1;i<gh.levelmax(); ++i)
			{	
				l = pow(2.0,i+1);
				yc = ((double)gh.offset_abs(i+1,1)+0.5*(double)gh.size(i+1,1))/l;
				ofst << ","<< std::setw(12) << yc;
			}
			ofst << "\nz_refine="<< std::setw(12) << zc;
			for( unsigned i=gh.levelmin()+1;i<gh.levelmax(); ++i)
			{	
				l = pow(2.0,i+1);
				zc = ((double)gh.offset_abs(i+1,2)+0.5*(double)gh.size(i+1,2))/l;
				ofst << ","<< std::setw(12) << zc;
			}
			
			ofst << "\nr_refine=";
			for(unsigned i=gh.levelmin();i<gh.levelmax(); ++i )
			{
				double r = (gh.size(i+1,0)-4.0)/pow(2.0,i+1);
				if( i==gh.levelmin() )
					ofst << std::setw(12) << r;
				else
					ofst << "," << std::setw(12) << r;
			}
			ofst << "\nexp_refine=10.0";
			for( unsigned i=gh.levelmin()+1;i<gh.levelmax(); ++i)
				ofst << ",10.0";
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

