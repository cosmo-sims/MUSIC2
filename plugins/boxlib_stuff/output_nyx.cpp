/*
 
 output_nyx.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 Copyright (C) 2012  Jan Frederik Engels
 
 */

#ifdef HAVE_BOXLIB


#include "../../output.hh"

#include <VisMF.H>
#include <Box.H>
#include <RealBox.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <PArray.H>

#define MAX_GRID_SIZE	32
#define BL_SPACEDIM 3


class nyx_output_plugin : public output_plugin
{
protected:
	
	struct patch_header{
		int component_rank;
		size_t component_size;
		std::vector<int> dimensions;
		int rank;
		std::vector<int> top_grid_dims;
		std::vector<int> top_grid_end;
		std::vector<int> top_grid_start;
	};
	
	struct sim_header{
		std::vector<int> dimensions;
		std::vector<int> offset;
		float a_start;
		float dx;
		float h0;
		float omega_b;
		float omega_m;
		float omega_v;
		float vfact;
		float boxlength;
		int   particle_idx;
	};

//	struct grid_on_one_level{
//		IntVect lo;
//		IntVect hi;
//	};

	int n_data_items;
	std::vector<std::string> field_name;
	int f_lev;
	int gridp;
	
	PArray<MultiFab> mfs;

//	std::vector<grid_on_one_level> grids;

	std::vector<BoxArray> boxarrays;

	sim_header the_sim_header;

	
	
	void dump_grid_data(int comp, std::string fieldname, const grid_hierarchy& gh, double factor = 1.0, double add = 0.0 )
	{
		std::cout << fieldname << " is dumped... to mf index " << comp << std::endl;

		//FIXME adapt for multiple levels!
		for(int mlevel=levelmin_; mlevel<=levelmax_; ++mlevel )
		{
			int blevel = mlevel-levelmin_;

			std::vector<int> ng;
			ng.push_back( gh.get_grid(mlevel)->size(0) );
			ng.push_back( gh.get_grid(mlevel)->size(1) );
			ng.push_back( gh.get_grid(mlevel)->size(2) );
			
			std::cout << ng[0] << " " << ng[1] << " " << ng[2] << std::endl;
	
			//write data to mf
			for(MFIter mfi(mfs[blevel]); mfi.isValid(); ++mfi) {
			    FArrayBox &myFab = mfs[blevel][mfi];	
			    const int  *fab_lo = mfi.validbox().loVect();
			    const int  *fab_hi = mfi.validbox().hiVect();
	
//			    int mk = 0;
	#ifdef OMP		    
			    #pragma omp parallel for default(shared)
	#endif
			    for (int k = fab_lo[2], mk=0; k <= fab_hi[2]; k++, mk++) {
//			      int mj = 0;
			      for (int j = fab_lo[1], mj=0; j <= fab_hi[1]; j++, mj++) {
//				int mi = 0;
			    	for (int i = fab_lo[0], mi=0; i <= fab_hi[0]; i++, mi++) {
				  if (mi>=ng[0])
					  std::cout << "mi (" << mi << ") too large " << ng[0] << std::endl;
				  if (mj>=ng[1])
					  std::cout << "mj (" << mj << ") too large " << ng[1] << std::endl;
				  if (mk>=ng[2])
					  std::cout << "mk (" << mk << ") too large " << ng[2] << std::endl;
	
				  IntVect iv(i,j,k);
				  double data = ( add + (*gh.get_grid(mlevel))(mi,mj,mk) )*factor;
				  int idx = myFab.box().index(iv);
				  myFab.dataPtr(comp)[idx] = data;

//				  mi++;
			    	}
//				mj++;
			      }
//			      mk++;
			    }
			} // MFI
		}
		
	//	char nyxname[256], filename[256];
	//	
	//	for(unsigned ilevel=levelmin_; ilevel<=levelmax_; ++ilevel )
	//	{
	}
	
public:
	
	nyx_output_plugin( config_file& cf )
	: output_plugin( cf )
	{ 
		int argc=1;
		char **argv;
		BoxLib::Initialize(argc,argv);

		bool bhave_hydro = cf_.getValue<bool>("setup","baryons");

		if (bhave_hydro)
			n_data_items = 10;
		else
			n_data_items = 6;

		field_name.resize(n_data_items);
		if (bhave_hydro)
		{
			field_name[0] = "baryon_density";
			field_name[1] = "baryon_vel_x";
			field_name[2] = "baryon_vel_y";
			field_name[3] = "baryon_vel_z";
			field_name[4] = "dm_pos_x";
			field_name[5] = "dm_pos_y";
			field_name[6] = "dm_pos_z";
			field_name[7] = "dm_vel_x";
			field_name[8] = "dm_vel_y";
			field_name[9] = "dm_vel_z";
			the_sim_header.particle_idx = 4;
		}
		else
		{
			field_name[0] = "dm_pos_x";
			field_name[1] = "dm_pos_y";
			field_name[2] = "dm_pos_z";
			field_name[3] = "dm_vel_x";
			field_name[4] = "dm_vel_y";
			field_name[5] = "dm_vel_z";
			the_sim_header.particle_idx = 0;
		}

		f_lev = levelmax_-levelmin_;
		std::cout << f_lev+1 << " level" << std::endl;

		mfs.resize(f_lev+1);

		Array<int> pmap(2);
		pmap[0]=0;
		pmap[1]=0;

		gridp = 1<<levelmin_;

		double off[] = {0, 0, 0};

		//at first we do this only for the topgrid...
		for(int lev = 0; lev <= f_lev; lev++)
		{
			BoxArray   domainBoxArray(1);
			
			int mlev = lev+levelmin_;
			int fac  = (1<<lev);

			off[0] += fac*offx_[lev];
			off[1] += fac*offy_[lev];
			off[2] += fac*offz_[lev];

			IntVect    pdLo(off[0], 
					off[1], 
					off[2]);
			IntVect    pdHi(off[0]+sizex_[lev]-1, 
					off[1]+sizey_[lev]-1, 
					off[2]+sizez_[lev]-1);
				
//			pdLo *= (1<<lev);
//			pdHi *= (1<<lev);

			std::cout << pdLo << std::endl;
			std::cout << pdHi << std::endl;

			// Start with a probDomain
//			IntVect    pdLo(0,0,0);
//			IntVect    pdHi(gridp-1,gridp-1,gridp-1);
			Box        probDomain(pdLo,pdHi);
			
			// We just have one box since we don't use mpi.
		   	domainBoxArray.set(0, probDomain);
			DistributionMapping domainDistMap(pmap);
	
			boxarrays.push_back(domainBoxArray);
	
			int ngrow(0);
			MultiFab *mf = new MultiFab;
			mfs.set(lev,mf);
			mfs[lev].define(domainBoxArray, n_data_items, ngrow, domainDistMap, Fab_allocate);
		}

//		if( mkdir( fname_.c_str(), 0777 ) )
//		{
//			perror( fname_.c_str() );
//			throw std::runtime_error("Error in nyx_output_plugin!");
//		}
		
		bool align_top			= cf.getValueSafe<bool>( "setup", "align_top", true );
		
		if( !align_top )
		{
            LOGERR("nyx output plug-in requires that \'align_top=true\'!");
            throw std::runtime_error("nyx output plug-in requires that \'align_top=true\'!");
		}
        
		the_sim_header.dimensions.push_back( 1<<levelmin_ );
		the_sim_header.dimensions.push_back( 1<<levelmin_ );
		the_sim_header.dimensions.push_back( 1<<levelmin_ );
		
		the_sim_header.offset.push_back( 0 );
		the_sim_header.offset.push_back( 0 );
		the_sim_header.offset.push_back( 0 );
		
		the_sim_header.a_start		= 1.0/(1.0+cf.getValue<double>("setup","zstart"));
		the_sim_header.dx			= cf.getValue<double>("setup","boxlength")/the_sim_header.dimensions[0]/(cf.getValue<double>("cosmology","H0")*0.01); // not sure?!?
		the_sim_header.boxlength=cf.getValue<double>("setup","boxlength");
		the_sim_header.h0			= cf.getValue<double>("cosmology","H0")*0.01;
		
		if( bhave_hydro )
			the_sim_header.omega_b		= cf.getValue<double>("cosmology","Omega_b");
		else
			the_sim_header.omega_b		= 0.0;
		
		the_sim_header.omega_m		= cf.getValue<double>("cosmology","Omega_m");
		the_sim_header.omega_v		= cf.getValue<double>("cosmology","Omega_L");
		the_sim_header.vfact		= cf.getValue<double>("cosmology","vfact")*the_sim_header.h0;   //.. need to multiply by h, nyx wants this factor for non h-1 units
		
		std::cout << "creating output object"  << std::endl;
	}
	
	~nyx_output_plugin()
	{
		std::string FullPath = fname_;
		if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
			BoxLib::CreateDirectoryFailed(FullPath);
		if (!FullPath.empty() && FullPath[FullPath.size()-1] != '/')
			FullPath += '/';
		FullPath += "Header";
		std::ofstream Header(FullPath.c_str());

		for(int lev=0; lev <= f_lev; lev++)
		{
			writeLevelPlotFile (	fname_,
						Header,
						VisMF::OneFilePerCPU,
						lev);
			//FIXME I would prefer VisMF::NFiles
		}
		Header.close();

		writeGridsFile(fname_);
		std::cout << "destroying output object"  << std::endl;
	}
	
	void write_dm_mass( const grid_hierarchy& gh )
	{	/* do nothing, not needed */	}
	
	
	void write_dm_density( const grid_hierarchy& gh )
	{	/* write the parameter file data */	
		// It's very useful to write a parameter file, but WHY here?	
        
	}
	
	
	void write_dm_velocity( int coord, const grid_hierarchy& gh )
	{
		char nyxname[256];
		sprintf( nyxname, "ParticleVelocities_%c", (char)('x'+coord) );
		
		double vunit = 1.0/(1.225e2*sqrt(the_sim_header.omega_m/the_sim_header.a_start));
		
		dump_grid_data(the_sim_header.particle_idx+3+coord, nyxname, gh);
	}
	
	
	void write_dm_position( int coord, const grid_hierarchy& gh )
	{
		char nyxname[256];
		sprintf( nyxname, "ParticleDisplacements_%c", (char)('x'+coord) );
        
		//dump_grid_data( nyxname, gh );
		dump_grid_data(the_sim_header.particle_idx+coord, nyxname, gh);
	}
	
	void write_dm_potential( const grid_hierarchy& gh )
	{ }
	
	void write_gas_potential( const grid_hierarchy& gh )
	{ }
	
	
	void write_gas_velocity( int coord, const grid_hierarchy& gh )
	{
		double vunit = 1.0/(1.225e2*sqrt(the_sim_header.omega_m/the_sim_header.a_start));
		
		char nyxname[256];
		sprintf( nyxname, "GridVelocities_%c", (char)('x'+coord) );
		dump_grid_data(coord+1, nyxname, gh);
	}
	
	
	void write_gas_position( int coord, const grid_hierarchy& gh )
	{	
		/* do nothing, not needed */	
	}
	
	
	void write_gas_density( const grid_hierarchy& gh )
	{
		char nyxname[256];
		sprintf( nyxname, "density" );
		//FIXME factor and add have to be adjusted to the
		//corresponding nyx units...
		dump_grid_data(0, nyxname, gh);
	}
	
	
	void finalize( void )
	{
		//
		//before finalizing we write out an inputs and a probin file for Nyx.
		//
		std::ofstream inputs("inputs");
		std::ofstream probin("probin");
		
		//at first the fortran stuff...
		probin << "&fortin" << std::endl;	
		probin << "  comoving_OmM = " << the_sim_header.omega_m << "d0" << std::endl;
		probin << "  comoving_OmB = " << the_sim_header.omega_b << "d0" << std::endl;
		probin << "  comoving_OmL = " << the_sim_header.omega_v << "d0" << std::endl;
		probin << "  comoving_h   = " << the_sim_header.h0      << "d0" << std::endl;
		probin << "/" << std::endl;
		probin << std::endl;

		//afterwards the cpp stuff...(for which we will need a template, which is read in by the code...)
		inputs << "nyx.final_a = 1.0 " << std::endl;
		inputs << "max_step = 100000 " << std::endl;
		inputs << "nyx.small_dens = 1e-4" << std::endl;
		inputs << "nyx.small_temp = 10" << std::endl;
		inputs << "nyx.cfl            = 0.9     # cfl number for hyperbolic system" << std::endl;
		inputs << "nyx.init_shrink    = 1.0     # scale back initial timestep" << std::endl;
		inputs << "nyx.change_max     = 1.05    # scale back initial timestep" << std::endl;
		inputs << "nyx.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt" << std::endl;
		inputs << "nyx.sum_interval   = 1      # timesteps between computing mass" << std::endl;
		inputs << "nyx.v              = 1       # verbosity in Castro.cpp" << std::endl;
		inputs << "gravity.v             = 1       # verbosity in Gravity.cpp" << std::endl;
		inputs << "amr.v                 = 1       # verbosity in Amr.cpp" << std::endl;
		inputs << "mg.v                  = 0       # verbosity in Amr.cpp" << std::endl;
		inputs << "particles.v           = 1       # verbosity in Particle class" << std::endl;
		inputs << "amr.ref_ratio       = 2 2 2 2 2 2 2 2 " << std::endl;
		inputs << "amr.regrid_int      = 2 2 2 2 2 2 2 2 " << std::endl;
		inputs << "amr.initial_grid_file = init/grids_file" << std::endl;
		inputs << "amr.useFixedCoarseGrids = 1" << std::endl;
		inputs << "amr.blocking_factor =8" << std::endl;
		inputs << "amr.check_file      = chk " << std::endl;
		inputs << "amr.check_int       = 10 " << std::endl;
		inputs << "amr.plot_file       = plt " << std::endl;
		inputs << "amr.plot_int        = 10 " << std::endl;
		inputs << "amr.derive_plot_vars = particle_count particle_mass_density pressure" << std::endl;
		inputs << "amr.plot_vars = ALL" << std::endl;
		inputs << "nyx.add_ext_src = 0" << std::endl;
		inputs << "gravity.gravity_type = PoissonGrav    " << std::endl;
		inputs << "gravity.no_sync      = 1              " << std::endl;
		inputs << "gravity.no_composite = 1              " << std::endl;
		inputs << "mg.bottom_solver = 1                  " << std::endl;
		inputs << "geometry.is_periodic =  1     1     1 " << std::endl;
		inputs << "geometry.coord_sys   =  0             " << std::endl;
		inputs << "amr.max_grid_size    = 32             " << std::endl;
		inputs << "nyx.lo_bc       =  0   0   0          " << std::endl;
		inputs << "nyx.hi_bc       =  0   0   0          " << std::endl;
		inputs << "nyx.do_grav  = 1                      " << std::endl;
		inputs << "nyx.do_dm_particles = 1               " << std::endl;
		inputs << "nyx.particle_init_type = Cosmological " << std::endl;
                inputs << "nyx.print_fortran_warnings = 0" << std::endl;
		inputs << "cosmo.initDirName  = init             " << std::endl;
		inputs << "nyx.particle_move_type = Gravitational" << std::endl;
		inputs << "amr.probin_file = probin              " << std::endl;
		inputs << "cosmo.ic-source = MUSIC               " << std::endl;
		

		inputs << "nyx.do_hydro = "<< (the_sim_header.omega_b>0?1:0) << std::endl;
		inputs << "amr.max_level       = " << levelmax_-levelmin_ << std::endl;
		inputs << "nyx.initial_z = " << 1/the_sim_header.a_start-1 << std::endl;
		inputs << "amr.n_cell           = " << sizex_[0] << " " << sizey_[0] << " " << sizez_[0] << std::endl;
		inputs << "nyx.n_particles      = " << sizex_[0] << " " << sizey_[0] << " " << sizez_[0] << std::endl;
		inputs << "geometry.prob_lo     = 0 0 0" << std::endl;  

		//double dx = the_sim_header.dx/the_sim_header.h0;
		double bl = the_sim_header.boxlength/the_sim_header.h0;
		inputs << "geometry.prob_hi     = " << bl << " " << bl << " " << bl << std::endl;



		probin.close();
		inputs.close();
		std::cout << "finalizing..." << std::endl;
		
	}
		
	



	void writeLevelPlotFile (const	std::string&	dir,
					std::ostream&	os,
		                        VisMF::How	how,
					int		level)
	{
		int i, n;
		
		const Real cur_time = 0.0;

		std::cout << "in writeLevelPlotFile" << std::endl;
		double h0 = cf_.getValue<double>("cosmology", "H0")*0.01;

//		for (MFIter mfi(mf); mfi.isValid(); ++mfi)
//		{ 
//			std::cout << "bla" << std::endl;
//			std::cout << mf[mfi] << std::endl;
//		}
	
		if (level == 0)
		{
			//
			// The first thing we write out is the plotfile type.
			//
			os << "MUSIC_for_Nyx_v0.1" << '\n';
			
			os << n_data_items << '\n';
			
			for (i = 0; i < n_data_items; i++)
				os << field_name[i] << '\n';
			
			os << 3 << '\n';
			os << 0 << '\n';
		
			os << f_lev << '\n';
		
			for (i = 0; i < BL_SPACEDIM; i++)
				os << 0 << ' '; //ProbLo
			os << '\n';
			double boxlength = cf_.getValue<double>("setup","boxlength");
			for (i = 0; i < BL_SPACEDIM; i++)
				os << boxlength/h0 << ' '; //ProbHi
			os << '\n';
		
			for (i = 0; i < f_lev; i++)
				os << 2 << ' '; //refinement factor
			os << '\n';
		
			IntVect    pdLo(0,0,0);
			IntVect    pdHi(gridp-1,gridp-1,gridp-1);
//			Box        probDomain(pdLo,pdHi);
			for (i = 0; i <= f_lev; i++) //Geom(i).Domain()
			{
//				IntVect    pdLo(offx_[i], offy_[i], offz_[i]);
//				IntVect    pdHi(offx_[i]+sizex_[i], offy_[i]+sizey_[i], offz_[i]+sizez_[i]);
				Box        probDomain(pdLo,pdHi);
				os << probDomain << ' ';
				pdHi *= 2;
				pdHi += 1;
			}
			os << '\n';
		
			for (i = 0; i <= f_lev; i++) //level steps
				os << 0 << ' ';
			os << '\n';
		
			double dx = cf_.getValue<double>("setup","boxlength")/gridp/h0;
			for (i = 0; i <= f_lev; i++)
			{
				for (int k = 0; k < BL_SPACEDIM; k++)
					os << dx << ' ';
				os << '\n';
				dx = dx/2.;
			}
			os << 0 << '\n';
			os << "0\n"; // Write bndry data.
		}
	
		//
		// Build the directory to hold the MultiFab at this level.
		// The name is relative to the directory containing the Header file.
		//
		static const std::string BaseName = "/Cell";
		
		std::string Level = BoxLib::Concatenate("Level_", level, 1);
		//
		// Now for the full pathname of that directory.
		//
		std::string FullPath = dir;
		if (!FullPath.empty() && FullPath[FullPath.size()-1] != '/')
			FullPath += '/';
		FullPath += Level;
		//
		// Only the I/O processor makes the directory if it doesn't already exist.
		//
		if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
			BoxLib::CreateDirectoryFailed(FullPath);
		
		os << level << ' ' << boxarrays[level].size() << ' ' << 0 << '\n';
		os << 0 << '\n';
	
		double cellsize[3];
		double dx = cf_.getValue<double>("setup","boxlength")/gridp/h0;
		for (n = 0; n < BL_SPACEDIM; n++)
		{
			cellsize[n] = dx;
		}
		for (i = 0; i < level; i++)
		{
			for (n = 0; n < BL_SPACEDIM; n++)
			{
				cellsize[n] /= 2.;
			}
		}
		std::cout << cellsize[0] << std::endl;
		for (i = 0; i < boxarrays[level].size(); ++i)
		{
			double problo[] = {0,0,0};
			RealBox gridloc = RealBox(boxarrays[level][i], cellsize, problo);
			for (n = 0; n < BL_SPACEDIM; n++)
				os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
		}
		//
		// The full relative pathname of the MultiFabs at this level.
		// The name is relative to the Header file containing this name.
		// It's the name that gets written into the Header.
		//
		std::string PathNameInHeader = Level;
		PathNameInHeader += BaseName;
		os << PathNameInHeader << '\n';
	
		//
		// Use the Full pathname when naming the MultiFab.
		//
		std::string TheFullPath = FullPath;
		TheFullPath += BaseName;
		VisMF::Write(mfs[level],TheFullPath,how,true);
	}

	void writeGridsFile (const	std::string&	dir)
	{
		int i, n;

		std::cout << "in writeGridsFile" << std::endl;
		
		std::string myFname = dir;
		if (!myFname.empty() && myFname[myFname.size()-1] != '/')
			myFname += '/';
		myFname += "grids_file";

		std::ofstream os(myFname.c_str());
		
		os << f_lev << '\n';
		
		for (int lev = 1; lev <= f_lev; lev++)
		{
			os << boxarrays[lev].size() << '\n';
			boxarrays[lev].coarsen(2);
			for (i=0; i < boxarrays[lev].size(); i++)
				os << boxarrays[lev][i] << "\n";
		}
		os.close();
	}
//	void get_grids (const grid_hierarchy &gh)
//	{
//		for(unsigned ilevel=levelmin_; ilevel<=levelmax_; ++ilevel )
//		{
//			grid_on_one_level gool;
//			
//			int xs = gh.get_grid(ilevel)->size(0);
//			int ys = gh.get_grid(ilevel)->size(1);
//			int zs = gh.get_grid(ilevel)->size(2);
//
//			int xo = gh.get_grid(ilevel)->offset(0);
//			int yo = gh.get_grid(ilevel)->offset(1);
//			int zo = gh.get_grid(ilevel)->offset(2);
//
//			IntVect gridLo(xo,yo,zo);
//			gool.lo = gridLo;
//
//			IntVect gridHi(xo+xs,yo+zs,zo+zs);
//			gool.hi = gridHi;
//
//			grids.push_back(gool);
//		}
//		
//	}
};

namespace{
	output_plugin_creator_concrete<nyx_output_plugin> creator("nyx");
}

#endif //HAVE_BOXLIB

