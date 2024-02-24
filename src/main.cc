// This file is part of MUSIC
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2010-2024 by Oliver Hahn
// 
// monofonIC is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// monofonIC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <math.h>

#include <thread>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>

#if defined(CMAKE_BUILD)
extern "C"
{
	extern const char *GIT_TAG;
	extern const char *GIT_REV;
	extern const char *GIT_BRANCH;
}
#endif

#include <exception>
#include <cfenv>

#include <general.hh>
#include <defaults.hh>
#include <output.hh>

#include <config_file.hh>

#include <poisson.hh>
#include <mg_solver.hh>
#include <fd_schemes.hh>
#include <random.hh>
#include <densities.hh>

#include <convolution_kernel.hh>
#include <perturbation_theory.hh>
#include <cosmology_parameters.hh>
#include <cosmology_calculator.hh>
#include <transfer_function.hh>

#define THE_CODE_NAME "music!"
#define THE_CODE_VERSION "2.0a"

// initialise with "default" values
namespace CONFIG{
// int  MPI_thread_support = -1;
int  MPI_task_rank = 0;
int  MPI_task_size = 1;
bool MPI_ok = false;
// bool MPI_threads_ok = false;
bool FFTW_threads_ok = false;
int  num_threads = 1;
}


namespace music
{

	struct framework
	{
		// transfer_function *the_transfer_function;
		// poisson_solver *the_poisson_solver;
		config_file *the_config_file;
		refinement_hierarchy *the_refinement_hierarchy;
	};

}

//... declare static class members here
// transfer_function *TransferFunction_real::ptf_ = NULL;
transfer_function *TransferFunction_k::ptf_ = NULL;
tf_type TransferFunction_k::type_;
// tf_type TransferFunction_real::type_;
// real_t TransferFunction_real::nspec_ = -1.0;
real_t TransferFunction_k::nspec_ = -1.0;

std::unique_ptr<cosmology::calculator>  the_cosmo_calc;

//... prototypes for routines used in main driver routine
void splash(void);
void modify_grid_for_TF(const refinement_hierarchy &rh_full, refinement_hierarchy &rh_TF, config_file &cf);
void print_hierarchy_stats(config_file &cf, const refinement_hierarchy &rh);
void store_grid_structure(config_file &cf, const refinement_hierarchy &rh);
double compute_finest_mean(grid_hierarchy &u);
double compute_finest_sigma(grid_hierarchy &u);

void splash(void)
{

	music::ilog << std::endl
			<< "           __    __     __  __     ______     __     ______      " << std::endl
			<< "          /\\ \"-./  \\   /\\ \\/\\ \\   /\\  ___\\   /\\ \\   /\\  ___\\  " << std::endl
			<< "          \\ \\ \\-./\\ \\  \\ \\ \\_\\ \\  \\ \\___  \\  \\ \\ \\  \\ \\ \\____ " << std::endl
			<< "           \\ \\_\\ \\ \\_\\  \\ \\_____\\  \\/\\_____\\  \\ \\_\\  \\ \\_____\\ " << std::endl
			<< "            \\/_/  \\/_/   \\/_____/   \\/_____/   \\/_/   \\/_____/ " << std::endl << std::endl
			<< "                                   this is " << THE_CODE_NAME << " version " << THE_CODE_VERSION << std::endl << std::endl;

// git and versioning info:
    music::ilog << "Version: git rev.: " << GIT_REV << ", tag: " << GIT_TAG << ", branch: " << GIT_BRANCH << std::endl;
    
    // Compilation CMake configuration, time etc info:
    music::ilog << "This " << CMAKE_BUILDTYPE_STR << " build was compiled at " << __TIME__ << " on " <<  __DATE__ << std::endl;

#ifdef __GNUC__
    music::ilog << "Compiled with GNU C++ version " << __VERSION__ <<std::endl;
#else
    music::ilog << "Compiled with " << __VERSION__ << std::endl;
#endif
}

void modify_grid_for_TF(const refinement_hierarchy &rh_full, refinement_hierarchy &rh_TF, config_file &cf)
{
	unsigned lbase, lbaseTF, lmax, overlap;

	lbase = cf.get_value<unsigned>("setup", "levelmin");
	lmax = cf.get_value<unsigned>("setup", "levelmax");
	lbaseTF = cf.get_value_safe<unsigned>("setup", "levelmin_TF", lbase);
	overlap = cf.get_value_safe<unsigned>("setup", "overlap", 4);
	rh_TF = rh_full;

	unsigned pad = overlap;

	for (unsigned i = lbase + 1; i <= lmax; ++i)
	{
		int x0[3], lx[3], lxmax = 0;

		for (int j = 0; j < 3; ++j)
		{
			lx[j] = rh_TF.size(i, j) + 2 * pad;
			x0[j] = rh_TF.offset_abs(i, j) - pad;

			if (lx[j] > lxmax)
				lxmax = lx[j];
		}

		//... make sure that grids are divisible by 4 for convolution.
		lxmax += lxmax % 4;

		for (int j = 0; j < 3; ++j)
		{
			double dl = 0.5 * ((double)(lxmax - lx[j]));
			int add_left = (int)ceil(dl);

			lx[j] = lxmax;
			x0[j] -= add_left;
			x0[j] += x0[j] % 2;
		}

		rh_TF.adjust_level(i, lx[0], lx[1], lx[2], x0[0], x0[1], x0[2]);
	}

	if (lbaseTF > lbase)
	{
		music::ilog << "- Will use levelmin = " << lbaseTF << " to compute density field...\n";

		for (unsigned i = lbase; i <= lbaseTF; ++i)
		{
			unsigned nfull = (unsigned)pow(2, i);
			rh_TF.adjust_level(i, nfull, nfull, nfull, 0, 0, 0);
		}
	}
}

void print_hierarchy_stats(config_file &cf, const refinement_hierarchy &rh)
{
	double omegam = cf.get_value<double>("cosmology", "Omega_m");
	double omegab = cf.get_value<double>("cosmology", "Omega_b");
	bool bbaryons = cf.get_value<bool>("setup", "baryons");
	double boxlength = cf.get_value<double>("setup", "boxlength");

	unsigned levelmin = rh.levelmin();
	double dx = boxlength / (double)(1 << levelmin), dx3 = dx * dx * dx;
	double rhom = 2.77519737e11; // h-1 M_o / (h-1 Mpc)**3
	double cmass, bmass(0.0), mtotgrid;
	if (bbaryons)
	{
		cmass = (omegam - omegab) * rhom * dx3;
		bmass = omegab * rhom * dx3;
	}
	else
		cmass = omegam * rhom * dx3;

	music::ilog << "-------------------------------------------------------------------------------" << std::endl;

	if (rh.get_shift(0) != 0 || rh.get_shift(1) != 0 || rh.get_shift(2) != 0)
		music::ilog << "- Domain will be shifted by (" << rh.get_shift(0) << ", " << rh.get_shift(1) << ", " << rh.get_shift(2) << ")\n"
							<< std::endl;

	music::ilog << "- Grid structure:\n";

	for (unsigned ilevel = rh.levelmin(); ilevel <= rh.levelmax(); ++ilevel)
	{
		double rfac = 1.0 / (1 << (ilevel - rh.levelmin())), rfac3 = rfac * rfac * rfac;

		mtotgrid = omegam * rhom * dx3 * rfac3 * rh.size(ilevel, 0) * rh.size(ilevel, 1) * rh.size(ilevel, 2);
		music::ilog
				<< "     Level " << std::setw(3) << ilevel << " :   offset = (" << std::setw(5) << rh.offset(ilevel, 0) << ", " << std::setw(5) << rh.offset(ilevel, 1) << ", " << std::setw(5) << rh.offset(ilevel, 2) << ")\n"
				<< "                     size = (" << std::setw(5) << rh.size(ilevel, 0) << ", " << std::setw(5) << rh.size(ilevel, 1) << ", " << std::setw(5) << rh.size(ilevel, 2) << ")\n";

		if (ilevel == rh.levelmax())
		{
			music::ilog << "-------------------------------------------------------------------------------" << std::endl;
			music::ilog << "- Finest level :\n";

			if (dx * rfac > 0.1)
				music::ilog << "                   extent =  " << dx * rfac * rh.size(ilevel, 0) << " x " << dx * rfac * rh.size(ilevel, 1) << " x " << dx * rfac * rh.size(ilevel, 2) << " h-3 Mpc**3\n";
			else if (dx * rfac > 1e-4)
				music::ilog << "                   extent =  " << dx * rfac * 1000.0 * rh.size(ilevel, 0) << " x " << dx * rfac * 1000.0 * rh.size(ilevel, 1) << " x " << dx * rfac * 1000.0 * rh.size(ilevel, 2) << " h-3 kpc**3\n";
			else
				music::ilog << "                   extent =  " << dx * rfac * 1.e6 * rh.size(ilevel, 0) << " x " << dx * rfac * 1.e6 * rh.size(ilevel, 1) << " x " << dx * rfac * 1.e6 * rh.size(ilevel, 2) << " h-3 pc**3\n";

			music::ilog << "                 mtotgrid =  " << mtotgrid << " h-1 M_o\n";
			music::ilog << "            particle mass =  " << cmass * rfac3 << " h-1 M_o\n";
			if (bbaryons)
				music::ilog << "         baryon mass/cell =  " << bmass * rfac3 << " h-1 M_o\n";
			if (dx * rfac > 0.1)
				music::ilog << "                       dx =  " << dx * rfac << " h-1 Mpc\n";
			else if (dx * rfac > 1e-4)
				music::ilog << "                       dx =  " << dx * rfac * 1000.0 << " h-1 kpc\n";
			else
				music::ilog << "                       dx =  " << dx * rfac * 1.e6 << " h-1 pc\n";
		}
	}
	music::ilog << "-------------------------------------------------------------------------------" << std::endl;
}

void store_grid_structure(config_file &cf, const refinement_hierarchy &rh)
{
	char str1[128], str2[128];
	for (unsigned i = rh.levelmin(); i <= rh.levelmax(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			snprintf(str1, 128, "offset(%d,%d)", i, j);
			snprintf(str2, 128, "%ld", rh.offset(i, j));
			cf.insert_value("setup", str1, str2);

			snprintf(str1, 128, "size(%d,%d)", i, j);
			snprintf(str2, 128, "%ld", rh.size(i, j));
			cf.insert_value("setup", str1, str2);
		}
	}
}

double compute_finest_sigma(grid_hierarchy &u)
{
	double sum = 0.0, sum2 = 0.0;
	#pragma omp parallel for reduction(+:sum,sum2)
	for (int ix = 0; ix < (int)(*u.get_grid(u.levelmax())).size(0); ++ix)
		for (int iy = 0; iy < (int)(*u.get_grid(u.levelmax())).size(1); ++iy)
			for (int iz = 0; iz < (int)(*u.get_grid(u.levelmax())).size(2); ++iz)
			{
				sum += (*u.get_grid(u.levelmax()))(ix, iy, iz);
				sum2 += (*u.get_grid(u.levelmax()))(ix, iy, iz) * (*u.get_grid(u.levelmax()))(ix, iy, iz);
			}

	size_t N = (size_t)(*u.get_grid(u.levelmax())).size(0) * (size_t)(*u.get_grid(u.levelmax())).size(1) * (size_t)(*u.get_grid(u.levelmax())).size(2);
	sum /= N;
	sum2 /= N;

	return sqrt(sum2 - sum * sum);
}

double compute_finest_absmax(grid_hierarchy &u)
{
	double valmax = 0.0;
	#pragma omp parallel for reduction(max:valmax)
	for (int ix = 0; ix < (int)(*u.get_grid(u.levelmax())).size(0); ++ix)
		for (int iy = 0; iy < (int)(*u.get_grid(u.levelmax())).size(1); ++iy)
			for (int iz = 0; iz < (int)(*u.get_grid(u.levelmax())).size(2); ++iz)
			{
				if (std::fabs((*u.get_grid(u.levelmax()))(ix, iy, iz)) > valmax)
					valmax = std::fabs((*u.get_grid(u.levelmax()))(ix, iy, iz));
			}

	return valmax;
}

double compute_finest_mean(grid_hierarchy &u)
{
	double sum = 0.0;
	size_t count = 0;
	#pragma omp parallel for reduction(+:sum,count)
	for (int ix = 0; ix < (int)(*u.get_grid(u.levelmax())).size(0); ++ix)
		for (int iy = 0; iy < (int)(*u.get_grid(u.levelmax())).size(1); ++iy)
			for (int iz = 0; iz < (int)(*u.get_grid(u.levelmax())).size(2); ++iz)
			{
				sum += (*u.get_grid(u.levelmax()))(ix, iy, iz);
				count++;
			}

	return sum/count;
}

void add_constant_value( grid_hierarchy &u, const double val )
{
	for( unsigned ilvl = u.levelmin(); ilvl<=u.levelmax(); ++ilvl )
	{
		#pragma omp parallel for
		for (int ix = 0; ix < (int)(*u.get_grid(ilvl)).size(0); ++ix)
			for (int iy = 0; iy < (int)(*u.get_grid(ilvl)).size(1); ++iy)
				for (int iz = 0; iz < (int)(*u.get_grid(ilvl)).size(2); ++iz)
				{
					(*u.get_grid(ilvl))(ix, iy, iz) += val;
				}
	}
}

#include <system_stat.hh>
void output_system_info()
{
	std::feclearexcept(FE_ALL_EXCEPT);

	//------------------------------------------------------------------------------
	// Write code configuration to screen
	//------------------------------------------------------------------------------
	// hardware related infos
	music::ilog << std::setw(32) << std::left << "CPU vendor string" << " : " << SystemStat::Cpu().get_CPUstring() << std::endl;
	
	// multi-threading related infos
	music::ilog << std::setw(32) << std::left << "Available HW threads / task" << " : " << std::thread::hardware_concurrency() << " (" << CONFIG::num_threads << " used)" << std::endl;

	// memory related infos
	SystemStat::Memory mem;

	unsigned availpmem = mem.get_AvailMem()/1024/1024;
	unsigned usedpmem = mem.get_UsedMem()/1024/1024;
	unsigned maxpmem = availpmem, minpmem = availpmem;
	unsigned maxupmem = usedpmem, minupmem = usedpmem;
	
	music::ilog << std::setw(32) << std::left << "Total system memory (phys)" << " : " << mem.get_TotalMem()/1024/1024 << " Mb" << std::endl;
	music::ilog << std::setw(32) << std::left << "Used system memory (phys)" << " : " << "Max: " << maxupmem << " Mb, Min: " << minupmem << " Mb" << std::endl;
	music::ilog << std::setw(32) << std::left << "Available system memory (phys)" << " : " <<  "Max: " << maxpmem << " Mb, Min: " << minpmem << " Mb" << std::endl;
			
	// Kernel related infos
	SystemStat::Kernel kern;
	auto kinfo = kern.get_kernel_info();
	music::ilog << std::setw(32) << std::left << "OS/Kernel version" << " : " << kinfo.kernel << " version " << kinfo.major << "." << kinfo.minor << " build " << kinfo.build_number << std::endl;

	// FFTW related infos
	music::ilog << std::setw(32) << std::left << "FFTW version" << " : " << FFTW_API(version) << std::endl;
	music::ilog << std::setw(32) << std::left << "FFTW supports multi-threading" << " : " << (CONFIG::FFTW_threads_ok? "yes" : "no") << std::endl;
	music::ilog << std::setw(32) << std::left << "FFTW mode" << " : ";
#if defined(FFTW_MODE_PATIENT)
	music::ilog << "FFTW_PATIENT" << std::endl;
#elif defined(FFTW_MODE_MEASURE)
    music::ilog << "FFTW_MEASURE" << std::endl;
#else
	music::ilog << "FFTW_ESTIMATE" << std::endl;
#endif
}

/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/

region_generator_plugin *the_region_generator;

int main(int argc, const char *argv[])
{
	const unsigned nbnd = 4;

	unsigned lbase, lmax, lbaseTF;

#if defined(NDEBUG)
	music::logger::set_level(music::log_level::info);
#else
	music::logger::set_level(music::log_level::debug);
#endif

	//------------------------------------------------------------------------------
	//... parse command line options
	//------------------------------------------------------------------------------

	
	if (argc != 2)
	{
		splash();
		std::cout << " This version is compiled with the following plug-ins:\n";

		cosmology::print_ParameterSets();
		print_region_generator_plugins();
		print_transfer_function_plugins();
		print_RNG_plugins();
		print_output_plugins();

		std::cerr << "\n In order to run, you need to specify a parameter file!\n\n";
		exit(0);
	}

	//------------------------------------------------------------------------------
	//... open log file
	//------------------------------------------------------------------------------

	char logfname[128];
	snprintf(logfname, 128, "%s_log.txt", argv[1]);
	music::logger::set_output(logfname);
	time_t ltime = time(NULL);

	splash();
	music::ilog.Print("Opening log file \'%s\'.", logfname);
	music::ulog.Print("Running %s, version %s", THE_CODE_NAME, THE_CODE_VERSION);
	music::ulog.Print("Log is for run started %s", asctime(localtime(&ltime)));

	//------------------------------------------------------------------------------
	//... read and interpret config file
	//------------------------------------------------------------------------------
	config_file cf(argv[1]);
	std::string tfname, randfname, temp;
	bool force_shift(false);

	//------------------------------------------------------------------------------
	//... init multi-threading
	//------------------------------------------------------------------------------
	CONFIG::FFTW_threads_ok = FFTW_API(init_threads)();
	CONFIG::num_threads = cf.get_value_safe<unsigned>("execution", "NumThreads",std::thread::hardware_concurrency());

	music::ilog << "-------------------------------------------------------------------------------" << std::endl;
	output_system_info();
	music::ilog << "-------------------------------------------------------------------------------" << std::endl;
  
	//------------------------------------------------------------------------------
	//... initialize some parameters about grid set-up
	//------------------------------------------------------------------------------

	lbase = cf.get_value<unsigned>("setup", "levelmin");
	lmax = cf.get_value<unsigned>("setup", "levelmax");
	lbaseTF = cf.get_value_safe<unsigned>("setup", "levelmin_TF", lbase);

	if (lbase == lmax && !force_shift)
		cf.insert_value("setup", "no_shift", "yes");

	if (lbaseTF < lbase)
	{
		music::ilog << " - WARNING: levelminTF < levelmin. This is not good!" << std::endl
							  << "            I will set levelminTF = levelmin." << std::endl;

		music::ulog.Print("levelminTF < levelmin. set levelminTF = levelmin.");

		lbaseTF = lbase;
		cf.insert_value("setup", "levelmin_TF", cf.get_value<std::string>("setup", "levelmin"));
	}


	//------------------------------------------------------------------------------
	//... initialize cosmology
	//------------------------------------------------------------------------------
	bool
			do_baryons = cf.get_value<bool>("setup", "baryons"),
			do_2LPT = cf.get_value_safe<bool>("setup", "use_2LPT", false),
			do_LLA = cf.get_value_safe<bool>("setup", "use_LLA", false),
			do_counter_mode = cf.get_value_safe<bool>("setup", "zero_zoom_velocity", false);

	the_cosmo_calc              = std::make_unique<cosmology::calculator>(cf);

	bool tf_has_velocities = the_cosmo_calc.get()->transfer_function_.get()->tf_has_velocities();
	//--------------------------------------------------------------------------------------------------------
	//! starting redshift
	const real_t zstart = cf.get_value<double>("setup", "zstart");
	const real_t astart = 1.0/(1.0+zstart);
	
	music::ilog << "- starting at a=" << 1.0/(1.0+zstart) << std::endl;

	double cosmo_dplus = the_cosmo_calc->get_growth_factor(astart) / the_cosmo_calc->get_growth_factor(1.0);
	double cosmo_vfact = the_cosmo_calc->get_vfact(astart);

	// if (!the_cosmo_calc.get()->transfer_function_.get()->tf_has_total0()){
	// 	the_cosmo_calc->cosmo_param_["pnorm"] *= cosmo_dplus * cosmo_dplus;
	// }
	double cosmo_pnorm = the_cosmo_calc->cosmo_param_["pnorm"];
	//... directly use the normalisation via a parameter rather than the calculated one
	cosmo_pnorm = cf.get_value_safe<double>("setup", "force_pnorm", cosmo_pnorm);

	double vfac2lpt = 1.0;

	if (the_cosmo_calc->transfer_function_->tf_velocity_units() && do_baryons)
	{
		vfac2lpt = cosmo_vfact; // if the velocities are in velocity units, we need to divide by vfact for the 2lPT term
		cosmo_vfact = 1.0;
	}

	
	{
		char tmpstr[128];
		snprintf(tmpstr, 128, "%.12g", cosmo_pnorm);
		cf.insert_value("cosmology", "pnorm", tmpstr);
		snprintf(tmpstr, 128, "%.12g", cosmo_dplus);
		cf.insert_value("cosmology", "dplus", tmpstr);
		snprintf(tmpstr, 128, "%.12g", cosmo_vfact);
		cf.insert_value("cosmology", "vfact", tmpstr);
	}

	the_region_generator = select_region_generator_plugin(cf);

	//------------------------------------------------------------------------------
	//... determine run parameters
	//------------------------------------------------------------------------------

	


	//------------------------------------------------------------------------------
	//... start up the random number generator plugin
	//... see if we need to set some grid building constraints
	noise_generator rand( cf );

	//------------------------------------------------------------------------------
	//... determine the refinement hierarchy
	//------------------------------------------------------------------------------

	refinement_hierarchy rh_Poisson(cf);
	store_grid_structure(cf, rh_Poisson);
	// rh_Poisson.output();
	print_hierarchy_stats(cf, rh_Poisson);

	refinement_hierarchy rh_TF(rh_Poisson);
	modify_grid_for_TF(rh_Poisson, rh_TF, cf);
	// rh_TF.output();

	music::ulog.Print("Grid structure for Poisson solver:");
	rh_Poisson.output_log();
	music::ulog.Print("Grid structure for density convolution:");
	rh_TF.output_log();

	//------------------------------------------------------------------------------
	//... initialize the output plug-in
	//------------------------------------------------------------------------------
	std::string outformat, outfname;
	outformat = cf.get_value<std::string>("output", "format");
	outfname = cf.get_value<std::string>("output", "filename");
	output_plugin *the_output_plugin = select_output_plugin(cf);

	//------------------------------------------------------------------------------
	//... initialize the random numbers
	//------------------------------------------------------------------------------
	music::ilog << "===============================================================================" << std::endl;
	music::ilog << "   GENERATING WHITE NOISE\n";
	music::ilog << "-------------------------------------------------------------------------------" << std::endl;
	music::ilog << "Computing white noise..." << std::endl;
	rand.initialize_for_grid_structure( rh_TF );

	//------------------------------------------------------------------------------
	//... initialize the Poisson solver
	//------------------------------------------------------------------------------
	bool bdefd = true; // we set this by default and don't allow it to be changed outside any more
	bool bglass = cf.get_value_safe<bool>("output", "glass", false);
	bool bsph = cf.get_value_safe<bool>("setup", "do_SPH", false) && do_baryons;
	bool bbshift = bsph && !bglass;

	bool kspace = cf.get_value_safe<bool>("poisson", "kspace", false);
	bool kspace2LPT = kspace;

	bool decic_DM = cf.get_value_safe<bool>("output", "glass_cicdeconvolve", false);
	bool decic_baryons = cf.get_value_safe<bool>("output", "glass_cicdeconvolve", false) & bsph;

	std::array<double,3> counter_mode_amp;

	//... if in unigrid mode, use k-space instead of hybrid
	if (bdefd && (lbase == lmax))
	{
		kspace = true;
		bdefd = false;
		kspace2LPT = false;
	}

	std::string poisson_solver_name;
	if (kspace)
		poisson_solver_name = std::string("fft_poisson");
	else
		poisson_solver_name = std::string("mg_poisson");

	unsigned grad_order = cf.get_value_safe<unsigned>("poisson", "grad_order", 4);

	//... switch off if using kspace anyway
	// bdefd &= !kspace;

	poisson_plugin_creator *the_poisson_plugin_creator = get_poisson_plugin_map()[poisson_solver_name];
	poisson_plugin *the_poisson_solver = the_poisson_plugin_creator->create(cf);

	// .. this parameter needs to be read after the random module is initialised as it will be overwritten by it
	const bool use_fourier_coarsening = cf.get_value_safe<bool>("setup", "fourier_splicing", true);

	//---------------------------------------------------------------------------------
	//... THIS IS THE MAIN DRIVER BRANCHING TREE RUNNING THE VARIOUS PARTS OF THE CODE
	//---------------------------------------------------------------------------------
	bool bfatal = false;
	try
	{
		if (!do_2LPT)
		{
			music::ulog.Print("Entering 1LPT branch");

			//------------------------------------------------------------------------------
			//... cdm density and displacements
			//------------------------------------------------------------------------------
			music::ilog << "===============================================================================" << std::endl;
			music::ilog << "   COMPUTING DARK MATTER DISPLACEMENTS\n";
			music::ilog << "-------------------------------------------------------------------------------" << std::endl;
			music::ulog.Print("Computing dark matter displacements...");

			grid_hierarchy f(nbnd); //, u(nbnd);
			tf_type my_tf_type = delta_cdm;
			if (!do_baryons)
				my_tf_type = delta_matter;

			GenerateDensityHierarchy(cf, the_cosmo_calc.get(), my_tf_type, rh_TF, rand, f, false, false);
			coarsen_density(rh_Poisson, f, use_fourier_coarsening);
			f.add_refinement_mask(rh_Poisson.get_coord_shift());

			normalize_density(f);

			music::ulog.Print("Writing CDM data");
			the_output_plugin->write_dm_mass(f);
			the_output_plugin->write_dm_density(f);

			grid_hierarchy u(f);
			u.zero();
			the_poisson_solver->solve(f, u);

			if (!bdefd)
				f.deallocate();

			music::ulog.Print("Writing CDM potential");
			the_output_plugin->write_dm_potential(u);

			//------------------------------------------------------------------------------
			//... DM displacements
			//------------------------------------------------------------------------------
			{
				grid_hierarchy data_forIO(u);
				for (int icoord = 0; icoord < 3; ++icoord)
				{
					if (bdefd)
					{
						data_forIO.zero();
						*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
						poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order,
													 data_forIO.levelmin() == data_forIO.levelmax(), decic_DM);
						*data_forIO.get_grid(data_forIO.levelmax()) /= 1 << f.levelmax();
						the_poisson_solver->gradient_add(icoord, u, data_forIO);
					}
					else
						//... displacement
						the_poisson_solver->gradient(icoord, u, data_forIO);
					double dispmax = compute_finest_absmax(data_forIO);
					music::ilog.Print("\t - max. %c-displacement of HR particles is %f [mean dx]", 'x' + icoord, dispmax * (double)(1ll << data_forIO.levelmax()));
					coarsen_density(rh_Poisson, data_forIO, false);
					
					//... compute counter-mode to minimize advection errors
					counter_mode_amp[icoord] = compute_finest_mean(data_forIO); 
					if( do_counter_mode ) add_constant_value( data_forIO, -counter_mode_amp[icoord] );
					
					music::ulog.Print("Writing CDM displacements");
					the_output_plugin->write_dm_position(icoord, data_forIO);
				}
				if (do_baryons)
					u.deallocate();
				data_forIO.deallocate();
			}

			//------------------------------------------------------------------------------
			//... gas density
			//------------------------------------------------------------------------------
			if (do_baryons)
			{
				music::ilog << "===============================================================================" << std::endl;
				music::ilog << "   COMPUTING BARYON DENSITY\n";
				music::ilog << "-------------------------------------------------------------------------------" << std::endl;
				music::ulog.Print("Computing baryon density...");
				GenerateDensityHierarchy(cf, the_cosmo_calc.get(), delta_baryon, rh_TF, rand, f, false, bbshift);
				coarsen_density(rh_Poisson, f, use_fourier_coarsening);
				f.add_refinement_mask(rh_Poisson.get_coord_shift());
				normalize_density(f);

				if (!do_LLA)
				{
					music::ulog.Print("Writing baryon density");
					the_output_plugin->write_gas_density(f);
				}

				if (bsph)
				{
					u = f;
					u.zero();
					the_poisson_solver->solve(f, u);

					if (!bdefd)
						f.deallocate();

					grid_hierarchy data_forIO(u);
					for (int icoord = 0; icoord < 3; ++icoord)
					{
						if (bdefd)
						{
							data_forIO.zero();
							*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
							poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order,
														 data_forIO.levelmin() == data_forIO.levelmax(), decic_baryons);
							*data_forIO.get_grid(data_forIO.levelmax()) /= 1 << f.levelmax();
							the_poisson_solver->gradient_add(icoord, u, data_forIO);
						}
						else
							//... displacement
							the_poisson_solver->gradient(icoord, u, data_forIO);

						coarsen_density(rh_Poisson, data_forIO, false);
						music::ulog.Print("Writing baryon displacements");
						the_output_plugin->write_gas_position(icoord, data_forIO);
					}
					u.deallocate();
					data_forIO.deallocate();
					if (bdefd)
						f.deallocate();
				}
				else if (do_LLA)
				{
					u = f;
					u.zero();
					the_poisson_solver->solve(f, u);
					compute_LLA_density(u, f, grad_order);
					u.deallocate();
					normalize_density(f);
					music::ulog.Print("Writing baryon density");
					the_output_plugin->write_gas_density(f);
				}

				f.deallocate();
			}

			//------------------------------------------------------------------------------
			//... velocities
			//------------------------------------------------------------------------------
			if ((!tf_has_velocities || !do_baryons) && !bsph)
			{
				music::ilog << "===============================================================================" << std::endl;
				music::ilog << "   COMPUTING VELOCITIES\n";
				music::ilog << "-------------------------------------------------------------------------------" << std::endl;
				music::ulog.Print("Computing velocitites...");

				if (do_baryons || tf_has_velocities)
				{
					music::ulog.Print("Generating velocity perturbations...");
					GenerateDensityHierarchy(cf, the_cosmo_calc.get(), theta_cdm, rh_TF, rand, f, false, false);
					coarsen_density(rh_Poisson, f, use_fourier_coarsening);
					f.add_refinement_mask(rh_Poisson.get_coord_shift());
					normalize_density(f);
					u = f;
					u.zero();
					the_poisson_solver->solve(f, u);

					if (!bdefd)
						f.deallocate();
				}
				grid_hierarchy data_forIO(u);
				for (int icoord = 0; icoord < 3; ++icoord)
				{
					//... displacement
					if (bdefd)
					{
						data_forIO.zero();
						*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
						poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order,
													 data_forIO.levelmin() == data_forIO.levelmax(), decic_baryons);
						*data_forIO.get_grid(data_forIO.levelmax()) /= 1 << f.levelmax();
						the_poisson_solver->gradient_add(icoord, u, data_forIO);
					}
					else
						the_poisson_solver->gradient(icoord, u, data_forIO);

					//... multiply to get velocity
					data_forIO *= cosmo_vfact;

					//... velocity kick to keep refined region centered?

					double sigv = compute_finest_sigma(data_forIO);
					music::ulog.Print("sigma of %c-velocity of high-res particles is %f", 'x' + icoord, sigv);

					double meanv = compute_finest_mean(data_forIO);
					music::ulog.Print("mean of %c-velocity of high-res particles is %f", 'x' + icoord, meanv);

					double maxv = compute_finest_absmax(data_forIO);
					music::ulog.Print("max of abs of %c-velocity of high-res particles is %f", 'x' + icoord, maxv);

					coarsen_density(rh_Poisson, data_forIO, false);

					// add counter velocity-mode
					if( do_counter_mode ) add_constant_value( data_forIO, -counter_mode_amp[icoord]*cosmo_vfact );

					music::ulog.Print("Writing CDM velocities");
					the_output_plugin->write_dm_velocity(icoord, data_forIO);

					if (do_baryons)
					{
						music::ulog.Print("Writing baryon velocities");
						the_output_plugin->write_gas_velocity(icoord, data_forIO);
					}
				}

				u.deallocate();
				data_forIO.deallocate();
			}
			else
			{
				music::ilog.Print("Computing separate velocities for CDM and baryons:");
				music::ilog << "===============================================================================" << std::endl;
				music::ilog << "   COMPUTING DARK MATTER VELOCITIES" << std::endl;
				music::ilog << "-------------------------------------------------------------------------------" << std::endl;
				music::ulog.Print("Computing dark matter velocitites...");

				//... we do baryons and have velocity transfer functions, or we do SPH and not to shift
				//... do DM first
				GenerateDensityHierarchy(cf, the_cosmo_calc.get(), theta_cdm, rh_TF, rand, f, false, false);
				coarsen_density(rh_Poisson, f, use_fourier_coarsening);
				f.add_refinement_mask(rh_Poisson.get_coord_shift());
				normalize_density(f);

				u = f;
				u.zero();

				the_poisson_solver->solve(f, u);

				if (!bdefd)
					f.deallocate();

				grid_hierarchy data_forIO(u);
				for (int icoord = 0; icoord < 3; ++icoord)
				{
					//... displacement
					if (bdefd)
					{
						data_forIO.zero();
						*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
						poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order,
													 data_forIO.levelmin() == data_forIO.levelmax(), decic_DM);
						*data_forIO.get_grid(data_forIO.levelmax()) /= 1 << f.levelmax();
						the_poisson_solver->gradient_add(icoord, u, data_forIO);
					}
					else
						the_poisson_solver->gradient(icoord, u, data_forIO);

					//... multiply to get velocity
					data_forIO *= cosmo_vfact;

					double sigv = compute_finest_sigma(data_forIO);
					music::ulog.Print("sigma of %c-velocity of high-res DM is %f", 'x' + icoord, sigv);

					double meanv = compute_finest_mean(data_forIO);
					music::ulog.Print("mean of %c-velocity of high-res particles is %f", 'x' + icoord, meanv);

					double maxv = compute_finest_absmax(data_forIO);
					music::ulog.Print("max of abs of %c-velocity of high-res particles is %f", 'x' + icoord, maxv);

					coarsen_density(rh_Poisson, data_forIO, false);

					// add counter velocity mode
					if( do_counter_mode ) add_constant_value( data_forIO, -counter_mode_amp[icoord]*cosmo_vfact );

					music::ulog.Print("Writing CDM velocities");
					the_output_plugin->write_dm_velocity(icoord, data_forIO);
				}
				u.deallocate();
				data_forIO.deallocate();
				f.deallocate();

				music::ilog << "===============================================================================" << std::endl;
				music::ilog << "   COMPUTING BARYON VELOCITIES" << std::endl;
				music::ilog << "-------------------------------------------------------------------------------" << std::endl;
				music::ulog.Print("Computing baryon velocitites...");
				//... do baryons
				GenerateDensityHierarchy(cf, the_cosmo_calc.get(), theta_baryon, rh_TF, rand, f, false, bbshift);
				coarsen_density(rh_Poisson, f, use_fourier_coarsening);
				f.add_refinement_mask(rh_Poisson.get_coord_shift());
				normalize_density(f);

				u = f;
				u.zero();

				the_poisson_solver->solve(f, u);

				if (!bdefd)
					f.deallocate();

				data_forIO = u;
				for (int icoord = 0; icoord < 3; ++icoord)
				{
					//... displacement
					if (bdefd)
					{
						data_forIO.zero();
						*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
						poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order,
													 data_forIO.levelmin() == data_forIO.levelmax(), decic_baryons);
						*data_forIO.get_grid(data_forIO.levelmax()) /= 1 << f.levelmax();
						the_poisson_solver->gradient_add(icoord, u, data_forIO);
					}
					else
						the_poisson_solver->gradient(icoord, u, data_forIO);

					//... multiply to get velocity
					data_forIO *= cosmo_vfact;

					double sigv = compute_finest_sigma(data_forIO);
					music::ulog.Print("sigma of %c-velocity of high-res baryons is %f", 'x' + icoord, sigv);

					double meanv = compute_finest_mean(data_forIO);
					music::ulog.Print("mean of %c-velocity of high-res baryons is %f", 'x' + icoord, meanv);

					double maxv = compute_finest_absmax(data_forIO);
					music::ulog.Print("max of abs of %c-velocity of high-res baryons is %f", 'x' + icoord, maxv);

					coarsen_density(rh_Poisson, data_forIO, false);

					// add counter velocity mode
					if( do_counter_mode ) add_constant_value( data_forIO, -counter_mode_amp[icoord]*cosmo_vfact );

					music::ulog.Print("Writing baryon velocities");
					the_output_plugin->write_gas_velocity(icoord, data_forIO);
				}
				u.deallocate();
				f.deallocate();
				data_forIO.deallocate();
			}
			/*********************************************************************************************/
			/*********************************************************************************************/
			/*** 2LPT ************************************************************************************/
			/*********************************************************************************************/
		}
		else
		{
			//.. use 2LPT ...
			music::ulog.Print("Entering 2LPT branch");

			grid_hierarchy f(nbnd), u1(nbnd), u2LPT(nbnd), f2LPT(nbnd);

			tf_type my_tf_type = theta_cdm;
			bool dm_only = !do_baryons;
			if (!do_baryons || !tf_has_velocities)
				my_tf_type = theta_matter;

			music::ilog << "===============================================================================" << std::endl;
			if (my_tf_type == theta_matter)
			{
				music::ilog << "   COMPUTING VELOCITIES" << std::endl;
			}
			else
			{
				music::ilog << "   COMPUTING DARK MATTER VELOCITIES" << std::endl;
			}
			music::ilog << "-------------------------------------------------------------------------------" << std::endl;

			GenerateDensityHierarchy(cf, the_cosmo_calc.get(), my_tf_type, rh_TF, rand, f, false, false);
			coarsen_density(rh_Poisson, f, use_fourier_coarsening);
			f.add_refinement_mask(rh_Poisson.get_coord_shift());
			normalize_density(f);

			if (dm_only)
			{
				the_output_plugin->write_dm_density(f);
				the_output_plugin->write_dm_mass(f);
			}

			u1 = f;
			u1.zero();

			//... compute 1LPT term
			the_poisson_solver->solve(f, u1);

			//... compute 2LPT term
			if (bdefd)
				f2LPT = f;
			else
				f.deallocate();

			music::ilog.Print("- Computing 2LPT term....");
			if (!kspace2LPT)
				compute_2LPT_source(u1, f2LPT, grad_order);
			else
			{
				music::ulog.Print("  computing term using FFT");
				compute_2LPT_source_FFT(cf, u1, f2LPT);
			}

			music::ilog.Print("- Solving 2LPT Poisson equation");
			u2LPT = u1;
			u2LPT.zero();
			the_poisson_solver->solve(f2LPT, u2LPT);

			//... if doing the hybrid step, we need a combined source term
			if (bdefd)
			{
				f2LPT *= 6.0 / 7.0 / vfac2lpt;
				f += f2LPT;

				if (!dm_only)
					f2LPT.deallocate();
			}

			//... add the 2LPT contribution
			u2LPT *= 6.0 / 7.0 / vfac2lpt;
			u1 += u2LPT;

			grid_hierarchy data_forIO(u1);
			for (int icoord = 0; icoord < 3; ++icoord)
			{
				if (bdefd)
				{
					data_forIO.zero();
					*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
					poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order,
												 data_forIO.levelmin() == data_forIO.levelmax(), decic_DM);
					*data_forIO.get_grid(data_forIO.levelmax()) /= (1 << f.levelmax());
					the_poisson_solver->gradient_add(icoord, u1, data_forIO);
				}
				else
					the_poisson_solver->gradient(icoord, u1, data_forIO);

				data_forIO *= cosmo_vfact;

				double sigv = compute_finest_sigma(data_forIO);

				double meanv = compute_finest_mean(data_forIO);
				music::ulog.Print("mean of %c-velocity of high-res particles is %f", 'x' + icoord, meanv);

				double maxv = compute_finest_absmax(data_forIO);
				music::ulog.Print("max of abs of %c-velocity of high-res particles is %f", 'x' + icoord, maxv);

				music::ilog << "\t - velocity component " << icoord << " : sigma = " << sigv << std::endl;
				music::ilog << "\t - velocity component " << icoord << " : mean = " << meanv << std::endl;

				coarsen_density(rh_Poisson, data_forIO, false);

				//... compute counter-mode to minimize advection errors
				counter_mode_amp[icoord] = compute_finest_mean(data_forIO); 
				if( do_counter_mode ) add_constant_value( data_forIO, -counter_mode_amp[icoord] );

				music::ulog.Print("Writing CDM velocities");
				the_output_plugin->write_dm_velocity(icoord, data_forIO);

				if (do_baryons && !tf_has_velocities && !bsph)
				{
					music::ulog.Print("Writing baryon velocities");
					the_output_plugin->write_gas_velocity(icoord, data_forIO);
				}
			}
			data_forIO.deallocate();
			if (!dm_only)
				u1.deallocate();

			if (do_baryons && (tf_has_velocities || bsph))
			{
				music::ilog << "===============================================================================" << std::endl;
				music::ilog << "   COMPUTING BARYON VELOCITIES" << std::endl;
				music::ilog << "-------------------------------------------------------------------------------" << std::endl;
				music::ulog.Print("Computing baryon displacements...");

				GenerateDensityHierarchy(cf, the_cosmo_calc.get(), theta_baryon, rh_TF, rand, f, false, bbshift);
				coarsen_density(rh_Poisson, f, use_fourier_coarsening);
				f.add_refinement_mask(rh_Poisson.get_coord_shift());
				normalize_density(f);

				u1 = f;
				u1.zero();

				if (bdefd)
					f2LPT = f;

				//... compute 1LPT term
				the_poisson_solver->solve(f, u1);

				music::ilog.Print("Writing baryon potential");
				the_output_plugin->write_gas_potential(u1);

				//... compute 2LPT term
				u2LPT = f;
				u2LPT.zero();

				if (!kspace2LPT)
					compute_2LPT_source(u1, f2LPT, grad_order);
				else
					compute_2LPT_source_FFT(cf, u1, f2LPT);

				the_poisson_solver->solve(f2LPT, u2LPT);

				//... if doing the hybrid step, we need a combined source term
				if (bdefd)
				{
					f2LPT *= 6.0 / 7.0 / vfac2lpt;
					f += f2LPT;

					f2LPT.deallocate();
				}

				//... add the 2LPT contribution
				u2LPT *= 6.0 / 7.0 / vfac2lpt;
				u1 += u2LPT;
				u2LPT.deallocate();

				// grid_hierarchy data_forIO(u1);
				data_forIO = u1;
				for (int icoord = 0; icoord < 3; ++icoord)
				{
					if (bdefd)
					{
						data_forIO.zero();
						*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
						poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order,
													 data_forIO.levelmin() == data_forIO.levelmax(), decic_baryons);
						*data_forIO.get_grid(data_forIO.levelmax()) /= (1 << f.levelmax());
						the_poisson_solver->gradient_add(icoord, u1, data_forIO);
					}
					else
						the_poisson_solver->gradient(icoord, u1, data_forIO);

					data_forIO *= cosmo_vfact;

					double sigv = compute_finest_sigma(data_forIO);

					double meanv = compute_finest_mean(data_forIO);
					music::ulog.Print("mean of %c-velocity of high-res baryons is %f", 'x' + icoord, meanv);

					double maxv = compute_finest_absmax(data_forIO);
					music::ulog.Print("max of abs of %c-velocity of high-res baryons is %f", 'x' + icoord, maxv);

					music::ilog << "\t - velocity component " << icoord << " : sigma = " << sigv << std::endl;
					music::ilog << "\t - velocity component " << icoord << " : mean = " << meanv << std::endl;

					coarsen_density(rh_Poisson, data_forIO, false);

					// add counter velocity mode
					if( do_counter_mode ) add_constant_value( data_forIO, -counter_mode_amp[icoord] );

					music::ulog.Print("Writing baryon velocities");
					the_output_plugin->write_gas_velocity(icoord, data_forIO);
				}
				data_forIO.deallocate();
				u1.deallocate();
			}

			music::ilog << "===============================================================================" << std::endl;
			music::ilog << "   COMPUTING DARK MATTER DISPLACEMENTS" << std::endl;
			music::ilog << "-------------------------------------------------------------------------------" << std::endl;
			music::ulog.Print("Computing dark matter displacements...");

			//... if baryons are enabled, the displacements have to be recomputed
			//... otherwise we can compute them directly from the velocities
			if (!dm_only)
			{
				// my_tf_type is cdm if do_baryons==true, total otherwise
				my_tf_type = delta_cdm;
				if (!do_baryons || !the_cosmo_calc->transfer_function_->tf_is_distinct())
					my_tf_type = delta_matter;

				GenerateDensityHierarchy(cf, the_cosmo_calc.get(), my_tf_type, rh_TF, rand, f, false, false);
				coarsen_density(rh_Poisson, f, use_fourier_coarsening);
				f.add_refinement_mask(rh_Poisson.get_coord_shift());
				normalize_density(f);

				music::ulog.Print("Writing CDM data");
				the_output_plugin->write_dm_density(f);
				the_output_plugin->write_dm_mass(f);
				u1 = f;
				u1.zero();

				if (bdefd)
					f2LPT = f;

				//... compute 1LPT term
				the_poisson_solver->solve(f, u1);

				//... compute 2LPT term
				u2LPT = f;
				u2LPT.zero();

				if (!kspace2LPT)
					compute_2LPT_source(u1, f2LPT, grad_order);
				else
					compute_2LPT_source_FFT(cf, u1, f2LPT);

				the_poisson_solver->solve(f2LPT, u2LPT);

				if (bdefd)
				{
					f2LPT *= 3.0 / 7.0;
					f += f2LPT;
					f2LPT.deallocate();
				}

				u2LPT *= 3.0 / 7.0;
				u1 += u2LPT;
				u2LPT.deallocate();
			}
			else
			{
				//... reuse prior data
				/*f-=f2LPT;
				the_output_plugin->write_dm_density(f);
				the_output_plugin->write_dm_mass(f);
				f+=f2LPT;*/

				u2LPT *= 0.5;
				u1 -= u2LPT;
				u2LPT.deallocate();

				if (bdefd)
				{
					f2LPT *= 0.5;
					f -= f2LPT;
					f2LPT.deallocate();
				}
			}

			data_forIO = u1;

			for (int icoord = 0; icoord < 3; ++icoord)
			{
				//... displacement
				if (bdefd)
				{
					data_forIO.zero();
					*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
					poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order,
												 data_forIO.levelmin() == data_forIO.levelmax(), decic_DM);
					*data_forIO.get_grid(data_forIO.levelmax()) /= 1 << f.levelmax();
					the_poisson_solver->gradient_add(icoord, u1, data_forIO);
				}
				else
					the_poisson_solver->gradient(icoord, u1, data_forIO);

				double dispmax = compute_finest_absmax(data_forIO);
				music::ilog.Print("\t - max. %c-displacement of HR particles is %f [mean dx]", 'x' + icoord, dispmax * (double)(1ll << data_forIO.levelmax()));

				coarsen_density(rh_Poisson, data_forIO, false);

				// add counter mode
				if( do_counter_mode ) add_constant_value( data_forIO, -counter_mode_amp[icoord]/cosmo_vfact );

				music::ulog.Print("Writing CDM displacements");
				the_output_plugin->write_dm_position(icoord, data_forIO);
			}

			data_forIO.deallocate();
			u1.deallocate();

			if (do_baryons && !bsph)
			{
				music::ilog << "===============================================================================" << std::endl;
				music::ilog << "   COMPUTING BARYON DENSITY" << std::endl;
				music::ilog << "-------------------------------------------------------------------------------" << std::endl;
				music::ulog.Print("Computing baryon density...");

				GenerateDensityHierarchy(cf, the_cosmo_calc.get(), delta_baryon, rh_TF, rand, f, true, false);
				coarsen_density(rh_Poisson, f, use_fourier_coarsening);
				f.add_refinement_mask(rh_Poisson.get_coord_shift());
				normalize_density(f);

				if (!do_LLA)
					the_output_plugin->write_gas_density(f);
				else
				{
					u1 = f;
					u1.zero();

					//... compute 1LPT term
					the_poisson_solver->solve(f, u1);

					//... compute 2LPT term
					u2LPT = f;
					u2LPT.zero();

					if (!kspace2LPT)
						compute_2LPT_source(u1, f2LPT, grad_order);
					else
						compute_2LPT_source_FFT(cf, u1, f2LPT);

					the_poisson_solver->solve(f2LPT, u2LPT);
					u2LPT *= 3.0 / 7.0;
					u1 += u2LPT;
					u2LPT.deallocate();

					compute_LLA_density(u1, f, grad_order);
					normalize_density(f);

					music::ulog.Print("Writing baryon density");
					the_output_plugin->write_gas_density(f);
				}
			}
			else if (do_baryons && bsph)
			{
				music::ilog << "===============================================================================" << std::endl;
				music::ilog << "   COMPUTING BARYON DISPLACEMENTS" << std::endl;
				music::ilog << "-------------------------------------------------------------------------------" << std::endl;
				music::ulog.Print("Computing baryon displacements...");

				GenerateDensityHierarchy(cf, the_cosmo_calc.get(), delta_baryon, rh_TF, rand, f, false, bbshift);
				coarsen_density(rh_Poisson, f, use_fourier_coarsening);
				f.add_refinement_mask(rh_Poisson.get_coord_shift());
				normalize_density(f);

				music::ulog.Print("Writing baryon density");
				the_output_plugin->write_gas_density(f);
				u1 = f;
				u1.zero();

				if (bdefd)
					f2LPT = f;

				//... compute 1LPT term
				the_poisson_solver->solve(f, u1);

				//... compute 2LPT term
				u2LPT = f;
				u2LPT.zero();

				if (!kspace2LPT)
					compute_2LPT_source(u1, f2LPT, grad_order);
				else
					compute_2LPT_source_FFT(cf, u1, f2LPT);

				the_poisson_solver->solve(f2LPT, u2LPT);

				if (bdefd)
				{
					f2LPT *= 3.0 / 7.0;
					f += f2LPT;
					f2LPT.deallocate();
				}

				u2LPT *= 3.0 / 7.0;
				u1 += u2LPT;
				u2LPT.deallocate();

				data_forIO = u1;

				for (int icoord = 0; icoord < 3; ++icoord)
				{
					//... displacement
					if (bdefd)
					{
						data_forIO.zero();
						*data_forIO.get_grid(data_forIO.levelmax()) = *f.get_grid(f.levelmax());
						poisson_hybrid(*data_forIO.get_grid(data_forIO.levelmax()), icoord, grad_order,
													 data_forIO.levelmin() == data_forIO.levelmax(), decic_baryons);
						*data_forIO.get_grid(data_forIO.levelmax()) /= 1 << f.levelmax();
						the_poisson_solver->gradient_add(icoord, u1, data_forIO);
					}
					else
						the_poisson_solver->gradient(icoord, u1, data_forIO);

					coarsen_density(rh_Poisson, data_forIO, false);

					// add counter mode
					if( do_counter_mode ) add_constant_value( data_forIO, -counter_mode_amp[icoord]/cosmo_vfact );


					music::ulog.Print("Writing baryon displacements");
					the_output_plugin->write_gas_position(icoord, data_forIO);
				}
			}
		}

		//------------------------------------------------------------------------------
		//... finish output
		//------------------------------------------------------------------------------

		the_output_plugin->finalize();
		delete the_output_plugin;
	}
	catch (std::runtime_error &excp)
	{
		music::elog.Print("Fatal error occured. Code will exit:");
		music::elog.Print("Exception: %s", excp.what());
		std::cerr << " - " << excp.what() << std::endl;
		std::cerr << " - A fatal error occured. We need to exit...\n";
		bfatal = true;
	}

	music::ilog << "===============================================================================" << std::endl;
	if (!bfatal)
	{
		music::ilog << " - Wrote output file \'" << outfname << "\'\n     using plugin \'" << outformat << "\'...\n";
		music::ulog.Print("Wrote output file \'%s\'.", outfname.c_str());
	}



	//------------------------------------------------------------------------------
	//... clean up
	//------------------------------------------------------------------------------
	// delete the_transfer_function_plugin;
	delete the_poisson_solver;

	if( CONFIG::FFTW_threads_ok )
		FFTW_API(cleanup_threads)();

	//------------------------------------------------------------------------------
	//... we are done !
	//------------------------------------------------------------------------------
	music::ilog << " - Done!" << std::endl << std::endl;

	ltime = time(NULL);

	music::ulog.Print("Run finished succesfully on %s", asctime(localtime(&ltime)));

	cf.dump_to_log();

	return 0;
}
