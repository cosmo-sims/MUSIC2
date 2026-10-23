/*

output_pasta.cc - Native Pasta output plugin for MUSIC
Generates initial conditions for the Pasta astrophysics simulation code.

Stores grid dimensions, cosmology parameters, scale factor a_i, box size,
and particle/gas datasets without legacy Enzo header requirements.
Automatically writes a Pasta-formatted parameter file (pasta.txt) appropriate
for either pure N-body (particles only) or hybrid (particles + baryons) runs.

Complies with Phase 1.4 of the Cosmological Galaxy Simulations Action Plan.

*/

#ifdef HAVE_HDF5

#include <sys/types.h>
#include <sys/stat.h>
#include <cerrno>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "output.hh"
#include "HDF_IO.hh"

#define MAX_SLAB_SIZE 268435456 // = 256 MBytes

class pasta_output_plugin : public output_plugin
{
protected:
	struct patch_header
	{
		std::vector<int> dimensions;
		std::vector<int> top_grid_dims;
		std::vector<int> top_grid_end;
		std::vector<int> top_grid_start;
	};

	struct sim_header
	{
		std::vector<int> dimensions;
		float a_start;
		float box_size;
		float dx;
		float h0;
		float omega_b;
		float omega_m;
		float omega_v;
	};

	sim_header the_sim_header;

	void write_sim_header(std::string fname, const sim_header &h)
	{
		HDFWriteGroupAttribute(fname, "/", "Dimensions", h.dimensions);
		HDFWriteGroupAttribute(fname, "/", "a_start", h.a_start);
		HDFWriteGroupAttribute(fname, "/", "box_size", h.box_size);
		HDFWriteGroupAttribute(fname, "/", "dx", h.dx);
		HDFWriteGroupAttribute(fname, "/", "h0", h.h0);
		HDFWriteGroupAttribute(fname, "/", "omega_b", h.omega_b);
		HDFWriteGroupAttribute(fname, "/", "omega_m", h.omega_m);
		HDFWriteGroupAttribute(fname, "/", "omega_v", h.omega_v);
	}

	void write_patch_header(std::string fname, std::string dsetname, const patch_header &h)
	{
		HDFWriteDatasetAttribute(fname, dsetname, "Dimensions", h.dimensions);
		HDFWriteDatasetAttribute(fname, dsetname, "TopGridDims", h.top_grid_dims);
		HDFWriteDatasetAttribute(fname, dsetname, "TopGridEnd", h.top_grid_end);
		HDFWriteDatasetAttribute(fname, dsetname, "TopGridStart", h.top_grid_start);
	}

	void dump_mask(const grid_hierarchy &gh)
	{
		char pastaname[256], filename[512];
		std::string fieldname("RefinementMask");

		for (unsigned ilevel = levelmin_; ilevel <= levelmax_; ++ilevel)
		{
			std::vector<int> ng;
			ng.push_back(gh.get_grid(ilevel)->size(0));
			ng.push_back(gh.get_grid(ilevel)->size(1));
			ng.push_back(gh.get_grid(ilevel)->size(2));

			size_t all_data_size = (size_t)ng[0] * (size_t)ng[1] * (size_t)ng[2];
			size_t max_slab_size = std::min((size_t)MAX_SLAB_SIZE / sizeof(double), all_data_size);
			max_slab_size = std::max(((size_t)ng[0] * (size_t)ng[1]), max_slab_size);

			size_t slices_in_slab = (size_t)((double)max_slab_size / ((size_t)ng[0] * (size_t)ng[1]));
			size_t nsz[3] = {size_t(ng[2]), size_t(ng[1]), size_t(ng[0])};

			if (levelmin_ != levelmax_)
				snprintf(pastaname, 256, "%s.%d", fieldname.c_str(), ilevel - levelmin_);
			else
				snprintf(pastaname, 256, "%s", fieldname.c_str());

			snprintf(filename, 512, "%s/%s", fname_.c_str(), pastaname);

			HDFCreateFile(filename);
			write_sim_header(filename, the_sim_header);

			HDFHyperslabWriter3Ds<int> *slab_writer = new HDFHyperslabWriter3Ds<int>(filename, pastaname, nsz);
			int *data_buf = new int[slices_in_slab * (size_t)ng[0] * (size_t)ng[1]];

			size_t slices_written = 0;
			while (slices_written < (size_t)ng[2])
			{
				slices_in_slab = std::min((size_t)ng[2] - slices_written, slices_in_slab);

#pragma omp parallel for
				for (int k = 0; k < (int)slices_in_slab; ++k)
					for (int j = 0; j < ng[1]; ++j)
						for (int i = 0; i < ng[0]; ++i)
						{
							int mask_val = -1;
							if (gh.is_in_mask(ilevel, i, j, k + slices_written))
							{
								if (gh.is_refined(ilevel, i, j, k + slices_written))
									mask_val = 1;
								else
									mask_val = 0;
							}
							data_buf[(size_t)(k * ng[1] + j) * (size_t)ng[0] + (size_t)i] = mask_val;
						}

				size_t count[3], offset[3];
				count[0] = slices_in_slab;
				count[1] = ng[1];
				count[2] = ng[0];

				offset[0] = slices_written;
				offset[1] = 0;
				offset[2] = 0;

				slab_writer->write_slab(data_buf, count, offset);
				slices_written += slices_in_slab;
			}

			delete[] data_buf;
			delete slab_writer;

			patch_header ph;
			ph.dimensions = ng;
			ph.top_grid_dims.assign(3, 1 << levelmin_);

			double rfac = 1.0 / (1 << (ilevel - levelmin_));
			ph.top_grid_start.push_back((int)(gh.offset_abs(ilevel, 0) * rfac));
			ph.top_grid_start.push_back((int)(gh.offset_abs(ilevel, 1) * rfac));
			ph.top_grid_start.push_back((int)(gh.offset_abs(ilevel, 2) * rfac));

			ph.top_grid_end.push_back(ph.top_grid_start[0] + (int)(ng[0] * rfac));
			ph.top_grid_end.push_back(ph.top_grid_start[1] + (int)(ng[1] * rfac));
			ph.top_grid_end.push_back(ph.top_grid_start[2] + (int)(ng[2] * rfac));

			write_patch_header(filename, pastaname, ph);
		}
	}

	void dump_grid_data(std::string fieldname, const grid_hierarchy &gh, double factor = 1.0, double add = 0.0)
	{
		char pastaname[256], filename[512];

		for (unsigned ilevel = levelmin_; ilevel <= levelmax_; ++ilevel)
		{
			std::vector<int> ng;
			ng.push_back(gh.get_grid(ilevel)->size(0));
			ng.push_back(gh.get_grid(ilevel)->size(1));
			ng.push_back(gh.get_grid(ilevel)->size(2));

			size_t all_data_size = (size_t)ng[0] * (size_t)ng[1] * (size_t)ng[2];
			size_t max_slab_size = std::min((size_t)MAX_SLAB_SIZE / sizeof(double), all_data_size);
			max_slab_size = std::max(((size_t)ng[0] * (size_t)ng[1]), max_slab_size);

			size_t slices_in_slab = (size_t)((double)max_slab_size / ((size_t)ng[0] * (size_t)ng[1]));
			size_t nsz[3] = {size_t(ng[2]), size_t(ng[1]), size_t(ng[0])};

			if (levelmin_ != levelmax_)
				snprintf(pastaname, 256, "%s.%d", fieldname.c_str(), ilevel - levelmin_);
			else
				snprintf(pastaname, 256, "%s", fieldname.c_str());

			snprintf(filename, 512, "%s/%s", fname_.c_str(), pastaname);

			HDFCreateFile(filename);
			write_sim_header(filename, the_sim_header);

			HDFHyperslabWriter3Ds<real_t> *slab_writer = new HDFHyperslabWriter3Ds<real_t>(filename, pastaname, nsz);
			real_t *data_buf = new real_t[slices_in_slab * (size_t)ng[0] * (size_t)ng[1]];

			size_t slices_written = 0;
			while (slices_written < (size_t)ng[2])
			{
				slices_in_slab = std::min((size_t)ng[2] - slices_written, slices_in_slab);

#pragma omp parallel for
				for (int k = 0; k < (int)slices_in_slab; ++k)
					for (int j = 0; j < ng[1]; ++j)
						for (int i = 0; i < ng[0]; ++i)
							data_buf[(size_t)(k * ng[1] + j) * (size_t)ng[0] + (size_t)i] =
									(add + (*gh.get_grid(ilevel))(i, j, k + slices_written)) * factor;

				size_t count[3], offset[3];
				count[0] = slices_in_slab;
				count[1] = ng[1];
				count[2] = ng[0];

				offset[0] = slices_written;
				offset[1] = 0;
				offset[2] = 0;

				slab_writer->write_slab(data_buf, count, offset);
				slices_written += slices_in_slab;
			}

			delete[] data_buf;
			delete slab_writer;

			patch_header ph;
			ph.dimensions = ng;
			ph.top_grid_dims.assign(3, 1 << levelmin_);

			double rfac = 1.0 / (1 << (ilevel - levelmin_));
			ph.top_grid_start.push_back((int)(gh.offset_abs(ilevel, 0) * rfac));
			ph.top_grid_start.push_back((int)(gh.offset_abs(ilevel, 1) * rfac));
			ph.top_grid_start.push_back((int)(gh.offset_abs(ilevel, 2) * rfac));

			ph.top_grid_end.push_back(ph.top_grid_start[0] + (int)(ng[0] * rfac));
			ph.top_grid_end.push_back(ph.top_grid_start[1] + (int)(ng[1] * rfac));
			ph.top_grid_end.push_back(ph.top_grid_start[2] + (int)(ng[2] * rfac));

			write_patch_header(filename, pastaname, ph);
		}
	}

	void write_pasta_parameter_file(const grid_hierarchy &gh)
	{
		bool bhave_hydro = cf_.get_value<bool>("setup", "baryons");
		bool is_zoom = (levelmin_ != levelmax_);
		unsigned nbase = 1 << levelmin_;
		int default_block = std::min(16, (int)nbase);
		int block_size = cf_.get_value_safe<int>("output", "pasta_block_size", default_block);
		std::string param_file_name = cf_.get_value_safe<std::string>("output", "pasta_param_file", "pasta.txt");

		char filename[512];
		snprintf(filename, 512, "%s/%s", fname_.c_str(), param_file_name.c_str());

		std::ofstream ofs(filename, std::ios::trunc);
		if (!ofs.is_open())
		{
			throw std::runtime_error("Error in pasta_output_plugin: cannot create parameter file " + std::string(filename));
		}

		ofs << "#======================================================================\n";
		ofs << "# Pasta Cosmological Simulation Parameter File\n";
		ofs << "# Generated automatically by MUSIC (pasta output plugin)\n";
		if (bhave_hydro)
			ofs << "# Type: Hybrid (Baryons + Dark Matter Particles)\n";
		else
			ofs << "# Type: Pure N-body (Dark Matter Only)\n";
		if (is_zoom)
			ofs << "# Mesh: Nested Zoom Hierarchy (levels " << levelmin_ << " to " << levelmax_ << ")\n";
		else
			ofs << "# Mesh: Unigrid (level " << levelmin_ << ")\n";
		ofs << "#======================================================================\n\n";

		ofs << "[problem]\n";
		ofs << "user/problem_name = cosmology_simulation\n\n";

		ofs << "[cosmology_simulation]\n";
		ofs << "cosmology_simulation/music_dir             = " << fname_ << "\n";
		if (bhave_hydro)
			ofs << "cosmology_simulation/initial_temperature   = -1.0\n";
		if (is_zoom)
		{
			ofs << "cosmology_simulation/mark_must_refine      = true\n";
			ofs << "cosmology_simulation/must_refine_min_level = 1\n";
		}
		else
		{
			ofs << "cosmology_simulation/mark_must_refine      = false\n";
		}
		ofs << "\n";

		ofs << "[mesh domain and resolution]\n";
		ofs << "mesh/nx1   = " << nbase << "\n";
		ofs << "mesh/x1min = 0.0\n";
		ofs << "mesh/x1max = 1.0\n\n";
		ofs << "mesh/nx2   = " << nbase << "\n";
		ofs << "mesh/x2min = 0.0\n";
		ofs << "mesh/x2max = 1.0\n\n";
		ofs << "mesh/nx3   = " << nbase << "\n";
		ofs << "mesh/x3min = 0.0\n";
		ofs << "mesh/x3max = 1.0\n\n";
		ofs << "mesh/nbd   = 2\n\n";

		ofs << "[size of each block]\n";
		ofs << "block/nx1 = " << block_size << "\n";
		ofs << "block/nx2 = " << block_size << "\n";
		ofs << "block/nx3 = " << block_size << "\n\n";

		if (is_zoom)
		{
			ofs << "[AMR settings]\n";
			ofs << "mesh/amr              = true\n";
			ofs << "mesh/amr_interval     = 4\n";
			ofs << "mesh/amr_max_level    = " << (levelmax_ - levelmin_) << "\n";
			ofs << "mesh/deref_threshold  = 3\n\n";
		}

		ofs << "[integrators]\n";
		ofs << "mesh/fluid_data           = " << (bhave_hydro ? "true" : "false") << "\n";
		ofs << "mesh/fluid_integrator     = " << (bhave_hydro ? "true" : "false") << "\n";
		ofs << "mesh/particle_data        = true\n";
		ofs << "mesh/particle_integrator  = true\n";
		ofs << "mesh/gravity_integrator   = true\n\n";

		ofs << "[particle]\n";
		ofs << "particle/type             = dark_matter\n";
		ofs << "particle/deposit          = cic\n";
		ofs << "particle/pusher           = kdk\n";
		ofs << "particle/self_gravity     = true\n";
		ofs << "particle/cfl              = 0.4\n\n";

		if (bhave_hydro)
		{
			ofs << "[EOS]\n";
			ofs << "eos/gamma = 1.666666666666667\n\n";
		}

		double omega_m = the_sim_header.omega_m;
		double omega_b = the_sim_header.omega_b;
		double omega_cdm = bhave_hydro ? (omega_m - omega_b) : omega_m;

		ofs << "[cosmology]\n";
		ofs << "cosmology/comoving_coordinates = true\n";
		ofs << "cosmology/hubble_constant_now  = " << the_sim_header.h0 << "\n";
		ofs << "cosmology/omega_matter_now     = " << omega_m << "\n";
		ofs << "cosmology/omega_lambda_now     = " << the_sim_header.omega_v << "\n";
		ofs << "cosmology/omega_radiation_now  = 0.0\n";
		ofs << "cosmology/omega_baryon_now     = " << omega_b << "\n";
		ofs << "cosmology/omega_cdm_now        = " << omega_cdm << "\n";
		ofs << "cosmology/comoving_box_size    = " << the_sim_header.box_size << "\n";
		ofs << "cosmology/initial_redshift     = " << cf_.get_value<double>("setup", "zstart") << "\n";
		ofs << "cosmology/final_redshift       = 0.0\n";
		ofs << "cosmology/max_expansion_rate   = 0.01\n\n";

		ofs << "[time related]\n";
		ofs << "time/cfl_num       = 0.3\n";
		ofs << "time/explicit_type = rk2\n";
		ofs << "time/stop_cycle    = 10\n\n";

		if (bhave_hydro)
		{
			ofs << "[fluid solver]\n";
			ofs << "fluid/solver = hllc\n";
			ofs << "fluid/order  = plm\n\n";
		}

		ofs << "[gravity]\n";
		ofs << "gravity/solver                = fft\n";
		ofs << "gravity/grav_bc               = periodic\n";
		if (bhave_hydro)
			ofs << "gravity/source                = fluid\n";
		ofs << "gravity/subtract_mean_density = true\n";
		ofs << "gravity/solve_per_substage    = true\n\n";

		ofs << "[boundary condition]\n";
		if (bhave_hydro)
		{
			ofs << "fluid/left_x1_bd    = periodic\n";
			ofs << "fluid/right_x1_bd   = periodic\n";
			ofs << "fluid/left_x2_bd    = periodic\n";
			ofs << "fluid/right_x2_bd   = periodic\n";
			ofs << "fluid/left_x3_bd    = periodic\n";
			ofs << "fluid/right_x3_bd   = periodic\n\n";
		}
		ofs << "gravity/left_x1_bd  = periodic\n";
		ofs << "gravity/right_x1_bd = periodic\n";
		ofs << "gravity/left_x2_bd  = periodic\n";
		ofs << "gravity/right_x2_bd = periodic\n";
		ofs << "gravity/left_x3_bd  = periodic\n";
		ofs << "gravity/right_x3_bd = periodic\n\n";

		ofs << "[output]\n";
		ofs << "output1/type = hst\n";
		ofs << "output1/dt   = 0.001\n";

		ofs.close();

		music::ilog.Print("Wrote Pasta parameter file '%s'", filename);
	}

public:
	pasta_output_plugin(config_file &cf)
			: output_plugin(cf)
	{
		if (mkdir(fname_.c_str(), 0777))
		{
			if (errno != EEXIST)
			{
				perror(fname_.c_str());
				throw std::runtime_error("Error in pasta_output_plugin: cannot create output directory!");
			}
		}

		bool bhave_hydro = cf_.get_value<bool>("setup", "baryons");

		the_sim_header.dimensions.push_back(1 << levelmin_);
		the_sim_header.dimensions.push_back(1 << levelmin_);
		the_sim_header.dimensions.push_back(1 << levelmin_);

		the_sim_header.a_start = 1.0 / (1.0 + cf.get_value<double>("setup", "zstart"));
		the_sim_header.box_size = cf.get_value<double>("setup", "boxlength");
		the_sim_header.h0 = cf.get_value<double>("cosmology", "H0") * 0.01;
		the_sim_header.dx = the_sim_header.box_size / the_sim_header.dimensions[0] / the_sim_header.h0;

		if (bhave_hydro)
			the_sim_header.omega_b = cf.get_value<double>("cosmology", "Omega_b");
		else
			the_sim_header.omega_b = 0.0;

		the_sim_header.omega_m = cf.get_value<double>("cosmology", "Omega_m");
		the_sim_header.omega_v = cf.get_value<double>("cosmology", "Omega_L");
	}

	~pasta_output_plugin()
	{
	}

	void write_dm_mass(const grid_hierarchy &gh)
	{
		/* DM mass assigned from grid and cosmological parameters */
	}

	void write_dm_density(const grid_hierarchy &gh)
	{
		dump_mask(gh);
		write_pasta_parameter_file(gh);

		// Determine density maximum and minimum location
		real_t rhomax = -1e30, rhomin = 1e30;
		double loc_rhomax[3] = {0.0, 0.0, 0.0}, loc_rhomin[3] = {0.0, 0.0, 0.0};
		int lvl_rhomax = 0, lvl_rhomin = 0;
		real_t rhomax_lm = -1e30, rhomin_lm = 1e30;
		double loc_rhomax_lm[3] = {0.0, 0.0, 0.0}, loc_rhomin_lm[3] = {0.0, 0.0, 0.0};

		for (int ilevel = gh.levelmax(); ilevel >= (int)gh.levelmin(); --ilevel)
			for (unsigned i = 0; i < gh.get_grid(ilevel)->size(0); ++i)
				for (unsigned j = 0; j < gh.get_grid(ilevel)->size(1); ++j)
					for (unsigned k = 0; k < gh.get_grid(ilevel)->size(2); ++k)
						if (!gh.is_refined(ilevel, i, j, k))
						{
							real_t rho = (*gh.get_grid(ilevel))(i, j, k);

							if (rho > rhomax)
							{
								rhomax = rho;
								lvl_rhomax = ilevel;
								gh.cell_pos(ilevel, i, j, k, loc_rhomax);
							}

							if (rho < rhomin)
							{
								rhomin = rho;
								lvl_rhomin = ilevel;
								gh.cell_pos(ilevel, i, j, k, loc_rhomin);
							}

							if (ilevel == (int)gh.levelmax())
							{
								if (rho > rhomax_lm)
								{
									rhomax_lm = rho;
									gh.cell_pos(ilevel, i, j, k, loc_rhomax_lm);
								}

								if (rho < rhomin_lm)
								{
									rhomin_lm = rho;
									gh.cell_pos(ilevel, i, j, k, loc_rhomin_lm);
								}
							}
						}

		double h = 1.0 / (1 << levelmin_);
		double shift[3];
		shift[0] = -(double)cf_.get_value<int>("setup", "shift_x") * h;
		shift[1] = -(double)cf_.get_value<int>("setup", "shift_y") * h;
		shift[2] = -(double)cf_.get_value<int>("setup", "shift_z") * h;

		if (gh.levelmin() != gh.levelmax())
		{
			music::ilog.Print("Global density extrema: ");
			music::ilog.Print("  minimum: delta=%f at (%f,%f,%f) (level=%d)", rhomin, loc_rhomin[0], loc_rhomin[1], loc_rhomin[2], lvl_rhomin);
			music::ilog.Print("       shifted back at (%f,%f,%f)", loc_rhomin[0] + shift[0], loc_rhomin[1] + shift[1], loc_rhomin[2] + shift[2]);
			music::ilog.Print("  maximum: delta=%f at (%f,%f,%f) (level=%d)", rhomax, loc_rhomax[0], loc_rhomax[1], loc_rhomax[2], lvl_rhomax);
			music::ilog.Print("       shifted back at (%f,%f,%f)", loc_rhomax[0] + shift[0], loc_rhomax[1] + shift[1], loc_rhomax[2] + shift[2]);

			music::ilog.Print("Density extrema on finest level: ");
			music::ilog.Print("  minimum: delta=%f at (%f,%f,%f)", rhomin_lm, loc_rhomin_lm[0], loc_rhomin_lm[1], loc_rhomin_lm[2]);
			music::ilog.Print("       shifted back at (%f,%f,%f)", loc_rhomin_lm[0] + shift[0], loc_rhomin_lm[1] + shift[1], loc_rhomin_lm[2] + shift[2]);
			music::ilog.Print("  maximum: delta=%f at (%f,%f,%f)", rhomax_lm, loc_rhomax_lm[0], loc_rhomax_lm[1], loc_rhomax_lm[2]);
			music::ilog.Print("       shifted back at (%f,%f,%f)", loc_rhomax_lm[0] + shift[0], loc_rhomax_lm[1] + shift[1], loc_rhomax_lm[2] + shift[2]);
		}
		else
		{
			music::ilog.Print("Global density extrema: ");
			music::ilog.Print("  minimum: delta=%f at (%f,%f,%f)", rhomin, loc_rhomin[0], loc_rhomin[1], loc_rhomin[2]);
			music::ilog.Print("       shifted back at (%f,%f,%f)", loc_rhomin[0] + shift[0], loc_rhomin[1] + shift[1], loc_rhomin[2] + shift[2]);
			music::ilog.Print("  maximum: delta=%f at (%f,%f,%f)", rhomax, loc_rhomax[0], loc_rhomax[1], loc_rhomax[2]);
			music::ilog.Print("       shifted back at (%f,%f,%f)", loc_rhomax[0] + shift[0], loc_rhomax[1] + shift[1], loc_rhomax[2] + shift[2]);
		}
	}

	void write_dm_velocity(int coord, const grid_hierarchy &gh)
	{
		char pastaname[256];
		snprintf(pastaname, 256, "ParticleVelocities_%c", (char)('x' + coord));

		double vunit = 1.0 / (1.225e2 * sqrt(the_sim_header.omega_m / the_sim_header.a_start));
		dump_grid_data(pastaname, gh, vunit);
	}

	void write_dm_position(int coord, const grid_hierarchy &gh)
	{
		char pastaname[256];
		snprintf(pastaname, 256, "ParticleDisplacements_%c", (char)('x' + coord));

		dump_grid_data(pastaname, gh);
	}

	void write_dm_potential(const grid_hierarchy &gh)
	{
	}

	void write_gas_potential(const grid_hierarchy &gh)
	{
	}

	void write_gas_velocity(int coord, const grid_hierarchy &gh)
	{
		double vunit = 1.0 / (1.225e2 * sqrt(the_sim_header.omega_m / the_sim_header.a_start));

		char pastaname[256];
		snprintf(pastaname, 256, "GridVelocities_%c", (char)('x' + coord));
		dump_grid_data(pastaname, gh, vunit);
	}

	void write_gas_position(int coord, const grid_hierarchy &gh)
	{
	}

	void write_gas_density(const grid_hierarchy &gh)
	{
		char pastaname[256];
		snprintf(pastaname, 256, "GridDensity");
		dump_grid_data(pastaname, gh, the_sim_header.omega_b / the_sim_header.omega_m, 1.0);
	}

	void finalize(void)
	{
	}
};

namespace
{
output_plugin_creator_concrete<pasta_output_plugin> creator("pasta");
}

#endif
