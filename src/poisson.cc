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

/****** ABSTRACT FACTORY PATTERN IMPLEMENTATION *******/

#include <poisson.hh>
#include <Numerics.hh>

std::map<std::string, poisson_plugin_creator *> &
get_poisson_plugin_map()
{
	static std::map<std::string, poisson_plugin_creator *> poisson_plugin_map;
	return poisson_plugin_map;
}

void print_poisson_plugins()
{
	std::map<std::string, poisson_plugin_creator *> &m = get_poisson_plugin_map();

	std::map<std::string, poisson_plugin_creator *>::iterator it;
	it = m.begin();
	std::cout << " - Available poisson solver plug-ins:\n";
	while (it != m.end())
	{
		if ((*it).second)
			std::cout << "\t\'" << (*it).first << "\'\n";

		music::ilog.Print("Poisson plug-in :: %s", std::string((*it).first).c_str());

		++it;
	}
}

/****** CALL IMPLEMENTATIONS OF POISSON SOLVER CLASSES ******/

#include <mg_solver.hh>
#include <fd_schemes.hh>


typedef multigrid::solver<stencil_7P, interp_O3_fluxcorr, mg_straight> poisson_solver_O2;
typedef multigrid::solver<stencil_13P, interp_O5_fluxcorr, mg_straight> poisson_solver_O4;
typedef multigrid::solver<stencil_19P, interp_O7_fluxcorr, mg_straight> poisson_solver_O6;

/**************************************************************************************/
/**************************************************************************************/

real_t multigrid_poisson_plugin::solve(grid_hierarchy &f, grid_hierarchy &u)
{
	music::ulog.Print("Initializing multi-grid Poisson solver...");

	unsigned verbosity = cf_.get_value_safe<unsigned>("setup", "verbosity", 2);

	if (verbosity > 0)
	{
		music::ilog << "-------------------------------------------------------------------------------" << std::endl;
		music::ilog << "- Invoking multi-grid Poisson solver..." << std::endl;
	}

	real_t acc = 1e-5, err;
	std::string ps_smoother_name;
	unsigned ps_presmooth, ps_postsmooth, order;

	acc = cf_.get_value_safe<real_t>("poisson", "accuracy", acc);
	ps_presmooth = cf_.get_value_safe<unsigned>("poisson", "pre_smooth", 3);
	ps_postsmooth = cf_.get_value_safe<unsigned>("poisson", "post_smooth", 3);
	ps_smoother_name = cf_.get_value_safe<std::string>("poisson", "smoother", "gs");
	order = cf_.get_value_safe<unsigned>("poisson", "laplace_order", 4);

	multigrid::opt::smtype ps_smtype = multigrid::opt::sm_gauss_seidel;

	if (ps_smoother_name == std::string("gs"))
	{
		ps_smtype = multigrid::opt::sm_gauss_seidel;
		music::ulog.Print("Selected Gauss-Seidel multigrid smoother");
	}
	else if (ps_smoother_name == std::string("jacobi"))
	{
		ps_smtype = multigrid::opt::sm_jacobi;
		music::ulog.Print("Selected Jacobi multigrid smoother");
	}
	else if (ps_smoother_name == std::string("sor"))
	{
		ps_smtype = multigrid::opt::sm_sor;
		music::ulog.Print("Selected SOR multigrid smoother");
	}
	else
	{
		music::wlog.Print("Unknown multigrid smoother \'%s\' specified. Reverting to Gauss-Seidel.", ps_smoother_name.c_str());
		std::cerr << " - Warning: unknown smoother \'" << ps_smoother_name << "\' for multigrid solver!\n"
							<< "            reverting to \'gs\' (Gauss-Seidel)" << std::endl;
	}

	real_t tstart, tend;

#if defined(_OPENMP)
	tstart = omp_get_wtime();
#else
	tstart = (real_t)clock() / CLOCKS_PER_SEC;
#endif

	//----- run Poisson solver -----//
	if (order == 2)
	{
		music::ulog.Print("Running multigrid solver with 2nd order Laplacian...");
		poisson_solver_O2 ps(f, ps_smtype, ps_presmooth, ps_postsmooth);
		err = ps.solve(u, acc, true);
	}
	else if (order == 4)
	{
		music::ulog.Print("Running multigrid solver with 4th order Laplacian...");
		poisson_solver_O4 ps(f, ps_smtype, ps_presmooth, ps_postsmooth);
		err = ps.solve(u, acc, true);
	}
	else if (order == 6)
	{
		music::ulog.Print("Running multigrid solver with 6th order Laplacian..");
		poisson_solver_O6 ps(f, ps_smtype, ps_presmooth, ps_postsmooth);
		err = ps.solve(u, acc, true);
	}
	else
	{
		music::elog.Print("Invalid order specified for Laplace operator");
		throw std::runtime_error("Invalid order specified for Laplace operator");
	}

	//------------------------------//

#if defined(_OPENMP)
	tend = omp_get_wtime();
	if (verbosity > 1)
		std::cout << " - Poisson solver took " << tend - tstart << "s with " << omp_get_max_threads() << " threads." << std::endl;
#else
	tend = (real_t)clock() / CLOCKS_PER_SEC;
	if (verbosity > 1)
		std::cout << " - Poisson solver took " << tend - tstart << "s." << std::endl;

#endif

	return err;
}

real_t multigrid_poisson_plugin::gradient(int dir, grid_hierarchy &u, grid_hierarchy &Du)
{
	Du = u;

	unsigned order = cf_.get_value_safe<unsigned>("poisson", "grad_order", 4);

	switch (order)
	{
	case 2:
		implementation().gradient_O2(dir, u, Du);
		break;
	case 4:
		implementation().gradient_O4(dir, u, Du);
		break;
	case 6:
		implementation().gradient_O6(dir, u, Du);
		break;
	default:
		music::elog.Print("Invalid order %d specified for gradient operator", order);
		throw std::runtime_error("Invalid order specified for gradient operator!");
	}

	return 0.0;
}

real_t multigrid_poisson_plugin::gradient_add(int dir, grid_hierarchy &u, grid_hierarchy &Du)
{
	// Du = u;

	unsigned order = cf_.get_value_safe<unsigned>("poisson", "grad_order", 4);

	switch (order)
	{
	case 2:
		implementation().gradient_add_O2(dir, u, Du);
		break;
	case 4:
		implementation().gradient_add_O4(dir, u, Du);
		break;
	case 6:
		implementation().gradient_add_O6(dir, u, Du);
		break;
	default:
		music::elog.Print("Invalid order %d specified for gradient operator!", order);
		throw std::runtime_error("Invalid order specified for gradient operator!");
	}

	return 0.0;
}

void multigrid_poisson_plugin::implementation::gradient_O2(int dir, grid_hierarchy &u, grid_hierarchy &Du)
{
	music::ulog.Print("Computing a 2nd order finite difference gradient...");

	for (unsigned ilevel = u.levelmin(); ilevel <= u.levelmax(); ++ilevel)
	{
		real_t h = pow(2.0, ilevel);
		meshvar_bnd *pvar = Du.get_grid(ilevel);

		if (dir == 0)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) = 0.5 * ((*u.get_grid(ilevel))(ix + 1, iy, iz) - (*u.get_grid(ilevel))(ix - 1, iy, iz)) * h;

		else if (dir == 1)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) = 0.5 * ((*u.get_grid(ilevel))(ix, iy + 1, iz) - (*u.get_grid(ilevel))(ix, iy - 1, iz)) * h;

		else if (dir == 2)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) = 0.5 * ((*u.get_grid(ilevel))(ix, iy, iz + 1) - (*u.get_grid(ilevel))(ix, iy, iz - 1)) * h;
	}

	music::ulog.Print("Done computing a 2nd order finite difference gradient.");
}

void multigrid_poisson_plugin::implementation::gradient_add_O2(int dir, grid_hierarchy &u, grid_hierarchy &Du)
{
	music::ulog.Print("Computing a 2nd order finite difference gradient...");

	for (unsigned ilevel = u.levelmin(); ilevel <= u.levelmax(); ++ilevel)
	{
		real_t h = pow(2.0, ilevel);
		meshvar_bnd *pvar = Du.get_grid(ilevel);

		if (dir == 0)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*Du.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*Du.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*Du.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) += 0.5 * ((*u.get_grid(ilevel))(ix + 1, iy, iz) - (*u.get_grid(ilevel))(ix - 1, iy, iz)) * h;

		else if (dir == 1)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*Du.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*Du.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*Du.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) += 0.5 * ((*u.get_grid(ilevel))(ix, iy + 1, iz) - (*u.get_grid(ilevel))(ix, iy - 1, iz)) * h;

		else if (dir == 2)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*Du.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*Du.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*Du.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) += 0.5 * ((*u.get_grid(ilevel))(ix, iy, iz + 1) - (*u.get_grid(ilevel))(ix, iy, iz - 1)) * h;
	}

	music::ulog.Print("Done computing a 4th order finite difference gradient.");
}

void multigrid_poisson_plugin::implementation::gradient_O4(int dir, grid_hierarchy &u, grid_hierarchy &Du)
{
	music::ulog.Print("Computing a 4th order finite difference gradient...");

	for (unsigned ilevel = u.levelmin(); ilevel <= u.levelmax(); ++ilevel)
	{
		real_t h = pow(2.0, ilevel);
		meshvar_bnd *pvar = Du.get_grid(ilevel);

		h /= 12.0;

		if (dir == 0)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) = ((*u.get_grid(ilevel))(ix - 2, iy, iz) - 8.0 * (*u.get_grid(ilevel))(ix - 1, iy, iz) + 8.0 * (*u.get_grid(ilevel))(ix + 1, iy, iz) - (*u.get_grid(ilevel))(ix + 2, iy, iz)) * h;

		else if (dir == 1)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) = ((*u.get_grid(ilevel))(ix, iy - 2, iz) - 8.0 * (*u.get_grid(ilevel))(ix, iy - 1, iz) + 8.0 * (*u.get_grid(ilevel))(ix, iy + 1, iz) - (*u.get_grid(ilevel))(ix, iy + 2, iz)) * h;

		else if (dir == 2)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) = ((*u.get_grid(ilevel))(ix, iy, iz - 2) - 8.0 * (*u.get_grid(ilevel))(ix, iy, iz - 1) + 8.0 * (*u.get_grid(ilevel))(ix, iy, iz + 1) - (*u.get_grid(ilevel))(ix, iy, iz + 2)) * h;
	}

	music::ulog.Print("Done computing a 4th order finite difference gradient.");
}

void multigrid_poisson_plugin::implementation::gradient_add_O4(int dir, grid_hierarchy &u, grid_hierarchy &Du)
{
	music::ulog.Print("Computing a 4th order finite difference gradient...");

	for (unsigned ilevel = u.levelmin(); ilevel <= u.levelmax(); ++ilevel)
	{
		real_t h = pow(2.0, ilevel);
		meshvar_bnd *pvar = Du.get_grid(ilevel);

		h /= 12.0;

		if (dir == 0)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) += ((*u.get_grid(ilevel))(ix - 2, iy, iz) - 8.0 * (*u.get_grid(ilevel))(ix - 1, iy, iz) + 8.0 * (*u.get_grid(ilevel))(ix + 1, iy, iz) - (*u.get_grid(ilevel))(ix + 2, iy, iz)) * h;

		else if (dir == 1)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) += ((*u.get_grid(ilevel))(ix, iy - 2, iz) - 8.0 * (*u.get_grid(ilevel))(ix, iy - 1, iz) + 8.0 * (*u.get_grid(ilevel))(ix, iy + 1, iz) - (*u.get_grid(ilevel))(ix, iy + 2, iz)) * h;

		else if (dir == 2)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) += ((*u.get_grid(ilevel))(ix, iy, iz - 2) - 8.0 * (*u.get_grid(ilevel))(ix, iy, iz - 1) + 8.0 * (*u.get_grid(ilevel))(ix, iy, iz + 1) - (*u.get_grid(ilevel))(ix, iy, iz + 2)) * h;
	}

	music::ulog.Print("Done computing a 4th order finite difference gradient.");
}

void multigrid_poisson_plugin::implementation::gradient_O6(int dir, grid_hierarchy &u, grid_hierarchy &Du)
{
	music::ulog.Print("Computing a 6th order finite difference gradient...");

	for (unsigned ilevel = u.levelmin(); ilevel <= u.levelmax(); ++ilevel)
	{
		real_t h = pow(2.0, ilevel);
		meshvar_bnd *pvar = Du.get_grid(ilevel);

		h /= 60.;
		if (dir == 0)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) =
								(-(*u.get_grid(ilevel))(ix - 3, iy, iz) + 9.0 * (*u.get_grid(ilevel))(ix - 2, iy, iz) - 45.0 * (*u.get_grid(ilevel))(ix - 1, iy, iz) + 45.0 * (*u.get_grid(ilevel))(ix + 1, iy, iz) - 9.0 * (*u.get_grid(ilevel))(ix + 2, iy, iz) + (*u.get_grid(ilevel))(ix + 3, iy, iz)) * h;

		else if (dir == 1)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) =
								(-(*u.get_grid(ilevel))(ix, iy - 3, iz) + 9.0 * (*u.get_grid(ilevel))(ix, iy - 2, iz) - 45.0 * (*u.get_grid(ilevel))(ix, iy - 1, iz) + 45.0 * (*u.get_grid(ilevel))(ix, iy + 1, iz) - 9.0 * (*u.get_grid(ilevel))(ix, iy + 2, iz) + (*u.get_grid(ilevel))(ix, iy + 3, iz)) * h;

		else if (dir == 2)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) =
								(-(*u.get_grid(ilevel))(ix, iy, iz - 3) + 9.0 * (*u.get_grid(ilevel))(ix, iy, iz - 2) - 45.0 * (*u.get_grid(ilevel))(ix, iy, iz - 1) + 45.0 * (*u.get_grid(ilevel))(ix, iy, iz + 1) - 9.0 * (*u.get_grid(ilevel))(ix, iy, iz + 2) + (*u.get_grid(ilevel))(ix, iy, iz + 3)) * h;
	}

	music::ulog.Print("Done computing a 6th order finite difference gradient.");
}

void multigrid_poisson_plugin::implementation::gradient_add_O6(int dir, grid_hierarchy &u, grid_hierarchy &Du)
{
	music::ulog.Print("Computing a 6th order finite difference gradient...");

	for (unsigned ilevel = u.levelmin(); ilevel <= u.levelmax(); ++ilevel)
	{
		real_t h = pow(2.0, ilevel);
		meshvar_bnd *pvar = Du.get_grid(ilevel);

		h /= 60.;
		if (dir == 0)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) +=
								(-(*u.get_grid(ilevel))(ix - 3, iy, iz) + 9.0 * (*u.get_grid(ilevel))(ix - 2, iy, iz) - 45.0 * (*u.get_grid(ilevel))(ix - 1, iy, iz) + 45.0 * (*u.get_grid(ilevel))(ix + 1, iy, iz) - 9.0 * (*u.get_grid(ilevel))(ix + 2, iy, iz) + (*u.get_grid(ilevel))(ix + 3, iy, iz)) * h;

		else if (dir == 1)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) +=
								(-(*u.get_grid(ilevel))(ix, iy - 3, iz) + 9.0 * (*u.get_grid(ilevel))(ix, iy - 2, iz) - 45.0 * (*u.get_grid(ilevel))(ix, iy - 1, iz) + 45.0 * (*u.get_grid(ilevel))(ix, iy + 1, iz) - 9.0 * (*u.get_grid(ilevel))(ix, iy + 2, iz) + (*u.get_grid(ilevel))(ix, iy + 3, iz)) * h;

		else if (dir == 2)
#pragma omp parallel for
			for (int ix = 0; ix < (int)(*u.get_grid(ilevel)).size(0); ++ix)
				for (int iy = 0; iy < (int)(*u.get_grid(ilevel)).size(1); ++iy)
					for (int iz = 0; iz < (int)(*u.get_grid(ilevel)).size(2); ++iz)
						(*pvar)(ix, iy, iz) +=
								(-(*u.get_grid(ilevel))(ix, iy, iz - 3) + 9.0 * (*u.get_grid(ilevel))(ix, iy, iz - 2) - 45.0 * (*u.get_grid(ilevel))(ix, iy, iz - 1) + 45.0 * (*u.get_grid(ilevel))(ix, iy, iz + 1) - 9.0 * (*u.get_grid(ilevel))(ix, iy, iz + 2) + (*u.get_grid(ilevel))(ix, iy, iz + 3)) * h;
	}

	music::ulog.Print("Done computing a 6th order finite difference gradient.");
}

/**************************************************************************************/
/**************************************************************************************/
#include "general.hh"

real_t fft_poisson_plugin::solve(grid_hierarchy &f, grid_hierarchy &u)
{
	music::ulog.Print("Entering k-space Poisson solver...");

	unsigned verbosity = cf_.get_value_safe<unsigned>("setup", "verbosity", 2);

	if (f.levelmin() != f.levelmax())
	{
		music::elog.Print("Attempt to run k-space Poisson solver on non unigrid mesh.");
		throw std::runtime_error("fft_poisson_plugin::solve : k-space method can only be used in unigrid mode (levelmin=levelmax)");
	}

	if (verbosity > 0)
	{
		music::ilog << "-------------------------------------------------------------------------------" << std::endl;
		music::ilog << " - Invoking unigrid FFT Poisson solver..." << std::endl;
	}

	int nx, ny, nz, nzp;
	nx = f.get_grid(f.levelmax())->size(0);
	ny = f.get_grid(f.levelmax())->size(1);
	nz = f.get_grid(f.levelmax())->size(2);
	nzp = 2 * (nz / 2 + 1);

	//... copy data ..................................................
	real_t *data = new real_t[(size_t)nx * (size_t)ny * (size_t)nzp];
	complex_t *cdata = reinterpret_cast<complex_t *>(data);

#pragma omp parallel for
	for (int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j)
			for (int k = 0; k < nz; ++k)
			{
				size_t idx = (size_t)(i * ny + j) * (size_t)nzp + (size_t)k;
				data[idx] = (*f.get_grid(f.levelmax()))(i, j, k);
			}

	//... perform FFT and Poisson solve................................
	music::ulog.Print("Performing forward transform.");

	fftw_plan_t
			plan = FFTW_API(plan_dft_r2c_3d)(nx, ny, nz, data, cdata, FFTW_ESTIMATE),
			iplan = FFTW_API(plan_dft_c2r_3d)(nx, ny, nz, cdata, data, FFTW_ESTIMATE);

	FFTW_API(execute)(plan);

	real_t kfac = 2.0 * M_PI;
	real_t fac = -1.0 / (real_t)((size_t)nx * (size_t)ny * (size_t)nz);

#pragma omp parallel for
	for (int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j)
			for (int k = 0; k < nz / 2 + 1; ++k)
			{
				int ii = i;
				if (ii > nx / 2)
					ii -= nx;
				int jj = j;
				if (jj > ny / 2)
					jj -= ny;
				real_t ki = (real_t)ii;
				real_t kj = (real_t)jj;
				real_t kk = (real_t)k;

				real_t kk2 = kfac * kfac * (ki * ki + kj * kj + kk * kk);

				size_t idx = (size_t)(i * ny + j) * (size_t)(nzp / 2) + (size_t)k;

				RE(cdata[idx]) *= -1.0 / kk2 * fac;
				IM(cdata[idx]) *= -1.0 / kk2 * fac;
			}

	RE(cdata[0]) = 0.0;
	IM(cdata[0]) = 0.0;

	music::ulog.Print("Performing backward transform.");

	FFTW_API(execute)(iplan);
	FFTW_API(destroy_plan)(plan);
	FFTW_API(destroy_plan)(iplan);

//... copy data ..........................................
#pragma omp parallel for
	for (int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j)
			for (int k = 0; k < nz; ++k)
			{
				size_t idx = (size_t)(i * ny + j) * (size_t)nzp + (size_t)k;
				(*u.get_grid(u.levelmax()))(i, j, k) = data[idx];
			}

	delete[] data;

	//... set boundary values ................................
	int nb = u.get_grid(u.levelmax())->m_nbnd;
	for (int iy = -nb; iy < ny + nb; ++iy)
		for (int iz = -nb; iz < nz + nb; ++iz)
		{
			int iiy((iy + ny) % ny), iiz((iz + nz) % nz);

			for (int i = -nb; i < 0; ++i)
			{
				(*u.get_grid(u.levelmax()))(i, iy, iz) = (*u.get_grid(u.levelmax()))(nx + i, iiy, iiz);
				(*u.get_grid(u.levelmax()))(nx - 1 - i, iy, iz) = (*u.get_grid(u.levelmax()))(-1 - i, iiy, iiz);
			}
		}

	for (int ix = -nb; ix < nx + nb; ++ix)
		for (int iz = -nb; iz < nz + nb; ++iz)
		{
			int iix((ix + nx) % nx), iiz((iz + nz) % nz);

			for (int i = -nb; i < 0; ++i)
			{
				(*u.get_grid(u.levelmax()))(ix, i, iz) = (*u.get_grid(u.levelmax()))(iix, ny + i, iiz);
				(*u.get_grid(u.levelmax()))(ix, ny - 1 - i, iz) = (*u.get_grid(u.levelmax()))(iix, -1 - i, iiz);
			}
		}

	for (int ix = -nb; ix < nx + nb; ++ix)
		for (int iy = -nb; iy < ny + nb; ++iy)
		{
			int iix((ix + nx) % nx), iiy((iy + ny) % ny);

			for (int i = -nb; i < 0; ++i)
			{
				(*u.get_grid(u.levelmax()))(ix, iy, i) = (*u.get_grid(u.levelmax()))(iix, iiy, nz + i);
				(*u.get_grid(u.levelmax()))(ix, iy, nz - 1 - i) = (*u.get_grid(u.levelmax()))(iix, iiy, -1 - i);
			}
		}

	music::ulog.Print("Done with k-space Poisson solver.");
	return 0.0;
}

real_t fft_poisson_plugin::gradient(int dir, grid_hierarchy &u, grid_hierarchy &Du)
{

	music::ulog.Print("Computing a gradient in k-space...\n");

	if (u.levelmin() != u.levelmax())
		throw std::runtime_error("fft_poisson_plugin::gradient : k-space method can only be used in unigrid mode (levelmin=levelmax)");

	Du = u;
	int nx, ny, nz, nzp;
	nx = u.get_grid(u.levelmax())->size(0);
	ny = u.get_grid(u.levelmax())->size(1);
	nz = u.get_grid(u.levelmax())->size(2);
	nzp = 2 * (nz / 2 + 1);

	//... copy data ..................................................
	real_t *data = new real_t[(size_t)nx * (size_t)ny * (size_t)nzp];
	complex_t *cdata = reinterpret_cast<complex_t *>(data);

#pragma omp parallel for
	for (int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j)
			for (int k = 0; k < nz; ++k)
			{
				size_t idx = (size_t)(i * ny + j) * (size_t)nzp + (size_t)k;
				data[idx] = (*u.get_grid(u.levelmax()))(i, j, k);
			}

			//... perform FFT and Poisson solve................................
	fftw_plan_t
			plan = FFTW_API(plan_dft_r2c_3d)(nx, ny, nz, data, cdata, FFTW_ESTIMATE),
			iplan = FFTW_API(plan_dft_c2r_3d)(nx, ny, nz, cdata, data, FFTW_ESTIMATE);

	FFTW_API(execute)(plan);

	real_t fac = -1.0 / (real_t)((size_t)nx * (size_t)ny * (size_t)nz);
	real_t kfac = 2.0 * M_PI;

	bool do_glass = cf_.get_value_safe<bool>("output", "glass", false);
	bool deconvolve_cic = do_glass | cf_.get_value_safe<bool>("output", "glass_cicdeconvolve", false);

	if (deconvolve_cic)
		music::ilog.Print("CIC deconvolution is enabled for kernel!");

#pragma omp parallel for
	for (int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j)
			for (int k = 0; k < nz / 2 + 1; ++k)
			{
				size_t idx = (size_t)(i * ny + j) * (size_t)(nzp / 2) + (size_t)k;
				int ii = i;
				if (ii > nx / 2) ii -= nx;
				int jj = j;
				if (jj > ny / 2) jj -= ny;

				const real_t ki{(real_t)ii};
				const real_t kj{(real_t)jj};
				const real_t kk{(real_t)k};
				const real_t kkdir[3] = {kfac * ki, kfac * kj, kfac * kk};
				const real_t kdir = kkdir[dir];

				real_t re = RE(cdata[idx]);
				real_t im = IM(cdata[idx]);

				RE(cdata[idx]) = fac * im * kdir;
				IM(cdata[idx]) = -fac * re * kdir;

				if (deconvolve_cic)
				{
					real_t dfx, dfy, dfz;
					dfx = M_PI * ki / (real_t)nx;
					dfx = (i != 0) ? std::sin(dfx) / dfx : 1.0;
					dfy = M_PI * kj / (real_t)ny;
					dfy = (j != 0) ? std::sin(dfy) / dfy : 1.0;
					dfz = M_PI * kk / (real_t)nz;
					dfz = (k != 0) ? std::sin(dfz) / dfz : 1.0;

					dfx = 1.0 / (dfx * dfy * dfz);
					dfx = dfx * dfx;
					RE(cdata[idx]) *= dfx;
					IM(cdata[idx]) *= dfx;
				}

				if( (dir == 0 && i==nx/2) || (dir == 1 && j==ny/2) || (dir == 2 && k==nz/2) )
				{
					RE(cdata[idx]) = 0.0;
					IM(cdata[idx]) = 0.0;
				}
			}

	RE(cdata[0]) = 0.0;
	IM(cdata[0]) = 0.0;

	FFTW_API(execute)(iplan);
	FFTW_API(destroy_plan)(plan);
	FFTW_API(destroy_plan)(iplan);

	//... copy data ..........................................
	real_t dmax = 0.0;
	for (int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j)
			for (int k = 0; k < nz; ++k)
			{
				size_t idx = ((size_t)i * ny + (size_t)j) * nzp + (size_t)k;
				(*Du.get_grid(u.levelmax()))(i, j, k) = data[idx];
				if (fabs(data[idx]) > dmax)
					dmax = fabs(data[idx]);
			}

	delete[] data;

	music::ulog.Print("Done with k-space gradient.\n");

	return 0.0;
}

/**************************************************************************************/
/**************************************************************************************/

template <int order>
real_t poisson_hybrid_kernel(int idir, int i, int j, int k, int n)
{
	return 1.0;
}

template <>
inline real_t poisson_hybrid_kernel<2>(int idir, int i, int j, int k, int n)
{
	if (i == 0 && j == 0 && k == 0)
		return 0.0;

	real_t
			ki(M_PI * (real_t)i / (real_t)n),
			kj(M_PI * (real_t)j / (real_t)n),
			kk(M_PI * (real_t)k / (real_t)n),
			kr(sqrt(ki * ki + kj * kj + kk * kk));

	real_t grad = 1.0, laplace = 1.0;

	if (idir == 0)
		grad = std::sin(ki);
	else if (idir == 1)
		grad = std::sin(kj);
	else
		grad = std::sin(kk);

	laplace = 2.0 * ((-std::cos(ki) + 1.0) + (-std::cos(kj) + 1.0) + (-std::cos(kk) + 1.0));

	real_t kgrad = 1.0;
	if (idir == 0)
		kgrad = ki;
	else if (idir == 1)
		kgrad = kj;
	else if (idir == 2)
		kgrad = kk;

	return kgrad / kr / kr - grad / laplace;
}

template <>
inline real_t poisson_hybrid_kernel<4>(int idir, int i, int j, int k, int n)
{

	if (i == 0 && j == 0 && k == 0)
		return 0.0;

	real_t
			ki(M_PI * (real_t)i / (real_t)n),
			kj(M_PI * (real_t)j / (real_t)n),
			kk(M_PI * (real_t)k / (real_t)n),
			kr(sqrt(ki * ki + kj * kj + kk * kk));

	real_t grad = 1.0, laplace = 1.0;

	if (idir == 0)
		grad = 0.166666666667 * (-std::sin(2. * ki) + 8. * std::sin(ki));
	else if (idir == 1)
		grad = 0.166666666667 * (-std::sin(2. * kj) + 8. * std::sin(kj));
	else if (idir == 2)
		grad = 0.166666666667 * (-std::sin(2. * kk) + 8. * std::sin(kk));

	laplace = 0.1666666667 * ((std::cos(2 * ki) - 16. * std::cos(ki) + 15.) + (std::cos(2 * kj) - 16. * std::cos(kj) + 15.) + (std::cos(2 * kk) - 16. * std::cos(kk) + 15.));

	real_t kgrad = 1.0;
	if (idir == 0)
		kgrad = ki;
	else if (idir == 1)
		kgrad = kj;
	else if (idir == 2)
		kgrad = kk;

	return kgrad / kr / kr - grad / laplace;
}

template <>
inline real_t poisson_hybrid_kernel<6>(int idir, int i, int j, int k, int n)
{
	real_t
			ki(M_PI * (real_t)i / (real_t)n),
			kj(M_PI * (real_t)j / (real_t)n),
			kk(M_PI * (real_t)k / (real_t)n),
			kr(sqrt(ki * ki + kj * kj + kk * kk));

	if (i == 0 && j == 0 && k == 0)
		return 0.0;

	real_t grad = 1.0, laplace = 1.0;

	if (idir == 0)
		grad = 0.0333333333333 * (std::sin(3. * ki) - 9. * std::sin(2. * ki) + 45. * std::sin(ki));
	else if (idir == 1)
		grad = 0.0333333333333 * (std::sin(3. * kj) - 9. * std::sin(2. * kj) + 45. * std::sin(kj));
	else if (idir == 2)
		grad = 0.0333333333333 * (std::sin(3. * kk) - 9. * std::sin(2. * kk) + 45. * std::sin(kk));

	laplace = 0.01111111111111 * ((-2. * std::cos(3.0 * ki) + 27. * std::cos(2. * ki) - 270. * std::cos(ki) + 245.) + (-2. * std::cos(3.0 * kj) + 27. * std::cos(2. * kj) - 270. * std::cos(kj) + 245.) + (-2. * std::cos(3.0 * kk) + 27. * std::cos(2. * kk) - 270. * std::cos(kk) + 245.));

	real_t kgrad = 1.0;
	if (idir == 0)
		kgrad = ki;
	else if (idir == 1)
		kgrad = kj;
	else if (idir == 2)
		kgrad = kk;

	//	if( i*i+j*j+k*k >= n*n )
	//		kgrad = 0.0;

	return kgrad / kr / kr - grad / laplace;
}

template <int order>
void do_poisson_hybrid(real_t *data, int idir, int nxp, int nyp, int nzp, bool periodic, bool deconvolve_cic)
{
	real_t fftnorm = 1.0 / ((real_t)nxp * (real_t)nyp * (real_t)nzp);

	complex_t *cdata = reinterpret_cast<complex_t *>(data);

	if (deconvolve_cic)
		music::ilog.Print("CIC deconvolution step is enabled.");

	fftw_plan_t iplan, plan;
	plan = FFTW_API(plan_dft_r2c_3d)(nxp, nyp, nzp, data, cdata, FFTW_ESTIMATE);
	iplan = FFTW_API(plan_dft_c2r_3d)(nxp, nyp, nzp, cdata, data, FFTW_ESTIMATE);
	FFTW_API(execute)(plan);

#pragma omp parallel for
	for (int i = 0; i < nxp; ++i)
		for (int j = 0; j < nyp; ++j)
			for (int k = 0; k < nzp / 2 + 1; ++k)
			{

				size_t ii = (size_t)(i * nyp + j) * (size_t)(nzp / 2 + 1) + (size_t)k;

				int ki(i), kj(j), kk(k);
				if (ki > nxp / 2)
					ki -= nxp;
				if (kj > nyp / 2)
					kj -= nyp;

				//... apply hybrid correction
				real_t dk = poisson_hybrid_kernel<order>(idir, ki, kj, k, nxp / 2);

				real_t re = RE(cdata[ii]), im = IM(cdata[ii]);

				RE(cdata[ii]) = -im * dk * fftnorm;
				IM(cdata[ii]) = re * dk * fftnorm;

				if (deconvolve_cic)
				{
					real_t dfx, dfy, dfz;
					dfx = M_PI * ki / (real_t)nxp;
					dfx = (i != 0) ? std::sin(dfx) / dfx : 1.0;
					dfy = M_PI * kj / (real_t)nyp;
					dfy = (j != 0) ? std::sin(dfy) / dfy : 1.0;
					dfz = M_PI * kk / (real_t)nzp;
					dfz = (k != 0) ? std::sin(dfz) / dfz : 1.0;

					dfx = 1.0 / (dfx * dfy * dfz);
					dfx = dfx * dfx;
					RE(cdata[ii]) *= dfx;
					IM(cdata[ii]) *= dfx;
				}

				if ((idir == 0 && i == nxp / 2) || (idir == 1 && j == nyp / 2) || (idir == 2 && k == nzp / 2))
				{
					RE(cdata[ii]) = 0.0;
					IM(cdata[ii]) = 0.0;
				}
			}

	RE(cdata[0]) = 0.0;
	IM(cdata[0]) = 0.0;

	FFTW_API(execute)(iplan);
	FFTW_API(destroy_plan)(plan);
	FFTW_API(destroy_plan)(iplan);
}

template <typename T>
void poisson_hybrid(T &f, int idir, int order, bool periodic, bool deconvolve_cic)

{
	int nx = f.size(0), ny = f.size(1), nz = f.size(2), nxp, nyp, nzp;
	real_t *data;
	int xo = 0, yo = 0, zo = 0;
	int nmax = std::max(nx, std::max(ny, nz));

	music::ulog.Print("Entering hybrid Poisson solver...");

	const int boundary = 32;

	if (!periodic)
	{
		nxp = nmax + 2 * boundary; 
		nyp = nmax + 2 * boundary; 
		nzp = nmax + 2 * boundary; 
		xo = boundary;						 
		yo = boundary;						 
		zo = boundary;						 
	}
	else
	{
		nxp = nmax;
		nyp = nmax;
		nzp = nmax;
	}

	data = new real_t[(size_t)nxp * (size_t)nyp * (size_t)(nzp + 2)];

	if (idir == 0)
		music::ilog << "   - Performing hybrid Poisson step... (" << nxp << ", " << nyp << ", " << nzp << ")" << std::endl;

#pragma omp parallel for
	for (int i = 0; i < nxp; ++i)
		for (int j = 0; j < nyp; ++j)
			for (int k = 0; k <= nzp; ++k)
			{
				size_t idx = ((size_t)i * (size_t)nxp + (size_t)j) * (size_t)(nzp + 2) + (size_t)k;
				data[idx] = 0.0;
			}

#pragma omp parallel for
	for (int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j)
			for (int k = 0; k < nz; ++k)
			{
				size_t idx = (size_t)((i + xo) * nyp + j + yo) * (size_t)(nzp + 2) + (size_t)(k + zo);
				data[idx] = f(i, j, k);
			}

	switch (order)
	{
	case 2:
		do_poisson_hybrid<2>(data, idir, nxp, nyp, nzp, periodic, deconvolve_cic);
		break;
	case 4:
		do_poisson_hybrid<4>(data, idir, nxp, nyp, nzp, periodic, deconvolve_cic);
		break;
	case 6:
		do_poisson_hybrid<6>(data, idir, nxp, nyp, nzp, periodic, deconvolve_cic);
		break;
	default:
		std::cerr << " - ERROR: invalid operator order specified in deconvolution.";
		music::elog.Print("Invalid operator order specified in deconvolution.");
		break;
	}

	music::ulog.Print("Copying hybrid correction factor...");

#pragma omp parallel for
	for (int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j)
			for (int k = 0; k < nz; ++k)
			{
				size_t idx = ((size_t)(i + xo) * nyp + (size_t)(j + yo)) * (size_t)(nzp + 2) + (size_t)(k + zo);
				f(i, j, k) = data[idx];
			}

	delete[] data;

	music::ulog.Print("Done with hybrid Poisson solve.");
}

/**************************************************************************************/
/**************************************************************************************/

template void poisson_hybrid<MeshvarBnd<real_t>>(MeshvarBnd<real_t> &f, int idir, int order, bool periodic, bool deconvolve_cic);
template void poisson_hybrid<MeshvarBnd<float>>(MeshvarBnd<float> &f, int idir, int order, bool periodic, bool deconvolve_cic);

namespace
{
	poisson_plugin_creator_concrete<multigrid_poisson_plugin> multigrid_poisson_creator("mg_poisson");
	poisson_plugin_creator_concrete<fft_poisson_plugin> fft_poisson_creator("fft_poisson");
}
