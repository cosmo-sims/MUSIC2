// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2024 by Oliver Hahn
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
#pragma once

#include <vector>
#include <complex>

#include <gsl/gsl_linalg.h>

#include <general.hh>
#include <config_file.hh>
#include <transfer_function.hh>
#include <cosmology_calculator.hh>

//! matrix class serving as a gsl wrapper
class matrix
{
protected:
	gsl_matrix *m_;
	// double *data_;
	size_t M_, N_;

public:
	matrix(size_t M, size_t N)
			: M_(M), N_(N)
	{
		m_ = gsl_matrix_alloc(M_, N_);
	}

	matrix(size_t N)
			: M_(N), N_(N)
	{
		m_ = gsl_matrix_alloc(M_, N_);
	}

	matrix(const matrix &o)
	{
		M_ = o.M_;
		N_ = o.N_;
		m_ = gsl_matrix_alloc(M_, N_);
		gsl_matrix_memcpy(m_, o.m_);
	}

	~matrix()
	{
		gsl_matrix_free(m_);
	}

	double &operator()(size_t i, size_t j)
	{
		return *gsl_matrix_ptr(m_, i, j);
	}

	const double &operator()(size_t i, size_t j) const
	{
		return *gsl_matrix_const_ptr(m_, i, j);
	}

	matrix &operator=(const matrix &o)
	{
		gsl_matrix_free(m_);

		M_ = o.M_;
		N_ = o.N_;
		m_ = gsl_matrix_alloc(M_, N_);
		gsl_matrix_memcpy(m_, o.m_);
		return *this;
	}

	matrix &invert()
	{
		if (M_ != N_)
			throw std::runtime_error("Attempt to invert a non-square matrix!");

		int s;
		gsl_matrix *im = gsl_matrix_alloc(M_, N_);

		gsl_permutation *p = gsl_permutation_alloc(M_);
		gsl_linalg_LU_decomp(m_, p, &s);
		gsl_linalg_LU_invert(m_, p, im);

		gsl_matrix_memcpy(m_, im);

		gsl_permutation_free(p);
		gsl_matrix_free(im);
		return *this;
	}
};

//! class to impose constraints on the white noise field (van de Weygaert & Bertschinger 1996)
class constraint_set
{

public:
	enum constr_type
	{
		halo,
		peak
	};

protected:
	struct constraint
	{
		constr_type type;
		double x, y, z;
		double gx, gy, gz;
		double Rg, Rg2;
		double gRg, gRg2;
		double sigma;
	};

	config_file *pcf_;
	std::vector<constraint> cset_;
	transfer_function *ptf_;
	const cosmology::calculator *pccalc_;
	const cosmology::parameters *pcosmo_;
	double dplus0_;
	unsigned constr_level_;

	inline std::complex<double> eval_constr(size_t icon, double kx, double ky, double kz)
	{
		double re, im, kdotx, k2;

		kdotx = cset_[icon].gx * kx + cset_[icon].gy * ky + cset_[icon].gz * kz;
		k2 = kx * kx + ky * ky + kz * kz;

		re = im = exp(-k2 * cset_[icon].gRg2 / 2.0);
		re *= cos(kdotx);
		im *= sin(kdotx);

		return std::complex<double>(re, im);
	}
	//! apply constraints to the white noise
	void wnoise_constr_corr(double dx, size_t nx, size_t ny, size_t nz, std::vector<double> &g0, matrix &cinv, complex_t *cw);

	//! measure sigma for each constraint in the unconstrained noise
	void wnoise_constr_corr(double dx, complex_t *cw, size_t nx, size_t ny, size_t nz, std::vector<double> &g0);

	//! compute the covariance between the constraints
	void icov_constr(double dx, size_t nx, size_t ny, size_t nz, matrix &cij);

public:
	//! constructor
	constraint_set(config_file &cf, transfer_function *ptf);

	//! destructor
	~constraint_set()
	{
		delete pccalc_;
		delete pcosmo_;
	}

	template <typename rng>
	void apply(unsigned ilevel, int x0[], int lx[], rng *wnoise)
	{
		if (cset_.size() == 0 || constr_level_ != ilevel)
			return;

		unsigned nlvl = 1 << ilevel;
		double boxlength = pcf_->get_value<double>("setup", "boxlength");

		//... compute constraint coordinates for grid
		for (size_t i = 0; i < cset_.size(); ++i)
		{
			cset_[i].gx = cset_[i].x * (double)nlvl;
			cset_[i].gy = cset_[i].y * (double)nlvl;
			cset_[i].gz = cset_[i].z * (double)nlvl;
			cset_[i].gRg = cset_[i].Rg / boxlength * (double)nlvl;
			cset_[i].gRg2 = cset_[i].gRg * cset_[i].gRg;

			if (cset_[i].gRg > 0.5 * lx[0])
				music::wlog.Print("Constraint %d appears to be too large scale", i);
		}

		std::vector<double> g0;

		//		unsigned levelmax = pcf_->get_value<unsigned>("setup","levelmax");
		unsigned levelmin = pcf_->get_value<unsigned>("setup", "levelmin_TF");

		bool bperiodic = ilevel == levelmin;
		double dx = pcf_->get_value<double>("setup", "boxlength") / (1 << ilevel);

		music::ilog.Print("Computing constrained realization...");

		if (bperiodic)
		{
			//... we are operating on the periodic coarse grid
			size_t nx = lx[0], ny = lx[1], nz = lx[2], nzp = nz + 2;
			real_t *w = new real_t[nx * ny * nzp];

			complex_t *cw = reinterpret_cast<complex_t *>(w);
			fftw_plan_t p = FFTW_API(plan_dft_r2c_3d)(nx, ny, nz, w, cw, FFTW_ESTIMATE),
									ip = FFTW_API(plan_dft_c2r_3d)(nx, ny, nz, cw, w, FFTW_ESTIMATE);

			double fftnorm = 1.0 / sqrt(nx * ny * nz);

#pragma omp parallel for
			for (int i = 0; i < (int)nx; i++)
				for (int j = 0; j < (int)ny; j++)
					for (int k = 0; k < (int)nz; k++)
					{
						size_t q = ((size_t)i * ny + (size_t)j) * nzp + (size_t)k;
						w[q] = (*wnoise)((x0[0] + i) % nx, (x0[1] + j) % ny, (x0[2] + k) % nz) * fftnorm;
					}

			FFTW_API(execute)(p);
			wnoise_constr_corr(dx, cw, nx, ny, nz, g0);

			matrix c(2, 2);
			icov_constr(dx, nx, ny, nz, c);

			wnoise_constr_corr(dx, nx, ny, nz, g0, c, cw);

			FFTW_API(execute)(ip);

#pragma omp parallel for
			for (int i = 0; i < (int)nx; i++)
				for (int j = 0; j < (int)ny; j++)
					for (int k = 0; k < (int)nz; k++)
					{
						size_t q = ((size_t)i * ny + (size_t)j) * nzp + (size_t)k;
						(*wnoise)((x0[0] + i), (x0[1] + j), (x0[2] + k)) = w[q] * fftnorm;
					}

			music::ilog.Print("Applied constraints to level %d.", ilevel);

			delete[] w;

			FFTW_API(destroy_plan)(p);
			FFTW_API(destroy_plan)(ip);
		}
		else
		{

			//... we are operating on a refinement grid, not necessarily the finest

			size_t nx = lx[0], ny = lx[1], nz = lx[2], nzp = nz + 2;
			real_t *w = new real_t[nx * ny * nzp];

			complex_t *cw = reinterpret_cast<complex_t *>(w);
			fftw_plan_t p = FFTW_API(plan_dft_r2c_3d)(nx, ny, nz, w, cw, FFTW_ESTIMATE),
									ip = FFTW_API(plan_dft_c2r_3d)(nx, ny, nz, cw, w, FFTW_ESTIMATE);

			double fftnorm = 1.0 / sqrt(nx * ny * nz);

			int il = nx / 4, ir = 3 * nx / 4, jl = ny / 4, jr = 3 * ny / 4, kl = nz / 4, kr = 3 * nz / 4;

#pragma omp parallel for
			for (int i = 0; i < (int)nx; i++)
				for (int j = 0; j < (int)ny; j++)
					for (int k = 0; k < (int)nz; k++)
					{
						size_t q = ((size_t)i * ny + (size_t)j) * nzp + (size_t)k;

						if (i >= il && i < ir && j >= jl && j < jr && k >= kl && k < kr)
							w[q] = (*wnoise)((x0[0] + i), (x0[1] + j), (x0[2] + k)) * fftnorm;
						else
							w[q] = 0.0;
					}

			int nlvl05 = 1 << (ilevel - 1);
			int xs = nlvl05 - x0[0], ys = nlvl05 - x0[1], zs = nlvl05 - x0[2];

			for (size_t i = 0; i < cset_.size(); ++i)
			{
				cset_[i].gx -= xs;
				cset_[i].gy -= ys;
				cset_[i].gz -= zs;
			}

			FFTW_API(execute)(p);

			wnoise_constr_corr(dx, cw, nx, ny, nz, g0);

			matrix c(2, 2);
			icov_constr(dx, nx, ny, nz, c);

			wnoise_constr_corr(dx, nx, ny, nz, g0, c, cw);

			FFTW_API(execute)(ip);

#pragma omp parallel for
			for (int i = 0; i < (int)nx; i++)
				for (int j = 0; j < (int)ny; j++)
					for (int k = 0; k < (int)nz; k++)
					{
						size_t q = ((size_t)i * ny + (size_t)j) * nzp + (size_t)k;
						if (i >= il && i < ir && j >= jl && j < jr && k >= kl && k < kr)
							(*wnoise)((x0[0] + i), (x0[1] + j), (x0[2] + k)) = w[q] * fftnorm;
					}

			music::ilog.Print("Applied constraints to level %d.", ilevel);

			delete[] w;

			FFTW_API(destroy_plan)(p);
			FFTW_API(destroy_plan)(ip);
		}
	}
};
