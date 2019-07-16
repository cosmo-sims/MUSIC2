/*
 
 densities.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#include <cstring>

#include "densities.hh"
#include "convolution_kernel.hh"

//TODO: this should be a larger number by default, just to maintain consistency with old default
#define DEF_RAN_CUBE_SIZE 32

double blend_sharpness = 0.5;

double Blend_Function(double k, double kmax)
{
#if 0
    const static double pihalf = 0.5*M_PI;
    if( fabs(k)>kmax ) return 0.0;
    return cos(k/kmax * pihalf );
#else
	float kabs = fabs(k);
	double const eps = blend_sharpness;
	float kp = (1.0f - 2.0f * eps) * kmax;

	if (kabs >= kmax)
		return 0.;
	if (kabs > kp)
		return 1.0f / (expf((kp - kmax) / (k - kp) + (kp - kmax) / (k - kmax)) + 1.0f);
	return 1.0f;
#endif
}

template <typename m1, typename m2>
void fft_coarsen(m1 &v, m2 &V)
{
	size_t nxf = v.size(0), nyf = v.size(1), nzf = v.size(2), nzfp = nzf + 2;
	size_t nxF = V.size(0), nyF = V.size(1), nzF = V.size(2), nzFp = nzF + 2;

	fftw_real *rcoarse = new fftw_real[nxF * nyF * nzFp];
	fftw_complex *ccoarse = reinterpret_cast<fftw_complex *>(rcoarse);

	fftw_real *rfine = new fftw_real[nxf * nyf * nzfp];
	fftw_complex *cfine = reinterpret_cast<fftw_complex *>(rfine);

#ifdef FFTW3
#ifdef SINGLE_PRECISION
	fftwf_plan
		pf = fftwf_plan_dft_r2c_3d(nxf, nyf, nzf, rfine, cfine, FFTW_ESTIMATE),
		ipc = fftwf_plan_dft_c2r_3d(nxF, nyF, nzF, ccoarse, rcoarse, FFTW_ESTIMATE);
#else
	fftw_plan
		pf = fftw_plan_dft_r2c_3d(nxf, nyf, nzf, rfine, cfine, FFTW_ESTIMATE),
		ipc = fftw_plan_dft_c2r_3d(nxF, nyF, nzF, ccoarse, rcoarse, FFTW_ESTIMATE);
#endif

#else
	rfftwnd_plan
		pf = rfftw3d_create_plan(nxf, nyf, nzf, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE),
		ipc = rfftw3d_create_plan(nxF, nyF, nzF, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
#endif

#pragma omp parallel for
	for (int i = 0; i < (int)nxf; i++)
		for (int j = 0; j < (int)nyf; j++)
			for (int k = 0; k < (int)nzf; k++)
			{
				size_t q = ((size_t)i * nyf + (size_t)j) * nzfp + (size_t)k;
				rfine[q] = v(i, j, k);
			}

#ifdef FFTW3
#ifdef SINGLE_PRECISION
	fftwf_execute(pf);
#else
	fftw_execute(pf);
#endif
#else
#ifndef SINGLETHREAD_FFTW
	rfftwnd_threads_one_real_to_complex(omp_get_max_threads(), pf, rfine, NULL);
#else
	rfftwnd_one_real_to_complex(pf, rfine, NULL);
#endif
#endif

	double fftnorm = 1.0 / ((double)nxF * (double)nyF * (double)nzF);

#pragma omp parallel for
	for (int i = 0; i < (int)nxF; i++)
		for (int j = 0; j < (int)nyF; j++)
			for (int k = 0; k < (int)nzF / 2 + 1; k++)
			{
				int ii(i), jj(j), kk(k);

				if (i > (int)nxF / 2)
					ii += (int)nxf / 2;
				if (j > (int)nyF / 2)
					jj += (int)nyf / 2;

				size_t qc, qf;

				double kx = (i <= (int)nxF / 2) ? (double)i : (double)(i - (int)nxF);
				double ky = (j <= (int)nyF / 2) ? (double)j : (double)(j - (int)nyF);
				double kz = (k <= (int)nzF / 2) ? (double)k : (double)(k - (int)nzF);

				qc = ((size_t)i * nyF + (size_t)j) * (nzF / 2 + 1) + (size_t)k;
				qf = ((size_t)ii * nyf + (size_t)jj) * (nzf / 2 + 1) + (size_t)kk;

				std::complex<double> val_fine(RE(cfine[qf]), IM(cfine[qf]));
				double phase = (kx / nxF + ky / nyF + kz / nzF) * 0.5 * M_PI;

				std::complex<double> val_phas(cos(phase), sin(phase));

				val_fine *= val_phas * fftnorm / 8.0; //sqrt(8.0);

				RE(ccoarse[qc]) = val_fine.real();
				IM(ccoarse[qc]) = val_fine.imag();
			}

	delete[] rfine;

#ifdef FFTW3
#ifdef SINGLE_PRECISION
	fftwf_execute(ipc);
#else
	fftw_execute(ipc);
#endif
#else
#ifndef SINGLETHREAD_FFTW
	rfftwnd_threads_one_complex_to_real(omp_get_max_threads(), ipc, ccoarse, NULL);
#else
	rfftwnd_one_complex_to_real(ipc, ccoarse, NULL);
#endif
#endif

#pragma omp parallel for
	for (int i = 0; i < (int)nxF; i++)
		for (int j = 0; j < (int)nyF; j++)
			for (int k = 0; k < (int)nzF; k++)
			{
				size_t q = ((size_t)i * nyF + (size_t)j) * nzFp + (size_t)k;
				V(i, j, k) = rcoarse[q];
			}

	delete[] rcoarse;

#ifdef FFTW3
#ifdef SINGLE_PRECISION
	fftwf_destroy_plan(pf);
	fftwf_destroy_plan(ipc);
#else
	fftw_destroy_plan(pf);
	fftw_destroy_plan(ipc);
#endif
#else
	rfftwnd_destroy_plan(pf);
	rfftwnd_destroy_plan(ipc);
#endif
}

//#define NO_COARSE_OVERLAP

template <typename m1, typename m2>
void fft_interpolate(m1 &V, m2 &v, bool from_basegrid = false)
{
	int oxf = v.offset(0), oyf = v.offset(1), ozf = v.offset(2);
	size_t nxf = v.size(0), nyf = v.size(1), nzf = v.size(2), nzfp = nzf + 2;
	size_t nxF = V.size(0), nyF = V.size(1), nzF = V.size(2);

	if (!from_basegrid)
	{
#ifdef NO_COARSE_OVERLAP
		oxf += nxF / 4;
		oyf += nyF / 4;
		ozf += nzF / 4;
#else
		oxf += nxF / 4 - nxf / 8;
		oyf += nyF / 4 - nyf / 8;
		ozf += nzF / 4 - nzf / 8;
	}
	else
	{
		oxf -= nxf / 8;
		oyf -= nyf / 8;
		ozf -= nzf / 8;
#endif
	}

	LOGUSER("FFT interpolate: offset=%d,%d,%d size=%d,%d,%d", oxf, oyf, ozf, nxf, nyf, nzf);

	// cut out piece of coarse grid that overlaps the fine:
	assert(nxf % 2 == 0 && nyf % 2 == 0 && nzf % 2 == 0);

	size_t nxc = nxf / 2, nyc = nyf / 2, nzc = nzf / 2, nzcp = nzf / 2 + 2;

	fftw_real *rcoarse = new fftw_real[nxc * nyc * nzcp];
	fftw_complex *ccoarse = reinterpret_cast<fftw_complex *>(rcoarse);

	fftw_real *rfine = new fftw_real[nxf * nyf * nzfp];
	fftw_complex *cfine = reinterpret_cast<fftw_complex *>(rfine);

	// copy coarse data to rcoarse[.]
	memset(rcoarse, 0, sizeof(fftw_real) * nxc * nyc * nzcp);

#ifdef NO_COARSE_OVERLAP
#pragma omp parallel for
	for (int i = 0; i < (int)nxc / 2; ++i)
		for (int j = 0; j < (int)nyc / 2; ++j)
			for (int k = 0; k < (int)nzc / 2; ++k)
			{
				int ii(i + nxc / 4);
				int jj(j + nyc / 4);
				int kk(k + nzc / 4);
				size_t q = ((size_t)ii * nyc + (size_t)jj) * nzcp + (size_t)kk;
				rcoarse[q] = V(oxf + i, oyf + j, ozf + k);
			}
#else
#pragma omp parallel for
	for (int i = 0; i < (int)nxc; ++i)
		for (int j = 0; j < (int)nyc; ++j)
			for (int k = 0; k < (int)nzc; ++k)
			{
				int ii(i);
				int jj(j);
				int kk(k);
				size_t q = ((size_t)ii * nyc + (size_t)jj) * nzcp + (size_t)kk;
				rcoarse[q] = V(oxf + i, oyf + j, ozf + k);
			}
#endif

#pragma omp parallel for
	for (int i = 0; i < (int)nxf; ++i)
		for (int j = 0; j < (int)nyf; ++j)
			for (int k = 0; k < (int)nzf; ++k)
			{
				size_t q = ((size_t)i * nyf + (size_t)j) * nzfp + (size_t)k;
				rfine[q] = v(i, j, k);
			}

#ifdef FFTW3
#ifdef SINGLE_PRECISION
	fftwf_plan
		pc = fftwf_plan_dft_r2c_3d(nxc, nyc, nzc, rcoarse, ccoarse, FFTW_ESTIMATE),
		pf = fftwf_plan_dft_r2c_3d(nxf, nyf, nzf, rfine, cfine, FFTW_ESTIMATE),
		ipf = fftwf_plan_dft_c2r_3d(nxf, nyf, nzf, cfine, rfine, FFTW_ESTIMATE);
	fftwf_execute(pc);
	fftwf_execute(pf);
#else
	fftw_plan
		pc = fftw_plan_dft_r2c_3d(nxc, nyc, nzc, rcoarse, ccoarse, FFTW_ESTIMATE),
		pf = fftw_plan_dft_r2c_3d(nxf, nyf, nzf, rfine, cfine, FFTW_ESTIMATE),
		ipf = fftw_plan_dft_c2r_3d(nxf, nyf, nzf, cfine, rfine, FFTW_ESTIMATE);
	fftw_execute(pc);
	fftw_execute(pf);
#endif
#else
	rfftwnd_plan
		pc = rfftw3d_create_plan(nxc, nyc, nzc, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE),
		pf = rfftw3d_create_plan(nxf, nyf, nzf, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE),
		ipf = rfftw3d_create_plan(nxf, nyf, nzf, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

#ifndef SINGLETHREAD_FFTW
	rfftwnd_threads_one_real_to_complex(omp_get_max_threads(), pc, rcoarse, NULL);
	rfftwnd_threads_one_real_to_complex(omp_get_max_threads(), pf, rfine, NULL);
#else
	rfftwnd_one_real_to_complex(pc, rcoarse, NULL);
	rfftwnd_one_real_to_complex(pf, rfine, NULL);
#endif
#endif

	/*************************************************/
	//.. perform actual interpolation
	double fftnorm = 1.0 / ((double)nxf * (double)nyf * (double)nzf);
	double sqrt8 = 8.0; //sqrt(8.0);
	double phasefac = -0.5;

	// this enables filtered splicing of coarse and fine modes
#if 1
	for (int i = 0; i < (int)nxc; i++)
		for (int j = 0; j < (int)nyc; j++)
			for (int k = 0; k < (int)nzc / 2 + 1; k++)
			{
				int ii(i), jj(j), kk(k);

				if (i > (int)nxc / 2)
					ii += (int)nxf / 2;
				if (j > (int)nyc / 2)
					jj += (int)nyf / 2;
				if (k > (int)nzc / 2)
					kk += (int)nzf / 2;

				size_t qc, qf;
				qc = ((size_t)i * (size_t)nyc + (size_t)j) * (nzc / 2 + 1) + (size_t)k;
				qf = ((size_t)ii * (size_t)nyf + (size_t)jj) * (nzf / 2 + 1) + (size_t)kk;

				double kx = (i <= (int)nxc / 2) ? (double)i : (double)(i - (int)nxc);
				double ky = (j <= (int)nyc / 2) ? (double)j : (double)(j - (int)nyc);
				double kz = (k <= (int)nzc / 2) ? (double)k : (double)(k - (int)nzc);

				double phase = phasefac * (kx / nxc + ky / nyc + kz / nzc) * M_PI;

				std::complex<double> val_phas(cos(phase), sin(phase));

				std::complex<double> val(RE(ccoarse[qc]), IM(ccoarse[qc]));
				val *= sqrt8 * val_phas;

				double blend_coarse = Blend_Function(sqrt(kx * kx + ky * ky + kz * kz), nxc / 2);
				double blend_fine = 1.0 - blend_coarse;

				RE(cfine[qf]) = blend_fine * RE(cfine[qf]) + blend_coarse * val.real();
				IM(cfine[qf]) = blend_fine * IM(cfine[qf]) + blend_coarse * val.imag();
			}

#else

	// 0 0
#pragma omp parallel for
	for (int i = 0; i < (int)nxc / 2 + 1; i++)
		for (int j = 0; j < (int)nyc / 2 + 1; j++)
			for (int k = 0; k < (int)nzc / 2 + 1; k++)
			{
				int ii(i), jj(j), kk(k);
				size_t qc, qf;
				qc = ((size_t)i * (size_t)nyc + (size_t)j) * (nzc / 2 + 1) + (size_t)k;
				qf = ((size_t)ii * (size_t)nyf + (size_t)jj) * (nzf / 2 + 1) + (size_t)kk;

				double kx = (i <= (int)nxc / 2) ? (double)i : (double)(i - (int)nxc);
				double ky = (j <= (int)nyc / 2) ? (double)j : (double)(j - (int)nyc);
				double kz = (k <= (int)nzc / 2) ? (double)k : (double)(k - (int)nzc);

				double phase = phasefac * (kx / nxc + ky / nyc + kz / nzc) * M_PI;
				std::complex<double> val_phas(cos(phase), sin(phase));

				std::complex<double> val(RE(ccoarse[qc]), IM(ccoarse[qc]));
				val *= sqrt8 * val_phas;

				RE(cfine[qf]) = val.real();
				IM(cfine[qf]) = val.imag();
			}

// 1 0
#pragma omp parallel for
	for (int i = nxc / 2; i < (int)nxc; i++)
		for (int j = 0; j < (int)nyc / 2 + 1; j++)
			for (int k = 0; k < (int)nzc / 2 + 1; k++)
			{
				int ii(i + nxf / 2), jj(j), kk(k);
				size_t qc, qf;
				qc = ((size_t)i * (size_t)nyc + (size_t)j) * (nzc / 2 + 1) + (size_t)k;
				qf = ((size_t)ii * (size_t)nyf + (size_t)jj) * (nzf / 2 + 1) + (size_t)kk;

				double kx = (i <= (int)nxc / 2) ? (double)i : (double)(i - (int)nxc);
				double ky = (j <= (int)nyc / 2) ? (double)j : (double)(j - (int)nyc);
				double kz = (k <= (int)nzc / 2) ? (double)k : (double)(k - (int)nzc);

				double phase = phasefac * (kx / nxc + ky / nyc + kz / nzc) * M_PI;
				std::complex<double> val_phas(cos(phase), sin(phase));

				std::complex<double> val(RE(ccoarse[qc]), IM(ccoarse[qc]));
				val *= sqrt8 * val_phas;

				RE(cfine[qf]) = val.real();
				IM(cfine[qf]) = val.imag();
			}

// 0 1
#pragma omp parallel for
	for (int i = 0; i < (int)nxc / 2 + 1; i++)
		for (int j = nyc / 2; j < (int)nyc; j++)
			for (int k = 0; k < (int)nzc / 2 + 1; k++)
			{
				int ii(i), jj(j + nyf / 2), kk(k);
				size_t qc, qf;
				qc = ((size_t)i * (size_t)nyc + (size_t)j) * (nzc / 2 + 1) + (size_t)k;
				qf = ((size_t)ii * (size_t)nyf + (size_t)jj) * (nzf / 2 + 1) + (size_t)kk;

				double kx = (i <= (int)nxc / 2) ? (double)i : (double)(i - (int)nxc);
				double ky = (j <= (int)nyc / 2) ? (double)j : (double)(j - (int)nyc);
				double kz = (k <= (int)nzc / 2) ? (double)k : (double)(k - (int)nzc);

				double phase = phasefac * (kx / nxc + ky / nyc + kz / nzc) * M_PI;
				std::complex<double> val_phas(cos(phase), sin(phase));

				std::complex<double> val(RE(ccoarse[qc]), IM(ccoarse[qc]));
				val *= sqrt8 * val_phas;

				RE(cfine[qf]) = val.real();
				IM(cfine[qf]) = val.imag();
			}

// 1 1
#pragma omp parallel for
	for (int i = nxc / 2; i < (int)nxc; i++)
		for (int j = nyc / 2; j < (int)nyc; j++)
			for (int k = 0; k < (int)nzc / 2 + 1; k++)
			{
				int ii(i + nxf / 2), jj(j + nyf / 2), kk(k);
				size_t qc, qf;
				qc = ((size_t)i * (size_t)nyc + (size_t)j) * (nzc / 2 + 1) + (size_t)k;
				qf = ((size_t)ii * (size_t)nyf + (size_t)jj) * (nzf / 2 + 1) + (size_t)kk;

				double kx = (i <= (int)nxc / 2) ? (double)i : (double)(i - (int)nxc);
				double ky = (j <= (int)nyc / 2) ? (double)j : (double)(j - (int)nyc);
				double kz = (k <= (int)nzc / 2) ? (double)k : (double)(k - (int)nzc);

				double phase = phasefac * (kx / nxc + ky / nyc + kz / nzc) * M_PI;
				std::complex<double> val_phas(cos(phase), sin(phase));

				std::complex<double> val(RE(ccoarse[qc]), IM(ccoarse[qc]));
				val *= sqrt8 * val_phas;

				RE(cfine[qf]) = val.real();
				IM(cfine[qf]) = val.imag();
			}
#endif

	delete[] rcoarse;

	/*************************************************/

#ifdef FFTW3
#ifdef SINGLE_PRECISION
	fftwf_execute(ipf);
	fftwf_destroy_plan(pf);
	fftwf_destroy_plan(pc);
	fftwf_destroy_plan(ipf);
#else
	fftw_execute(ipf);
	fftw_destroy_plan(pf);
	fftw_destroy_plan(pc);
	fftw_destroy_plan(ipf);
#endif
#else
#ifndef SINGLETHREAD_FFTW
	rfftwnd_threads_one_complex_to_real(omp_get_max_threads(), ipf, cfine, NULL);
#else
	rfftwnd_one_complex_to_real(ipf, cfine, NULL);
#endif
	fftwnd_destroy_plan(pf);
	fftwnd_destroy_plan(pc);
	fftwnd_destroy_plan(ipf);
#endif

// copy back and normalize
#pragma omp parallel for
	for (int i = 0; i < (int)nxf; ++i)
		for (int j = 0; j < (int)nyf; ++j)
			for (int k = 0; k < (int)nzf; ++k)
			{
				size_t q = ((size_t)i * nyf + (size_t)j) * nzfp + (size_t)k;
				v(i, j, k) = rfine[q] * fftnorm;
			}

	delete[] rfine;
}

/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/

void GenerateDensityUnigrid(config_file &cf, transfer_function *ptf, tf_type type,
							refinement_hierarchy &refh, rand_gen &rand, grid_hierarchy &delta, bool smooth, bool shift)
{
	unsigned levelmin, levelmax, levelminPoisson;

	levelminPoisson = cf.getValue<unsigned>("setup", "levelmin");
	levelmin = cf.getValueSafe<unsigned>("setup", "levelmin_TF", levelminPoisson);
	levelmax = cf.getValue<unsigned>("setup", "levelmax");

	bool kspace = cf.getValue<bool>("setup", "kspace_TF");

	bool fix  = cf.getValueSafe<bool>("setup","fix_mode_amplitude",false);
	bool flip = cf.getValueSafe<bool>("setup","flip_mode_amplitude",false);

	unsigned nbase = 1 << levelmin;

	std::cerr << " - Running unigrid version\n";
	LOGUSER("Running unigrid density convolution...");

	//... select the transfer function to be used
	convolution::kernel_creator *the_kernel_creator;

	if (kspace)
	{
		std::cout << " - Using k-space transfer function kernel.\n";
		LOGUSER("Using k-space transfer function kernel.");

#ifdef SINGLE_PRECISION
		the_kernel_creator = convolution::get_kernel_map()["tf_kernel_k_float"];
#else
		the_kernel_creator = convolution::get_kernel_map()["tf_kernel_k_double"];
#endif
	}
	else
	{
		std::cout << " - Using real-space transfer function kernel.\n";
		LOGUSER("Using real-space transfer function kernel.");

#ifdef SINGLE_PRECISION
		the_kernel_creator = convolution::get_kernel_map()["tf_kernel_real_float"];
#else
		the_kernel_creator = convolution::get_kernel_map()["tf_kernel_real_double"];
#endif
	}

	//... initialize convolution kernel
	convolution::kernel *the_tf_kernel = the_kernel_creator->create(cf, ptf, refh, type);

	//...
	std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmax << " ..." << std::endl;
	LOGUSER("Performing noise convolution on level %3d", levelmax);

	//... create convolution mesh
	DensityGrid<real_t> *top = new DensityGrid<real_t>(nbase, nbase, nbase);

	//... fill with random numbers
	rand.load(*top, levelmin);

	//... load convolution kernel
	the_tf_kernel->fetch_kernel(levelmin, false);

	//... perform convolution
	convolution::perform<real_t>(the_tf_kernel, reinterpret_cast<void *>(top->get_data_ptr()), shift, fix, flip);

	//... clean up kernel
	delete the_tf_kernel;

	//... create multi-grid hierarchy
	delta.create_base_hierarchy(levelmin);

	//... copy convolved field to multi-grid hierarchy
	top->copy(*delta.get_grid(levelmin));

	//... delete convolution grid
	delete top;
}

/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/

void GenerateDensityHierarchy(config_file &cf, transfer_function *ptf, tf_type type,
							  refinement_hierarchy &refh, rand_gen &rand,
							  grid_hierarchy &delta, bool smooth, bool shift)
{
	unsigned levelmin, levelmax, levelminPoisson;
	std::vector<long> rngseeds;
	std::vector<std::string> rngfnames;
	bool kspaceTF;

	double tstart, tend;

#ifndef SINGLETHREAD_FFTW
	tstart = omp_get_wtime();
#else
	tstart = (double)clock() / CLOCKS_PER_SEC;
#endif

	levelminPoisson = cf.getValue<unsigned>("setup", "levelmin");
	levelmin = cf.getValueSafe<unsigned>("setup", "levelmin_TF", levelminPoisson);
	levelmax = cf.getValue<unsigned>("setup", "levelmax");
	kspaceTF = cf.getValue<bool>("setup", "kspace_TF");

	bool fix  = cf.getValueSafe<bool>("setup","fix_mode_amplitude",false);
	bool flip = cf.getValueSafe<bool>("setup","flip_mode_amplitude",false);

	if( fix && levelmin != levelmax ){
		LOGWARN("You have chosen mode fixing for a zoom. This is not well tested,\n please proceed at your own risk...");
	}

	blend_sharpness = cf.getValueSafe<double>("setup", "kspace_filter", blend_sharpness);

	unsigned nbase = 1 << levelmin;

	convolution::kernel_creator *the_kernel_creator;

	if (kspaceTF)
	{
		std::cout << " - Using k-space transfer function kernel.\n";
		LOGUSER("Using k-space transfer function kernel.");

#ifdef SINGLE_PRECISION
		the_kernel_creator = convolution::get_kernel_map()["tf_kernel_k_float"];
#else
		the_kernel_creator = convolution::get_kernel_map()["tf_kernel_k_double"];
#endif
	}
	else
	{
		std::cout << " - Using real-space transfer function kernel.\n";
		LOGUSER("Using real-space transfer function kernel.");
#ifdef SINGLE_PRECISION
		the_kernel_creator = convolution::get_kernel_map()["tf_kernel_real_float"];
#else
		the_kernel_creator = convolution::get_kernel_map()["tf_kernel_real_double"];
#endif
	}

	convolution::kernel *the_tf_kernel = the_kernel_creator->create(cf, ptf, refh, type);

	/***** PERFORM CONVOLUTIONS *****/
	if (kspaceTF)
	{

		//... create and initialize density grids with white noise
		DensityGrid<real_t> *top(NULL);
		PaddedDensitySubGrid<real_t> *coarse(NULL), *fine(NULL);
		int nlevels = (int)levelmax - (int)levelmin + 1;

		// do coarse level
		top = new DensityGrid<real_t>(nbase, nbase, nbase);
		LOGINFO("Performing noise convolution on level %3d", levelmin);
		rand.load(*top, levelmin);
		convolution::perform<real_t>(the_tf_kernel->fetch_kernel(levelmin, false), reinterpret_cast<void *>(top->get_data_ptr()), shift, fix, flip);

		delta.create_base_hierarchy(levelmin);
		top->copy(*delta.get_grid(levelmin));

		for (int i = 1; i < nlevels; ++i)
		{
			LOGINFO("Performing noise convolution on level %3d...", levelmin + i);
			/////////////////////////////////////////////////////////////////////////
			//... add new refinement patch
			LOGUSER("Allocating refinement patch");
			LOGUSER("   offset=(%5d,%5d,%5d)", refh.offset(levelmin + i, 0),
					refh.offset(levelmin + i, 1), refh.offset(levelmin + i, 2));
			LOGUSER("   size  =(%5d,%5d,%5d)", refh.size(levelmin + i, 0),
					refh.size(levelmin + i, 1), refh.size(levelmin + i, 2));

			fine = new PaddedDensitySubGrid<real_t>(refh.offset(levelmin + i, 0),
													refh.offset(levelmin + i, 1),
													refh.offset(levelmin + i, 2),
													refh.size(levelmin + i, 0),
													refh.size(levelmin + i, 1),
													refh.size(levelmin + i, 2));
			/////////////////////////////////////////////////////////////////////////

			// load white noise for patch
			rand.load(*fine, levelmin + i);

			convolution::perform<real_t>(the_tf_kernel->fetch_kernel(levelmin + i, true),
										 reinterpret_cast<void *>(fine->get_data_ptr()), shift, fix, flip);

			if (i == 1)
				fft_interpolate(*top, *fine, true);
			else
				fft_interpolate(*coarse, *fine, false);

			delta.add_patch(refh.offset(levelmin + i, 0),
							refh.offset(levelmin + i, 1),
							refh.offset(levelmin + i, 2),
							refh.size(levelmin + i, 0),
							refh.size(levelmin + i, 1),
							refh.size(levelmin + i, 2));

			fine->copy_unpad(*delta.get_grid(levelmin + i));

			if (i == 1)
				delete top;
			else
				delete coarse;

			coarse = fine;
		}

		delete coarse;
	}
	else
	{

		//... create and initialize density grids with white noise
		PaddedDensitySubGrid<real_t> *coarse(NULL), *fine(NULL);
		DensityGrid<real_t> *top(NULL);

		if (levelmax == levelmin)
		{
			std::cout << " - Performing noise convolution on level "
					  << std::setw(2) << levelmax << " ..." << std::endl;
			LOGUSER("Performing noise convolution on level %3d...", levelmax);

			top = new DensityGrid<real_t>(nbase, nbase, nbase);
			//rand_gen.load( *top, levelmin );
			rand.load(*top, levelmin);

			convolution::perform<real_t>(the_tf_kernel->fetch_kernel(levelmax),
										 reinterpret_cast<void *>(top->get_data_ptr()),
										 shift, fix, flip);
			the_tf_kernel->deallocate();

			delta.create_base_hierarchy(levelmin);
			top->copy(*delta.get_grid(levelmin));
			delete top;
		}

		for (int i = 0; i < (int)levelmax - (int)levelmin; ++i)
		{
			//.......................................................................................................//
			//... GENERATE/FILL WITH RANDOM NUMBERS .................................................................//
			//.......................................................................................................//

			if (i == 0)
			{
				top = new DensityGrid<real_t>(nbase, nbase, nbase);
				rand.load(*top, levelmin);
			}

			fine = new PaddedDensitySubGrid<real_t>(refh.offset(levelmin + i + 1, 0),
													refh.offset(levelmin + i + 1, 1),
													refh.offset(levelmin + i + 1, 2),
													refh.size(levelmin + i + 1, 0),
													refh.size(levelmin + i + 1, 1),
													refh.size(levelmin + i + 1, 2));

			rand.load(*fine, levelmin + i + 1);

			//.......................................................................................................//
			//... PERFORM CONVOLUTIONS ..............................................................................//
			//.......................................................................................................//
			if (i == 0)
			{
				/**********************************************************************************************************\
			 *	multi-grid: top-level grid grids .....
			 \**********************************************************************************************************/
				std::cout << " - Performing noise convolution on level "
						  << std::setw(2) << levelmin + i << " ..." << std::endl;
				LOGUSER("Performing noise convolution on level %3d", levelmin + i);

				LOGUSER("Creating base hierarchy...");
				delta.create_base_hierarchy(levelmin);

				DensityGrid<real_t> top_save(*top);

				the_tf_kernel->fetch_kernel(levelmin);

				//... 1) compute standard convolution for levelmin
				LOGUSER("Computing density self-contribution");
				convolution::perform<real_t>(the_tf_kernel, reinterpret_cast<void *>(top->get_data_ptr()), shift, fix, flip);
				top->copy(*delta.get_grid(levelmin));

				//... 2) compute contribution to finer grids from non-refined region
				LOGUSER("Computing long-range component for finer grid.");
				*top = top_save;
				top_save.clear();
				top->zero_subgrid(refh.offset(levelmin + i + 1, 0), refh.offset(levelmin + i + 1, 1), refh.offset(levelmin + i + 1, 2),
								  refh.size(levelmin + i + 1, 0) / 2, refh.size(levelmin + i + 1, 1) / 2, refh.size(levelmin + i + 1, 2) / 2);

				convolution::perform<real_t>(the_tf_kernel, reinterpret_cast<void *>(top->get_data_ptr()), shift, fix, flip);
				the_tf_kernel->deallocate();

				meshvar_bnd delta_longrange(*delta.get_grid(levelmin));
				top->copy(delta_longrange);
				delete top;

				//... inject these contributions to the next level
				LOGUSER("Allocating refinement patch");
				LOGUSER("   offset=(%5d,%5d,%5d)", refh.offset(levelmin + 1, 0), refh.offset(levelmin + 1, 1), refh.offset(levelmin + 1, 2));
				LOGUSER("   size  =(%5d,%5d,%5d)", refh.size(levelmin + 1, 0), refh.size(levelmin + 1, 1), refh.size(levelmin + 1, 2));

				delta.add_patch(refh.offset(levelmin + 1, 0), refh.offset(levelmin + 1, 1), refh.offset(levelmin + 1, 2),
								refh.size(levelmin + 1, 0), refh.size(levelmin + 1, 1), refh.size(levelmin + 1, 2));

				LOGUSER("Injecting long range component");
				//mg_straight().prolong( delta_longrange, *delta.get_grid(levelmin+1) );
				//mg_cubic_mult().prolong( delta_longrange, *delta.get_grid(levelmin+1) );

				mg_cubic().prolong(delta_longrange, *delta.get_grid(levelmin + 1));
			}
			else
			{
				/**********************************************************************************************************\
			 *	multi-grid: intermediate sub-grids .....
			 \**********************************************************************************************************/

				std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmin + i << " ..." << std::endl;
				LOGUSER("Performing noise convolution on level %3d", levelmin + i);

				//... add new refinement patch
				LOGUSER("Allocating refinement patch");
				LOGUSER("   offset=(%5d,%5d,%5d)", refh.offset(levelmin + i + 1, 0), refh.offset(levelmin + i + 1, 1), refh.offset(levelmin + i + 1, 2));
				LOGUSER("   size  =(%5d,%5d,%5d)", refh.size(levelmin + i + 1, 0), refh.size(levelmin + i + 1, 1), refh.size(levelmin + i + 1, 2));

				delta.add_patch(refh.offset(levelmin + i + 1, 0), refh.offset(levelmin + i + 1, 1), refh.offset(levelmin + i + 1, 2),
								refh.size(levelmin + i + 1, 0), refh.size(levelmin + i + 1, 1), refh.size(levelmin + i + 1, 2));

				//... copy coarse grid long-range component to fine grid
				LOGUSER("Injecting long range component");
				//mg_straight().prolong( *delta.get_grid(levelmin+i), *delta.get_grid(levelmin+i+1) );
				mg_cubic().prolong(*delta.get_grid(levelmin + i), *delta.get_grid(levelmin + i + 1));

				PaddedDensitySubGrid<real_t> coarse_save(*coarse);
				the_tf_kernel->fetch_kernel(levelmin + i);

				//... 1) the inner region
				LOGUSER("Computing density self-contribution");
				coarse->subtract_boundary_oct_mean();
				convolution::perform<real_t>(the_tf_kernel, reinterpret_cast<void *>(coarse->get_data_ptr()), shift, fix, flip);
				coarse->copy_add_unpad(*delta.get_grid(levelmin + i));

				//... 2) the 'BC' for the next finer grid
				LOGUSER("Computing long-range component for finer grid.");
				*coarse = coarse_save;
				coarse->subtract_boundary_oct_mean();
				coarse->zero_subgrid(refh.offset(levelmin + i + 1, 0), refh.offset(levelmin + i + 1, 1), refh.offset(levelmin + i + 1, 2),
									 refh.size(levelmin + i + 1, 0) / 2, refh.size(levelmin + i + 1, 1) / 2, refh.size(levelmin + i + 1, 2) / 2);

				convolution::perform<real_t>(the_tf_kernel, reinterpret_cast<void *>(coarse->get_data_ptr()), shift, fix, flip);

				//... interpolate to finer grid(s)
				meshvar_bnd delta_longrange(*delta.get_grid(levelmin + i));
				coarse->copy_unpad(delta_longrange);

				LOGUSER("Injecting long range component");
				//mg_straight().prolong_add( delta_longrange, *delta.get_grid(levelmin+i+1) );

				mg_cubic().prolong_add(delta_longrange, *delta.get_grid(levelmin + i + 1));

				//... 3) the coarse-grid correction
				LOGUSER("Computing coarse grid correction");
				*coarse = coarse_save;
				coarse->subtract_oct_mean();
				convolution::perform<real_t>(the_tf_kernel, reinterpret_cast<void *>(coarse->get_data_ptr()), shift, fix, flip);
				coarse->subtract_mean();
				coarse->upload_bnd_add(*delta.get_grid(levelmin + i - 1));

				//... clean up
				the_tf_kernel->deallocate();
				delete coarse;
			}

			coarse = fine;
		}

		//... and convolution for finest grid (outside loop)
		if (levelmax > levelmin)
		{
			/**********************************************************************************************************\
			 *	multi-grid: finest sub-grid .....
			\**********************************************************************************************************/
			std::cout << " - Performing noise convolution on level " << std::setw(2) << levelmax << " ..." << std::endl;
			LOGUSER("Performing noise convolution on level %3d", levelmax);

			//... 1) grid self-contribution
			LOGUSER("Computing density self-contribution");
			PaddedDensitySubGrid<real_t> coarse_save(*coarse);

			//... create convolution kernel
			the_tf_kernel->fetch_kernel(levelmax);

			//... subtract oct mean on boundary but not in interior
			coarse->subtract_boundary_oct_mean();

			//... perform convolution
			convolution::perform<real_t>(the_tf_kernel, reinterpret_cast<void *>(coarse->get_data_ptr()), shift, fix, flip);

			//... copy to grid hierarchy
			coarse->copy_add_unpad(*delta.get_grid(levelmax));

			//... 2) boundary correction to top grid
			LOGUSER("Computing coarse grid correction");
			*coarse = coarse_save;

			//... subtract oct mean
			coarse->subtract_oct_mean();

			//... perform convolution
			convolution::perform<real_t>(the_tf_kernel, reinterpret_cast<void *>(coarse->get_data_ptr()), shift, fix, flip);

			the_tf_kernel->deallocate();

			coarse->subtract_mean();

			//... upload data to coarser grid
			coarse->upload_bnd_add(*delta.get_grid(levelmax - 1));

			delete coarse;
		}
	}

	delete the_tf_kernel;

#ifndef SINGLETHREAD_FFTW
	tend = omp_get_wtime();
	if (true) //verbosity > 1 )
		std::cout << " - Density calculation took " << tend - tstart << "s with " << omp_get_max_threads() << " threads." << std::endl;
#else
	tend = (double)clock() / CLOCKS_PER_SEC;
	if (true) //verbosity > 1 )
		std::cout << " - Density calculation took " << tend - tstart << "s." << std::endl;
#endif

	LOGUSER("Finished computing the density field in %fs", tend - tstart);
}

/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/

void normalize_density(grid_hierarchy &delta)
{
	//return;

	long double sum = 0.0;
	unsigned levelmin = delta.levelmin(), levelmax = delta.levelmax();

	{
		size_t nx, ny, nz;

		nx = delta.get_grid(levelmin)->size(0);
		ny = delta.get_grid(levelmin)->size(1);
		nz = delta.get_grid(levelmin)->size(2);

#pragma omp parallel for reduction(+ \
								   : sum)
		for (int ix = 0; ix < (int)nx; ++ix)
			for (size_t iy = 0; iy < ny; ++iy)
				for (size_t iz = 0; iz < nz; ++iz)
					sum += (*delta.get_grid(levelmin))(ix, iy, iz);

		sum /= (double)(nx * ny * nz);
	}

	std::cout << " - Top grid mean density is off by " << sum << ", correcting..." << std::endl;
	LOGUSER("Grid mean density is %g. Correcting...", sum);

	for (unsigned i = levelmin; i <= levelmax; ++i)
	{
		size_t nx, ny, nz;
		nx = delta.get_grid(i)->size(0);
		ny = delta.get_grid(i)->size(1);
		nz = delta.get_grid(i)->size(2);

#pragma omp parallel for
		for (int ix = 0; ix < (int)nx; ++ix)
			for (size_t iy = 0; iy < ny; ++iy)
				for (size_t iz = 0; iz < nz; ++iz)
					(*delta.get_grid(i))(ix, iy, iz) -= sum;
	}
}

void coarsen_density(const refinement_hierarchy &rh, GridHierarchy<real_t> &u, bool kspace)
{
	unsigned levelmin_TF = u.levelmin();

	/*for( int i=rh.levelmax(); i>0; --i )
        mg_straight().restrict( *(u.get_grid(i)), *(u.get_grid(i-1)) );*/

	if (kspace)
	{
		for (int i = levelmin_TF; i >= (int)rh.levelmin(); --i)
			fft_coarsen(*(u.get_grid(i)), *(u.get_grid(i - 1)));
	}
	else
	{
		for (int i = levelmin_TF; i >= (int)rh.levelmin(); --i)
			mg_straight().restrict(*(u.get_grid(i)), *(u.get_grid(i - 1)));
	}

	for (unsigned i = 1; i <= rh.levelmax(); ++i)
	{
		if (rh.offset(i, 0) != u.get_grid(i)->offset(0) || rh.offset(i, 1) != u.get_grid(i)->offset(1) || rh.offset(i, 2) != u.get_grid(i)->offset(2) || rh.size(i, 0) != u.get_grid(i)->size(0) || rh.size(i, 1) != u.get_grid(i)->size(1) || rh.size(i, 2) != u.get_grid(i)->size(2))
		{
			u.cut_patch(i, rh.offset_abs(i, 0), rh.offset_abs(i, 1), rh.offset_abs(i, 2),
						rh.size(i, 0), rh.size(i, 1), rh.size(i, 2));
		}
	}

	for (int i = rh.levelmax(); i > 0; --i)
		mg_straight().restrict(*(u.get_grid(i)), *(u.get_grid(i - 1)));
}
