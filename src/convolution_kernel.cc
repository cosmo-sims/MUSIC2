/*
 
 convolution_kernel.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010-19  Oliver Hahn
 
*/

#include "general.hh"
#include "densities.hh"
#include "convolution_kernel.hh"

#if defined(FFTW3) && defined(SINGLE_PRECISION)
//#define fftw_complex fftwf_complex
typedef fftw_complex fftwf_complex;
#endif

double T0 = 1.0;

namespace convolution
{

std::map<std::string, kernel_creator *> &
get_kernel_map()
{
	static std::map<std::string, kernel_creator *> kernel_map;
	return kernel_map;
}

template <typename real_t>
void perform(kernel *pk, void *pd, bool shift, bool fix, bool flip)
{
	//return;

	parameters cparam_ = pk->cparam_;
	double fftnormp = 1.0/sqrt((double)cparam_.nx * (double)cparam_.ny * (double)cparam_.nz);
	double fftnorm = pow(2.0 * M_PI, 1.5) / sqrt(cparam_.lx * cparam_.ly * cparam_.lz) * fftnormp;

	fftw_complex *cdata, *ckernel;
	fftw_real *data;

	data = reinterpret_cast<fftw_real *>(pd);
	cdata = reinterpret_cast<fftw_complex *>(data);
	ckernel = reinterpret_cast<fftw_complex *>(pk->get_ptr());

	std::cout << "   - Performing density convolution... ("
			  << cparam_.nx << ", " << cparam_.ny << ", " << cparam_.nz << ")\n";

	LOGUSER("Performing kernel convolution on (%5d,%5d,%5d) grid", cparam_.nx, cparam_.ny, cparam_.nz);
	LOGUSER("Performing forward FFT...");
#ifdef FFTW3
#ifdef SINGLE_PRECISION
	fftwf_plan plan, iplan;
	plan = fftwf_plan_dft_r2c_3d(cparam_.nx, cparam_.ny, cparam_.nz, data, cdata, FFTW_ESTIMATE);
	iplan = fftwf_plan_dft_c2r_3d(cparam_.nx, cparam_.ny, cparam_.nz, cdata, data, FFTW_ESTIMATE);

	fftwf_execute(plan);
#else
	fftw_plan plan, iplan;
	plan = fftw_plan_dft_r2c_3d(cparam_.nx, cparam_.ny, cparam_.nz, data, cdata, FFTW_ESTIMATE);
	iplan = fftw_plan_dft_c2r_3d(cparam_.nx, cparam_.ny, cparam_.nz, cdata, data, FFTW_ESTIMATE);

	fftw_execute(plan);
#endif
#else
	rfftwnd_plan iplan, plan;

	plan = rfftw3d_create_plan(cparam_.nx, cparam_.ny, cparam_.nz,
							   FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);

	iplan = rfftw3d_create_plan(cparam_.nx, cparam_.ny, cparam_.nz,
								FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

#ifndef SINGLETHREAD_FFTW
	rfftwnd_threads_one_real_to_complex(omp_get_max_threads(), plan, data, NULL);
#else
	rfftwnd_one_real_to_complex(plan, data, NULL);
#endif

#endif
	//..... need a phase shift for baryons for SPH
	double dstag = 0.0;

	if (shift)
	{
		double boxlength = pk->pcf_->getValue<double>("setup", "boxlength");
		double stagfact = pk->pcf_->getValueSafe<double>("setup", "baryon_staggering", 0.5);
		int lmax = pk->pcf_->getValue<int>("setup", "levelmax");
		double dxmax = boxlength / (1 << lmax);
		double dxcur = cparam_.lx / cparam_.nx;
		//std::cerr << "Performing staggering shift for SPH\n";
		LOGUSER("Performing staggering shift for SPH");
		dstag = stagfact * 2.0 * M_PI / cparam_.nx * dxmax / dxcur;
	}

	//.............................................

	std::complex<double> dcmode(RE(cdata[0]), IM(cdata[0]));

	if (!pk->is_ksampled())
	{

#pragma omp parallel for
		for (int i = 0; i < cparam_.nx; ++i)
			for (int j = 0; j < cparam_.ny; ++j)
				for (int k = 0; k < cparam_.nz / 2 + 1; ++k)
				{
					size_t ii = (size_t)(i * cparam_.ny + j) * (size_t)(cparam_.nz / 2 + 1) + (size_t)k;

					double kx, ky, kz;

					kx = (double)i;
					ky = (double)j;
					kz = (double)k;

					if (kx > cparam_.nx / 2)
						kx -= cparam_.nx;
					if (ky > cparam_.ny / 2)
						ky -= cparam_.ny;

					double arg = (kx + ky + kz) * dstag;
					std::complex<double> carg(cos(arg), sin(arg));

					std::complex<double>
						ccdata(RE(cdata[ii]), IM(cdata[ii])),
						cckernel(RE(ckernel[ii]), IM(ckernel[ii]));

					if( fix ){
						ccdata = ccdata / std::abs(ccdata);
					}
					if( flip ){
						ccdata = -ccdata;
					}

					ccdata = ccdata * cckernel * fftnorm * carg;

					RE(cdata[ii]) = ccdata.real();
					IM(cdata[ii]) = ccdata.imag();
				}
	}
	else
	{

#pragma omp parallel
		{

			const size_t veclen = cparam_.nz / 2 + 1;

			double *kvec = new double[veclen];
			double *Tkvec = new double[veclen];
			double *argvec = new double[veclen];

#pragma omp for
			for (int i = 0; i < cparam_.nx; ++i)
				for (int j = 0; j < cparam_.ny; ++j)
				{

					for (int k = 0; k < cparam_.nz / 2 + 1; ++k)
					{
						double kx, ky, kz;

						kx = (double)i;
						ky = (double)j;
						kz = (double)k;

						if (kx > cparam_.nx / 2)
							kx -= cparam_.nx;
						if (ky > cparam_.ny / 2)
							ky -= cparam_.ny;

						kvec[k] = sqrt(kx * kx + ky * ky + kz * kz);
						argvec[k] = (kx + ky + kz) * dstag;
					}

					pk->at_k(veclen, kvec, Tkvec);

					for (int k = 0; k < cparam_.nz / 2 + 1; ++k)
					{
						size_t ii = (size_t)(i * cparam_.ny + j) * (size_t)(cparam_.nz / 2 + 1) + (size_t)k;
						std::complex<double> carg(cos(argvec[k]), sin(argvec[k]));

						std::complex<double> ccdata(RE(cdata[ii]), IM(cdata[ii]));

						if( fix ){
							ccdata = ccdata / std::abs(ccdata) / fftnormp;
						}
						if( flip ){
							ccdata = -ccdata;
						}

						ccdata = ccdata * Tkvec[k] * fftnorm * carg;

						RE(cdata[ii]) = ccdata.real();
						IM(cdata[ii]) = ccdata.imag();
					}
				}

			delete[] kvec;
			delete[] Tkvec;
			delete[] argvec;
		}

		// we now set the correct DC mode below...
		RE(cdata[0]) = 0.0;
		IM(cdata[0]) = 0.0;
	}

	LOGUSER("Performing backward FFT...");

#ifdef FFTW3
#ifdef SINGLE_PRECISION
	fftwf_execute(iplan);
	fftwf_destroy_plan(plan);
	fftwf_destroy_plan(iplan);
#else
	fftw_execute(iplan);
	fftw_destroy_plan(plan);
	fftw_destroy_plan(iplan);

#endif
#else
#ifndef SINGLETHREAD_FFTW
	rfftwnd_threads_one_complex_to_real(omp_get_max_threads(), iplan, cdata, NULL);
#else
	rfftwnd_one_complex_to_real(iplan, cdata, NULL);
#endif

	rfftwnd_destroy_plan(plan);
	rfftwnd_destroy_plan(iplan);
#endif

	// set the DC mode here to avoid a possible truncation error in single precision
	if (pk->is_ksampled())
	{
		size_t nelem = (size_t)cparam_.nx * (size_t)cparam_.ny * (size_t)cparam_.nz;
		real_t mean = dcmode.real() * fftnorm / (real_t)nelem;

#pragma omp parallel for
		for (size_t i = 0; i < nelem; ++i)
			data[i] += mean;
	}
}

template void perform<double>(kernel *pk, void *pd, bool shift, bool fix, bool flip);
template void perform<float>(kernel *pk, void *pd, bool shift, bool fix, bool flip);

/*****************************************************************************************/
/***    SPECIFIC KERNEL IMPLEMENTATIONS      *********************************************/
/*****************************************************************************************/

template <typename real_t>
class kernel_k : public kernel
{
protected:
	/**/
	double boxlength_, patchlength_, nspec_, pnorm_, volfac_, kfac_, kmax_;
	TransferFunction_k *tfk_;

public:
	kernel_k(config_file &cf, transfer_function *ptf, refinement_hierarchy &refh, tf_type type)
		: kernel(cf, ptf, refh, type)
	{
		boxlength_ = pcf_->getValue<double>("setup", "boxlength");
		nspec_ = pcf_->getValue<double>("cosmology", "nspec");
		pnorm_ = pcf_->getValue<double>("cosmology", "pnorm");
		volfac_ = 1.0; //pow(boxlength,3)/pow(2.0*M_PI,3);
		kfac_ = 2.0 * M_PI / boxlength_;
		kmax_ = kfac_ / 2;
		tfk_ = new TransferFunction_k(type_, ptf_, nspec_, pnorm_);

		cparam_.nx = 1;
		cparam_.ny = 1;
		cparam_.nz = 1;
		cparam_.lx = boxlength_;
		cparam_.ly = boxlength_;
		cparam_.lz = boxlength_;
		cparam_.pcf = pcf_;
		patchlength_ = boxlength_;
	}

	kernel *fetch_kernel(int ilevel, bool isolated = false)
	{
		if (!isolated)
		{
			cparam_.nx = prefh_->size(ilevel, 0);
			cparam_.ny = prefh_->size(ilevel, 1);
			cparam_.nz = prefh_->size(ilevel, 2);

			cparam_.lx = (double)cparam_.nx / (double)(1 << ilevel) * boxlength_;
			cparam_.ly = (double)cparam_.ny / (double)(1 << ilevel) * boxlength_;
			cparam_.lz = (double)cparam_.nz / (double)(1 << ilevel) * boxlength_;

			patchlength_ = cparam_.lx;
			kfac_ = 2.0 * M_PI / patchlength_;
			kmax_ = kfac_ * cparam_.nx / 2;
		}
		else
		{
			cparam_.nx = 2 * prefh_->size(ilevel, 0);
			cparam_.ny = 2 * prefh_->size(ilevel, 1);
			cparam_.nz = 2 * prefh_->size(ilevel, 2);

			cparam_.lx = (double)cparam_.nx / (double)(1 << ilevel) * boxlength_;
			cparam_.ly = (double)cparam_.ny / (double)(1 << ilevel) * boxlength_;
			cparam_.lz = (double)cparam_.nz / (double)(1 << ilevel) * boxlength_;

			patchlength_ = cparam_.lx;
			kfac_ = 2.0 * M_PI / patchlength_;
			kmax_ = kfac_ * cparam_.nx / 2;
		}

		return this;
	}

	void *get_ptr() { return NULL; }

	bool is_ksampled() { return true; }

	void at_k(size_t len, const double *in_k, double *out_Tk)
	{
		for (size_t i = 0; i < len; ++i)
		{
			double kk = kfac_ * in_k[i];
			out_Tk[i] = volfac_ * tfk_->compute(kk);
		}
	}

	~kernel_k() { delete tfk_; }

	void deallocate() {}
};

////////////////////////////////////////////////////////////////////////////

template <typename real_t>
class kernel_real_cached : public kernel
{
protected:
	std::vector<real_t> kdata_;
	void precompute_kernel(transfer_function *ptf, tf_type type, const refinement_hierarchy &refh);

public:
	kernel_real_cached(config_file &cf, transfer_function *ptf, refinement_hierarchy &refh, tf_type type)
		: kernel(cf, ptf, refh, type)
	{
		precompute_kernel(ptf, type, refh);
	}

	kernel *fetch_kernel(int ilevel, bool isolated = false);

	void *get_ptr()
	{
		return reinterpret_cast<void *>(&kdata_[0]);
	}

	bool is_ksampled()
	{
		return false;
	}

	void at_k(size_t, const double *, double *) {}

	~kernel_real_cached()
	{
		deallocate();
	}

	void deallocate()
	{
		std::vector<real_t>().swap(kdata_);
	}
};

template <typename real_t>
kernel *kernel_real_cached<real_t>::fetch_kernel(int ilevel, bool isolated)
{
	char cachefname[128];
	sprintf(cachefname, "temp_kernel_level%03d.tmp", ilevel);
	FILE *fp = fopen(cachefname, "r");

	std::cout << " - Fetching kernel for level " << ilevel << std::endl;

	LOGUSER("Loading kernel for level %3d from file \'%s\'...", ilevel, cachefname);

	if (fp == NULL)
	{
		LOGERR("Could not open kernel file \'%s\'.", cachefname);
		throw std::runtime_error("Internal error: cached convolution kernel does not exist on disk!");
	}

	unsigned nx, ny, nz;
	size_t nread = 0;

	nread = fread(reinterpret_cast<void *>(&nx), sizeof(unsigned), 1, fp);
	nread = fread(reinterpret_cast<void *>(&ny), sizeof(unsigned), 1, fp);
	nread = fread(reinterpret_cast<void *>(&nz), sizeof(unsigned), 1, fp);

	kdata_.assign((size_t)nx * (size_t)ny * (size_t)nz, 0.0);

	for (size_t ix = 0; ix < nx; ++ix)
	{
		const size_t sz = ny * nz;
		nread = fread(reinterpret_cast<void *>(&kdata_[(size_t)ix * sz]), sizeof(fftw_real), sz, fp);
		assert(nread == sz);
	}

	fclose(fp);

	//... set parameters

	double boxlength = pcf_->getValue<double>("setup", "boxlength");
	double dx = boxlength / (1 << ilevel);

	cparam_.nx = nx;
	cparam_.ny = ny;
	cparam_.nz = nz - 2;
	cparam_.lx = dx * cparam_.nx;
	cparam_.ly = dx * cparam_.ny;
	cparam_.lz = dx * cparam_.nz;
	cparam_.pcf = pcf_;

	fftw_real *rkernel = reinterpret_cast<fftw_real *>(&kdata_[0]);

#ifdef FFTW3
#ifdef SINGLE_PRECISION
	fftwf_complex *kkernel = reinterpret_cast<fftwf_complex *>(&rkernel[0]);
	fftwf_plan plan = fftwf_plan_dft_r2c_3d(cparam_.nx, cparam_.ny, cparam_.nz, rkernel, kkernel, FFTW_ESTIMATE);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
#else
	fftw_complex *kkernel = reinterpret_cast<fftw_complex *>(&rkernel[0]);
	fftw_plan plan = fftw_plan_dft_r2c_3d(cparam_.nx, cparam_.ny, cparam_.nz, rkernel, kkernel, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
#endif
#else
	rfftwnd_plan plan = rfftw3d_create_plan(cparam_.nx, cparam_.ny, cparam_.nz,
											FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);

#ifndef SINGLETHREAD_FFTW
	rfftwnd_threads_one_real_to_complex(omp_get_max_threads(), plan, rkernel, NULL);
#else
	rfftwnd_one_real_to_complex(plan, rkernel, NULL);
#endif

	rfftwnd_destroy_plan(plan);
#endif
	return this;
}

template <typename real_t>
inline real_t sqr(real_t x)
{
	return x * x;
}

template <typename real_t>
inline real_t eval_split_recurse(const TransferFunction_real *tfr, real_t *xmid, real_t dx, real_t prevval, int nsplit)
{
	const real_t abs_err = 1e-8, rel_err = 1e-6;
	const int nmaxsplits = 12;

	real_t dxnew = dx / 2, dxnew2 = dx / 4;
	real_t dV = dxnew * dxnew * dxnew;

	real_t xl[3] = {xmid[0] - dxnew2, xmid[1] - dxnew2, xmid[2] - dxnew2};
	real_t xr[3] = {xmid[0] + dxnew2, xmid[1] + dxnew2, xmid[2] + dxnew2};

	real_t xc[8][3] =
		{
			{xl[0], xl[1], xl[2]},
			{xl[0], xl[1], xr[2]},
			{xl[0], xr[1], xl[2]},
			{xl[0], xr[1], xr[2]},
			{xr[0], xl[1], xl[2]},
			{xr[0], xl[1], xr[2]},
			{xr[0], xr[1], xl[2]},
			{xr[0], xr[1], xr[2]},
		};

	real_t rr2, res[8], ressum = 0.;

	for (int i = 0; i < 8; ++i)
	{
		rr2 = sqr(xc[i][0]) + sqr(xc[i][1]) + sqr(xc[i][2]);
		res[i] = tfr->compute_real(rr2) * dV;
		if (res[i] != res[i])
		{
			LOGERR("NaN encountered at r=%f, dx=%f, dV=%f : TF = %f", sqrt(rr2), dx, dV, res[i]);
			abort();
		}
		ressum += res[i];
	}

	real_t ae = fabs((prevval - ressum));
	real_t re = fabs(ae / ressum);

	if (ae < abs_err || re < rel_err)
		return ressum;

	if (nsplit > nmaxsplits)
	{
		//LOGWARN("reached maximum number of supdivisions in eval_split_recurse. Ending recursion... : abs. err.=%f, rel. err.=%f",ae, re);
		return ressum;
	}

	//otherwise keep splitting
	ressum = 0;
	for (int i = 0; i < 8; ++i)
		ressum += eval_split_recurse(tfr, xc[i], dxnew, res[i], nsplit + 1);

	return ressum;
}

template <typename real_t>
inline real_t eval_split_recurse(const TransferFunction_real *tfr, real_t *xmid, real_t dx, int nsplit = 0)
{
	//sqr(xmid[0])+sqr(xmid[1])+sqr(xmid[2])

	real_t rr2 = sqr(xmid[0]) + sqr(xmid[1]) + sqr(xmid[2]);
	real_t prevval = tfr->compute_real(rr2) * dx * dx * dx;
	return eval_split_recurse(tfr, xmid, dx, prevval, nsplit);
}

#define OLD_KERNEL_SAMPLING

template <typename real_t>
void kernel_real_cached<real_t>::precompute_kernel(transfer_function *ptf, tf_type type, const refinement_hierarchy &refh)
{
	//... compute kernel for finest level
	int nx, ny, nz;
	real_t dx, lx, ly, lz;

	real_t
		boxlength = pcf_->getValue<double>("setup", "boxlength"),
		boxlength2 = 0.5 * boxlength;

	int
		levelmax = refh.levelmax(),
		levelmin = refh.levelmin();

	LOGUSER("Precomputing transfer function kernels...");

	nx = refh.size(refh.levelmax(), 0);
	ny = refh.size(refh.levelmax(), 1);
	nz = refh.size(refh.levelmax(), 2);

	if (levelmax != levelmin)
	{
		nx *= 2;
		ny *= 2;
		nz *= 2;
	}

	dx = boxlength / (1 << refh.levelmax());
	lx = dx * nx;
	ly = dx * ny;
	lz = dx * nz;

	real_t
		kny = M_PI / dx,
		fac = lx * ly * lz / pow(2.0 * M_PI, 3) / ((double)nx * (double)ny * (double)nz),
		nspec = pcf_->getValue<double>("cosmology", "nspec"),
		pnorm = pcf_->getValue<double>("cosmology", "pnorm");

	bool
		bperiodic = pcf_->getValueSafe<bool>("setup", "periodic_TF", true),
		deconv = pcf_->getValueSafe<bool>("setup", "deconvolve", true);
	//		bool deconv_baryons = true;//pcf_->getValueSafe<bool>("setup","deconvolve_baryons",false) || do_SPH;
	bool bsmooth_baryons = false; //type==baryon && !deconv_baryons;
	//bool bbaryons = type==baryon | type==vbaryon;
	bool kspacepoisson = ((pcf_->getValueSafe<bool>("poisson", "fft_fine", true) |
						   pcf_->getValueSafe<bool>("poisson", "kspace", false))); // & !(type==baryon&!do_SPH);//&!baryons ;

	std::cout << "   - Computing transfer function kernel...\n";

	TransferFunction_real *tfr = new TransferFunction_real(boxlength, 1 << levelmax, type, ptf, nspec, pnorm,
														   0.25 * dx, 2.0 * boxlength, kny, (int)pow(2, levelmax + 2));

	fftw_real *rkernel = new fftw_real[(size_t)nx * (size_t)ny * ((size_t)nz + 2)], *rkernel_coarse;

#pragma omp parallel for
	for (int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j)
			for (int k = 0; k < nz; ++k)
			{
				size_t q = ((size_t)(i)*ny + (size_t)(j)) * (size_t)(nz + 2) + (size_t)(k);
				rkernel[q] = 0.0;
			}

	LOGUSER("Computing fine kernel (level %d)...", levelmax);

#ifdef OLD_KERNEL_SAMPLING
	int ref_fac = (deconv && kspacepoisson) ? 2 : 0;
	const int ql = -ref_fac / 2 + 1, qr = ql + ref_fac;
	const double rf8 = pow(ref_fac, 3);
	const double dx05 = 0.5 * dx, dx025 = 0.25 * dx;
#endif

	if (bperiodic)
	{
#pragma omp parallel for
		for (int i = 0; i <= nx / 2; ++i)
			for (int j = 0; j <= ny / 2; ++j)
				for (int k = 0; k <= nz / 2; ++k)
				{
					int iix(i), iiy(j), iiz(k);
					real_t rr[3];

					if (iix > (int)nx / 2)
						iix -= nx;
					if (iiy > (int)ny / 2)
						iiy -= ny;
					if (iiz > (int)nz / 2)
						iiz -= nz;

					//... speed up 8x by copying data to other octants
					size_t idx[8];

					idx[0] = ((size_t)(i)*ny + (size_t)(j)) * 2 * (nz / 2 + 1) + (size_t)(k);
					idx[1] = ((size_t)(nx - i) * ny + (size_t)(j)) * 2 * (nz / 2 + 1) + (size_t)(k);
					idx[2] = ((size_t)(i)*ny + (size_t)(ny - j)) * 2 * (nz / 2 + 1) + (size_t)(k);
					idx[3] = ((size_t)(nx - i) * ny + (size_t)(ny - j)) * 2 * (nz / 2 + 1) + (size_t)(k);
					idx[4] = ((size_t)(i)*ny + (size_t)(j)) * 2 * (nz / 2 + 1) + (size_t)(nz - k);
					idx[5] = ((size_t)(nx - i) * ny + (size_t)(j)) * 2 * (nz / 2 + 1) + (size_t)(nz - k);
					idx[6] = ((size_t)(i)*ny + (size_t)(ny - j)) * 2 * (nz / 2 + 1) + (size_t)(nz - k);
					idx[7] = ((size_t)(nx - i) * ny + (size_t)(ny - j)) * 2 * (nz / 2 + 1) + (size_t)(nz - k);

					if (i == 0 || i == nx / 2)
					{
						idx[1] = idx[3] = idx[5] = idx[7] = (size_t)-1;
					}
					if (j == 0 || j == ny / 2)
					{
						idx[2] = idx[3] = idx[6] = idx[7] = (size_t)-1;
					}
					if (k == 0 || k == nz / 2)
					{
						idx[4] = idx[5] = idx[6] = idx[7] = (size_t)-1;
					}

					double val = 0.0;

					for (int ii = -1; ii <= 1; ++ii)
						for (int jj = -1; jj <= 1; ++jj)
							for (int kk = -1; kk <= 1; ++kk)
							{
								rr[0] = ((double)iix) * dx + ii * boxlength;
								rr[1] = ((double)iiy) * dx + jj * boxlength;
								rr[2] = ((double)iiz) * dx + kk * boxlength;

								if (rr[0] > -boxlength && rr[0] <= boxlength && rr[1] > -boxlength && rr[1] <= boxlength && rr[2] > -boxlength && rr[2] <= boxlength)
								{
#ifdef OLD_KERNEL_SAMPLING
									if (ref_fac > 0)
									{
										double rrr[3];
										register double rrr2[3];
										for (int iii = ql; iii < qr; ++iii)
										{
											rrr[0] = rr[0] + (double)iii * dx05 - dx025;
											rrr2[0] = rrr[0] * rrr[0];
											for (int jjj = ql; jjj < qr; ++jjj)
											{
												rrr[1] = rr[1] + (double)jjj * dx05 - dx025;
												rrr2[1] = rrr[1] * rrr[1];
												for (int kkk = ql; kkk < qr; ++kkk)
												{
													rrr[2] = rr[2] + (double)kkk * dx05 - dx025;
													rrr2[2] = rrr[2] * rrr[2];
													val += tfr->compute_real(rrr2[0] + rrr2[1] + rrr2[2]) / rf8;
												}
											}
										}
									}
									else
									{
										val += tfr->compute_real(rr[0] * rr[0] + rr[1] * rr[1] + rr[2] * rr[2]);
									}

#else // !OLD_KERNEL_SAMPLING
									val += eval_split_recurse(tfr, rr, dx) / (dx * dx * dx);
#endif
								}
							}

					val *= fac;

					for (int q = 0; q < 8; ++q)
						if (idx[q] != (size_t)-1)
							rkernel[idx[q]] = val;
				}
	}
	else
	{
#pragma omp parallel for
		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j)
				for (int k = 0; k < nz; ++k)
				{
					int iix(i), iiy(j), iiz(k);
					real_t rr[3];

					if (iix > (int)nx / 2)
						iix -= nx;
					if (iiy > (int)ny / 2)
						iiy -= ny;
					if (iiz > (int)nz / 2)
						iiz -= nz;

					//size_t idx = ((size_t)i*ny + (size_t)j) * 2*(nz/2+1) + (size_t)k;

					rr[0] = ((double)iix) * dx;
					rr[1] = ((double)iiy) * dx;
					rr[2] = ((double)iiz) * dx;

					//rkernel[idx] = 0.0;

					//rr2 = rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2];

					//... speed up 8x by copying data to other octants
					size_t idx[8];

					idx[0] = ((size_t)(i)*ny + (size_t)(j)) * 2 * (nz / 2 + 1) + (size_t)(k);
					idx[1] = ((size_t)(nx - i) * ny + (size_t)(j)) * 2 * (nz / 2 + 1) + (size_t)(k);
					idx[2] = ((size_t)(i)*ny + (size_t)(ny - j)) * 2 * (nz / 2 + 1) + (size_t)(k);
					idx[3] = ((size_t)(nx - i) * ny + (size_t)(ny - j)) * 2 * (nz / 2 + 1) + (size_t)(k);
					idx[4] = ((size_t)(i)*ny + (size_t)(j)) * 2 * (nz / 2 + 1) + (size_t)(nz - k);
					idx[5] = ((size_t)(nx - i) * ny + (size_t)(j)) * 2 * (nz / 2 + 1) + (size_t)(nz - k);
					idx[6] = ((size_t)(i)*ny + (size_t)(ny - j)) * 2 * (nz / 2 + 1) + (size_t)(nz - k);
					idx[7] = ((size_t)(nx - i) * ny + (size_t)(ny - j)) * 2 * (nz / 2 + 1) + (size_t)(nz - k);

					if (i == 0 || i == nx / 2)
					{
						idx[1] = idx[3] = idx[5] = idx[7] = (size_t)-1;
					}
					if (j == 0 || j == ny / 2)
					{
						idx[2] = idx[3] = idx[6] = idx[7] = (size_t)-1;
					}
					if (k == 0 || k == nz / 2)
					{
						idx[4] = idx[5] = idx[6] = idx[7] = (size_t)-1;
					}

					double val = 0.0; //(fftw_real)tfr->compute_real(rr2)*fac;

#ifdef OLD_KERNEL_SAMPLING

					if (ref_fac > 0)
					{
						double rrr[3];
						register double rrr2[3];
						for (int iii = ql; iii < qr; ++iii)
						{
							rrr[0] = rr[0] + (double)iii * dx05 - dx025;
							rrr2[0] = rrr[0] * rrr[0];
							for (int jjj = ql; jjj < qr; ++jjj)
							{
								rrr[1] = rr[1] + (double)jjj * dx05 - dx025;
								rrr2[1] = rrr[1] * rrr[1];
								for (int kkk = ql; kkk < qr; ++kkk)
								{
									rrr[2] = rr[2] + (double)kkk * dx05 - dx025;
									rrr2[2] = rrr[2] * rrr[2];
									val += tfr->compute_real(rrr2[0] + rrr2[1] + rrr2[2]) / rf8;
								}
							}
						}
					}
					else
					{
						val = tfr->compute_real(rr[0] * rr[0] + rr[1] * rr[1] + rr[2] * rr[2]);
					}

#else
					if (i == 0 && j == 0 && k == 0)
						continue;

					// use new exact volume integration scheme
					val = eval_split_recurse(tfr, rr, dx) / (dx * dx * dx);

#endif

					//if( rr2 <= boxlength2*boxlength2 )
					//rkernel[idx] += (fftw_real)tfr->compute_real(rr2)*fac;
					val *= fac;

					for (int q = 0; q < 8; ++q)
						if (idx[q] != (size_t)-1)
							rkernel[idx[q]] = val;
				}
	}
	{
#ifdef OLD_KERNEL_SAMPLING
		rkernel[0] = tfr->compute_real(0.0) * fac;
#else
		real_t xmid[3] = {0.0, 0.0, 0.0};
		rkernel[0] = fac * eval_split_recurse(tfr, xmid, dx) / (dx * dx * dx);
#endif
	}
	/*************************************************************************************/
	/*************************************************************************************/
	/******* perform deconvolution *******************************************************/

	//bool baryons = type==baryon||type==vbaryon;
	if (deconv)
	{

		LOGUSER("Deconvolving fine kernel...");
		std::cout << " - Deconvolving density kernel...\n";

		double fftnorm = 1.0 / ((size_t)nx * (size_t)ny * (size_t)nz);
		double k0 = rkernel[0];

		fftw_complex *kkernel = reinterpret_cast<fftw_complex *>(&rkernel[0]);

		//... subtract white noise component before deconvolution
		if (!bsmooth_baryons)
			rkernel[0] = 0.0;

#ifdef FFTW3
#ifdef SINGLE_PRECISION
		fftwf_plan
			plan = fftwf_plan_dft_r2c_3d(nx, ny, nz, rkernel, kkernel, FFTW_ESTIMATE),
			iplan = fftwf_plan_dft_c2r_3d(nx, ny, nz, kkernel, rkernel, FFTW_ESTIMATE);

		fftwf_execute(plan);
#else
		fftw_plan
			plan = fftw_plan_dft_r2c_3d(nx, ny, nz, rkernel, kkernel, FFTW_ESTIMATE),
			iplan = fftw_plan_dft_c2r_3d(nx, ny, nz, kkernel, rkernel, FFTW_ESTIMATE);

		fftw_execute(plan);
#endif
#else
		rfftwnd_plan plan = rfftw3d_create_plan(nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE),
					 iplan = rfftw3d_create_plan(nx, ny, nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

#ifndef SINGLETHREAD_FFTW
		rfftwnd_threads_one_real_to_complex(omp_get_max_threads(), plan, rkernel, NULL);
#else
		rfftwnd_one_real_to_complex(plan, rkernel, NULL);
#endif
#endif

		if (deconv)
		{

			double ksum = 0.0;
			size_t kcount = 0;
			double kmax = 0.5 * M_PI / std::max(nx, std::max(ny, nz));

#pragma omp parallel for reduction(+ \
								   : ksum, kcount)
			for (int i = 0; i < nx; ++i)
				for (int j = 0; j < ny; ++j)
					for (int k = 0; k < nz / 2 + 1; ++k)
					{
						double kx, ky, kz;

						kx = (double)i;
						ky = (double)j;
						kz = (double)k;

						if (kx > nx / 2)
							kx -= nx;
						if (ky > ny / 2)
							ky -= ny;

						double kkmax = kmax;
						size_t q = ((size_t)i * ny + (size_t)j) * (size_t)(nz / 2 + 1) + (size_t)k;

						if (!bsmooth_baryons)
						{
							if (kspacepoisson)
							{
								//... Use child average response function to emulate sub-sampling
								double ipix = cos(kx * kkmax) * cos(ky * kkmax) * cos(kz * kkmax);

								RE(kkernel[q]) /= ipix;
								IM(kkernel[q]) /= ipix;
							}
							else
							{

								//... Use piecewise constant response function (NGP-kernel)
								//... for finite difference methods
								kkmax = kmax;
								double ipix = 1.0;
								if (i > 0)
									ipix /= sin(kx * 2.0 * kkmax) / (kx * 2.0 * kkmax);
								if (j > 0)
									ipix /= sin(ky * 2.0 * kkmax) / (ky * 2.0 * kkmax);
								if (k > 0)
									ipix /= sin(kz * 2.0 * kkmax) / (kz * 2.0 * kkmax);

								RE(kkernel[q]) *= ipix;
								IM(kkernel[q]) *= ipix;
							}
						}
#if 1
						else
						{
							//... if smooth==true, convolve with
							//... NGP kernel to get CIC smoothness

							//kkmax *= 2.0;

							double ipix = 1.0;
							if (i > 0)
								ipix /= sin(kx * 2.0 * kkmax) / (kx * 2.0 * kkmax);
							if (j > 0)
								ipix /= sin(ky * 2.0 * kkmax) / (ky * 2.0 * kkmax);
							if (k > 0)
								ipix /= sin(kz * 2.0 * kkmax) / (kz * 2.0 * kkmax);

							RE(kkernel[q]) /= ipix;
							IM(kkernel[q]) /= ipix;
						}
#endif

						//... store k-space average
						if (k == 0 || k == nz / 2)
						{
							ksum += RE(kkernel[q]);
							kcount++;
						}
						else
						{
							ksum += 2.0 * (RE(kkernel[q]));
							kcount += 2;
						}
					}

			double dk;

			//... re-add white noise component for finest grid
			dk = k0 - ksum / kcount;

			//... set white noise component to zero if smoothing is enabled
			//if( false )//cparam_.smooth )
			if (bsmooth_baryons)
				dk = 0.0;

				//... enforce the r=0 component by adjusting the k-space mean
#pragma omp parallel for reduction(+ \
								   : ksum, kcount)
			for (int i = 0; i < nx; ++i)
				for (int j = 0; j < ny; ++j)
					for (int k = 0; k < (nz / 2 + 1); ++k)
					{
						size_t q = ((size_t)i * ny + (size_t)j) * (nz / 2 + 1) + (size_t)k;

						RE(kkernel[q]) += dk;

						RE(kkernel[q]) *= fftnorm;
						IM(kkernel[q]) *= fftnorm;
					}
		}

#ifdef FFTW3
#ifdef SINGLE_PRECISION
		fftwf_execute(iplan);
		fftwf_destroy_plan(plan);
		fftwf_destroy_plan(iplan);
#else
		fftw_execute(iplan);
		fftw_destroy_plan(plan);
		fftw_destroy_plan(iplan);
#endif
#else
#ifndef SINGLETHREAD_FFTW
		rfftwnd_threads_one_complex_to_real(omp_get_max_threads(), iplan, kkernel, NULL);
#else
		rfftwnd_one_complex_to_real(iplan, kkernel, NULL);
#endif
		rfftwnd_destroy_plan(plan);
		rfftwnd_destroy_plan(iplan);
#endif
	}

	/*************************************************************************************/
	/*************************************************************************************/
	/*************************************************************************************/

	char cachefname[128];
	sprintf(cachefname, "temp_kernel_level%03d.tmp", levelmax);
	LOGUSER("Storing kernel in temp file \'%s\'.", cachefname);

	FILE *fp = fopen(cachefname, "w+");
	unsigned q = nx;
	fwrite(reinterpret_cast<void *>(&q), sizeof(unsigned), 1, fp);
	q = ny;
	fwrite(reinterpret_cast<void *>(&q), sizeof(unsigned), 1, fp);
	q = 2 * (nz / 2 + 1);
	fwrite(reinterpret_cast<void *>(&q), sizeof(unsigned), 1, fp);

	for (int ix = 0; ix < nx; ++ix)
	{
		size_t sz = ny * 2 * (nz / 2 + 1);
		//fwrite( reinterpret_cast<void*>(&rkernel[0]), sizeof(fftw_real), nx*ny*2*(nz/2+1), fp );
		fwrite(reinterpret_cast<void *>(&rkernel[(size_t)ix * sz]), sizeof(fftw_real), sz, fp);
	}

	fclose(fp);

	//... average and fill for other levels
	for (int ilevel = levelmax - 1; ilevel >= levelmin; ilevel--)
	{
		LOGUSER("Computing coarse kernel (level %d)...", ilevel);

		int nxc, nyc, nzc;
		real_t dxc, lxc, lyc, lzc;

		nxc = refh.size(ilevel, 0);
		nyc = refh.size(ilevel, 1);
		nzc = refh.size(ilevel, 2);

		if (ilevel != levelmin)
		{
			nxc *= 2;
			nyc *= 2;
			nzc *= 2;
		}

		dxc = boxlength / (1 << ilevel);
		lxc = dxc * nxc;
		lyc = dxc * nyc;
		lzc = dxc * nzc;

		rkernel_coarse = new fftw_real[(size_t)nxc * (size_t)nyc * 2 * ((size_t)nzc / 2 + 1)];
		fac = lxc * lyc * lzc / pow(2.0 * M_PI, 3) / ((double)nxc * (double)nyc * (double)nzc);

		if (bperiodic)
		{
#pragma omp parallel for
			for (int i = 0; i <= nxc / 2; ++i)
				for (int j = 0; j <= nyc / 2; ++j)
					for (int k = 0; k <= nzc / 2; ++k)
					{
						int iix(i), iiy(j), iiz(k);
						real_t rr[3], rr2;

						if (iix > (int)nxc / 2)
							iix -= nxc;
						if (iiy > (int)nyc / 2)
							iiy -= nyc;
						if (iiz > (int)nzc / 2)
							iiz -= nzc;

						//... speed up 8x by copying data to other octants
						size_t idx[8];

						idx[0] = ((size_t)(i)*nyc + (size_t)(j)) * 2 * (nzc / 2 + 1) + (size_t)(k);
						idx[1] = ((size_t)(nxc - i) * nyc + (size_t)(j)) * 2 * (nzc / 2 + 1) + (size_t)(k);
						idx[2] = ((size_t)(i)*nyc + (size_t)(nyc - j)) * 2 * (nzc / 2 + 1) + (size_t)(k);
						idx[3] = ((size_t)(nxc - i) * nyc + (size_t)(nyc - j)) * 2 * (nzc / 2 + 1) + (size_t)(k);
						idx[4] = ((size_t)(i)*nyc + (size_t)(j)) * 2 * (nzc / 2 + 1) + (size_t)(nzc - k);
						idx[5] = ((size_t)(nxc - i) * nyc + (size_t)(j)) * 2 * (nzc / 2 + 1) + (size_t)(nzc - k);
						idx[6] = ((size_t)(i)*nyc + (size_t)(nyc - j)) * 2 * (nzc / 2 + 1) + (size_t)(nzc - k);
						idx[7] = ((size_t)(nxc - i) * nyc + (size_t)(nyc - j)) * 2 * (nzc / 2 + 1) + (size_t)(nzc - k);

						if (i == 0 || i == nxc / 2)
						{
							idx[1] = idx[3] = idx[5] = idx[7] = (size_t)-1;
						}
						if (j == 0 || j == nyc / 2)
						{
							idx[2] = idx[3] = idx[6] = idx[7] = (size_t)-1;
						}
						if (k == 0 || k == nzc / 2)
						{
							idx[4] = idx[5] = idx[6] = idx[7] = (size_t)-1;
						}

						double val = 0.0;

						for (int ii = -1; ii <= 1; ++ii)
							for (int jj = -1; jj <= 1; ++jj)
								for (int kk = -1; kk <= 1; ++kk)
								{
									rr[0] = ((double)iix) * dxc + ii * boxlength;
									rr[1] = ((double)iiy) * dxc + jj * boxlength;
									rr[2] = ((double)iiz) * dxc + kk * boxlength;

									if (rr[0] > -boxlength && rr[0] < boxlength && rr[1] > -boxlength && rr[1] < boxlength && rr[2] > -boxlength && rr[2] < boxlength)
									{
#ifdef OLD_KERNEL_SAMPLING
										rr2 = rr[0] * rr[0] + rr[1] * rr[1] + rr[2] * rr[2];
										val += tfr->compute_real(rr2);
#else // ! OLD_KERNEL_SAMPLING
										val += eval_split_recurse(tfr, rr, dxc) / (dxc * dxc * dxc);
#endif
									}
								}

						val *= fac;

						for (int qq = 0; qq < 8; ++qq)
							if (idx[qq] != (size_t)-1)
								rkernel_coarse[idx[qq]] = val;
					}
		}
		else
		{
#pragma omp parallel for
			for (int i = 0; i < nxc; ++i)
				for (int j = 0; j < nyc; ++j)
					for (int k = 0; k < nzc; ++k)
					{
						int iix(i), iiy(j), iiz(k);
						real_t rr[3];

						if (iix > (int)nxc / 2)
							iix -= nxc;
						if (iiy > (int)nyc / 2)
							iiy -= nyc;
						if (iiz > (int)nzc / 2)
							iiz -= nzc;

						size_t idx = ((size_t)i * nyc + (size_t)j) * 2 * (nzc / 2 + 1) + (size_t)k;

						rr[0] = ((double)iix) * dxc;
						rr[1] = ((double)iiy) * dxc;
						rr[2] = ((double)iiz) * dxc;

#ifdef OLD_KERNEL_SAMPLING
						rkernel_coarse[idx] = 0.0;

						real_t rr2 = rr[0] * rr[0] + rr[1] * rr[1] + rr[2] * rr[2];
						if (fabs(rr[0]) <= boxlength2 || fabs(rr[1]) <= boxlength2 || fabs(rr[2]) <= boxlength2)
							rkernel_coarse[idx] += (fftw_real)tfr->compute_real(rr2) * fac;
#else

						rkernel_coarse[idx] = 0.0;

						//if( i==0 && j==0 && k==0 ) continue;
						real_t val = eval_split_recurse(tfr, rr, dxc) / (dxc * dxc * dxc);

						if (fabs(rr[0]) <= boxlength2 || fabs(rr[1]) <= boxlength2 || fabs(rr[2]) <= boxlength2)
							rkernel_coarse[idx] += val * fac;

#endif
					}
		}

#ifdef OLD_KERNEL_SAMPLING
		LOGUSER("Averaging fine kernel to coarse kernel...");

//... copy averaged and convolved fine kernel to coarse kernel
#pragma omp parallel for
		for (int ix = 0; ix < nx; ix += 2)
			for (int iy = 0; iy < ny; iy += 2)
				for (int iz = 0; iz < nz; iz += 2)
				{
					int iix(ix / 2), iiy(iy / 2), iiz(iz / 2);
					if (ix > nx / 2)
						iix += nxc - nx / 2;
					if (iy > ny / 2)
						iiy += nyc - ny / 2;
					if (iz > nz / 2)
						iiz += nzc - nz / 2;

					if (ix == nx / 2 || iy == ny / 2 || iz == nz / 2)
						continue;

					for (int i = 0; i <= 1; ++i)
						for (int j = 0; j <= 1; ++j)
							for (int k = 0; k <= 1; ++k)
								if (i == 0 && k == 0 && j == 0)
									rkernel_coarse[ACC_RC(iix, iiy, iiz)] =
										0.125 * (rkernel[ACC_RF(ix - i, iy - j, iz - k)] + rkernel[ACC_RF(ix - i + 1, iy - j, iz - k)] + rkernel[ACC_RF(ix - i, iy - j + 1, iz - k)] + rkernel[ACC_RF(ix - i, iy - j, iz - k + 1)] + rkernel[ACC_RF(ix - i + 1, iy - j + 1, iz - k)] + rkernel[ACC_RF(ix - i + 1, iy - j, iz - k + 1)] + rkernel[ACC_RF(ix - i, iy - j + 1, iz - k + 1)] + rkernel[ACC_RF(ix - i + 1, iy - j + 1, iz - k + 1)]);

								else
								{

									rkernel_coarse[ACC_RC(iix, iiy, iiz)] +=
										0.125 * (rkernel[ACC_RF(ix - i, iy - j, iz - k)] + rkernel[ACC_RF(ix - i + 1, iy - j, iz - k)] + rkernel[ACC_RF(ix - i, iy - j + 1, iz - k)] + rkernel[ACC_RF(ix - i, iy - j, iz - k + 1)] + rkernel[ACC_RF(ix - i + 1, iy - j + 1, iz - k)] + rkernel[ACC_RF(ix - i + 1, iy - j, iz - k + 1)] + rkernel[ACC_RF(ix - i, iy - j + 1, iz - k + 1)] + rkernel[ACC_RF(ix - i + 1, iy - j + 1, iz - k + 1)]);
								}
				}

#endif // #OLD_KERNEL_SAMPLING
		sprintf(cachefname, "temp_kernel_level%03d.tmp", ilevel);
		LOGUSER("Storing kernel in temp file \'%s\'.", cachefname);
		fp = fopen(cachefname, "w+");
		q = nxc;
		fwrite(reinterpret_cast<void *>(&q), sizeof(unsigned), 1, fp);
		q = nyc;
		fwrite(reinterpret_cast<void *>(&q), sizeof(unsigned), 1, fp);
		q = 2 * (nzc / 2 + 1);
		fwrite(reinterpret_cast<void *>(&q), sizeof(unsigned), 1, fp);

		for (int ix = 0; ix < nxc; ++ix)
		{
			size_t sz = nyc * 2 * (nzc / 2 + 1);
			//fwrite( reinterpret_cast<void*>(&rkernel_coarse[0]), sizeof(fftw_real), nxc*nyc*2*(nzc/2+1), fp );
			fwrite(reinterpret_cast<void *>(&rkernel_coarse[ix * sz]), sizeof(fftw_real), sz, fp);
		}

		fclose(fp);

		delete[] rkernel;

		//... prepare for next iteration
		nx = nxc;
		ny = nyc;
		nz = nzc;
		lx = lxc;
		ly = lyc;
		lz = lzc;
		dx = dxc;
		rkernel = rkernel_coarse;
	}

	//... clean up
	delete[] rkernel;
}

} // namespace convolution

/**************************************************************************************/
/**************************************************************************************/

namespace
{
convolution::kernel_creator_concrete<convolution::kernel_real_cached<double>> creator_d("tf_kernel_real_double");
convolution::kernel_creator_concrete<convolution::kernel_real_cached<float>> creator_f("tf_kernel_real_float");

convolution::kernel_creator_concrete<convolution::kernel_k<double>> creator_kd("tf_kernel_k_double");
convolution::kernel_creator_concrete<convolution::kernel_k<float>> creator_kf("tf_kernel_k_float");
} // namespace
