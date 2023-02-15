/*
 
 convolution_kernel.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions for cosmological simulations 
 
 Copyright (C) 2010-23  Oliver Hahn
 
*/

#include <general.hh>
#include <densities.hh>
#include <convolution_kernel.hh>

namespace convolution
{

std::map<std::string, kernel_creator *> &
get_kernel_map()
{
	static std::map<std::string, kernel_creator *> kernel_map;
	return kernel_map;
}

void perform(kernel *pk, void *pd, bool shift, bool fix, bool flip)
{
	//return;

	parameters cparam_ = pk->cparam_;
	double fftnormp = 1.0/sqrt((double)cparam_.nx * (double)cparam_.ny * (double)cparam_.nz);
	double fftnorm = pow(2.0 * M_PI, 1.5) / sqrt(cparam_.lx * cparam_.ly * cparam_.lz) * fftnormp;

	complex_t *cdata;
	[[maybe_unused]] complex_t *ckernel;
	real_t *data;

	data = reinterpret_cast<real_t *>(pd);
	cdata = reinterpret_cast<complex_t *>(data);
	ckernel = reinterpret_cast<complex_t *>(pk->get_ptr());

	std::cout << "   - Performing density convolution... ("
			  << cparam_.nx << ", " << cparam_.ny << ", " << cparam_.nz << ")\n";

	music::ulog.Print("Performing kernel convolution on (%5d,%5d,%5d) grid", cparam_.nx, cparam_.ny, cparam_.nz);
	music::ulog.Print("Performing forward FFT...");

	fftw_plan_t plan, iplan;
	plan = FFTW_API(plan_dft_r2c_3d)(cparam_.nx, cparam_.ny, cparam_.nz, data, cdata, FFTW_ESTIMATE);
	iplan = FFTW_API(plan_dft_c2r_3d)(cparam_.nx, cparam_.ny, cparam_.nz, cdata, data, FFTW_ESTIMATE);

	FFTW_API(execute)(plan);

	//..... need a phase shift for baryons for SPH
	double dstag = 0.0;

	if (shift)
	{
		double boxlength = pk->pcf_->get_value<double>("setup", "boxlength");
		double stagfact = pk->pcf_->get_value_safe<double>("setup", "baryon_staggering", 0.5);
		int lmax = pk->pcf_->get_value<int>("setup", "levelmax");
		double dxmax = boxlength / (1 << lmax);
		double dxcur = cparam_.lx / cparam_.nx;
		//std::cerr << "Performing staggering shift for SPH\n";
		music::ulog.Print("Performing staggering shift for SPH");
		dstag = stagfact * 2.0 * M_PI / cparam_.nx * dxmax / dxcur;
	}

	//.............................................

	std::complex<double> dcmode(RE(cdata[0]), IM(cdata[0]));

	
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


	music::ulog.Print("Performing backward FFT...");

	FFTW_API(execute)(iplan);
	FFTW_API(destroy_plan)(plan);
	FFTW_API(destroy_plan)(iplan);

	// set the DC mode here to avoid a possible truncation error in single precision
	{
		size_t nelem = (size_t)cparam_.nx * (size_t)cparam_.ny * (size_t)cparam_.nz;
		real_t mean = dcmode.real() * fftnorm / (real_t)nelem;

#pragma omp parallel for
		for (size_t i = 0; i < nelem; ++i)
			data[i] += mean;
	}
}

void perform(kernel *pk, void *pd, bool shift, bool fix, bool flip);

/*****************************************************************************************/
/***    SPECIFIC KERNEL IMPLEMENTATIONS      *********************************************/
/*****************************************************************************************/

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
		boxlength_ = pcf_->get_value<double>("setup", "boxlength");
		nspec_ = pcf_->get_value<double>("cosmology", "nspec");
		pnorm_ = pcf_->get_value<double>("cosmology", "pnorm");
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
			if( prefh_->get_margin() < 0 ){
				cparam_.nx = 2 * prefh_->size(ilevel, 0);
				cparam_.ny = 2 * prefh_->size(ilevel, 1);
				cparam_.nz = 2 * prefh_->size(ilevel, 2);
			}else{
				cparam_.nx = prefh_->size(ilevel, 0) + 2*prefh_->get_margin();
				cparam_.ny = prefh_->size(ilevel, 1) + 2*prefh_->get_margin();
				cparam_.nz = prefh_->size(ilevel, 2) + 2*prefh_->get_margin();
			}

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

///////////////////////////////////////////////////////////////////////////

} // namespace convolution

/**************************************************************************************/
/**************************************************************************************/

convolution::kernel_creator_concrete<convolution::kernel_k> creator_kd("tf_kernel_k");
