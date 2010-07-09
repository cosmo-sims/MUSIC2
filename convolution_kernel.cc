/*
 
 convolution_kernel.cc - This file is part of MUSIC -
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

#ifdef SINGLE_PRECISION
	#ifdef SINGLETHREAD_FFTW
		#include <srfftw.h>
	#else
		#include <srfftw_threads.h>
	#endif
#else
	#ifdef SINGLETHREAD_FFTW
		#include <drfftw.h>
	#else
		#include <drfftw_threads.h>
	#endif
#endif


#include "densities.hh"
#include "convolution_kernel.hh"



namespace convolution{

	std::map< std::string, kernel_creator *>& 
	get_kernel_map()
	{
		static std::map< std::string, kernel_creator* > kernel_map;
		return kernel_map;
	}
	
	template< typename real_t >
	void perform( kernel * pk, void *pd )
	{
		parameters cparam_ = pk->cparam_;
		double fftnorm = pow(2.0*M_PI,1.5)/sqrt(cparam_.lx*cparam_.ly*cparam_.lz)/sqrt((double)(cparam_.nx*cparam_.ny*cparam_.nz));
		
		fftw_complex	*cdata,*ckernel;
		fftw_real		*data;
		
		data		= reinterpret_cast<fftw_real*>(pd);
		cdata		= reinterpret_cast<fftw_complex*>(data);
		ckernel		= reinterpret_cast<fftw_complex*>( pk->get_ptr() );
		
		rfftwnd_plan	iplan, plan;
		
		plan  = rfftw3d_create_plan( cparam_.nx, cparam_.ny, cparam_.nz,
									FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE);
		
		iplan = rfftw3d_create_plan( cparam_.nx, cparam_.ny, cparam_.nz,
									FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
		
		
		#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), plan, data, NULL );
		#else
		rfftwnd_one_real_to_complex( plan, data, NULL );
		#endif
		
		#pragma omp parallel for
		for( int i=0; i<cparam_.nx; ++i )
			for( int j=0; j<cparam_.ny; ++j )
				for( int k=0; k<cparam_.nz/2+1; ++k )
				{
					unsigned ii = (i*cparam_.ny + j) * (cparam_.nz/2+1) + k;
					
					complex ccdata(cdata[ii].re,cdata[ii].im), cckernel(ckernel[ii].re,ckernel[ii].im);
					ccdata = ccdata * cckernel *fftnorm;
					
					cdata[ii].re = ccdata.real();
					cdata[ii].im = ccdata.imag();
				}
		
		#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), iplan, cdata, NULL);
		#else		
		rfftwnd_one_complex_to_real(iplan, cdata, NULL);
		#endif
		
		rfftwnd_destroy_plan(plan);
		rfftwnd_destroy_plan(iplan);
	}
	
	
	template< typename real_t >
	void perform_filtered( kernel * pk, void *pd )
	{
		parameters cparam_ = pk->cparam_;
		double 
			kny, kmax,
			fftnorm = pow(2.0*M_PI,1.5)/sqrt(cparam_.lx*cparam_.ly*cparam_.lz)/sqrt((double)(cparam_.nx*cparam_.ny*cparam_.nz));
		
		fftw_complex	*cdata,*ckernel;
		fftw_real		*data;
		
		kny			= cparam_.nx*M_PI/cparam_.lx;
		kmax        = M_PI/2.0;
		
		data		= reinterpret_cast<fftw_real*>(pd);
		cdata		= reinterpret_cast<fftw_complex*>(data);
		ckernel		= reinterpret_cast<fftw_complex*>( pk->get_ptr() );
		
		rfftwnd_plan	iplan, plan;
		
		plan  = rfftw3d_create_plan( cparam_.nx, cparam_.ny, cparam_.nz,
									FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE);
		
		iplan = rfftw3d_create_plan( cparam_.nx, cparam_.ny, cparam_.nz,
									FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
		
		
		#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), plan, data, NULL );
		#else
		rfftwnd_one_real_to_complex( plan, data, NULL );
		#endif
		
		//double kmax2 = kmax*kmax;
		
		#pragma omp parallel for
		for( int i=0; i<cparam_.nx; ++i )
			for( int j=0; j<cparam_.ny; ++j )
				for( int k=0; k<cparam_.nz/2+1; ++k )
				{
					unsigned ii = (i*cparam_.ny + j) * (cparam_.nz/2+1) + k;
					
					complex ccdata(cdata[ii].re,cdata[ii].im), cckernel(ckernel[ii].re,ckernel[ii].im);
					
					int ik(i), jk(j);
					
					if( ik > cparam_.nx/2 )
						ik -= cparam_.nx;
					if( jk > cparam_.ny/2 )
						jk -= cparam_.ny;
					
					double 
						kx( M_PI*(double)ik/cparam_.nx ),
						ky( M_PI*(double)jk/cparam_.ny ),
					kz( M_PI*(double)k/cparam_.nz );//,
													//						kk( kx*kx+ky*ky+kz*kz );
					
					//... cos(k) is the Hanning filter function (cf. Bertschinger 2001)
					double filter = 0.0;
					if( true ){//kk <= kmax2 ){
						//filter = cos( kx )*cos( ky )*cos( kz );
						filter = cos( 0.5*kx )*cos( 0.5*ky )*cos( 0.5*kz );
							   //filter = 1.0;
						//filter *=filter;
					}
					//filter = 1.0;
					
					ccdata = ccdata * cckernel *fftnorm * filter;
					
					cdata[ii].re = ccdata.real();
					cdata[ii].im = ccdata.imag();
				}
		
		#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), iplan, cdata, NULL);
		#else		
		rfftwnd_one_complex_to_real(iplan, cdata, NULL);
		#endif
		
		rfftwnd_destroy_plan(plan);
		rfftwnd_destroy_plan(iplan);
	}
	
	
	template void perform<double>( kernel* pk, void *pd );
	template void perform<float>( kernel* pk, void *pd );
	template void perform_filtered<double>( kernel* pk, void *pd );
	template void perform_filtered<float>( kernel* pk, void *pd );
	
	
	void truncate( kernel* pk, double R, double alpha )
	{
		parameters cparam_ = pk->cparam_;

		double 
			dx		= cparam_.lx/cparam_.nx, 
			dy		= cparam_.ly/cparam_.ny, 
			dz		= cparam_.lz/cparam_.nz;
		
		double fftnorm = 1.0/(cparam_.nx * cparam_.ny * cparam_.nz);
		
		rfftwnd_plan	iplan, plan;
		
		plan  = rfftw3d_create_plan( cparam_.nx, cparam_.ny, cparam_.nz,
									FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE);
		
		iplan = rfftw3d_create_plan( cparam_.nx, cparam_.ny, cparam_.nz,
									FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
		
		fftw_real		*rkernel	= reinterpret_cast<fftw_real*>( pk->get_ptr() );
		fftw_complex	*ckernel	= reinterpret_cast<fftw_complex*>( pk->get_ptr() );
		
		#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), iplan, ckernel, NULL);
		#else		
		rfftwnd_one_complex_to_real(iplan, ckernel, NULL);
		#endif
		
		
		#pragma omp parallel for
		for( int i=0; i<cparam_.nx; ++i )
			for( int j=0; j<cparam_.ny; ++j )
				for( int k=0; k<cparam_.nz; ++k )
				{
					int iix(i), iiy(j), iiz(k);
					double rr[3], rr2;
					
					if( iix > (int)cparam_.nx/2 ) iix -= cparam_.nx;
					if( iiy > (int)cparam_.ny/2 ) iiy -= cparam_.ny;
					if( iiz > (int)cparam_.nz/2 ) iiz -= cparam_.nz;
					
					rr[0] = (double)iix * dx;
					rr[1] = (double)iiy * dy;
					rr[2] = (double)iiz * dz;
					
					rr2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2];
					
					unsigned idx = (i*cparam_.ny + j) * 2*(cparam_.nz/2+1) + k;
					rkernel[idx] *= 0.5*(erfc((sqrt(rr2)-R)*alpha))*fftnorm;
				}
		
		
		
		#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), plan, rkernel, NULL );
		#else
		rfftwnd_one_real_to_complex( plan, rkernel, NULL );
		#endif
		
	}
	
	void truncate_sharp( kernel* pk, double R )
	{
		parameters cparam_ = pk->cparam_;
		
		double 
		dx		= cparam_.lx/cparam_.nx, 
		dy		= cparam_.ly/cparam_.ny, 
		dz		= cparam_.lz/cparam_.nz,
		R2		= R*R;
		
		double fftnorm = 1.0/(cparam_.nx * cparam_.ny * cparam_.nz);
		
		rfftwnd_plan	iplan, plan;
		
		plan  = rfftw3d_create_plan( cparam_.nx, cparam_.ny, cparam_.nz,
									FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE);
		
		iplan = rfftw3d_create_plan( cparam_.nx, cparam_.ny, cparam_.nz,
									FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
		
		fftw_real		*rkernel	= reinterpret_cast<fftw_real*>( pk->get_ptr() );
		fftw_complex	*ckernel	= reinterpret_cast<fftw_complex*>( pk->get_ptr() );
		
#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), iplan, ckernel, NULL);
#else		
		rfftwnd_one_complex_to_real(iplan, ckernel, NULL);
#endif
		
		
#pragma omp parallel for
		for( int i=0; i<cparam_.nx; ++i )
			for( int j=0; j<cparam_.ny; ++j )
				for( int k=0; k<cparam_.nz; ++k )
				{
					int iix(i), iiy(j), iiz(k);
					double rr[3], rr2;
					
					if( iix > (int)cparam_.nx/2 ) iix -= cparam_.nx;
					if( iiy > (int)cparam_.ny/2 ) iiy -= cparam_.ny;
					if( iiz > (int)cparam_.nz/2 ) iiz -= cparam_.nz;
					
					rr[0] = (double)iix * dx;
					rr[1] = (double)iiy * dy;
					rr[2] = (double)iiz * dz;
					
					rr2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2];
					unsigned idx = (i*cparam_.ny + j) * 2*(cparam_.nz/2+1) + k;
					if( rr2 > R2 )
					{
						
						rkernel[idx] = 0.0;
					}else {
						rkernel[idx] *= fftnorm;
					}

				}
		
		
		
#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), plan, rkernel, NULL );
#else
		rfftwnd_one_real_to_complex( plan, rkernel, NULL );
#endif
		
	}


	/*****************************************************************************************\
	 ***    SPECIFIC KERNEL IMPLEMENTATIONS      *********************************************
	\*****************************************************************************************/	
	
	template< typename real_t >
	class kernel_real : public kernel
	{
	protected:
		std::vector<real_t> kdata_;

		void compute_kernel( void );
		
	public:
		kernel_real( const parameters& cp )
		: kernel( cp )
		{
			kdata_.assign( cparam_.nx*cparam_.ny*2*(cparam_.nz/2+1), 0.0 );
			compute_kernel();
		}
		
		void *get_ptr()
		{	return reinterpret_cast<void*> (&kdata_[0]);	}
		
		~kernel_real()
		{	std::vector<real_t>().swap( kdata_ );	}
		
	};


	
	
	
	template< typename real_t >
	void kernel_real<real_t>::compute_kernel( void )
	{
		double 
			kny			= cparam_.nx*M_PI/cparam_.lx,
			fac			= cparam_.lx*cparam_.ly*cparam_.lz/pow(2.0*M_PI,3)/(cparam_.nx*cparam_.ny*cparam_.nz),
			dx			= cparam_.lx/cparam_.nx, 
			dy			= cparam_.ly/cparam_.ny, 
			dz			= cparam_.lz/cparam_.nz,
			boxlength	= cparam_.pcf->getValue<double>("setup","boxlength"),
			nspec		= cparam_.pcf->getValue<double>("cosmology","nspec"),//cosmo.nspect,
			pnorm		= cparam_.pcf->getValue<double>("cosmology","pnorm"),//cosmo.pnorm,
		dplus		= cparam_.pcf->getValue<double>("cosmology","dplus");//,//cosmo.dplus;
																		 //			t00f        = cparam_.t0scale;
		unsigned 
		levelmax	= cparam_.pcf->getValue<unsigned>("setup","levelmax");
		
		bool
		bperiodic	= cparam_.pcf->getValueSafe<bool>("setup","periodic_TF",true);
		
		//if( cparam_.normalize )
		//	kny *= 2.0;
		
		TransferFunction_real *tfr = new TransferFunction_real(cparam_.ptf,nspec,pnorm,dplus,0.25*cparam_.lx/(double)cparam_.nx,2.0*boxlength,kny, (int)pow(2,levelmax+2));
		
		if( bperiodic )
		{
			
			#pragma omp parallel for
			for( int i=0; i<cparam_.nx; ++i )
				for( int j=0; j<cparam_.ny; ++j )
					for( int k=0; k<cparam_.nz; ++k )
					{
						int iix(i), iiy(j), iiz(k);
						double rr[3], rr2;
						
						if( iix > (int)cparam_.nx/2 ) iix -= cparam_.nx;
						if( iiy > (int)cparam_.ny/2 ) iiy -= cparam_.ny;
						if( iiz > (int)cparam_.nz/2 ) iiz -= cparam_.nz;
						
						unsigned idx = (i*cparam_.ny + j) * 2*(cparam_.nz/2+1) + k;
						
						for( int ii=-1; ii<=1; ++ii )
							for( int jj=-1; jj<=1; ++jj )
								for( int kk=-1; kk<=1; ++kk )
								{
									rr[0] = ((double)iix ) * dx + ii*boxlength;
									rr[1] = ((double)iiy ) * dy + jj*boxlength;
									rr[2] = ((double)iiz ) * dz + kk*boxlength;
									
									if( fabs(rr[0]) < boxlength
									 && fabs(rr[1]) < boxlength
									 && fabs(rr[2]) < boxlength )
									{
										rr2 = rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2];
										kdata_[idx] += tfr->compute_real(rr2);
									}
								}
						
						kdata_[idx] *= fac;
						
					}
		}else{
			#pragma omp parallel for
			for( int i=0; i<cparam_.nx; ++i )
				for( int j=0; j<cparam_.ny; ++j )
					for( int k=0; k<cparam_.nz; ++k )
					{
						int iix(i), iiy(j), iiz(k);
						double rr[3], rr2;
						
						if( iix > (int)cparam_.nx/2 ) iix -= cparam_.nx;
						if( iiy > (int)cparam_.ny/2 ) iiy -= cparam_.ny;
						if( iiz > (int)cparam_.nz/2 ) iiz -= cparam_.nz;
						
						unsigned idx = (i*cparam_.ny + j) * 2*(cparam_.nz/2+1) + k;
						
						rr[0] = ((double)iix ) * dx;
						rr[1] = ((double)iiy ) * dy;
						rr[2] = ((double)iiz ) * dz;
						
						rr2 = rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2];
						kdata_[idx] = tfr->compute_real(rr2)*fac;
						
						
					}
		}		
		//kdata_[0] = tfr->compute_real(dx*dx*0.25)*fac;
		
		//std::cerr << "T(r=0) = " << kdata_[0]/fac << std::endl;
		//if( cparam_.normalize )
		//	kdata_[0] *= 0.125;//= 0.0;//*= 0.125;
		
		//... determine normalization
		double sum = 0.0;
		#pragma omp parallel for reduction(+:sum)
		for( int i=0; i<cparam_.nx; ++i )
			for( int j=0; j<cparam_.ny; ++j )
				for( int k=0; k<cparam_.nz; ++k )
				{	
					unsigned idx = (i*cparam_.ny + j) * 2*(cparam_.nz/2+1) + k;
					//if( idx > 0 )
						sum += kdata_[idx];
				}
		
		
		//std::cerr << " - The convolution kernel has avg of " << (sum+kdata_[0])/(cparam_.nx*cparam_.ny*cparam_.nz) << std::endl;
		
		sum /= cparam_.nx*cparam_.ny*cparam_.nz;//-1;
		
		if( false )//cparam_.normalize )
		{
			
			#pragma omp parallel for
			for( int i=0; i<cparam_.nx; ++i )
				for( int j=0; j<cparam_.ny; ++j )
					for( int k=0; k<cparam_.nz; ++k )
					{
						unsigned idx = (i*cparam_.ny + j) * 2*(cparam_.nz/2+1) + k;
						//if( idx > 0 )
							kdata_[idx] -= sum;
					}
		}
		
		
		//if( t00f < 0.0 )
		//	kdata_[0] += (tfr->compute_real(0.125*dx*dx) - tfr->compute_real(0.0))*fac;
		
		fftw_real		*rkernel	= reinterpret_cast<fftw_real*>( &kdata_[0] );
		rfftwnd_plan	plan		= rfftw3d_create_plan( cparam_.nx, cparam_.ny, cparam_.nz,
												FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE);
		
		#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), plan, rkernel, NULL );
		#else
		rfftwnd_one_real_to_complex( plan, rkernel, NULL );
		#endif
	
		rfftwnd_destroy_plan(plan);
		
		delete tfr;
	}
	



	template< typename real_t >
	class kernel_k : public kernel
	{
	protected:
		std::vector<real_t> kdata_;
		
		void compute_kernel( void );
		
	public:
		kernel_k( const parameters& cp )
		: kernel( cp )
		{
			kdata_.assign( cparam_.nx*cparam_.ny*2*(cparam_.nz/2+1), 0.0 );
			compute_kernel();
		}
		
		void *get_ptr()
		{	return reinterpret_cast<void*> (&kdata_[0]);	}
		
		~kernel_k()
		{	std::vector<real_t>().swap( kdata_ );	}
		
	};
	
	template< typename real_t >
	void kernel_k<real_t>::compute_kernel( void )
	{
		double 
		//kny			= cparam_.nx*M_PI/cparam_.lx,
		fac			= cparam_.lx*cparam_.ly*cparam_.lz/pow(2.0*M_PI,3),// /(cparam_.nx*cparam_.ny*cparam_.nz),
																	   //	dx			= cparam_.lx/cparam_.nx, 
																	   //			dy			= cparam_.ly/cparam_.ny, 
																	   //			dz			= cparam_.lz/cparam_.nz,
			boxlength	= cparam_.pcf->getValue<double>("setup","boxlength"),
			nspec		= cparam_.pcf->getValue<double>("cosmology","nspec"),
			pnorm		= cparam_.pcf->getValue<double>("cosmology","pnorm"),
			dplus		= cparam_.pcf->getValue<double>("cosmology","dplus");
		
		TransferFunction_k *tfk = new TransferFunction_k(cparam_.ptf,nspec,pnorm,dplus);
		
		fftw_complex *kdata = reinterpret_cast<fftw_complex*> ( this->get_ptr() );
		
		unsigned nx = cparam_.nx, ny = cparam_.ny, nz = cparam_.nz, nzp = (nz/2+1);
		fac =1.0;//*= 1.0/sqrt(nx*ny*nz);
		
		double kfac = 2.0*M_PI/boxlength;
		unsigned q=0;
		for( int i=0; i<cparam_.nx; ++i )
			for( int j=0; j<cparam_.ny; ++j )
				for( int k=0; k<cparam_.nz/2+1; ++k )
				{
					double kx,ky,kz;
					
					kx = (double)i;
					ky = (double)j;
					kz = (double)k;
					
					if( kx > nx/2 ) kx -= nx;
					if( ky > ny/2 ) ky -= ny;
					
					q = (i*ny+j)*nzp+k;
					kdata[q].re = fac*tfk->compute(kfac*sqrt(kx*kx+ky*ky+kz*kz));
					kdata[q].im = 0.0;
					
				}
		
		
		delete tfk;
	}
};
	
namespace{
	convolution::kernel_creator_concrete< convolution::kernel_real<double> > creator_d("tf_kernel_real_double");
	convolution::kernel_creator_concrete< convolution::kernel_real<float> > creator_f("tf_kernel_real_float");
	convolution::kernel_creator_concrete< convolution::kernel_k<double> > creator_kd("tf_kernel_k_double");
	convolution::kernel_creator_concrete< convolution::kernel_k<float> > creator_kf("tf_kernel_k_float");

}



