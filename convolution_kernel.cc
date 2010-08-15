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

double T0 = 1.0;

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
		
		std::cout << "   - Performing density convolution... (" 
			<< cparam_.nx <<  ", " << cparam_.ny << ", " << cparam_.nz << ")\n";
		
		
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
	
	
	template void perform<double>( kernel* pk, void *pd );
	template void perform<float>( kernel* pk, void *pd );
	
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
			nspec		= cparam_.pcf->getValue<double>("cosmology","nspec"),
			pnorm		= cparam_.pcf->getValue<double>("cosmology","pnorm"),
			dplus		= cparam_.pcf->getValue<double>("cosmology","dplus");

		unsigned 
			levelmax	= cparam_.pcf->getValue<unsigned>("setup","levelmax");
		
		bool
			bperiodic	= cparam_.pcf->getValueSafe<bool>("setup","periodic_TF",true);
		
		const int ref_fac = (int)(pow(2,cparam_.coarse_fact)+0.5), ql = -ref_fac/2+1, qr=ql+ref_fac, rf8=(int)pow(ref_fac,3);
				
		std::cout << "   - Computing transfer function kernel...\n";
		
		if(! cparam_.is_finest )
			kny *= pow(2,cparam_.coarse_fact);
		
		TransferFunction_real *tfr = 
				new TransferFunction_real(cparam_.ptf,nspec,pnorm,dplus,
							0.25*cparam_.lx/(double)cparam_.nx,2.0*boxlength,kny, (int)pow(2,levelmax+2));
		
		fftw_real		*rkernel	= reinterpret_cast<fftw_real*>( &kdata_[0] );
		fftw_complex	*kkernel	= reinterpret_cast<fftw_complex*>( &kdata_[0] );
		
		rfftwnd_plan	plan		= rfftw3d_create_plan( cparam_.nx, cparam_.ny, cparam_.nz,
												FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE);
						
		int nzp = cparam_.nz/2+1;
		double kmax = 0.5*M_PI/std::max(cparam_.nx,std::max(cparam_.ny,cparam_.nz));
		
		T0 = tfr->compute_real(0.0);
		
		//... obtain a grid-version of the real space transfer function
		if( bperiodic  )
		{		
			#pragma omp parallel for //schedule(static,4)
			for( int i=0; i<=cparam_.nx/2; ++i )
				for( int j=0; j<=cparam_.ny/2; ++j )
					for( int k=0; k<=cparam_.nz/2; ++k )
					{
						int iix(i), iiy(j), iiz(k);
						double rr[3], rr2;
						
						if( iix > (int)cparam_.nx/2 ) iix -= cparam_.nx;
						if( iiy > (int)cparam_.ny/2 ) iiy -= cparam_.ny;
						if( iiz > (int)cparam_.nz/2 ) iiz -= cparam_.nz;
						
						
						//... speed up 8x by copying data to other octants
						int idx[8];
						int nx=cparam_.nx, ny = cparam_.ny, nz = cparam_.nz;
						
						idx[0] = ((i)*cparam_.ny + (j)) * 2*(cparam_.nz/2+1) + (k);
						idx[1] = ((nx-i)*cparam_.ny + (j)) * 2*(cparam_.nz/2+1) + (k);
						idx[2] = ((i)*cparam_.ny + (ny-j)) * 2*(cparam_.nz/2+1) + (k);
						idx[3] = ((nx-i)*cparam_.ny + (ny-j)) * 2*(cparam_.nz/2+1) + (k);
						idx[4] = ((i)*cparam_.ny + (j)) * 2*(cparam_.nz/2+1) + (nz-k);
						idx[5] = ((nx-i)*cparam_.ny + (j)) * 2*(cparam_.nz/2+1) + (nz-k);
						idx[6] = ((i)*cparam_.ny + (ny-j)) * 2*(cparam_.nz/2+1) + (nz-k);
						idx[7] = ((nx-i)*cparam_.ny + (ny-j)) * 2*(cparam_.nz/2+1) + (nz-k);
						
						if(i==0||i==nx/2){ idx[1]=idx[3]=idx[5]=idx[7]=-1;}
						if(j==0||j==ny/2){ idx[2]=idx[3]=idx[6]=idx[7]=-1;}
						if(k==0||k==nz/2){ idx[4]=idx[5]=idx[6]=idx[7]=-1;}
						
						double val = 0.0;
						
						for( int ii=-1; ii<=1; ++ii )
							for( int jj=-1; jj<=1; ++jj )
								for( int kk=-1; kk<=1; ++kk )
								{
									rr[0] = ((double)iix ) * dx + ii*boxlength;
									rr[1] = ((double)iiy ) * dy + jj*boxlength;
									rr[2] = ((double)iiz ) * dz + kk*boxlength;
									
									if( rr[0] > -boxlength && rr[0] < boxlength
									 && rr[1] > -boxlength && rr[1] < boxlength
									 && rr[2] > -boxlength && rr[2] < boxlength )
									{
										if( ref_fac > 1 )//&& fabs(tfr->get_grad(rr2)*dx/T0) > 1e-4 )
										{
											double rrr[3];
											register double rrr2[3];
											for( int iii=ql; iii<qr; ++iii )
											{	
												rrr[0] = rr[0]+(double)iii*0.5*dx - 0.25*dx;
												rrr2[0]= rrr[0]*rrr[0];
												for( int jjj=ql; jjj<qr; ++jjj )
												{
													rrr[1] = rr[1]+(double)jjj*0.5*dx - 0.25*dx;
													rrr2[1]= rrr[1]*rrr[1];
													for( int kkk=ql; kkk<qr; ++kkk )
													{
														rrr[2] = rr[2]+(double)kkk*0.5*dx - 0.25*dx;
														rrr2[2]= rrr[2]*rrr[2];
														rr2 = rrr2[0]+rrr2[1]+rrr2[2];
														val += tfr->compute_real(rr2)/rf8;
													}
												}
											}
											
										}else{
											rr2 = rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2];
											val += tfr->compute_real(rr2);	
										}
									}
								}
						
						val *= fac;
						
						for(int q=0;q<8;++q)
							if(idx[q]!=-1)  
								kdata_[idx[q]] = val;
						
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
						
						if( ref_fac > 1 )
						{
							double rrr[3];
							register double rrr2[3];
							for( int iii=ql; iii<qr; ++iii )
							{	
								rrr[0] = rr[0]+(double)iii*0.5*dx - 0.25*dx;
								rrr2[0]= rrr[0]*rrr[0];
								for( int jjj=ql; jjj<qr; ++jjj )
								{
									rrr[1] = rr[1]+(double)jjj*0.5*dx - 0.25*dx;
									rrr2[1]= rrr[1]*rrr[1];
									for( int kkk=ql; kkk<qr; ++kkk )
									{
										rrr[2] = rr[2]+(double)kkk*0.5*dx - 0.25*dx;
										rrr2[2]= rrr[2]*rrr[2];
										rr2 = rrr2[0]+rrr2[1]+rrr2[2];
										kdata_[idx] += (fftw_real)(tfr->compute_real(rr2)/rf8);	
									}
								}
							}
							
						}else{
							rr2 = rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2];
							kdata_[idx] += (fftw_real)tfr->compute_real(rr2);
						}
					}
		}
				
		T0 = kdata_[0];
		delete tfr;
		
		double k0 = kdata_[0];

		//... subtract white noise component before deconvolution
		if( cparam_.deconvolve )//&& cparam_.is_finest)
			kdata_[0] = 0.0;
		
		#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), plan, rkernel, NULL );
		#else
		rfftwnd_one_real_to_complex( plan, rkernel, NULL );
		#endif
		
		if( cparam_.deconvolve || cparam_.smooth )
		{
			
			double ksum = 0.0;
			unsigned kcount = 0;
			
			#pragma omp parallel for reduction(+:ksum,kcount)
			for( int i=0; i<cparam_.nx; ++i )
				for( int j=0; j<cparam_.ny; ++j )
					for( int k=0; k<nzp; ++k )
					{
						double kx,ky,kz;
						
						kx = (double)i;
						ky = (double)j;
						kz = (double)k;
						
						if( kx > cparam_.nx/2 ) kx -= cparam_.nx;
						if( ky > cparam_.ny/2 ) ky -= cparam_.ny;
						

						double ipix = 1.0;
						
						//... perform k-space 'averaging' for coarser grids
						/*if( !cparam_.is_finest )
						{
							for( unsigned c=1; c<= cparam_.coarse_fact; ++c )
							{
								double kkmax = 0.5*kmax/pow(2,c-1);
								ipix *= (cos(kx*kkmax)*cos(ky*kkmax)*cos(kz*kkmax));								
							}
						}*/
						
						double kkmax = kmax;// / pow(2,cparam_.coarse_fact);
							
						if( true )//cparam_.is_finest||cparam_.smooth )
						{
							
							//... deconvolve with grid-cell (NGP) kernel
							if( i > 0 )
								ipix /= sin(kx*2.0*kkmax)/(kx*2.0*kkmax);
							if( j > 0 )
								ipix /= sin(ky*2.0*kkmax)/(ky*2.0*kkmax);
							if( k > 0 )
								ipix /= sin(kz*2.0*kkmax)/(kz*2.0*kkmax);
						}
						
						unsigned q  = (i*cparam_.ny+j)*nzp+k;
						
						if( !cparam_.smooth )
						{
							kkernel[q].re *= ipix;
							kkernel[q].im *= ipix;	
						}else{
							//... if smooth==true, convolve with 
							//... NGP kernel to get CIC smoothness
							kkernel[q].re /= ipix;
							kkernel[q].im /= ipix;	
						}
						
						//... store k-space average
						if( k==0 || k==cparam_.nz/2 )
						{
							ksum  += kkernel[q].re;
							kcount++;
						}else{
							ksum  += 2.0*(kkernel[q].re);
							kcount+=2;
						}
					}
			
			double dk;

			//... re-add white noise component for finest grid
			if( cparam_.is_finest )
				dk = k0-ksum/kcount;
			else
				dk = k0-ksum/kcount;
			
			//... set white noise component to zero if smoothing is enabled
			if( cparam_.smooth )
				dk = 0.0;
			
			//... enforce the r=0 component by adjusting the k-space mean
			#pragma omp parallel for reduction(+:ksum,kcount)
			for( int i=0; i<cparam_.nx; ++i )
				for( int j=0; j<cparam_.ny; ++j )
					for( int k=0; k<nzp; ++k )
					{
						unsigned q  = (i*cparam_.ny+j)*nzp+k;
						kkernel[q].re += dk;
					}
		}		
		
		rfftwnd_destroy_plan(plan);
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
			fac			= cparam_.lx*cparam_.ly*cparam_.lz/pow(2.0*M_PI,3),
			boxlength	= cparam_.pcf->getValue<double>("setup","boxlength"),
			nspec		= cparam_.pcf->getValue<double>("cosmology","nspec"),
			pnorm		= cparam_.pcf->getValue<double>("cosmology","pnorm"),
			dplus		= cparam_.pcf->getValue<double>("cosmology","dplus");
		
		TransferFunction_k *tfk = new TransferFunction_k(cparam_.ptf,nspec,pnorm,dplus);
		
		fftw_complex *kdata = reinterpret_cast<fftw_complex*> ( this->get_ptr() );
		
		unsigned nx = cparam_.nx, ny = cparam_.ny, nz = cparam_.nz, nzp = (nz/2+1);
		fac =1.0;
		
		double kfac = 2.0*M_PI/boxlength, ksum = 0.0;
		unsigned q=0, kcount = 0;
		
		#pragma omp parallel for reduction(+:ksum,kcount)
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
					
					if( k==0 || k==cparam_.nz/2 )
					{
						ksum  += kdata[q].re;
						kcount++;
					}else{
						ksum  += 2.0*(kdata[q].re);
						kcount+=2;
					}
					
				}
		
		delete tfk;
	}
}

	
namespace{
	convolution::kernel_creator_concrete< convolution::kernel_real<double> > creator_d("tf_kernel_real_double");
	convolution::kernel_creator_concrete< convolution::kernel_real<float> > creator_f("tf_kernel_real_float");
	convolution::kernel_creator_concrete< convolution::kernel_k<double> > creator_kd("tf_kernel_k_double");
	convolution::kernel_creator_concrete< convolution::kernel_k<float> > creator_kf("tf_kernel_k_float");

}



