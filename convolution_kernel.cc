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
					
					std::complex<double> ccdata(cdata[ii].re,cdata[ii].im), cckernel(ckernel[ii].re,ckernel[ii].im);
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

		void compute_kernel( tf_type type );
		
	public:
		kernel_real( const parameters& cp, tf_type type )
		: kernel( cp, type )
		{
			kdata_.assign( cparam_.nx*cparam_.ny*2*(cparam_.nz/2+1), 0.0 );
			compute_kernel( type );
		}
		
		void *get_ptr()
		{	return reinterpret_cast<void*> (&kdata_[0]);	}
		
		~kernel_real()
		{	std::vector<real_t>().swap( kdata_ );	}
		
	};
	
	
	template< typename real_t >
	void kernel_real<real_t>::compute_kernel( tf_type type )
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
		
		bool
			bavgfine	= cparam_.pcf->getValueSafe<bool>("setup","avg_fine",false),
			bdefd		= cparam_.pcf->getValueSafe<bool> ( "poisson" , "deconvolve", false ),
			bperiodic	= cparam_.pcf->getValueSafe<bool>("setup","periodic_TF",true);
		
		unsigned 
			levelmax	= cparam_.pcf->getValue<unsigned>("setup","levelmax");
		
			
		if( bavgfine )
			cparam_.coarse_fact++;
		
		const int ref_fac = (int)(pow(2,cparam_.coarse_fact)+0.5), ql = -ref_fac/2+1, qr=ql+ref_fac, rf8=(int)pow(ref_fac,3);
				
		std::cout << "   - Computing transfer function kernel...\n";
		
		if(! cparam_.is_finest )
			kny *= pow(2,cparam_.coarse_fact);
		
		TransferFunction_real *tfr = 
				new TransferFunction_real(type, cparam_.ptf,nspec,pnorm,dplus,
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
		
		if( cparam_.is_finest )
			kdata_[0] = tfr->compute_real(0.0)*fac;
				
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
							
							if(!(bdefd&&cparam_.is_finest))
							{
								kkernel[q].re *= ipix;
								kkernel[q].im *= ipix;		
							}else{
								kkmax /= pow(2,cparam_.coarse_fact);
								kkernel[q].re /= cos(kx*kkmax)*cos(ky*kkmax)*cos(kz*kkmax);
								kkernel[q].im /= cos(kx*kkmax)*cos(ky*kkmax)*cos(kz*kkmax);
							}
							
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
		
		if(bavgfine)
			cparam_.coarse_fact--;
	}
	



	template< typename real_t >
	class kernel_k : public kernel
	{
	protected:
		std::vector<real_t> kdata_;
		
		void compute_kernel( tf_type type );
		
	public:
		kernel_k( const parameters& cp, tf_type type )
		: kernel( cp, type )
		{
			kdata_.assign( cparam_.nx*cparam_.ny*2*(cparam_.nz/2+1), 0.0 );
			compute_kernel(type);
		}
		
		void *get_ptr()
		{	return reinterpret_cast<void*> (&kdata_[0]);	}
		
		~kernel_k()
		{	std::vector<real_t>().swap( kdata_ );	}
		
	};
	
	template< typename real_t >
	void kernel_k<real_t>::compute_kernel( tf_type type )
	{
		double 
			fac			= cparam_.lx*cparam_.ly*cparam_.lz/pow(2.0*M_PI,3),
			boxlength	= cparam_.pcf->getValue<double>("setup","boxlength"),
			nspec		= cparam_.pcf->getValue<double>("cosmology","nspec"),
			pnorm		= cparam_.pcf->getValue<double>("cosmology","pnorm"),
			dplus		= cparam_.pcf->getValue<double>("cosmology","dplus");
		
		TransferFunction_k *tfk = new TransferFunction_k(type,cparam_.ptf,nspec,pnorm,dplus);
		
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

	

/**************************************************************************************/
/**************************************************************************************/
#pragma mark -


template<int order>
double deconv_kernel( int idir, int i, int j, int k, int n )
{
	return 1.0;
}

template<>
inline double deconv_kernel<2>(int idir, int i, int j, int k, int n )
{
	double 
		ki(M_PI*(double)i/(double)n), 
		kj(M_PI*(double)j/(double)n), 
		kk(M_PI*(double)k/(double)n), 
		kr(sqrt(ki*ki+kj*kj+kk*kk));
	
	double grad = 1.0, laplace = 1.0;

	if( idir==0&&i!=0 )
		grad = sin(ki)/ki;
	else if( idir==1&&j!=0 )
		grad = sin(kj)/kj;
	else if( idir==2&&k!=0 )
		grad = sin(kk)/kk;
	
	laplace = 2.0*((-cos(ki)+1.0)+(-cos(kj)+1.0)+(-cos(kk)+1.0));
	
	if(i!=0&&j!=0&&k!=0)
		laplace /= (kr*kr);
	else
		laplace=1.0;

	return laplace/grad;
	
}

template<>
inline double deconv_kernel<4>(int idir, int i, int j, int k, int n )
{
	double 
	ki(M_PI*(double)i/(double)n), 
	kj(M_PI*(double)j/(double)n), 
	kk(M_PI*(double)k/(double)n), 
	kr(sqrt(ki*ki+kj*kj+kk*kk));
	
	double grad = 1.0, laplace = 1.0;
	
	if( idir==0&&i!=0 )
		grad = 0.166666666667*(-sin(2.*ki)+8.*sin(ki))/ki;
	else if( idir==1&&j!=0 )
		grad = 0.166666666667*(-sin(2.*kj)+8.*sin(kj))/kj;
	else if( idir==2&&k!=0 )
		grad = 0.166666666667*(-sin(2.*kk)+8.*sin(kk))/kk;
	
	laplace = 0.1666666667*((cos(2*ki)-16.*cos(ki)+15.)
							+(cos(2*kj)-16.*cos(kj)+15.)
							+(cos(2*kk)-16.*cos(kk)+15.));
	
	if(i!=0&&j!=0&&k!=0)
		laplace /= (kr*kr);
	else
		laplace=1.0;
	
	return laplace/grad;
	
}

template<>
inline double deconv_kernel<6>(int idir, int i, int j, int k, int n )
{
	double 
	ki(M_PI*(double)i/(double)n), 
	kj(M_PI*(double)j/(double)n), 
	kk(M_PI*(double)k/(double)n), 
	//ki(M_PI*(1.0-(double)(n-i)/(double)n)),
	//kj(M_PI*(1.0-(double)(n-j)/(double)n)),
	//kk(M_PI*(1.0-(double)(n-k)/(double)n)),
	
	kr(sqrt(ki*ki+kj*kj+kk*kk));
	
	if(i==0&&j==0&&k==0)
		return 0.0;
		
	double grad = 1.0, laplace = 1.0;
	
	if( idir==0 )//&&i!=0&&i!=n )
		grad = 0.0333333333333*(sin(3.*ki)-9.*sin(2.*ki)+45.*sin(ki));
	else if( idir==1 )//&&j!=0&&j!=n )
		grad = 0.0333333333333*(sin(3.*kj)-9.*sin(2.*kj)+45.*sin(kj));
	else if( idir==2 ) //&&k!=0&&k!=n )
		grad = 0.0333333333333*(sin(3.*kk)-9.*sin(2.*kk)+45.*sin(kk));
	
	laplace = 0.01111111111111*(
		(-2.*cos(3.0*ki)+27.*cos(2.*ki)-270.*cos(ki)+245.)
		+(-2.*cos(3.0*kj)+27.*cos(2.*kj)-270.*cos(kj)+245.)
		+(-2.*cos(3.0*kk)+27.*cos(2.*kk)-270.*cos(kk)+245.));
	
	//if(i!=0&&j!=0&&k!=0)
/*	if(i==0&&j==0&&k==0)
		laplace=1.0;
	else*/
		//laplace -= (kr*kr);
	
		
	double kgrad = 1.0;
	if( idir==0 )
		kgrad = ki;
	else if( idir ==1)
		kgrad = kj;
	else if( idir ==2)
		kgrad = kk;
	
	return (kgrad/kr/kr-grad/laplace)*M_PI/n*M_PI/n;//laplace/grad;
	//return kgrad/kr/kr*M_PI/n*M_PI/n;
	//return grad/laplace*M_PI/n*M_PI/n;
	//return (kgrad/kr/kr-grad/laplace);
}


template<int order>
void do_deconvolve( fftw_real* data, int idir, int nxp, int nyp, int nzp, bool periodic )
{
	double fftnorm = 1.0/(nxp*nyp*nzp);
	
	
	if(periodic)
		fftnorm *= nxp/(4.0*M_PI*M_PI);//sqrt(8.0);
	
	
	/*if(periodic)
		fftnorm *= M_PI*M_PI/nxp/nyp/nzp;
	else
		fftnorm *= M_PI*M_PI/nxp/nyp/nzp;
	
	if( idir == 0 )
		fftnorm *= nxp;
	else if( idir == 1 )
		fftnorm *= nyp;
	else
		fftnorm *= nzp;*/
		
	fftw_complex	*cdata = reinterpret_cast<fftw_complex*>(data);
	rfftwnd_plan	iplan, plan;
	
	plan  = rfftw3d_create_plan( nxp, nyp, nzp,
								FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE);
	
	iplan = rfftw3d_create_plan( nxp, nyp, nzp, 
								FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
	
	#ifndef SINGLETHREAD_FFTW		
	rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), plan, data, NULL );
#else
	rfftwnd_one_real_to_complex( plan, data, NULL );
#endif
	
	double ksum = 0.0;
	unsigned kcount = 0;
	
	for( int i=0; i<nxp; ++i )
		for( int j=0; j<nyp; ++j )
			for( int k=0; k<nzp/2+1; ++k )
			{
				unsigned ii = (i*nyp + j) * (nzp/2+1) + k;
				
				if( k==0 || k==nzp/2 )
				{
					ksum  += cdata[ii].re;
					kcount++;
				}else{
					ksum  += 2.0*(cdata[ii].re);
					kcount+=2;
				}
			}
	
	ksum /= kcount;
	kcount = 0;
	
	cdata[0].re = 0.0;
	cdata[0].im = 0.0;
	
#pragma omp parallel for
	for( int i=0; i<nxp; ++i )
		for( int j=0; j<nyp; ++j )
			for( int k=0; k<nzp/2+1; ++k )
			{
				unsigned ii = (i*nyp + j) * (nzp/2+1) + k;
				
				int ki(i), kj(j);
				if( ki > nxp/2 ) ki-=nxp;
				if( kj > nyp/2 ) kj-=nyp;
				
				double dk = deconv_kernel<order>(idir, ki, kj, k, nxp/2 );
				//cdata[ii].re -= ksum;
				
				double re = cdata[ii].re, im = cdata[ii].im;
				
				cdata[ii].re = -im*dk*fftnorm;
				cdata[ii].im = re*dk*fftnorm;
				
				//cdata[ii].re += ksum*fftnorm;
				
			}
	cdata[0].re = 0.0;
	cdata[0].im = 0.0;
	
#ifndef SINGLETHREAD_FFTW		
	rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), iplan, cdata, NULL);
#else		
	rfftwnd_one_complex_to_real(iplan, cdata, NULL);
#endif
	
	rfftwnd_destroy_plan(plan);
	rfftwnd_destroy_plan(iplan);
	
}
void deconvolve( MeshvarBnd<double>& f, int idir, int order, bool periodic )
{
	int nx=f.size(0), ny=f.size(1), nz=f.size(2), nxp, nyp, nzp;
	fftw_real		*data;
	int xo=0,yo=0,zo=0;
	int nmax = std::max(nx,std::max(ny,nz));
	
	if(!periodic)
	{
		nxp = 2*nmax;
		nyp = 2*nmax;
		nzp = 2*nmax;
		xo  = nmax/2;
		yo  = nmax/2;
		zo  = nmax/2;
	}
	else
	{
		nxp = nmax;
		nyp = nmax;
		nzp = nmax;
	}
	

	
	data		= new fftw_real[nxp*nyp*2*(nzp/2+1)];
	
	if(idir==0)
		std::cout << "   - Performing de-convolution... (" << nxp <<  ", " << nyp << ", " << nzp << ")\n";
	
	
	for( int i=0; i<nxp*nyp*2*(nzp/2+1); ++i )
		data[i]=0.0;
	
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )
			for( int k=0; k<nz; ++k )
			{
				int idx = ((i+xo)*nyp + j+yo) * 2*(nzp/2+1) + k+zo;
				data[idx] = f(i,j,k);
			}
	
	switch (order) {
		case 2:
			do_deconvolve<2>(data, idir, nxp, nyp, nzp, periodic);
			break;
		case 4:
			do_deconvolve<4>(data, idir, nxp, nyp, nzp, periodic);
			break;
		case 6:
			do_deconvolve<6>(data, idir, nxp, nyp, nzp, periodic);
			break;
		default:
			std::cerr << " - ERROR: invalid operator order specified in deconvolution.";
			break;
	}
	
	
	for( int i=0; i<nx; ++i )
		for( int j=0; j<ny; ++j )
			for( int k=0; k<nz; ++k )
			{
				//int idx = (i*nyp + j) * 2*(nzp/2+1) + k;
				int idx = ((i+xo)*nyp + j+yo) * 2*(nzp/2+1) + k+zo;
				
				f(i,j,k) = data[idx];
			}
	
	delete[] data;
	
}








namespace{
	convolution::kernel_creator_concrete< convolution::kernel_real<double> > creator_d("tf_kernel_real_double");
	convolution::kernel_creator_concrete< convolution::kernel_real<float> > creator_f("tf_kernel_real_float");
	convolution::kernel_creator_concrete< convolution::kernel_k<double> > creator_kd("tf_kernel_k_double");
	convolution::kernel_creator_concrete< convolution::kernel_k<float> > creator_kf("tf_kernel_k_float");

}



