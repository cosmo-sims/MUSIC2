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

#include "general.hh"
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
		//return;
		parameters cparam_ = pk->cparam_;
		double fftnorm = pow(2.0*M_PI,1.5)/sqrt(cparam_.lx*cparam_.ly*cparam_.lz)/sqrt((double)(cparam_.nx*cparam_.ny*cparam_.nz));
		
		fftw_complex	*cdata,*ckernel;
		fftw_real		*data;
		
		data		= reinterpret_cast<fftw_real*>(pd);
		cdata		= reinterpret_cast<fftw_complex*>(data);
		ckernel		= reinterpret_cast<fftw_complex*>( pk->get_ptr() );
		
		
		
		std::cout << "   - Performing density convolution... (" 
			<< cparam_.nx <<  ", " << cparam_.ny << ", " << cparam_.nz << ")\n";
		
		
		LOGUSER("Performing kernel convolution on (%5d,%5d,%5d) grid",cparam_.nx ,cparam_.ny ,cparam_.nz );
		LOGUSER("Performing forward FFT...");
#ifdef FFTW3
		fftw_plan plan, iplan;
		plan = fftw_plan_dft_r2c_3d( cparam_.nx, cparam_.ny, cparam_.nz, data, cdata, FFTW_ESTIMATE);
		iplan = fftw_plan_dft_c2r_3d(cparam_.nx, cparam_.ny, cparam_.nz, cdata, data, FFTW_ESTIMATE);
		
		fftw_execute(plan);
#else
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
		
#endif
		
		#pragma omp parallel for
		for( int i=0; i<cparam_.nx; ++i )
			for( int j=0; j<cparam_.ny; ++j )
				for( int k=0; k<cparam_.nz/2+1; ++k )
				{
					unsigned ii = (i*cparam_.ny + j) * (cparam_.nz/2+1) + k;
#ifdef FFTW3
					std::complex<double> ccdata(cdata[ii][0],cdata[ii][1]), cckernel(ckernel[ii][0],ckernel[ii][1]);
					ccdata = ccdata * cckernel *fftnorm;
					
					cdata[ii][0] = ccdata.real();
					cdata[ii][1] = ccdata.imag();
#else
					std::complex<double> ccdata(cdata[ii].re,cdata[ii].im), cckernel(ckernel[ii].re,ckernel[ii].im);
					ccdata = ccdata * cckernel *fftnorm;
					
					cdata[ii].re = ccdata.real();
					cdata[ii].im = ccdata.imag();
#endif
				}

		LOGUSER("Performing backward FFT...");
		
#ifdef FFTW3
		fftw_execute(iplan);
		fftw_destroy_plan(plan);
		fftw_destroy_plan(iplan);
		
#else
		#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), iplan, cdata, NULL);
		#else		
		rfftwnd_one_complex_to_real(iplan, cdata, NULL);
		#endif
		
		rfftwnd_destroy_plan(plan);
		rfftwnd_destroy_plan(iplan);
#endif
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
		kernel_real( config_file& cf, transfer_function* ptf, refinement_hierarchy& refh, tf_type type )
		: kernel( cf, ptf, refh, type )
		{	}
		
		void *get_ptr()
		{	return reinterpret_cast<void*> (&kdata_[0]);	}
		
		kernel* fetch_kernel( int ilevel, bool isolated=false );
		
		~kernel_real()
		{	deallocate();	}
		
		void deallocate()
		{	std::vector<real_t>().swap( kdata_ );	}
		

		
	};
	
	
	
	
	template< typename real_t >
	kernel* kernel_real<real_t>::fetch_kernel( int ilevel, bool isolated )
	{
		double
			boxlength	= pcf_->getValue<double>("setup","boxlength"),
			boxlength2  = 0.5*boxlength;
		
		int
			levelmax	= prefh_->levelmax(),
			levelmin	= prefh_->levelmin(),
			nx = prefh_->size(ilevel,0),
			ny = prefh_->size(ilevel,1),
			nz = prefh_->size(ilevel,2);
		
		/*if( isolated )
		{	nx *= 2; ny *= 2; nz *= 2;	}*/
		if( ilevel != levelmin )
		{	nx*=2;	ny*=2;	nz*=2; }
		
		kdata_.assign( nx*ny*2*(nz/2+1), 0.0 );
				
		double
			dx = boxlength / (1<<ilevel),
			lx = dx * nx,
			ly = dx * ny,
			lz = dx * nz;
		
		double 
			kny			= nx*M_PI/lx,
			fac			= lx*ly*lz/pow(2.0*M_PI,3)/(nx*ny*nz),
			nspec		= pcf_->getValue<double>("cosmology","nspec"),
			pnorm		= pcf_->getValue<double>("cosmology","pnorm"),
			dplus		= pcf_->getValue<double>("cosmology","dplus");
		
		bool
			bavgfine	= pcf_->getValueSafe<bool>("setup","avg_fine",true),
			bperiodic	= pcf_->getValueSafe<bool>("setup","periodic_TF",true),
			kspacepoisson = (pcf_->getValueSafe<bool>("poisson","fft_fine",true)
							 |pcf_->getValueSafe<bool>("poisson","kspace",false));
		
		cparam_.nx = nx;
		cparam_.ny = ny;
		cparam_.nz = nz;
		cparam_.lx = dx * cparam_.nx;
		cparam_.ly = dx * cparam_.ny;
		cparam_.lz = dx * cparam_.nz;
		cparam_.pcf = pcf_;
		
		if( pcf_->getValueSafe<bool>("setup","deconvolve",true) )
			cparam_.deconvolve = true;
		else
			cparam_.deconvolve = false;
		
		cparam_.coarse_fact = (prefh_->levelmax()-ilevel);
				
		if( bavgfine )
			cparam_.coarse_fact++;

		const int ref_fac = 1<<cparam_.coarse_fact, ql = -ref_fac/2+1, qr=ql+ref_fac, rf8=(int)pow(ref_fac,3);
				
		std::cout << "   - Computing transfer function kernel...\n";
		
		if(! cparam_.is_finest )
				kny *= pow(2,cparam_.coarse_fact);

		
		TransferFunction_real *tfr = 
				new TransferFunction_real( boxlength, 1<<levelmax, type_, ptf_,nspec,pnorm,dplus,
							0.25*lx/(double)nx,2.0*boxlength,kny, (int)pow(2,levelmax+2));
		
		fftw_real		*rkernel	= reinterpret_cast<fftw_real*>( &kdata_[0] );
		fftw_complex	*kkernel	= reinterpret_cast<fftw_complex*>( &kdata_[0] );

#ifdef FFTW3
		fftw_plan		plan		= fftw_plan_dft_r2c_3d( cparam_.nx, cparam_.ny, cparam_.nz, rkernel, kkernel, FFTW_ESTIMATE);
#else
		rfftwnd_plan	plan		= rfftw3d_create_plan( cparam_.nx, cparam_.ny, cparam_.nz,
												FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE);
#endif
		int nzp = cparam_.nz/2+1;
		double kmax = 0.5*M_PI/std::max(cparam_.nx,std::max(cparam_.ny,cparam_.nz));
		
		T0 = tfr->compute_real(0.0);
		

		
		if(bavgfine)
			cparam_.coarse_fact--;
				
		//... obtain a grid-version of the real space transfer function
		if( bperiodic  )
		{		
			#pragma omp parallel for //schedule(static,4)
			for( int i=0; i<=nx/2; ++i )
				for( int j=0; j<=ny/2; ++j )
					for( int k=0; k<=nz/2; ++k )
					{
						int iix(i), iiy(j), iiz(k);
						double rr[3], rr2;
						
						if( iix > (int)nx/2 ) iix -= nx;
						if( iiy > (int)ny/2 ) iiy -= ny;
						if( iiz > (int)nz/2 ) iiz -= nz;
						
						
						//... speed up 8x by copying data to other octants
						size_t idx[8];
						//int nx=cparam_.nx, ny = cparam_.ny, nz = cparam_.nz;
						
						idx[0] = ((size_t)(i)*ny + (size_t)(j)) * 2*(size_t)(nz/2+1) + (size_t)(k);
						idx[1] = ((size_t)(nx-i)*ny + (size_t)(j)) * 2*(nz/2+1) + (size_t)(k);
						idx[2] = ((size_t)(i)*ny + (size_t)(ny-j)) * 2*(nz/2+1) + (size_t)(k);
						idx[3] = ((size_t)(nx-i)*ny + (size_t)(ny-j)) * 2*(nz/2+1) + (size_t)(k);
						idx[4] = ((size_t)(i)*ny + (size_t)(j)) * 2*(nz/2+1) + (size_t)(nz-k);
						idx[5] = ((size_t)(nx-i)*ny + (size_t)(j)) * 2*(nz/2+1) + (size_t)(nz-k);
						idx[6] = ((size_t)(i)*ny + (size_t)(ny-j)) * 2*(nz/2+1) + (size_t)(nz-k);
						idx[7] = ((size_t)(nx-i)*ny + (size_t)(ny-j)) * 2*(nz/2+1) + (size_t)(nz-k);
						
						if(i==0||i==nx/2){ idx[1]=idx[3]=idx[5]=idx[7]=-1;}
						if(j==0||j==ny/2){ idx[2]=idx[3]=idx[6]=idx[7]=-1;}
						if(k==0||k==nz/2){ idx[4]=idx[5]=idx[6]=idx[7]=-1;}
						
						double val = 0.0;
						
						for( int ii=-1; ii<=1; ++ii )
							for( int jj=-1; jj<=1; ++jj )
								for( int kk=-1; kk<=1; ++kk )
								{
									rr[0] = ((double)iix ) * dx + ii*boxlength;
									rr[1] = ((double)iiy ) * dx + jj*boxlength;
									rr[2] = ((double)iiz ) * dx + kk*boxlength;
									
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
			for( int i=0; i<nx; ++i )
				for( int j=0; j<ny; ++j )
					for( int k=0; k<nz; ++k )
					{
						int iix(i), iiy(j), iiz(k);
						double rr[3], rr2;
						
						if( iix > (int)nx/2 ) iix -= nx;
						if( iiy > (int)ny/2 ) iiy -= ny;
						if( iiz > (int)nz/2 ) iiz -= nz;
						
						size_t idx = ((size_t)i*ny + (size_t)j) * 2*(size_t)(nz/2+1) + (size_t)k;
						
						rr[0] = ((double)iix ) * dx;
						rr[1] = ((double)iiy ) * dx;
						rr[2] = ((double)iiz ) * dx;
						
						kdata_[idx] = 0.0;
						
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
										if( fabs(rr[0])<=boxlength2||fabs(rr[1])<=boxlength2||fabs(rr[2])<=boxlength2 )
											kdata_[idx] += (fftw_real)(tfr->compute_real(rr2)/rf8)*fac;	
									}
								}
							}
							if( false )//cparam_.deconvolve )//&& cparam_.is_finest  )
							{
								double rx((double)iix*dx/boxlength),ry((double)iiy*dx/boxlength),rz((double)iiz*dx/boxlength), ico;
								ico = 1.0;
								if( i>0 && fabs(dx*(double)iix) < 0.5*boxlength)//&& i<cparam_.nx/2)
									ico *= cos(M_PI*rx);//sin(M_PI*rx)/(M_PI*rx);
								if( j>0 && fabs(dx*(double)iiy) < 0.5*boxlength)//&& j<cparam_.ny/2 )
									ico *= cos(M_PI*ry);//sin(M_PI*ry)/(M_PI*ry);
								if( k>0 && fabs(dx*(double)iiz) < 0.5*boxlength)//&& k<cparam_.nz/2)
									ico *= cos(M_PI*rz);//sin(M_PI*rz)/(M_PI*rz);
								//kdata_[idx] /= ico;
								//if( fabs(ico)>1e-8 );
								kdata_[idx] /= ico;
								
							}
						}else{
							rr2 = rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2];
							if( fabs(rr[0])<=boxlength2||fabs(rr[1])<=boxlength2||fabs(rr[2])<=boxlength2 )
								kdata_[idx] += (fftw_real)tfr->compute_real(rr2)*fac;
							
						}
					}
		}
		
		
		if( ilevel == levelmax )
			kdata_[0] = tfr->compute_real(0.0)*fac;

		T0 = kdata_[0];
		delete tfr;
		
		double k0 = kdata_[0];
		
		//... subtract white noise component before deconvolution
		if( cparam_.deconvolve )//&& cparam_.is_finest)
			kdata_[0] = 0.0;
		
#ifdef FFTW3
		fftw_execute(plan);
#else
		#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), plan, rkernel, NULL );
		#else
		rfftwnd_one_real_to_complex( plan, rkernel, NULL );
		#endif
#endif
/*#define XI_SAMPLE
#ifdef XI_SAMPLE
		#pragma omp parallel for
		for( int i=0; i<cparam_.nx; ++i )
			for( int j=0; j<cparam_.ny; ++j )
				for( int k=0; k<nzp; ++k )	
				{
					unsigned q  = (i*cparam_.ny+j)*nzp+k;
					std::complex<double> z(kkernel[q].re,kkernel[q].im);
					z=sqrt(z);
					kkernel[q].re = z.real();
					kkernel[q].im = z.imag();
				}
		
#endif*/
		
		
		
		//... do some k-space filter correction stuff
		if( cparam_.deconvolve || cparam_.smooth )
		{
			
			double ksum = 0.0;
			size_t kcount = 0;
			
			#pragma omp parallel for reduction(+:ksum,kcount)
			for( int i=0; i<nx; ++i )
				for( int j=0; j<ny; ++j )
					for( int k=0; k<nzp; ++k )
					{
						double kx,ky,kz;
						
						kx = (double)i;
						ky = (double)j;
						kz = (double)k;
						
						if( kx > nx/2 ) kx -= nx;
						if( ky > ny/2 ) ky -= ny;
						

						double ipix = 1.0;
						
						
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
						
						size_t q  = ((size_t)i*ny+(size_t)j)*nzp+(size_t)k;
						
						if( !cparam_.smooth )
						{
							if( cparam_.is_finest )
							{
								if( kspacepoisson )
								{
#ifdef FFTW3
									kkernel[q][0] /= cos(kx*kkmax)*cos(ky*kkmax)*cos(kz*kkmax);
									kkernel[q][1] /= cos(kx*kkmax)*cos(ky*kkmax)*cos(kz*kkmax);
#else
									kkernel[q].re /= cos(kx*kkmax)*cos(ky*kkmax)*cos(kz*kkmax);
									kkernel[q].im /= cos(kx*kkmax)*cos(ky*kkmax)*cos(kz*kkmax);
#endif
									///////kkernel[q].re *= ipix;
									///////kkernel[q].im *= ipix;	
									//kkernel[q].re *= ipix;
									//kkernel[q].im *= ipix;
								}else{
									
									////////kkernel[q].re /= cos(kx*kkmax)*cos(ky*kkmax)*cos(kz*kkmax);
									////////kkernel[q].im /= cos(kx*kkmax)*cos(ky*kkmax)*cos(kz*kkmax);
#ifdef FFTW3
									kkernel[q][0] *= ipix;
									kkernel[q][1] *= ipix;	
#else
									kkernel[q].re *= ipix;
									kkernel[q].im *= ipix;	
#endif
								}
								
							}else{
#ifdef FFTW3
								kkernel[q][0] /= cos(kx*kkmax)*cos(ky*kkmax)*cos(kz*kkmax);
								kkernel[q][1] /= cos(kx*kkmax)*cos(ky*kkmax)*cos(kz*kkmax);
#else
								kkernel[q].re /= cos(kx*kkmax)*cos(ky*kkmax)*cos(kz*kkmax);
								kkernel[q].im /= cos(kx*kkmax)*cos(ky*kkmax)*cos(kz*kkmax);
#endif
								////////kkernel[q].re *= ipix;
								////////kkernel[q].im *= ipix;
							}
														
						}else{
							//... if smooth==true, convolve with 
							//... NGP kernel to get CIC smoothness
#ifdef FFTW3
							kkernel[q][0] /= ipix;
							kkernel[q][1] /= ipix;	
#else
							kkernel[q].re /= ipix;
							kkernel[q].im /= ipix;
#endif
						}
						
						//... store k-space average
#ifdef FFTW3
						if( k==0 || k==nz/2 )
						{
							ksum  += kkernel[q][0];
							kcount++;
						}else{
							ksum  += 2.0*(kkernel[q][0]);
							kcount+=2;
						}
#else
						if( k==0 || k==nz/2 )
						{
							ksum  += kkernel[q].re;
							kcount++;
						}else{
							ksum  += 2.0*(kkernel[q].re);
							kcount+=2;
						}
#endif
						
					}
			
			double dk;

			//... re-add white noise component for finest grid
			if( ilevel == levelmax )
				dk = k0-ksum/kcount;
			else
				dk = k0-ksum/kcount;
			
			//... set white noise component to zero if smoothing is enabled
			if( cparam_.smooth )
				dk = 0.0;
			
			//dk = 0.0;

			//... enforce the r=0 component by adjusting the k-space mean
			#pragma omp parallel for reduction(+:ksum,kcount)
			for( int i=0; i<nx; ++i )
				for( int j=0; j<ny; ++j )
					for( int k=0; k<nzp; ++k )
					{
						size_t q  = ((size_t)i*ny+(size_t)j)*nzp+(size_t)k;
#ifdef FFTW3
						kkernel[q][0] += dk;
#else
						kkernel[q].re += dk;
#endif
					}
		}		
		
#ifdef FFTW3
		fftw_destroy_plan(plan);
#else
		rfftwnd_destroy_plan(plan);
#endif
		
		return this;
	}
	



	template< typename real_t >
	class kernel_k : public kernel
	{
	protected:
		std::vector<real_t> kdata_;
		
		void compute_kernel( tf_type type );
		
	public:
		kernel_k( config_file& cf, transfer_function* ptf, refinement_hierarchy& refh, tf_type type )
		: kernel( cf, ptf, refh, type )
		{	}
		
		kernel* fetch_kernel( int ilevel, bool isolated=false );
		
		void *get_ptr()
		{	return reinterpret_cast<void*> (&kdata_[0]);	}
		
		~kernel_k()
		{	deallocate();	}
		
		void deallocate()
		{	std::vector<real_t>().swap( kdata_ );	}
		
	};
	
	template< typename real_t >
	kernel* kernel_k<real_t>::fetch_kernel( int ilevel, bool isolated )
	{
		double 
			boxlength	= pcf_->getValue<double>("setup","boxlength"),
			nspec		= pcf_->getValue<double>("cosmology","nspec"),
			pnorm		= pcf_->getValue<double>("cosmology","pnorm"),
			dplus		= pcf_->getValue<double>("cosmology","dplus"),
			fac			= pow(boxlength,3)/pow(2.0*M_PI,3);
		
		TransferFunction_k *tfk = new TransferFunction_k(type_,ptf_,nspec,pnorm,dplus);
		
		int 
			nx = prefh_->size(prefh_->levelmax(),0),
			ny = prefh_->size(prefh_->levelmax(),1),
			nz = prefh_->size(prefh_->levelmax(),2);
		
		cparam_.nx = nx;
		cparam_.ny = ny;
		cparam_.nz = nz;
		cparam_.lx = boxlength;
		cparam_.ly = boxlength;
		cparam_.lz = boxlength;
		cparam_.pcf = pcf_;
		
		kdata_.assign( nx*ny*2*(nz/2+1), 0.0 );
		
		fftw_complex *kdata = reinterpret_cast<fftw_complex*> ( this->get_ptr() );
		
		unsigned nzp = (nz/2+1);
		fac =1.0;
		
		double kfac = 2.0*M_PI/boxlength, ksum = 0.0;
		size_t kcount = 0;
		
		#pragma omp parallel for reduction(+:ksum,kcount)
		for( int i=0; i<nx; ++i )
			for( int j=0; j<ny; ++j )
				for( int k=0; k<nz/2+1; ++k )
				{
					double kx,ky,kz;
					
					kx = (double)i;
					ky = (double)j;
					kz = (double)k;
					
					if( kx > nx/2 ) kx -= nx;
					if( ky > ny/2 ) ky -= ny;
					
					size_t q = ((size_t)i*ny+(size_t)j)*nzp+(size_t)k;
#ifdef FFTW3
					kdata[q][0] = fac*tfk->compute(kfac*sqrt(kx*kx+ky*ky+kz*kz));
					kdata[q][1] = 0.0;
					
					if( k==0 || k==nz/2 )
					{
						ksum  += kdata[q][0];
						kcount++;
					}else{
						ksum  += 2.0*(kdata[q][0]);
						kcount+=2;
					}
#else
					kdata[q].re = fac*tfk->compute(kfac*sqrt(kx*kx+ky*ky+kz*kz));
					kdata[q].im = 0.0;
					
					if( k==0 || k==nz/2 )
					{
						ksum  += kdata[q].re;
						kcount++;
					}else{
						ksum  += 2.0*(kdata[q].re);
						kcount+=2;
					}
					
#endif
				}
		
		delete tfk;
		
		return this;
	}
	
	
	////////////////////////////////////////////////////////////////////////////
	
	template< typename real_t >
	class kernel_real_cached : public kernel
	{
	protected:
		std::vector<real_t> kdata_;
		void precompute_kernel( transfer_function* ptf, tf_type type, const refinement_hierarchy& refh );
		
	public:
		kernel_real_cached( config_file& cf, transfer_function* ptf, refinement_hierarchy& refh, tf_type type )
		: kernel( cf, ptf, refh, type )
		{
			precompute_kernel(ptf, type, refh);
		}
		
		kernel* fetch_kernel( int ilevel, bool isolated=false );
		
		void *get_ptr()
		{	return reinterpret_cast<void*> (&kdata_[0]);	}
		
		~kernel_real_cached()
		{	
			deallocate();
		}
		
		void deallocate()
		{
			std::vector<real_t>().swap( kdata_ );	
		}
		
	};
	
	template< typename real_t >
	kernel* kernel_real_cached<real_t>::fetch_kernel( int ilevel, bool isolated )
	{
		char cachefname[128];
		sprintf(cachefname,"temp_kernel_level%03d.tmp",ilevel);
		FILE *fp = fopen(cachefname,"r");
		
		
		std::cout << " - Fetching kernel for level " << ilevel << std::endl;
		
		LOGUSER("Loading kernel for level %3d from file \'%s\'...",ilevel,cachefname);
		
		if( fp == NULL )
		{	
			LOGERR("Could not open kernel file \'%s\'.",cachefname);
			throw std::runtime_error("Internal error: cached convolution kernel does not exist on disk!");
		}
		
		unsigned nx,ny,nz;
		
		fread(	reinterpret_cast<void*> (&nx), sizeof(unsigned), 1, fp );
		fread(	reinterpret_cast<void*> (&ny), sizeof(unsigned), 1, fp );
		fread(	reinterpret_cast<void*> (&nz), sizeof(unsigned), 1, fp );
		
		kdata_.assign((size_t)nx*(size_t)ny*(size_t)nz,0.0);
		//kdata_.assign((size_t)nx*(size_t)ny*2*((size_t)nz/2+1),0.0);
		
		for( size_t ix=0;ix<nx; ++ix )
		{
			//fread( reinterpret_cast<void*>(&kdata_[0]), sizeof(fftw_real), nx*ny*nz, fp );		
			size_t sz = ny*nz;
			fread( reinterpret_cast<void*>(&kdata_[(size_t)ix * sz]), sizeof(fftw_real), sz, fp );			
		}

		fclose(fp);
		
		//... set parameters
		
		double boxlength = pcf_->getValue<double>("setup","boxlength");
		double dx = boxlength / (1<<ilevel);
		
		cparam_.nx = nx;
		cparam_.ny = ny;
		cparam_.nz = nz-2;
		cparam_.lx = dx * cparam_.nx;
		cparam_.ly = dx * cparam_.ny;
		cparam_.lz = dx * cparam_.nz;
		cparam_.pcf = pcf_;
		
		fftw_real		*rkernel	= reinterpret_cast<fftw_real*>( &kdata_[0] );
		
#ifdef FFTW3
		fftw_complex		*kkernel	= reinterpret_cast<fftw_complex*> (&rkernel[0]);
		fftw_plan plan = fftw_plan_dft_r2c_3d( cparam_.nx, cparam_.ny, cparam_.nz, rkernel, kkernel, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
		
#else
		rfftwnd_plan	plan		= rfftw3d_create_plan( cparam_.nx, cparam_.ny, cparam_.nz,
														  FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE);
		
		
#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), plan, rkernel, NULL );
#else
		rfftwnd_one_real_to_complex( plan, rkernel, NULL );
#endif
		
		rfftwnd_destroy_plan( plan );
#endif
		return this;
	}
	
	
	template< typename real_t >
	void kernel_real_cached<real_t>::precompute_kernel( transfer_function* ptf, tf_type type, const refinement_hierarchy& refh )
	{
		//... compute kernel for finest level
		int nx,ny,nz;
		double dx,lx,ly,lz;
		
		double
			boxlength	= pcf_->getValue<double>("setup","boxlength"),
			boxlength2  = 0.5*boxlength;
		
		int
			levelmax	= refh.levelmax(),
			levelmin	= refh.levelmin();
		
		LOGUSER("Precomputing transfer function kernels...");
			
		nx = refh.size(refh.levelmax(),0);
		ny = refh.size(refh.levelmax(),1);
		nz = refh.size(refh.levelmax(),2);
		
		if( levelmax != levelmin )
		{
			nx *= 2;
			ny *= 2;
			nz *= 2;
		}
		
		dx = boxlength / (1<<refh.levelmax());
		lx = dx * nx;
		ly = dx * ny;
		lz = dx * nz;
		
		double 
		kny			= M_PI/dx,
		fac			= lx*ly*lz/pow(2.0*M_PI,3)/((double)nx*(double)ny*(double)nz),
		
		nspec		= pcf_->getValue<double>("cosmology","nspec"),
		pnorm		= pcf_->getValue<double>("cosmology","pnorm"),
		dplus		= pcf_->getValue<double>("cosmology","dplus");
		
		bool
		bperiodic	= pcf_->getValueSafe<bool>("setup","periodic_TF",true);
		
		std::cout << "   - Computing transfer function kernel...\n";
		
		TransferFunction_real *tfr = 
						new TransferFunction_real( boxlength, 1<<levelmax, type, ptf,nspec,pnorm,dplus,
								  0.25*dx,2.0*boxlength,kny, (int)pow(2,levelmax+2));
		
		
		fftw_real *rkernel = new fftw_real[(size_t)nx*(size_t)ny*2*((size_t)nz/2+1)], *rkernel_coarse;
		
		#pragma omp parallel for
		for( int i=0; i<nx; ++i )
			for( int j=0; j<ny; ++j )
				for( int k=0; k<nz; ++k )
				{
					size_t q=((size_t)(i)*ny + (size_t)(j)) * 2*(nz/2+1) + (size_t)(k);
					rkernel[q] = 0.0;
					
				}
		
		LOGUSER("Computing fine kernel (level %d)...", levelmax);
		
		if( bperiodic  )
		{		
			
			
			#pragma omp parallel for
			for( int i=0; i<=nx/2; ++i )
				for( int j=0; j<=ny/2; ++j )
					for( int k=0; k<=nz/2; ++k )
					{
						int iix(i), iiy(j), iiz(k);
						double rr[3], rr2;
						
						if( iix > (int)nx/2 ) iix -= nx;
						if( iiy > (int)ny/2 ) iiy -= ny;
						if( iiz > (int)nz/2 ) iiz -= nz;
						
						
						//... speed up 8x by copying data to other octants
						size_t idx[8];
						
						idx[0] = ((size_t)(i)*ny + (size_t)(j)) * 2*(nz/2+1) + (size_t)(k);
						idx[1] = ((size_t)(nx-i)*ny + (size_t)(j)) * 2*(nz/2+1) + (size_t)(k);
						idx[2] = ((size_t)(i)*ny + (size_t)(ny-j)) * 2*(nz/2+1) + (size_t)(k);
						idx[3] = ((size_t)(nx-i)*ny + (size_t)(ny-j)) * 2*(nz/2+1) + (size_t)(k);
						idx[4] = ((size_t)(i)*ny + (size_t)(j)) * 2*(nz/2+1) + (size_t)(nz-k);
						idx[5] = ((size_t)(nx-i)*ny + (size_t)(j)) * 2*(nz/2+1) + (size_t)(nz-k);
						idx[6] = ((size_t)(i)*ny + (size_t)(ny-j)) * 2*(nz/2+1) + (size_t)(nz-k);
						idx[7] = ((size_t)(nx-i)*ny + (size_t)(ny-j)) * 2*(nz/2+1) + (size_t)(nz-k);
						
						if(i==0||i==nx/2){ idx[1]=idx[3]=idx[5]=idx[7]=(size_t)-1;}
						if(j==0||j==ny/2){ idx[2]=idx[3]=idx[6]=idx[7]=(size_t)-1;}
						if(k==0||k==nz/2){ idx[4]=idx[5]=idx[6]=idx[7]=(size_t)-1;}
						
						double val = 0.0;
						
						for( int ii=-1; ii<=1; ++ii )
							for( int jj=-1; jj<=1; ++jj )
								for( int kk=-1; kk<=1; ++kk )
								{
									rr[0] = ((double)iix ) * dx + ii*boxlength;
									rr[1] = ((double)iiy ) * dx + jj*boxlength;
									rr[2] = ((double)iiz ) * dx + kk*boxlength;
									
									if( rr[0] > -boxlength && rr[0] <= boxlength
									   && rr[1] > -boxlength && rr[1] <= boxlength
									   && rr[2] > -boxlength && rr[2] <= boxlength )
									{
										rr2 = rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2];
										val += tfr->compute_real(rr2);	
									}
								}
						
						val *= fac;
						
						for(int q=0;q<8;++q)
							if(idx[q]!=(size_t)-1)  
								rkernel[idx[q]] = val;	
					}
			
		}else{
			#pragma omp parallel for
			for( int i=0; i<nx; ++i )
				for( int j=0; j<ny; ++j )
					for( int k=0; k<nz; ++k )
					{
						int iix(i), iiy(j), iiz(k);
						double rr[3], rr2;
						
						if( iix > (int)nx/2 ) iix -= nx;
						if( iiy > (int)ny/2 ) iiy -= ny;
						if( iiz > (int)nz/2 ) iiz -= nz;
						
						size_t idx = ((size_t)i*ny + (size_t)j) * 2*(nz/2+1) + (size_t)k;
						
						rr[0] = ((double)iix ) * dx;
						rr[1] = ((double)iiy ) * dx;
						rr[2] = ((double)iiz ) * dx;
						
						rkernel[idx] = 0.0;
						
						rr2 = rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2];
						if( fabs(rr[0])<=boxlength2||fabs(rr[1])<=boxlength2||fabs(rr[2])<=boxlength2 )
							//if( rr2 <= boxlength2*boxlength2 )
							rkernel[idx] += (fftw_real)tfr->compute_real(rr2)*fac;
						
					}
		}
		
		rkernel[0] = tfr->compute_real(0.0)*fac;
		
		/*************************************************************************************/
		/*************************************************************************************/
		/******* perform deconvolution *******************************************************/
		bool deconv = pcf_->getValueSafe<bool>("setup","deconvolve",true);
		bool deconv_baryons = pcf_->getValueSafe<bool>("setup","deconvolve_baryons",false);
		bool bsmooth_baryons = false;//type==baryon;// && !deconv_baryons;
		bool baryons = type==baryon||type==vbaryon;
		//if( deconv )
		if( deconv && !(type==baryon&&!deconv_baryons) )
		{
				
			LOGUSER("Deconvolving fine kernel...");
			std::cout << " - Deconvoling density kernel...\n";
			
			bool kspacepoisson = ((pcf_->getValueSafe<bool>("poisson","fft_fine",true)|
							 pcf_->getValueSafe<bool>("poisson","kspace",false)))&!baryons ;
			
			double fftnorm = 1.0/((size_t)nx*(size_t)ny*(size_t)nz);
			double k0 = rkernel[0];
			
			fftw_complex	*kkernel	= reinterpret_cast<fftw_complex*>( &rkernel[0] );
			
			//... subtract white noise component before deconvolution
			if(!bsmooth_baryons)
				rkernel[0] = 0.0;
			
#ifdef FFTW3
			fftw_plan 
				plan  = fftw_plan_dft_r2c_3d(nx,ny,nz, rkernel, kkernel, FFTW_ESTIMATE),
				iplan = fftw_plan_dft_c2r_3d(nx,ny,nz, kkernel, rkernel, FFTW_ESTIMATE);			
			
			fftw_execute(plan);
#else
			rfftwnd_plan	plan		= rfftw3d_create_plan( nx,ny,nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE),
							iplan		= rfftw3d_create_plan( nx,ny,nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
			
	#ifndef SINGLETHREAD_FFTW		
			rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), plan, rkernel, NULL );
	#else
			rfftwnd_one_real_to_complex( plan, rkernel, NULL );
	#endif
#endif		
			{
				
				double ksum = 0.0;
				size_t kcount = 0;
				double kmax = 0.5*M_PI/std::max(nx,std::max(ny,nz));
				
				#pragma omp parallel for reduction(+:ksum,kcount)
				for( int i=0; i<nx; ++i )
					for( int j=0; j<ny; ++j )
						for( int k=0; k<nz/2+1; ++k )
						{
							double kx,ky,kz;
							
							kx = (double)i;
							ky = (double)j;
							kz = (double)k;
							
							if( kx > nx/2 ) kx -= nx;
							if( ky > ny/2 ) ky -= ny;
							
							
							double kkmax = kmax;
							size_t q  = ((size_t)i*ny+(size_t)j)*(nz/2+1)+(size_t)k;
							
							if( !bsmooth_baryons )
							{
								if( kspacepoisson )
								{
									//... Use child average response function to emulate sub-sampling
									double ipix = cos(kx*kkmax)*cos(ky*kkmax)*cos(kz*kkmax);
#ifdef FFTW3
									kkernel[q][0] /= ipix;
									kkernel[q][1] /= ipix;
#else
									kkernel[q].re /= ipix;
									kkernel[q].im /= ipix;
#endif
								}else{
									
									//... Use piecewise constant response function (NGP-kernel)
									//... for finite difference methods
									kkmax = kmax;//*2.0;
									double ipix = 1.0;
									if( i > 0 )
										ipix /= sin(kx*2.0*kkmax)/(kx*2.0*kkmax);
									if( j > 0 )
										ipix /= sin(ky*2.0*kkmax)/(ky*2.0*kkmax);
									if( k > 0 )
										ipix /= sin(kz*2.0*kkmax)/(kz*2.0*kkmax);
									
#ifdef FFTW3
									kkernel[q][0] *= ipix;
									kkernel[q][1] *= ipix;
#else
									kkernel[q].re *= ipix;
									kkernel[q].im *= ipix;
#endif
								}
							}
#if 1					
							else{
								//... if smooth==true, convolve with 
								//... NGP kernel to get CIC smoothness
								
								//kkmax *= 2.0;
								
								double ipix = 1.0;
								if( i > 0 )
									ipix /= sin(kx*2.0*kkmax)/(kx*2.0*kkmax);
								if( j > 0 )
									ipix /= sin(ky*2.0*kkmax)/(ky*2.0*kkmax);
								if( k > 0 )
									ipix /= sin(kz*2.0*kkmax)/(kz*2.0*kkmax);
								
#ifdef FFTW3
								kkernel[q][0] /= ipix;
								kkernel[q][1] /= ipix;
#else
								kkernel[q].re /= ipix;
								kkernel[q].im /= ipix;
#endif	
							}
#endif
							//... store k-space average
#ifdef FFTW3
							if( k==0 || k==nz/2 )
							{
								ksum  += kkernel[q][0];
								kcount++;
							}else{
								ksum  += 2.0*(kkernel[q][1]);
								kcount+=2;
							}
#else
							if( k==0 || k==nz/2 )
							{
								ksum  += kkernel[q].re;
								kcount++;
							}else{
								ksum  += 2.0*(kkernel[q].re);
								kcount+=2;
							}
#endif
						}
				
				double dk;
				
				//... re-add white noise component for finest grid
				dk = k0-ksum/kcount;
				
				//... set white noise component to zero if smoothing is enabled
				//if( false )//cparam_.smooth )
				if( bsmooth_baryons )
					dk = 0.0;
				
				//... enforce the r=0 component by adjusting the k-space mean
	#pragma omp parallel for reduction(+:ksum,kcount)
				for( int i=0; i<nx; ++i )
					for( int j=0; j<ny; ++j )
						for( int k=0; k<(nz/2+1); ++k )
						{
							size_t q  = ((size_t)i*ny+(size_t)j)*(nz/2+1)+(size_t)k;
#ifdef FFTW3
							kkernel[q][0] += dk;
							
							kkernel[q][0] *= fftnorm;
							kkernel[q][1] *= fftnorm;
							
#else
							kkernel[q].re += dk;
							
							kkernel[q].re *= fftnorm;
							kkernel[q].im *= fftnorm;
#endif
						}
			}		
			
			
			
#ifdef FFTW3
			fftw_execute(iplan);
			fftw_destroy_plan(plan);
			fftw_destroy_plan(iplan);
#else
	#ifndef SINGLETHREAD_FFTW		
			rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), iplan, kkernel, NULL );
	#else
			rfftwnd_one_complex_to_real( iplan, kkernel, NULL );
	#endif
			rfftwnd_destroy_plan(plan);
			rfftwnd_destroy_plan(iplan);
#endif
		}
		
		/*************************************************************************************/
		/*************************************************************************************/
		/*************************************************************************************/
		
		
		char cachefname[128];
		sprintf(cachefname,"temp_kernel_level%03d.tmp",levelmax);
		LOGUSER("Storing kernel in temp file \'%s\'.",cachefname);

		FILE *fp = fopen(cachefname,"w+");
		unsigned q = nx;
		fwrite(	reinterpret_cast<void*> (&q), sizeof(unsigned), 1, fp );
		q = ny;
		fwrite(	reinterpret_cast<void*> (&q), sizeof(unsigned), 1, fp );
		q = 2*(nz/2+1);
		fwrite(	reinterpret_cast<void*> (&q), sizeof(unsigned), 1, fp );
		
		
		for( int ix=0; ix<nx; ++ix )
		{
			size_t sz = ny*2*(nz/2+1);
			//fwrite( reinterpret_cast<void*>(&rkernel[0]), sizeof(fftw_real), nx*ny*2*(nz/2+1), fp );			
			fwrite( reinterpret_cast<void*>(&rkernel[(size_t)ix * sz]), sizeof(fftw_real), sz, fp );			
		}

		
		
		
		fclose(fp);
		
		//... average and fill for other levels
		for( int ilevel=levelmax-1; ilevel>=levelmin; ilevel-- )
		{
			LOGUSER("Computing coarse kernel (level %d)...",ilevel);
			
			int nxc, nyc, nzc;
			double dxc,lxc,lyc,lzc;
			
			nxc = refh.size(ilevel,0);
			nyc = refh.size(ilevel,1);
			nzc = refh.size(ilevel,2);
			
			if( ilevel!=levelmin )
			{
				nxc *= 2;
				nyc *= 2;
				nzc *= 2;
			}
			
			
			dxc = boxlength / (1<<ilevel);
			lxc = dxc * nxc;
			lyc = dxc * nyc;
			lzc = dxc * nzc;
			
			rkernel_coarse = new fftw_real[(size_t)nxc*(size_t)nyc*2*((size_t)nzc/2+1)];
			fac			= lxc*lyc*lzc/pow(2.0*M_PI,3)/((double)nxc*(double)nyc*(double)nzc);
			
			if( bperiodic  )
			{		
				#pragma omp parallel for
				for( int i=0; i<=nxc/2; ++i )
					for( int j=0; j<=nyc/2; ++j )
						for( int k=0; k<=nzc/2; ++k )
						{
							int iix(i), iiy(j), iiz(k);
							double rr[3], rr2;
							
							if( iix > (int)nxc/2 ) iix -= nxc;
							if( iiy > (int)nyc/2 ) iiy -= nyc;
							if( iiz > (int)nzc/2 ) iiz -= nzc;
							
							//... speed up 8x by copying data to other octants
							size_t idx[8];
							
							idx[0] = ((size_t)(i)*nyc + (size_t)(j)) * 2*(nzc/2+1) + (size_t)(k);
							idx[1] = ((size_t)(nxc-i)*nyc + (size_t)(j)) * 2*(nzc/2+1) + (size_t)(k);
							idx[2] = ((size_t)(i)*nyc + (size_t)(nyc-j)) * 2*(nzc/2+1) + (size_t)(k);
							idx[3] = ((size_t)(nxc-i)*nyc + (size_t)(nyc-j)) * 2*(nzc/2+1) + (size_t)(k);
							idx[4] = ((size_t)(i)*nyc + (size_t)(j)) * 2*(nzc/2+1) + (size_t)(nzc-k);
							idx[5] = ((size_t)(nxc-i)*nyc + (size_t)(j)) * 2*(nzc/2+1) + (size_t)(nzc-k);
							idx[6] = ((size_t)(i)*nyc + (size_t)(nyc-j)) * 2*(nzc/2+1) + (size_t)(nzc-k);
							idx[7] = ((size_t)(nxc-i)*nyc + (size_t)(nyc-j)) * 2*(nzc/2+1) + (size_t)(nzc-k);
							
							if(i==0||i==nxc/2){ idx[1]=idx[3]=idx[5]=idx[7]=(size_t)-1;}
							if(j==0||j==nyc/2){ idx[2]=idx[3]=idx[6]=idx[7]=(size_t)-1;}
							if(k==0||k==nzc/2){ idx[4]=idx[5]=idx[6]=idx[7]=(size_t)-1;}
							
							double val = 0.0;
							
							for( int ii=-1; ii<=1; ++ii )
								for( int jj=-1; jj<=1; ++jj )
									for( int kk=-1; kk<=1; ++kk )
									{
										rr[0] = ((double)iix ) * dxc + ii*boxlength;
										rr[1] = ((double)iiy ) * dxc + jj*boxlength;
										rr[2] = ((double)iiz ) * dxc + kk*boxlength;
										
										if( rr[0] > -boxlength && rr[0] < boxlength
										   && rr[1] > -boxlength && rr[1] < boxlength
										   && rr[2] > -boxlength && rr[2] < boxlength )
										{
											rr2 = rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2];
											val += tfr->compute_real(rr2);	
										}
									}
							
							val *= fac;
							
							for(int qq=0;qq<8;++qq)
								if(idx[qq]!=(size_t)-1)  
									rkernel_coarse[idx[qq]] = val;	
						}
				
			}else{
				#pragma omp parallel for
				for( int i=0; i<nxc; ++i )
					for( int j=0; j<nyc; ++j )
						for( int k=0; k<nzc; ++k )
						{
							int iix(i), iiy(j), iiz(k);
							double rr[3], rr2;
							
							if( iix > (int)nxc/2 ) iix -= nxc;
							if( iiy > (int)nyc/2 ) iiy -= nyc;
							if( iiz > (int)nzc/2 ) iiz -= nzc;
							
							size_t idx = ((size_t)i*nyc + (size_t)j) * 2*(nzc/2+1) + (size_t)k;
							
							rr[0] = ((double)iix ) * dxc;
							rr[1] = ((double)iiy ) * dxc;
							rr[2] = ((double)iiz ) * dxc;
							
							rkernel_coarse[idx] = 0.0;
							
							rr2 = rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2];
							if( fabs(rr[0])<=boxlength2||fabs(rr[1])<=boxlength2||fabs(rr[2])<=boxlength2 )
								rkernel_coarse[idx] += (fftw_real)tfr->compute_real(rr2)*fac;
							
						}
			}

			LOGUSER("Averaging fine kernel to coarse kernel...");

			//... copy averaged and convolved fine kernel to coarse kernel
			#pragma omp parallel for
			for( int ix=0; ix<nx; ix+=2 )
				for( int iy=0; iy<ny; iy+=2 )
					for( int iz=0; iz<nz; iz+=2 )
					{
						int iix(ix/2),iiy(iy/2),iiz(iz/2);
						if( ix>nx/2 ) iix+=nxc-nx/2;
						if( iy>ny/2 ) iiy+=nyc-ny/2;
						if( iz>nz/2 ) iiz+=nzc-nz/2;
						
						if( ix==nx/2||iy==ny/2||iz==nz/2 ) continue;
						
						for( int i=0; i<=1; ++i )
							for( int j=0; j<=1; ++j )
								for( int k=0; k<=1; ++k )
									if(i==0&&k==0&&j==0 )
										rkernel_coarse[ ACC_RC(iix,iiy,iiz) ] =
										0.125 * ( rkernel[ACC_RF(ix-i,iy-j,iz-k)] + rkernel[ACC_RF(ix-i+1,iy-j,iz-k)] 
												 +rkernel[ACC_RF(ix-i,iy-j+1,iz-k)] + rkernel[ACC_RF(ix-i,iy-j,iz-k+1)]
												 +rkernel[ACC_RF(ix-i+1,iy-j+1,iz-k)] + rkernel[ACC_RF(ix-i+1,iy-j,iz-k+1)]
												 +rkernel[ACC_RF(ix-i,iy-j+1,iz-k+1)] + rkernel[ACC_RF(ix-i+1,iy-j+1,iz-k+1)]);
						
									else
									{	
										
										rkernel_coarse[ ACC_RC(iix,iiy,iiz) ] +=
										0.125 * ( rkernel[ACC_RF(ix-i,iy-j,iz-k)] + rkernel[ACC_RF(ix-i+1,iy-j,iz-k)] 
												 +rkernel[ACC_RF(ix-i,iy-j+1,iz-k)] + rkernel[ACC_RF(ix-i,iy-j,iz-k+1)]
												 +rkernel[ACC_RF(ix-i+1,iy-j+1,iz-k)] + rkernel[ACC_RF(ix-i+1,iy-j,iz-k+1)]
												 +rkernel[ACC_RF(ix-i,iy-j+1,iz-k+1)] + rkernel[ACC_RF(ix-i+1,iy-j+1,iz-k+1)]);
									}
						
					}

			
			sprintf(cachefname,"temp_kernel_level%03d.tmp",ilevel);
			LOGUSER("Storing kernel in temp file \'%s\'.",cachefname);
			fp = fopen(cachefname,"w+");
			q = nxc;
			fwrite(	reinterpret_cast<void*> (&q), sizeof(unsigned), 1, fp );
			q = nyc;
			fwrite(	reinterpret_cast<void*> (&q), sizeof(unsigned), 1, fp );
			q = 2*(nzc/2+1);
			fwrite(	reinterpret_cast<void*> (&q), sizeof(unsigned), 1, fp );
			
			for( int ix=0; ix<nxc; ++ix )
			{
				size_t sz = nyc*2*(nzc/2+1);
				//fwrite( reinterpret_cast<void*>(&rkernel_coarse[0]), sizeof(fftw_real), nxc*nyc*2*(nzc/2+1), fp );
				fwrite( reinterpret_cast<void*>(&rkernel_coarse[ix*sz]), sizeof(fftw_real), sz, fp );
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
	
}


/**************************************************************************************/
/**************************************************************************************/


namespace{
	//convolution::kernel_creator_concrete< convolution::kernel_real<double> > creator_d("tf_kernel_real_double");
	//convolution::kernel_creator_concrete< convolution::kernel_real<float> > creator_f("tf_kernel_real_float");
	convolution::kernel_creator_concrete< convolution::kernel_real_cached<double> > creator_d("tf_kernel_real_double");
	convolution::kernel_creator_concrete< convolution::kernel_real_cached<float> > creator_f("tf_kernel_real_float");
	convolution::kernel_creator_concrete< convolution::kernel_k<double> > creator_kd("tf_kernel_k_double");
	convolution::kernel_creator_concrete< convolution::kernel_k<float> > creator_kf("tf_kernel_k_float");
}



