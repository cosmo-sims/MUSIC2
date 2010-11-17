/*
 
 random.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#include "random.hh"


template< typename T >
random_numbers<T>::random_numbers( unsigned res, unsigned cubesize, long baseseed, int *x0, int *lx )
: res_( res ), cubesize_( cubesize ), ncubes_( 1 ), baseseed_( baseseed )
{
	std::cout << " - Generating random numbers (1) with seed " << baseseed << std::endl;
	initialize();
	fill_subvolume( x0, lx );
}

template< typename T >
random_numbers<T>::random_numbers( unsigned res, unsigned cubesize, long baseseed, bool zeromean )
: res_( res ), cubesize_( cubesize ), ncubes_( 1 ), baseseed_( baseseed )
{
	std::cout << " - Generating random numbers (2) with seed " << baseseed << std::endl;
	double mean = 0.0;
	initialize();
	mean = fill_all();
	
	if( zeromean )
	{
		mean = 0.0;
		
#pragma omp parallel for reduction(+:mean)
		for(int i=0; i<(int)res_; ++i )
			for( unsigned j=0; j<res_; ++j )
				for( unsigned k=0; k<res_; ++k )
					mean += (*this)(i,j,k);
		
		mean *= 1.0/(res_*res_*res_);
		
#pragma omp parallel for
		for(int i=0; i<(int)res_; ++i )
			for( unsigned j=0; j<res_; ++j )
				for( unsigned k=0; k<res_; ++k )
					(*this)(i,j,k) = (*this)(i,j,k) - mean;
	}
	
}

template< typename T >
random_numbers<T>::random_numbers( unsigned res, std::string randfname )
: res_( res ), cubesize_( res ), ncubes_(1)
{
	rnums_.push_back( new Meshvar<T>( res, 0, 0, 0 ) );
	
	std::ifstream ifs(randfname.c_str(), std::ios::binary);
	if( !ifs )
		throw std::runtime_error(std::string("Could not open random number file \'")+randfname+std::string("\'!"));
	
	unsigned vartype;
	unsigned nx,ny,nz,blksz;
	int iseed;
	long seed;
	//ifs.read( (char*)&vartype, sizeof(unsigned) );
	
	//... read header .../
	ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
	
	if( blksz != 4*sizeof(int) )
		throw std::runtime_error("corrupt random number file");
	
	ifs.read( reinterpret_cast<char*> (&nx), sizeof(unsigned) );
	ifs.read( reinterpret_cast<char*> (&ny), sizeof(unsigned) );
	ifs.read( reinterpret_cast<char*> (&nz), sizeof(unsigned) );
	ifs.read( reinterpret_cast<char*> (&iseed), sizeof(int) );
	seed = (long)iseed;
	
	if( nx!=res_ || ny!=res_ || nz!=res_ )
	{	
		char errmsg[128];
		sprintf(errmsg,"White noise file dimensions do not match level dimensions: %ux%ux%u vs. %u**3",nx,ny,nz,res_);
		throw std::runtime_error(errmsg);
		
	}
	
	ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
	
	//... read data ...//
	
	ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
	if( blksz == nx*ny*sizeof(float) )
		vartype = 4;
	else if( blksz == nx*ny*sizeof(double) )
		vartype = 8;
	else
		throw std::runtime_error("corrupt random number file");
	
	ifs.seekg(-sizeof(int),std::ios::cur);
	
	std::vector<float> in_float;
	std::vector<double> in_double;
	
	std::cout << " - Random number file \'" << randfname << "\'\n"
	<< "   contains " << nx*ny*nz << " numbers. Reading..." << std::endl;
	
	double sum = 0.0, sum2 = 0.0;
	unsigned count = 0;
	
	if( vartype == 4 )
	{
		for( int ii=0; ii<(int)nz; ++ii )
		{
			ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
			if( blksz != nx*ny*sizeof(float) )
				throw std::runtime_error("corrupt random number file");
			
			in_float.assign(nx*ny,0.0f);
			ifs.read( (char*)&in_float[0], nx*ny*sizeof(float) );
			
			for( int jj=0,q=0; jj<(int)ny; ++jj )
				for( int kk=0; kk<(int)nx; ++kk ){
					sum += in_float[q];
					sum2 += in_float[q]*in_float[q];
					++count;
					
					(*rnums_[0])(kk,jj,ii) = -in_float[q++];
				}
			ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
			if( blksz != nx*ny*sizeof(float) )
				throw std::runtime_error("corrupt random number file");
			
		}
	}
	else if( vartype == 8 )
	{
		for( int ii=0; ii<(int)nz; ++ii )
		{
			ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
			if( blksz != nx*ny*sizeof(double) )
				throw std::runtime_error("corrupt random number file");
			
			in_double.assign(nx*ny,0.0f);				
			ifs.read( (char*)&in_double[0], nx*ny*sizeof(double) );
			
			for( int jj=0,q=0; jj<(int)ny; ++jj )
				for( int kk=0; kk<(int)nx; ++kk )
				{
					sum += in_double[q];
					sum2 += in_double[q]*in_double[q];
					++count;
					(*rnums_[0])(kk,jj,ii) = -in_double[q++];
				}
			ifs.read( reinterpret_cast<char*> (&blksz), sizeof(int) );
			if( blksz != nx*ny*sizeof(double) )
				throw std::runtime_error("corrupt random number file");
			
		}
	}
	
	double mean, var;
	mean = sum/count;
	var = sum2/count-mean*mean;
	
	std::cout << " - Random numbers in file have \n"
	<< "     mean = " << mean << " and var = " << var << std::endl;
	
}

//... copy construct by averaging down
template< typename T >
random_numbers<T>::random_numbers( /*const*/ random_numbers <T>& rc, bool kdegrade )
{
	//if( res > rc.m_res || (res/rc.m_res)%2 != 0 )
	//			throw std::runtime_error("Invalid restriction in random number container copy constructor.");
	
	double sum = 0.0, sum2 = 0.0;
	unsigned count = 0;
	
	if( kdegrade )
	{
		std::cout << " - Generating a coarse white noise field by k-space degrading\n";
		//... initialize properties of container		
		res_		= rc.res_/2;
		cubesize_	= res_;
		ncubes_		= 1;
		baseseed_	= -2;
		
		if( sizeof(fftw_real)!=sizeof(T) )
			throw std::runtime_error("type mismatch with fftw_real in k-space averaging");
		
		fftw_real 
		*rfine = new fftw_real[rc.res_*rc.res_*2*(rc.res_/2+1)],
		*rcoarse = new fftw_real[res_*res_*2*(res_/2+1)];
		
		fftw_complex
		*ccoarse = reinterpret_cast<fftw_complex*> (rcoarse),
		*cfine = reinterpret_cast<fftw_complex*> (rfine);
		
		int nx(rc.res_), ny(rc.res_), nz(rc.res_), nxc(res_), nyc(res_), nzc(res_);
#ifdef FFTW3
		fftw_plan
		pf = fftw_plan_dft_r2c_3d(nx, ny, nz, rfine, cfine, FFTW_ESTIMATE),
		ipc= fftw_plan_dft_c2r_3d(nxc, nyc, nzc, ccoarse, rcoarse, FFTW_ESTIMATE);
#else
		rfftwnd_plan 
		pf	= rfftw3d_create_plan( nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE),
		ipc	= rfftw3d_create_plan( nxc, nyc, nzc, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
#endif
		
#pragma omp parallel for
		for( int i=0; i<nx; i++ )
			for( int j=0; j<ny; j++ )
				for( int k=0; k<nz; k++ )
				{
					size_t q = ((size_t)i*ny+(size_t)j)*(nz+2)+(size_t)k;
					rfine[q] = rc(i,j,k);
				}
		
#ifdef FFTW3
		fftw_execute( pf );
#else
#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), pf, rfine, NULL );
#else
		rfftwnd_one_real_to_complex( pf, rfine, NULL );
#endif
#endif
		
		double fftnorm = 1.0/(nxc*nyc*nzc);
		
#pragma omp parallel for
		for( int i=0; i<nxc; i++ )
			for( int j=0; j<nyc; j++ )
				for( int k=0; k<nzc/2+1; k++ )
				{
					int ii(i),jj(j),kk(k);
					
					if( i > nxc/2 ) ii += nx/2;
					if( j > nyc/2 ) jj += ny/2;
					
					size_t qc,qf;
					
					qc = ((size_t)i*nyc+(size_t)j)*(nzc/2+1)+(size_t)k;
					qf = ((size_t)ii*ny+(size_t)jj)*(nz/2+1)+(size_t)kk;
					
#ifdef FFTW3
					ccoarse[qc][0] = 1.0/sqrt(8.0)*cfine[qf][0]*fftnorm;
					ccoarse[qc][1] = 1.0/sqrt(8.0)*cfine[qf][1]*fftnorm;
#else
					ccoarse[qc].re = 1.0/sqrt(8.0)*cfine[qf].re*fftnorm;
					ccoarse[qc].im = 1.0/sqrt(8.0)*cfine[qf].im*fftnorm;
#endif
				}
		
		delete[] rfine;
#ifdef FFTW3
		fftw_execute( ipc );
#else
#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), ipc, ccoarse, NULL );
#else
		rfftwnd_one_complex_to_real( ipc, ccoarse, NULL );
#endif
#endif
		rnums_.push_back( new Meshvar<T>( res_, 0, 0, 0 ) );
		
#pragma omp parallel for reduction(+:sum,sum2,count)
		for( int i=0; i<nxc; i++ )
			for( int j=0; j<nyc; j++ )
				for( int k=0; k<nzc; k++ )
				{
					size_t q = ((size_t)i*nyc+(size_t)j)*(nzc+2)+(size_t)k;
					(*rnums_[0])(i,j,k) = rcoarse[q];
					sum += (*rnums_[0])(i,j,k);
					sum2+= (*rnums_[0])(i,j,k) * (*rnums_[0])(i,j,k);
					++count;
				}
		
		delete[] rcoarse;
		
#ifdef FFTW3
		fftw_destroy_plan(pf);
		fftw_destroy_plan(ipc);
#else
		rfftwnd_destroy_plan(pf);
		rfftwnd_destroy_plan(ipc);
#endif
		
	}
	else
	{
		std::cout << " - Generating a coarse white noise field by averaging\n";
		if( rc.rnums_.size() == 1 )
		{
			//... initialize properties of container
			res_		= rc.res_/2;
			cubesize_	= res_;
			ncubes_		= 1;
			baseseed_	= -2;
			
			//... use restriction to get consistent random numbers on coarser grid
			mg_straight gop;
			rnums_.push_back( new Meshvar<T>( res_, 0, 0, 0 ) );		
			gop.restrict( *rc.rnums_[0], *rnums_[0] );
			
#pragma omp parallel for reduction(+:sum,sum2,count)
			for( int i=0; i< (int)rnums_[0]->size(0); ++i )
				for( unsigned j=0; j< rnums_[0]->size(1); ++j )
					for( unsigned k=0; k< rnums_[0]->size(2); ++k )
					{
						(*rnums_[0])(i,j,k) *= sqrt(8); //.. maintain that var(delta)=1
						sum += (*rnums_[0])(i,j,k);
						sum2+= (*rnums_[0])(i,j,k) * (*rnums_[0])(i,j,k);
						++count;
					}
		}
		else 
		{
			//... initialize properties of container
			res_		= rc.res_/2;
			cubesize_	= res_;
			ncubes_		= 1;
			baseseed_	= -2;
			
			rnums_.push_back( new Meshvar<T>( res_, 0, 0, 0 ) );
			double fac = 1.0/sqrt(8);
			
#pragma omp parallel for reduction(+:sum,sum2,count)
			for( int ii=0; ii<(int)rc.res_/2; ++ii )
			{	
				unsigned i=2*ii;
				
				for( unsigned j=0,jj=0; j<rc.res_; j+=2,++jj )
					for( unsigned k=0,kk=0; k<rc.res_; k+=2,++kk )
					{
						(*rnums_[0])(ii,jj,kk) = fac * 
						( rc(i,j,k)+rc(i+1,j,k)+rc(i,j+1,k)+rc(i,j,k+1)+
						 rc(i+1,j+1,k)+rc(i+1,j,k+1)+rc(i,j+1,k+1)+rc(i+1,j+1,k+1));
						
						sum += (*rnums_[0])(ii,jj,kk);
						sum2+= (*rnums_[0])(ii,jj,kk) * (*rnums_[0])(ii,jj,kk);
						++count;
					}
			}
		}
	}
	
	double rmean, rvar;
	rmean = sum/count;
	rvar = sum2/count-rmean*rmean;
	
	std::cout << " - Restricted random numbers have\n"
	<< "       mean = " << rmean << ", var = " << rvar << std::endl;
}


template< typename T >
random_numbers<T>::random_numbers( random_numbers<T>& rc, unsigned cubesize, long baseseed, 
								   bool kspace, int *x0_, int *lx_, bool zeromean )
: res_( 2*rc.res_ ), cubesize_( cubesize ), ncubes_( 1 ), baseseed_( baseseed )
{
	
	
	//double mean = 0.0;
	
	initialize();
	
	int x0[3],lx[3];
	if( x0_==NULL || lx_==NULL )
	{	
		for(int i=0;i<3;++i ){
			x0[i]=0;
			lx[i]=res_;
		}
		fill_all();
		
	}
	else
	{
		for(int i=0;i<3;++i ){
			x0[i]=x0_[i];
			lx[i]=lx_[i];
		}
		fill_subvolume( x0, lx );
	}
	
	
	
	if( kspace )
	{
		
		std::cout << " - Generating a constrained random number set with seed " << baseseed << "\n"
		<< "    using coarse mode replacement...\n";
		
		size_t nx=lx[0], ny=lx[1], nz=lx[2],
		nxc=lx[0]/2, nyc=lx[1]/2, nzc=lx[2]/2;
		
		
		fftw_real 
		*rcoarse = new fftw_real[nxc*nyc*(nzc+2)], 
		*rfine = new fftw_real[nx*ny*(nz+2)];
		
		fftw_complex
		*ccoarse = reinterpret_cast<fftw_complex*> (rcoarse),
		*cfine = reinterpret_cast<fftw_complex*> (rfine);
#ifdef FFTW3
		fftw_plan
		pc  = fftw_plan_dft_r2c_3d( nxc, nyc, nzc, rcoarse, ccoarse, FFTW_ESTIMATE),
		pf  = fftw_plan_dft_r2c_3d( nx, ny, nz, rfine, cfine, FFTW_ESTIMATE),
		ipf	= fftw_plan_dft_c2r_3d( nx, ny, nz, cfine, rfine, FFTW_ESTIMATE);
#else
		rfftwnd_plan 
		pc	= rfftw3d_create_plan( nxc, nyc, nzc, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE),
		pf	= rfftw3d_create_plan( nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE),
		ipf	= rfftw3d_create_plan( nx, ny, nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
#endif
		
#pragma omp parallel for
		for( int i=0; i<(int)nx; i++ )
			for( int j=0; j<(int)ny; j++ )
				for( int k=0; k<(int)nz; k++ )
				{
					size_t q = ((size_t)i*(size_t)ny+(size_t)j)*(size_t)(nz+2)+(size_t)k;
					rfine[q] = (*this)(x0[0]+i,x0[1]+j,x0[2]+k);
				}
		
#pragma omp parallel for
		for( int i=0; i<(int)nxc; i++ )
			for( int j=0; j<(int)nyc; j++ )
				for( int k=0; k<(int)nzc; k++ )
				{
					size_t q = ((size_t)i*(size_t)nyc+(size_t)j)*(size_t)(nzc+2)+(size_t)k;
					rcoarse[q] = rc(x0[0]/2+i,x0[1]/2+j,x0[2]/2+k);
				}
#ifdef FFTW3
		fftw_execute( pc );
		fftw_execute( pf );
#else
#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), pc, rcoarse, NULL );
		rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), pf, rfine, NULL );
#else
		rfftwnd_one_real_to_complex( pc, rcoarse, NULL );
		rfftwnd_one_real_to_complex( pf, rfine, NULL );
#endif
#endif
		
		double fftnorm = 1.0/(nx*ny*nz);
		
#pragma omp parallel for
		for( int i=0; i<(int)nxc; i++ )
			for( int j=0; j<(int)nyc; j++ )
				for( int k=0; k<(int)nzc/2+1; k++ )
				{
					int ii(i),jj(j),kk(k);
					
					if( i > (int)nxc/2 ) ii += nx/2;
					if( j > (int)nyc/2 ) jj += ny/2;
					
					size_t qc,qf;
					qc = ((size_t)i*(size_t)nyc+(size_t)j)*(nzc/2+1)+(size_t)k;
					qf = ((size_t)ii*(size_t)ny+(size_t)jj)*(nz/2+1)+(size_t)kk;
#ifdef FFTW3
					cfine[qf][0] = sqrt(8.0)*ccoarse[qc][0];
					cfine[qf][1] = sqrt(8.0)*ccoarse[qc][1];
#else
					cfine[qf].re = sqrt(8.0)*ccoarse[qc].re;
					cfine[qf].im = sqrt(8.0)*ccoarse[qc].im;
#endif
				}
		
		delete[] rcoarse;
		
#pragma omp parallel for
		for( int i=0; i<(int)nx; i++ )
			for( int j=0; j<(int)ny; j++ )
				for( int k=0; k<(int)nz/2+1; k++ )
				{
					size_t q = ((size_t)i*ny+(size_t)j)*(nz/2+1)+(size_t)k;
#ifdef FFTW3
					cfine[q][0] *= fftnorm;
					cfine[q][1] *= fftnorm;
#else
					cfine[q].re *= fftnorm;
					cfine[q].im *= fftnorm;
#endif
				}
		
#ifdef FFTW3
		fftw_execute( ipf );
#else
#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), ipf, cfine, NULL );
#else
		rfftwnd_one_complex_to_real( ipf, cfine, NULL );
#endif
#endif
		
#pragma omp parallel for
		for( int i=0; i<(int)nx; i++ )
			for( int j=0; j<(int)ny; j++ )
				for( int k=0; k<(int)nz; k++ )
				{
					size_t q = ((size_t)i*ny+(size_t)j)*(nz+2)+(size_t)k;
					(*this)(x0[0]+i,x0[1]+j,x0[2]+k) = rfine[q];
				}
		
		delete[] rfine;
		
#ifdef FFTW3
		fftw_destroy_plan(pf);
		fftw_destroy_plan(pc);
		fftw_destroy_plan(ipf);
#else
		fftwnd_destroy_plan(pf);
		fftwnd_destroy_plan(pc);
		fftwnd_destroy_plan(ipf);
#endif
		
	}
	else
	{
		std::cout << " - Generating a constrained random number set with seed " << baseseed << "\n"
		<< "    using Hoffman-Ribak constraints...\n";
		
		double fac = 1./sqrt(8.0);
		
		for( int i=x0[0],ii=x0[0]/2; ii<lx[0]; i+=2,++ii )
			for( int j=x0[1],jj=x0[1]/2; jj<lx[1]; j+=2,++jj )
				for( int k=x0[2],kk=x0[2]/2; kk<lx[2]; k+=2,++kk )
				{
					double topval = rc(ii,jj,kk);
					double locmean = 0.125*((*this)(i,j,k)+(*this)(i+1,j,k)+(*this)(i,j+1,k)+(*this)(i,j,k+1)+
											(*this)(i+1,j+1,k)+(*this)(i+1,j,k+1)+(*this)(i,j+1,k+1)+(*this)(i+1,j+1,k+1));
					double dif = fac*topval-locmean;
					
					(*this)(i,j,k) += dif;
					(*this)(i+1,j,k) += dif;
					(*this)(i,j+1,k) += dif;
					(*this)(i,j,k+1) += dif;
					(*this)(i+1,j+1,k) += dif;
					(*this)(i+1,j,k+1) += dif;
					(*this)(i,j+1,k+1) += dif;
					(*this)(i+1,j+1,k+1) += dif;
					
				}			
	}
}



template< typename T >
double random_numbers<T>::fill_cube( int i, int j, int k)
{
	
	gsl_rng	*RNG = gsl_rng_alloc( gsl_rng_mt19937 );
	
	i = (i+ncubes_)%ncubes_;
	j = (j+ncubes_)%ncubes_;
	k = (k+ncubes_)%ncubes_;
	
	long icube = (i*ncubes_+j)*ncubes_+k;
	long cubeseed = baseseed_+icube; //... each cube gets its unique seed
	
	gsl_rng_set( RNG, cubeseed );
	
	if( rnums_[icube] == NULL )
		rnums_[icube] =  new Meshvar<T>( cubesize_, 0, 0, 0 );
	
	double mean = 0.0;
	
	for( int ii=0; ii<(int)cubesize_; ++ii )
		for( int jj=0; jj<(int)cubesize_; ++jj )
			for( int kk=0; kk<(int)cubesize_; ++kk )
			{	
				(*rnums_[icube])(ii,jj,kk) = gsl_ran_ugaussian_ratio_method( RNG );
				mean += (*rnums_[icube])(ii,jj,kk);
			}
	
	gsl_rng_free( RNG );
	
	return mean/(cubesize_*cubesize_*cubesize_);
}

template< typename T >
void random_numbers<T>::subtract_from_cube( int i, int j, int k, double val )
{
	i = (i+ncubes_)%ncubes_;
	j = (j+ncubes_)%ncubes_;
	k = (k+ncubes_)%ncubes_;
	
	long icube = (i*ncubes_+j)*ncubes_+k;
	
	for( int ii=0; ii<(int)cubesize_; ++ii )
		for( int jj=0; jj<(int)cubesize_; ++jj )
			for( int kk=0; kk<(int)cubesize_; ++kk )
				(*rnums_[icube])(ii,jj,kk) -= val;
	
}

template< typename T >
void random_numbers<T>::free_cube( int i, int j, int k ) 
{
	
	i = (i+ncubes_)%ncubes_;
	j = (j+ncubes_)%ncubes_;
	k = (k+ncubes_)%ncubes_;
	size_t icube = ((size_t)i*(size_t)ncubes_+(size_t)j)*(size_t)ncubes_+(size_t)k;
	
	if( rnums_[icube] != NULL )
	{
		delete rnums_[icube];
		rnums_[icube] = NULL;
	}
}

template< typename T >
void random_numbers<T>::initialize( void )
{
	
	ncubes_ = std::max((int)((double)res_/cubesize_),1);
	if( res_ < cubesize_ )
	{	
		ncubes_ = 1;
		cubesize_ = res_;
	}
	
	std::cout << " - Generating random numbers w/ sample cube size of " << cubesize_ << std::endl;
	
	rnums_.assign( ncubes_*ncubes_*ncubes_, NULL );
}

template< typename T >
double random_numbers<T>::fill_subvolume( int *i0, int *n )
{
	int i0cube[3], ncube[3];
	
	i0cube[0] = (int)((double)i0[0]/cubesize_);
	i0cube[1] = (int)((double)i0[1]/cubesize_);
	i0cube[2] = (int)((double)i0[2]/cubesize_);
	
	ncube[0] = (int)((double)(i0[0]+n[0])/cubesize_+1.0)-i0cube[0];
	ncube[1] = (int)((double)(i0[1]+n[1])/cubesize_+1.0)-i0cube[1];
	ncube[2] = (int)((double)(i0[2]+n[2])/cubesize_+1.0)-i0cube[2];
	
	double mean = 0.0;
	
#pragma omp parallel for reduction(+:mean)
	for( int i=i0cube[0]; i<i0cube[0]+ncube[0]; ++i )
		for( int j=i0cube[1]; j<i0cube[1]+ncube[1]; ++j )
			for( int k=i0cube[2]; k<i0cube[2]+ncube[2]; ++k )
			{
				int ii(i),jj(j),kk(k);
				
				ii = (ii+ncubes_)%ncubes_;
				jj = (jj+ncubes_)%ncubes_;
				kk = (kk+ncubes_)%ncubes_;
				
				mean += fill_cube(ii, jj, kk);
			}
	return mean/(ncube[0]*ncube[1]*ncube[2]);
}

template< typename T >
double random_numbers<T>::fill_all( void )
{
	double sum = 0.0;
	
#pragma omp parallel for reduction(+:sum)
	for( int i=0; i<(int)ncubes_; ++i )
		for( int j=0; j<(int)ncubes_; ++j )
			for( int k=0; k<(int)ncubes_; ++k )
			{
				int ii(i),jj(j),kk(k);
				
				ii = (ii+ncubes_)%ncubes_;
				jj = (jj+ncubes_)%ncubes_;
				kk = (kk+ncubes_)%ncubes_;
				
				sum+=fill_cube(ii, jj, kk);
			}
	
	//... subtract mean
#pragma omp parallel for reduction(+:sum)
	for( int i=0; i<(int)ncubes_; ++i )
		for( int j=0; j<(int)ncubes_; ++j )
			for( int k=0; k<(int)ncubes_; ++k )
			{
				int ii(i),jj(j),kk(k);
				
				ii = (ii+ncubes_)%ncubes_;
				jj = (jj+ncubes_)%ncubes_;
				kk = (kk+ncubes_)%ncubes_;
				subtract_from_cube(ii,jj,kk,sum/(ncubes_*ncubes_*ncubes_));
			}
	
	/////////////////////////////////////////////////////
	
#if defined(DEGRADE_RAND1)
	
	{
		std::cerr << " - degrading field for 1 level...(" << res_ << ")\n";
		//unsigned ixoff=51,iyoff=51,izoff=51;
		//unsigned nx=52, ny=52, nz=52;
		int ixoff=102, iyoff=102, izoff=102;
		int nx=104, ny=104, nz=104;
		
#pragma omp parallel for
		for( int ix=0; ix<(int)res_; ix+=2 )
			for( int iy=0; iy<(int)res_; iy+=2 )
				for( int iz=0; iz<(int)res_; iz+=2 )
				{
					if( ix>=2*ixoff && ix < 2*ixoff+nx 
					   && iy>=2*iyoff && iy < 2*iyoff+ny
					   && iz>=2*izoff && iz < 2*izoff+nz )
					{
						continue;
					}
					
					double avg = 0.125*((*this)(ix,iy,iz)+(*this)(ix+1,iy,iz)+(*this)(ix,iy+1,iz)+(*this)(ix,iy,iz+1)+
										(*this)(ix+1,iy+1,iz)+(*this)(ix+1,iy,iz+1)+(*this)(ix,iy+1,iz+1)+(*this)(ix+1,iy+1,iz+1));
					
					
					
					(*this)(ix,iy,iz) = avg;
					(*this)(ix+1,iy,iz) = avg;
					(*this)(ix,iy+1,iz) = avg;
					(*this)(ix,iy,iz+1) = avg;
					(*this)(ix+1,iy+1,iz) = avg;
					(*this)(ix+1,iy,iz+1) = avg;
					(*this)(ix,iy+1,iz+1) = avg;
					(*this)(ix+1,iy+1,iz+1) = avg;
					
				}
		
		
	}
	
#elif defined(DEGRADE_RAND2)
	
	{
		std::cerr << " - degrading field for 2 level...(" << res_ << ")\n";
		
		/*unsigned ixoff2=102,iyoff2=102,izoff2=102;
		 unsigned nx2=52, ny2=52, nz2=52;
		 
		 unsigned ixoff1=86,iyoff1=86,izoff1=86;
		 unsigned nx1=168,ny1=168,nz1=168;*/
		
		int ixoff2=86,iyoff2=86,izoff2=86;
		int nx2=84,ny2=84,nz2=84;
		
		int ixoff1=102,iyoff1=102,izoff1=102;
		int nx1=104,ny1=104,nz1=104;
		
		
		//unsigned ixoff2=51,iyoff2=51,izoff2=51;
		//unsigned nx2=52, ny2=52, nz2=52;
		
		//unsigned ixoff1=34,iyoff1=34,izoff1=34;
		//unsigned nx1=120,ny1=120,nz1=120;
		
		//unsigned ixoff2=84, iyoff2=84, izoff2=84;
		//unsigned nx2=90, ny2=90, nz2=90;
		
		//unsigned ixoff1=100, iyoff1=100, izoff1=100;
		//unsigned nx1=112, ny1=112, nz1=112;
		
		
		/*unsigned ixoff2=100, iyoff2=100, izoff2=100;
		 unsigned nx2=112, ny2=112, nz2=112;
		 
		 unsigned ixoff1=84, iyoff1=84, izoff1=84;
		 unsigned nx1=180, ny1=180, nz1=180;*/
		
#pragma omp parallel for
		for( int ix=0; ix<(int)res_; ix+=2 )
			for( int iy=0; iy<(int)res_; iy+=2 )
				for( int iz=0; iz<(int)res_; iz+=2 )
				{
					if( ix>=2*ixoff2 && ix < 2*ixoff2+nx2 
					   && iy>=2*iyoff2 && iy < 2*iyoff2+ny2
					   && iz>=2*izoff2 && iz < 2*izoff2+nz2 )
					{
						continue;
					}
					
					double avg = 0.125*((*this)(ix,iy,iz)+(*this)(ix+1,iy,iz)+(*this)(ix,iy+1,iz)+(*this)(ix,iy,iz+1)+
										(*this)(ix+1,iy+1,iz)+(*this)(ix+1,iy,iz+1)+(*this)(ix,iy+1,iz+1)+(*this)(ix+1,iy+1,iz+1));
					
					
					
					(*this)(ix,iy,iz) = avg;
					(*this)(ix+1,iy,iz) = avg;
					(*this)(ix,iy+1,iz) = avg;
					(*this)(ix,iy,iz+1) = avg;
					(*this)(ix+1,iy+1,iz) = avg;
					(*this)(ix+1,iy,iz+1) = avg;
					(*this)(ix,iy+1,iz+1) = avg;
					(*this)(ix+1,iy+1,iz+1) = avg;
					
				}
		
#pragma omp parallel for
		for( int ix=0; ix<(int)res_; ix+=4 )
			for( int iy=0; iy<(int)res_; iy+=4 )
				for( int iz=0; iz<(int)res_; iz+=4 )
				{
					if( ix>=2*ixoff1 && ix < 2*ixoff1+nx1
					   && iy>=2*iyoff1 && iy < 2*iyoff1+ny1
					   && iz>=2*izoff1 && iz < 2*izoff1+nz1 )
					{
						continue;
					}
					double avg = 0.0;//0.125*((*this)(ix,iy,iz)+(*this)(ix+1,iy,iz)+(*this)(ix,iy+1,iz)+(*this)(ix,iy,iz+1)+
									 //(*this)(ix+1,iy+1,iz)+(*this)(ix+1,iy,iz+1)+(*this)(ix,iy+1,iz+1)+(*this)(ix+1,iy+1,iz+1));
					for( int i=0; i<4; ++i )
						for( int j=0; j<4; ++j )
							for( int k=0; k<4; ++k )
								avg += (*this)(ix+i,iy+j,iz+k);
					avg /=4.*4.*4.;
					
					for( int i=0; i<4; ++i )
						for( int j=0; j<4; ++j )
							for( int k=0; k<4; ++k )
								(*this)(ix+i,iy+j,iz+k) = avg;
				}
		
		
	}
#endif
	
	
	/////////////////////////////////////////////////////
	
	return sum/(ncubes_*ncubes_*ncubes_);
}


template< typename T >
void random_numbers<T>:: print_allocated( void )
{
	unsigned ncount = 0, ntot = rnums_.size();
	for( size_t i=0; i<rnums_.size(); ++i )
		if( rnums_[i]!=NULL ) ncount++;
	
	std::cerr << " -> " << ncount << " of " << ntot << " random number cubes currently allocated\n";
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma mark -
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////


template< typename rng, typename T >
random_number_generator<rng,T>::random_number_generator( config_file& cf, refinement_hierarchy& refh, transfer_function *ptf )
: pcf_( &cf ), prefh_( &refh ), constraints( cf, ptf )
{
	levelmin_ = prefh_->levelmin();
	levelmax_ = prefh_->levelmax();
	
	ran_cube_size_	= pcf_->getValueSafe<unsigned>("random","cubesize",DEF_RAN_CUBE_SIZE);
	disk_cached_ = pcf_->getValueSafe<bool>("random","disk_cached",true);
	mem_cache_.assign(levelmax_-levelmin_+1,NULL);
	
	
	////disk_cached_ = false;
	
	//... determine seed/white noise file data to be applied
	parse_rand_parameters();
	
	//... compute the actual random numbers
	compute_random_numbers();
}



template< typename rng, typename T >
random_number_generator<rng,T>::~random_number_generator()
{  
	
	//... clear memory caches
	for( unsigned i=0; i<mem_cache_.size(); ++i )
		if( mem_cache_[i] != NULL )
			delete mem_cache_[i];
	
	
	//... clear disk caches
	if( disk_cached_ )
	{
		for( int ilevel=levelmin_; ilevel<=levelmax_; ++ilevel )
		{
			char fname[128];
			sprintf(fname,"wnoise_%04d.bin",ilevel);
			unlink(fname);
		}
	}
}




template< typename rng, typename T >
bool random_number_generator<rng,T>::is_number(const std::string& s)
{
	for (size_t i = 0; i < s.length(); i++)
		if (!std::isdigit(s[i])&&s[i]!='-' )
			return false;
	
	return true;
}

template< typename rng, typename T >
void random_number_generator<rng,T>::parse_rand_parameters( void )
{
	//... parse random number options
	for( int i=0; i<=100; ++i )
	{
		char seedstr[128];
		std::string tempstr;
		sprintf(seedstr,"seed[%d]",i);
		if( pcf_->containsKey( "random", seedstr ) )
			tempstr = pcf_->getValue<std::string>( "random", seedstr );
		else
			// "-2" means that no seed entry was found for that level
			tempstr = std::string("-2");
		
		if( is_number( tempstr ) )
		{	
			long ltemp;
			pcf_->convert( tempstr, ltemp );
			rngfnames_.push_back( "" );
			if( ltemp < 0 )
				//... generate some dummy seed which only depends on the level, negative so we know it's not
				//... an actual seed and thus should not be used as a constraint for coarse levels
				rngseeds_.push_back( -abs((unsigned)(ltemp-i)%123+(unsigned)(ltemp+827342523521*i)%123456789) );
			else
				rngseeds_.push_back( ltemp );
		}else{
			rngfnames_.push_back( tempstr );
			rngseeds_.push_back(-1);
			std::cout << " - Random numbers for level " << std::setw(3) << i << " will be read from file.\n";
		}
		
	}
	
	//.. determine for which levels random seeds/random number files are given
	levelmin_seed_ = -1;
	for( unsigned ilevel = 0; ilevel < rngseeds_.size(); ++ilevel )
	{	
		if( levelmin_seed_ < 0 && (rngfnames_[ilevel].size() > 0 || rngseeds_[ilevel] > 0) )
			levelmin_seed_ = ilevel;
	}
	
}


template< typename rng, typename T >
void random_number_generator<rng,T>::correct_avg( int icoarse, int ifine )
{
	int shift[3], levelmin_poisson;
	shift[0] = pcf_->getValue<int>("setup","shift_x");
	shift[1] = pcf_->getValue<int>("setup","shift_y");
	shift[2] = pcf_->getValue<int>("setup","shift_z");
	
	levelmin_poisson = pcf_->getValue<unsigned>("setup","levelmin");
	
	int lfacc = 1<<(icoarse-levelmin_poisson);
	
	
	
	
	int nc[3], i0c[3], nf[3], i0f[3];
	if( icoarse != levelmin_ )
	{
		nc[0]  = 2*prefh_->size(icoarse, 0);
		nc[1]  = 2*prefh_->size(icoarse, 1);
		nc[2]  = 2*prefh_->size(icoarse, 2);
		i0c[0] = prefh_->offset_abs(icoarse, 0) - lfacc*shift[0] - nc[0]/4;
		i0c[1] = prefh_->offset_abs(icoarse, 1) - lfacc*shift[1] - nc[1]/4;
		i0c[2] = prefh_->offset_abs(icoarse, 2) - lfacc*shift[2] - nc[2]/4;
		
	}
	else
	{
		nc[0]  = prefh_->size(icoarse, 0);
		nc[1]  = prefh_->size(icoarse, 1);
		nc[2]  = prefh_->size(icoarse, 2);
		i0c[0] = - lfacc*shift[0];
		i0c[1] = - lfacc*shift[1];
		i0c[2] = - lfacc*shift[2];
	}
	nf[0]  = 2*prefh_->size(ifine, 0);
	nf[1]  = 2*prefh_->size(ifine, 1);
	nf[2]  = 2*prefh_->size(ifine, 2);
	i0f[0] = prefh_->offset_abs(ifine, 0) - 2*lfacc*shift[0] - nf[0]/4;
	i0f[1] = prefh_->offset_abs(ifine, 1) - 2*lfacc*shift[1] - nf[1]/4;
	i0f[2] = prefh_->offset_abs(ifine, 2) - 2*lfacc*shift[2] - nf[2]/4;
	
	//.................................
	if( disk_cached_ )
	{
		char fncoarse[128], fnfine[128];
		sprintf(fncoarse,"wnoise_%04d.bin",icoarse);
		sprintf(fnfine,"wnoise_%04d.bin",ifine);
		
		std::ifstream 
		iffine( fnfine, std::ios::binary ), 
		ifcoarse( fncoarse, std::ios::binary );
		
		int nxc,nyc,nzc,nxf,nyf,nzf;
		iffine.read( reinterpret_cast<char*> (&nxf), sizeof(unsigned) );
		iffine.read( reinterpret_cast<char*> (&nyf), sizeof(unsigned) );
		iffine.read( reinterpret_cast<char*> (&nzf), sizeof(unsigned) );
		
		ifcoarse.read( reinterpret_cast<char*> (&nxc), sizeof(unsigned) );
		ifcoarse.read( reinterpret_cast<char*> (&nyc), sizeof(unsigned) );
		ifcoarse.read( reinterpret_cast<char*> (&nzc), sizeof(unsigned) );
		
		if( nxf!=nf[0] || nyf!=nf[1] || nzf!=nf[2] || nxc!=nc[0] || nyc!=nc[1] || nzc!=nc[2] )
			throw std::runtime_error("White noise file mismatch. This should not happen. Notify a developer!");
		
		int nxd(nxf/2),nyd(nyf/2),nzd(nzf/2);
		std::vector<T> deg_rand( nxd*nyd*nzd, 0.0 );
		double fac = 1.0/sqrt(8.0);
		
		for( int i=0; i<nxf; i+=2 )
		{	
			std::vector<T> fine_rand( 2*nyf*nzf, 0.0 );
			iffine.read( reinterpret_cast<char*> (&fine_rand[0]), 2*nyf*nzf*sizeof(T) );
			
#pragma omp parallel for
			for( int j=0; j<nyf; j+=2 )
				for( int k=0; k<nzf; k+=2 )
				{
					unsigned qc = ((i/2)*nyd+(j/2))*nzd+(k/2);
					unsigned qf[8];
					qf[0] = (0*nyf+j+0)*nzf+k+0;
					qf[1] = (0*nyf+j+0)*nzf+k+1;
					qf[2] = (0*nyf+j+1)*nzf+k+0;
					qf[3] = (0*nyf+j+1)*nzf+k+1;
					qf[4] = (1*nyf+j+0)*nzf+k+0;
					qf[5] = (1*nyf+j+0)*nzf+k+1;
					qf[6] = (1*nyf+j+1)*nzf+k+0;
					qf[7] = (1*nyf+j+1)*nzf+k+1;
					
					for( int q=0; q<8; ++q )
						deg_rand[qc] += fac*fine_rand[qf[q]];
					
				}
		}
		
		//... now deg_rand holds the oct-averaged fine field, store this in the coarse field
		std::vector<T> coarse_rand(nxc*nyc*nzc,0.0);
		ifcoarse.read( reinterpret_cast<char*> (&coarse_rand[0]), nxc*nyc*nzc*sizeof(T) );
		
		int di,dj,dk;
		
		di = i0f[0]/2-i0c[0];
		dj = i0f[1]/2-i0c[1];
		dk = i0f[2]/2-i0c[2];
		
#pragma omp parallel for
		for( int i=0; i<nxd; i++ )
			for( int j=0; j<nyd; j++ )
				for( int k=0; k<nzd; k++ )
				{
					//unsigned qc = (((i+di+nxc)%nxc)*nyc+(((j+dj+nyc)%nyc)))*nzc+((k+dk+nzc)%nzc);
					
					if( i+di < 0 || i+di >= nxc || j+dj < 0 || j+dj >= nyc || k+dk < 0 || k+dk >= nzc )
						continue;
					
					unsigned qc = ((i+di)*nyc+(j+dj))*nzc+(k+dk);
					unsigned qcd = (i*nyd+j)*nzd+k;
					
					coarse_rand[qc] = deg_rand[qcd];
				}
		
		deg_rand.clear();
		
		ifcoarse.close();
		std::ofstream ofcoarse( fncoarse, std::ios::binary|std::ios::trunc );
		ofcoarse.write( reinterpret_cast<char*> (&nxc), sizeof(unsigned) );
		ofcoarse.write( reinterpret_cast<char*> (&nyc), sizeof(unsigned) );
		ofcoarse.write( reinterpret_cast<char*> (&nzc), sizeof(unsigned) );
		ofcoarse.write( reinterpret_cast<char*> (&coarse_rand[0]), nxc*nyc*nzc*sizeof(T) );
		ofcoarse.close();	
	}
	else
	{
		int nxc,nyc,nzc,nxf,nyf,nzf;
		nxc = nc[0]; nyc = nc[1]; nzc = nc[2];
		nxf = nf[0]; nyf = nf[1]; nzf = nf[2];
		int nxd(nxf/2),nyd(nyf/2),nzd(nzf/2);
		
		int di,dj,dk;
		
		di = i0f[0]/2-i0c[0];
		dj = i0f[1]/2-i0c[1];
		dk = i0f[2]/2-i0c[2];
		
		double fac = 1.0/sqrt(8.0);
		
#pragma omp parallel for
		for( int i=0; i<nxd; i++ )
			for( int j=0; j<nyd; j++ )
				for( int k=0; k<nzd; k++ )
				{
					if( i+di < 0 || i+di >= nxc || j+dj < 0 || j+dj >= nyc || k+dk < 0 || k+dk >= nzc )
						continue;
					
					unsigned qf[8];
					qf[0] = ((2*i+0)*nyf+2*j+0)*nzf+2*k+0;
					qf[1] = ((2*i+0)*nyf+2*j+0)*nzf+2*k+1;
					qf[2] = ((2*i+0)*nyf+2*j+1)*nzf+2*k+0;
					qf[3] = ((2*i+0)*nyf+2*j+1)*nzf+2*k+1;
					qf[4] = ((2*i+1)*nyf+2*j+0)*nzf+2*k+0;
					qf[5] = ((2*i+1)*nyf+2*j+0)*nzf+2*k+1;
					qf[6] = ((2*i+1)*nyf+2*j+1)*nzf+2*k+0;
					qf[7] = ((2*i+1)*nyf+2*j+1)*nzf+2*k+1;
					
					double finesum = 0.0;
					for( int q=0; q<8; ++q )
						finesum += fac*(*mem_cache_[ifine-levelmin_])[qf[q]];
					
					size_t qc = ((size_t)(i+di)*nyc+(size_t)(j+dj))*nzc+(size_t)(k+dk);
					
					(*mem_cache_[icoarse-levelmin_])[qc] = finesum;
				}						
	}
	
	
}



template< typename rng, typename T >
void random_number_generator<rng,T>::compute_random_numbers( void )
{
	bool kavg = pcf_->getValueSafe<bool>("random","kaveraging",true);
	
	std::vector< rng* > randc(std::max(levelmax_,levelmin_seed_)+1,(rng*)NULL);
	
	//--- FILL ALL WHITE NOISE ARRAYS FOR WHICH WE NEED THE FULL FIELD ---//
	
	//... seeds are given for a level coarser than levelmin
	if( levelmin_seed_ < levelmin_ )
	{
		if( rngfnames_[levelmin_seed_].size() > 0 )
			randc[levelmin_seed_] 
			= new rng( 1<<levelmin_seed_, rngfnames_[levelmin_seed_] );
		else
			randc[levelmin_seed_]
			= new rng( 1<<levelmin_seed_, ran_cube_size_, rngseeds_[levelmin_seed_], true );
		
		for( int i=levelmin_seed_+1; i<=levelmin_; ++i )
		{
#warning add possibility to read noise from file also here!
			randc[i] = new rng( *randc[i-1], ran_cube_size_, rngseeds_[i], kavg );
			delete randc[i-1];
			randc[i-1] = NULL;
		}
	}
	
	//... seeds are given for a level finer than levelmin, obtain by averaging
	if( levelmin_seed_ > levelmin_ )
	{
		randc[levelmin_seed_] = new rng( 1<<levelmin_seed_, ran_cube_size_, rngseeds_[levelmin_seed_], true );//, x0, lx );
		
		for( int ilevel = levelmin_seed_-1; ilevel >= (int)levelmin_; --ilevel ){
			if( rngseeds_[ilevel-levelmin_] > 0 )
				std::cerr << " - Warning: random seed for level " << ilevel << " will be ignored.\n"
				<< "            consistency requires that it is obtained by restriction from level " << levelmin_seed_ << std::endl;
			
			
			
			//if( levelmin_ == levelmax_ )
			if( ilevel >= levelmax_ )
				randc[ilevel] = new rng( *randc[ilevel+1], kavg );
			else
				randc[ilevel] = new rng( *randc[ilevel+1], false );
			
			if( ilevel+1 > levelmax_ )
			{
				delete randc[ilevel+1];
				randc[ilevel+1] = NULL;
			}
		}
		
	}
	
	//--- GENERATE AND STORE ALL LEVELS, INCLUDING REFINEMENTS ---//
	
	//... levelmin
	if( randc[levelmin_] == NULL )
		if( rngfnames_[levelmin_].size() > 0 )
			randc[levelmin_] = new rng( 1<<levelmin_, rngfnames_[levelmin_] );
		else
			randc[levelmin_] = new rng( 1<<levelmin_, ran_cube_size_, rngseeds_[levelmin_], true );
	
	//if( levelmax_ == levelmin_ )
	{
		//... apply constraints to coarse grid
		//... if no constraints are specified, or not for this level
		//... constraints.apply will return without doing anything
		int x0[3] = { 0, 0, 0 };
		int lx[3] = { 1<<levelmin_, 1<<levelmin_, 1<<levelmin_ };
		constraints.apply( levelmin_, x0, lx, randc[levelmin_] );
		
	}
	
	store_rnd( levelmin_, randc[levelmin_] );
	
	
	
	//... refinement levels
	for( int ilevel=levelmin_+1; ilevel<=levelmax_; ++ilevel )
	{
		int lx[3], x0[3];
		int shift[3], levelmin_poisson;
		shift[0] = pcf_->getValue<int>("setup","shift_x");
		shift[1] = pcf_->getValue<int>("setup","shift_y");
		shift[2] = pcf_->getValue<int>("setup","shift_z");
		
		levelmin_poisson = pcf_->getValue<unsigned>("setup","levelmin");
		
		int lfac = 1<<(ilevel-levelmin_poisson);
		
		lx[0] = 2*prefh_->size(ilevel, 0);
		lx[1] = 2*prefh_->size(ilevel, 1);
		lx[2] = 2*prefh_->size(ilevel, 2);
		x0[0] = prefh_->offset_abs(ilevel, 0) - lfac*shift[0] - lx[0]/4;
		x0[1] = prefh_->offset_abs(ilevel, 1) - lfac*shift[1] - lx[1]/4;
		x0[2] = prefh_->offset_abs(ilevel, 2) - lfac*shift[2] - lx[2]/4;
		
		if( randc[ilevel] == NULL )
			randc[ilevel] = new rng( *randc[ilevel-1], ran_cube_size_, rngseeds_[ilevel], kavg, x0, lx );
		delete randc[ilevel-1];
		randc[ilevel-1] = NULL;
		
		//... apply constraints to this level, if any
		//if( ilevel == levelmax_ )
		constraints.apply( ilevel, x0, lx, randc[ilevel] );
		
		//... store numbers
		store_rnd( ilevel, randc[ilevel] );
	}
	
	delete randc[levelmax_];
	randc[levelmax_] = NULL;
	
	//... make sure that the coarse grid contains oct averages where it overlaps with a fine grid
	//... this also ensures that constraints enforced on fine grids are carried to the coarser grids
	for( int ilevel=levelmax_; ilevel>levelmin_; --ilevel )
		correct_avg( ilevel-1, ilevel );
	
	
	
	
	
	//.. we do not have random numbers for a coarse level, generate them
	/*if( levelmax_rand_ >= (int)levelmin_ )
	 {
	 std::cerr << "lmaxread >= (int)levelmin\n";
	 randc[levelmax_rand_] = new rng( (unsigned)pow(2,levelmax_rand_), rngfnames_[levelmax_rand_] );
	 for( int ilevel = levelmax_rand_-1; ilevel >= (int)levelmin_; --ilevel )
	 randc[ilevel] = new rng( *randc[ilevel+1] );
	 }*/
}


template< typename rng, typename T >
void random_number_generator<rng,T>:: store_rnd( int ilevel, rng* prng )
{
	int shift[3], levelmin_poisson;
	shift[0] = pcf_->getValue<int>("setup","shift_x");
	shift[1] = pcf_->getValue<int>("setup","shift_y");
	shift[2] = pcf_->getValue<int>("setup","shift_z");
	
	levelmin_poisson = pcf_->getValue<unsigned>("setup","levelmin");
	
	int lfac = 1<<(ilevel-levelmin_poisson);
	
	if( disk_cached_ )
	{
		std::vector<T> data;
		if( ilevel == levelmin_ )
		{
			int N = 1<<levelmin_;
			int i0,j0,k0;
			
			i0 = -lfac*shift[0];
			j0 = -lfac*shift[1];
			k0 = -lfac*shift[2];
			
			char fname[128];
			sprintf(fname,"wnoise_%04d.bin",ilevel);
			
			LOGUSER("Storing white noise field in file \'%s\'...", fname );
			
			std::ofstream ofs(fname,std::ios::binary|std::ios::trunc);
			
			ofs.write( reinterpret_cast<char*> (&N), sizeof(unsigned) );
			ofs.write( reinterpret_cast<char*> (&N), sizeof(unsigned) );
			ofs.write( reinterpret_cast<char*> (&N), sizeof(unsigned) );
			
			data.assign( N*N, 0.0 );
			for( int i=0; i<N; ++i )
			{	
#pragma omp parallel for
				for( int j=0; j<N; ++j )
					for( int k=0; k<N; ++k )
						data[j*N+k] = (*prng)(i+i0,j+j0,k+k0);
				
				ofs.write(reinterpret_cast<char*> (&data[0]), N*N*sizeof(T) );
			}
			ofs.close();
		}
		else
		{
			int nx,ny,nz;
			int i0,j0,k0;
			
			nx = 2*prefh_->size(ilevel, 0);
			ny = 2*prefh_->size(ilevel, 1);
			nz = 2*prefh_->size(ilevel, 2);
			i0 = prefh_->offset_abs(ilevel, 0) - lfac*shift[0] - nx/4;
			j0 = prefh_->offset_abs(ilevel, 1) - lfac*shift[1] - nx/4;
			k0 = prefh_->offset_abs(ilevel, 2) - lfac*shift[2] - nx/4;
			
			char fname[128];
			sprintf(fname,"wnoise_%04d.bin",ilevel);
			
			LOGUSER("Storing white noise field in file \'%s\'...", fname );
			
			std::ofstream ofs(fname,std::ios::binary|std::ios::trunc);
			
			ofs.write( reinterpret_cast<char*> (&nx), sizeof(unsigned) );
			ofs.write( reinterpret_cast<char*> (&ny), sizeof(unsigned) );
			ofs.write( reinterpret_cast<char*> (&nz), sizeof(unsigned) );
			
			data.assign( ny*nz, 0.0 );
			for( int i=0; i<nx; ++i )
			{	
#pragma omp parallel for
				for( int j=0; j<ny; ++j )
					for( int k=0; k<nz; ++k )
						data[j*nz+k] = (*prng)(i+i0,j+j0,k+k0);
				
				ofs.write(reinterpret_cast<char*> (&data[0]), ny*nz*sizeof(T) );
			}
			ofs.close();
			
		}
		
	}
	else 
	{
		int nx,ny,nz;
		int i0,j0,k0;
		
		if( ilevel == levelmin_ )
		{
			i0 = -lfac*shift[0];
			j0 = -lfac*shift[1];
			k0 = -lfac*shift[2];
			
			nx = ny = nz = 1<<levelmin_;
		}
		else
		{
			nx = 2*prefh_->size(ilevel, 0);
			ny = 2*prefh_->size(ilevel, 1);
			nz = 2*prefh_->size(ilevel, 2);
			i0 = prefh_->offset_abs(ilevel, 0) - lfac*shift[0] - nx/4;
			j0 = prefh_->offset_abs(ilevel, 1) - lfac*shift[1] - nx/4;
			k0 = prefh_->offset_abs(ilevel, 2) - lfac*shift[2] - nx/4;
		}
		
		mem_cache_[ilevel-levelmin_] = new std::vector<T>(nx*ny*nz,0.0);
		
		LOGUSER("Copying white noise to mem cache...");
		
#pragma omp parallel for
		for( int i=0; i<nx; ++i )
			for( int j=0; j<ny; ++j )
				for( int k=0; k<nz; ++k )
					(*mem_cache_[ilevel-levelmin_])[((size_t)i*ny+(size_t)j)*nz+(size_t)k] = (*prng)(i+i0,j+j0,k+k0);
		
	}		
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma mark -
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////


template class random_numbers<float>;
template class random_numbers<double>;
template class random_number_generator< random_numbers<float>, float >;
template class random_number_generator< random_numbers<double>, double >;
