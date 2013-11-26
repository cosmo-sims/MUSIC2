/*
 
 random.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#include "random.hh"

// TODO: move all this into a plugin!!!

std::map< std::string, RNG_plugin_creator *>& 
get_RNG_plugin_map()
{
  static std::map< std::string, RNG_plugin_creator* > RNG_plugin_map;
  return RNG_plugin_map;
}

void print_RNG_plugins()
{
  std::map< std::string, RNG_plugin_creator *>& m = get_RNG_plugin_map();
  std::map< std::string, RNG_plugin_creator *>::iterator it;
  it = m.begin();
  std::cout << " - Available random number generator plug-ins:\n";
  while( it!=m.end() )
    {
      if( (*it).second )
	std::cout << "\t\'" << (*it).first << "\'\n";
      ++it;
    }	
}

RNG_plugin *select_RNG_plugin( config_file& cf )
{
	std::string rngname = cf.getValueSafe<std::string>( "random", "generator", "MUSIC" );
	
	RNG_plugin_creator *the_RNG_plugin_creator 
	= get_RNG_plugin_map()[ rngname ];
	
	if( !the_RNG_plugin_creator )
	{	
		std::cerr << " - Error: random number generator plug-in \'" << rngname << "\' not found." << std::endl;
		LOGERR("Invalid/Unregistered random number generator plug-in encountered : %s",rngname.c_str() );
		print_RNG_plugins();
		throw std::runtime_error("Unknown random number generator plug-in");
		
	}else
	{	
		std::cout << " - Selecting random number generator plug-in \'" << rngname << "\'..." << std::endl;
		LOGUSER("Selecting random number generator plug-in  : %s",rngname.c_str() );
	}
	
	RNG_plugin *the_RNG_plugin 
	= the_RNG_plugin_creator->create( cf );
	
	return the_RNG_plugin;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////


#if defined(FFTW3) && defined( SINGLE_PRECISION)
//#define fftw_complex fftwf_complex
typedef fftw_complex fftwf_complex;
#endif

template< typename T >
void rapid_proto_ngenic_rng( size_t res, long baseseed, random_numbers<T>& R )
{
  LOGUSER("Invoking the N-GenIC random number generator");

  unsigned *seedtable = new unsigned[res*res];

  gsl_rng *random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, baseseed);

  for(size_t i = 0; i < res / 2; i++)
    {
      size_t j;
      for(j = 0; j < i; j++)
	seedtable[i * res + j] = 0x7fffffff * gsl_rng_uniform(random_generator);
      for(j = 0; j < i + 1; j++)
	seedtable[j * res + i] = 0x7fffffff * gsl_rng_uniform(random_generator);
      for(j = 0; j < i; j++)
	seedtable[(res - 1 - i) * res + j] = 0x7fffffff * gsl_rng_uniform(random_generator);
      for(j = 0; j < i + 1; j++)
	seedtable[(res - 1 - j) * res + i] = 0x7fffffff * gsl_rng_uniform(random_generator);
      for(j = 0; j < i; j++)
	seedtable[i * res + (res - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);
      for(j = 0; j < i + 1; j++)
	seedtable[j * res + (res - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
      for(j = 0; j < i; j++)
	seedtable[(res - 1 - i) * res + (res - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);
      for(j = 0; j < i + 1; j++)
	seedtable[(res - 1 - j) * res + (res - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
    }
  
  fftw_real *rnoise = new fftw_real[ res*res*(res+2) ];
  fftw_complex *knoise = reinterpret_cast< fftw_complex* >( rnoise );

  double fnorm = 1./sqrt(res*res*res);

#warning need to check for race conditions below
  //#pragma omp parallel for
  for( size_t i = 0; i < res; i++ )
    {
      int ii = (int)res - (int)i;
      if(ii == (int)res)
	ii = 0;
      
      for( size_t j = 0; j < res; j++ )
	{
	  gsl_rng_set(random_generator, seedtable[i * res + j]);
	  
	  for( size_t k = 0; k < res / 2; k++ )
	    {
	      double phase = gsl_rng_uniform(random_generator) * 2 * M_PI;
	      double ampl;
	      do
		ampl = gsl_rng_uniform(random_generator);
	      while(ampl == 0);
	      
	      if(i == res / 2 || j == res / 2 || k == res / 2)
		continue;
	      if(i == 0 && j == 0 && k == 0)
		continue;

	      T rp = -sqrt(-log(ampl))*cos(phase) * fnorm;
	      T ip = -sqrt(-log(ampl))*sin(phase) * fnorm;

	      if(k > 0)
		{
		  RE(knoise[(i * res + j) * (res / 2 + 1) + k]) = rp;
		  IM(knoise[(i * res + j) * (res / 2 + 1) + k]) = ip;
		}
	      else	/* k=0 plane needs special treatment */
		{
		  if(i == 0)
		    {
		      if(j >= res / 2)
			continue;
		      else
			{
			  int jj = (int)res - (int)j;	/* note: j!=0 surely holds at this point */
			  
			  RE(knoise[(i * res + j) * (res / 2 + 1) + k]) = rp;
			  IM(knoise[(i * res + j) * (res / 2 + 1) + k]) = ip;
			  
			  RE(knoise[(i * res + jj) * (res / 2 + 1) + k]) = rp;
			  IM(knoise[(i * res + jj) * (res / 2 + 1) + k]) = -ip;
			}
		    }
		  else				
		    {
		      if(i >= res / 2) continue;
		      else
			{
			  int ii = (int)res - (int)i;
			  if(ii == (int)res)
			    ii = 0;
			  int jj = (int)res - (int)j;
			  if(jj == (int)res)
			    jj = 0;
									
			  RE(knoise[(i * res + j) * (res / 2 + 1) + k]) = rp;
			  IM(knoise[(i * res + j) * (res / 2 + 1) + k]) = ip;
			
			  if(ii >= 0 && ii < (int)res)
			    {
			      RE(knoise[(ii * res + jj) * (res / 2 + 1) + k]) = rp;
			      IM(knoise[(ii * res + jj) * (res / 2 + 1) + k]) = -ip;
			    }
			}
		    } 
		}
	    }
	}
      

    }

  delete[] seedtable;


  //... perform FT to real space

#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_plan plan = fftwf_plan_dft_c2r_3d( res, res, res, knoise, rnoise, FFTW_ESTIMATE); 
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
#else
  fftw_plan plan = fftw_plan_dft_c2r_3d( res, res, res, knoise, rnoise, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
#endif
#else
  rfftwnd_plan	plan = rfftw3d_create_plan( res, res, res, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
#ifndef SINGLETHREAD_FFTW		
  rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), plan, knoise, NULL );
#else
  rfftwnd_one_complex_to_real( plan, knoise, NULL );
#endif		
  rfftwnd_destroy_plan(plan);
#endif

  // copy to array that holds the random numbers

  #pragma omp parallel for
  for( int i=0; i<(int)res; ++i )
    for( size_t j=0; j<res; ++j )
      for( size_t k=0; k<res; ++k )
	R(i,j,k) = rnoise[ ((size_t)i*res+j)*res+k ];

  delete[] rnoise;

}

 

template< typename T >
random_numbers<T>::random_numbers( unsigned res, unsigned cubesize, long baseseed, int *x0, int *lx )
: res_( res ), cubesize_( cubesize ), ncubes_( 1 ), baseseed_( baseseed )
{
	LOGINFO("Generating random numbers (1) with seed %ld", baseseed);
    
	initialize();
	fill_subvolume( x0, lx );
}

template< typename T >
random_numbers<T>::random_numbers( unsigned res, unsigned cubesize, long baseseed, bool zeromean )
: res_( res ), cubesize_( cubesize ), ncubes_( 1 ), baseseed_( baseseed )
{
	LOGINFO("Generating random numbers (2) with seed %ld", baseseed);
    
	double mean = 0.0;

	bool musicnoise = true;
	if( !musicnoise )
	  cubesize_ = res_;
    
    if( !musicnoise )
        LOGERR("This currently breaks compatibility. Need to disable by hand! Make sure to not check into repo");

	initialize();

	if( musicnoise )
	  mean = fill_all();
	else
	  {
	    rnums_.push_back( new Meshvar<T>( res, 0, 0, 0 ) );
	    cubemap_[0] = 0; // create dummy map index
	    register_cube(0,0,0);
	    rapid_proto_ngenic_rng( res_, baseseed_, *this );
	  }

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
random_numbers<T>::random_numbers( unsigned res, std::string randfname, bool randsign )
: res_( res ), cubesize_( res ), ncubes_(1)
{
	rnums_.push_back( new Meshvar<T>( res, 0, 0, 0 ) );
    cubemap_[0] = 0; // create dummy map index
	
	std::ifstream ifs(randfname.c_str(), std::ios::binary);
	if( !ifs )
	{
        LOGERR("Could not open random number file \'%s\'!",randfname.c_str());
        throw std::runtime_error(std::string("Could not open random number file \'")+randfname+std::string("\'!"));
	}
    
	unsigned vartype;
	unsigned nx,ny,nz,blksz32;
    size_t blksz64;
	int iseed;
	//long seed;
    
    float sign4 = -1.0f;
    double sign8 = -1.0;
    
    int addrtype = 32;
    
    if( randsign ) // use grafic2 sign convention
    {
        sign4 = 1.0f;
        sign8 = 1.0;
    }
    
	//... read header and check if 32bit or 64bit block size .../
    ifs.read( reinterpret_cast<char*> (&blksz32), sizeof(int) );
	ifs.read( reinterpret_cast<char*> (&nx), sizeof(unsigned) );
	if( blksz32 != 4*sizeof(int) || nx != res_ )
    {
        addrtype = 64;
        
        ifs.seekg( 0 );
        ifs.read( reinterpret_cast<char*> (&blksz64), sizeof(size_t) );
        ifs.read( reinterpret_cast<char*> (&nx), sizeof(unsigned) );
        
        if( blksz64 != 4*sizeof(int) || nx != res_ )
            addrtype = -1;
    }
    ifs.seekg( 0 );
    
    if( addrtype < 0 )
		throw std::runtime_error("corrupt random number file");
	
    if( addrtype == 32 )
        ifs.read( reinterpret_cast<char*> (&blksz32), sizeof(int) );
    else
        ifs.read( reinterpret_cast<char*> (&blksz64), sizeof(size_t) );
    
	ifs.read( reinterpret_cast<char*> (&nx), sizeof(unsigned) );
	ifs.read( reinterpret_cast<char*> (&ny), sizeof(unsigned) );
	ifs.read( reinterpret_cast<char*> (&nz), sizeof(unsigned) );
	ifs.read( reinterpret_cast<char*> (&iseed), sizeof(int) );
	//seed = (long)iseed;
	
	if( nx!=res_ || ny!=res_ || nz!=res_ )
	{	
		char errmsg[128];
		sprintf(errmsg,"White noise file dimensions do not match level dimensions: %ux%ux%u vs. %u**3",nx,ny,nz,res_);
		throw std::runtime_error(errmsg);
		
	}
    
    if( addrtype == 32 )
        ifs.read( reinterpret_cast<char*> (&blksz32), sizeof(int) );
    else
        ifs.read( reinterpret_cast<char*> (&blksz64), sizeof(size_t) );
	
	//... read data ...//
    //check whether random numbers are single or double precision numbers
    if( addrtype == 32 )
    {
        ifs.read( reinterpret_cast<char*> (&blksz32), sizeof(int) );
        if( blksz32 == nx*ny*sizeof(float) )
            vartype = 4;
        else if( blksz32 == nx*ny*sizeof(double) )
            vartype = 8;
        else
            throw std::runtime_error("corrupt random number file");
    }else{
    
        ifs.read( reinterpret_cast<char*> (&blksz64), sizeof(size_t) );
        if( blksz64 == nx*ny*sizeof(float) )
            vartype = 4;
        else if( blksz64 == nx*ny*sizeof(double) )
            vartype = 8;
        else
            throw std::runtime_error("corrupt random number file");
    }
    
	//rewind to beginning of block
    if( addrtype == 32 )
        ifs.seekg(-sizeof(int),std::ios::cur);
    else
        ifs.seekg(-sizeof(size_t),std::ios::cur);
	
	std::vector<float> in_float;
	std::vector<double> in_double;
	
	LOGINFO("Random number file \'%s\'\n   contains %ld numbers. Reading...",randfname.c_str(),nx*ny*nz);
	
	double sum = 0.0, sum2 = 0.0;
	unsigned count = 0;
	
    //perform actual reading
	if( vartype == 4 )
	{
		for( int ii=0; ii<(int)nz; ++ii )
		{

			if( addrtype == 32 )
            {
                ifs.read( reinterpret_cast<char*> (&blksz32), sizeof(int) );
                if( blksz32 != nx*ny*sizeof(float) )
                    throw std::runtime_error("corrupt random number file");
            }
            else
            {
                ifs.read( reinterpret_cast<char*> (&blksz64), sizeof(size_t) );
                if( blksz64 != nx*ny*sizeof(float) )
                    throw std::runtime_error("corrupt random number file");
            }
			
			in_float.assign(nx*ny,0.0f);
			ifs.read( (char*)&in_float[0], nx*ny*sizeof(float) );
			
			for( int jj=0,q=0; jj<(int)ny; ++jj )
				for( int kk=0; kk<(int)nx; ++kk ){
					sum += in_float[q];
					sum2 += in_float[q]*in_float[q];
					++count;
					
					(*rnums_[0])(kk,jj,ii) = sign4 * in_float[q++];
				}
			
            if( addrtype == 32 )
            {
                ifs.read( reinterpret_cast<char*> (&blksz32), sizeof(int) );
                if( blksz32 != nx*ny*sizeof(float) )
                    throw std::runtime_error("corrupt random number file");
            }
            else
            {
                ifs.read( reinterpret_cast<char*> (&blksz64), sizeof(size_t) );
                if( blksz64 != nx*ny*sizeof(float) )
                    throw std::runtime_error("corrupt random number file");
            }
		}
	}
	else if( vartype == 8 )
	{
		for( int ii=0; ii<(int)nz; ++ii )
		{
			if( addrtype == 32 )
            {
                ifs.read( reinterpret_cast<char*> (&blksz32), sizeof(int) );
                if( blksz32 != nx*ny*sizeof(double) )
                    throw std::runtime_error("corrupt random number file");
            }
            else
            {
                ifs.read( reinterpret_cast<char*> (&blksz64), sizeof(size_t) );
                if( blksz64 != nx*ny*sizeof(double) )
                    throw std::runtime_error("corrupt random number file");
            }

			in_double.assign(nx*ny,0.0f);				
			ifs.read( (char*)&in_double[0], nx*ny*sizeof(double) );
			
			for( int jj=0,q=0; jj<(int)ny; ++jj )
				for( int kk=0; kk<(int)nx; ++kk )
				{
					sum += in_double[q];
					sum2 += in_double[q]*in_double[q];
					++count;
					(*rnums_[0])(kk,jj,ii) = sign8 * in_double[q++];
				}

            if( addrtype == 32 )
            {
                ifs.read( reinterpret_cast<char*> (&blksz32), sizeof(int) );
                if( blksz32 != nx*ny*sizeof(double) )
                    throw std::runtime_error("corrupt random number file");
            }
            else
            {
                ifs.read( reinterpret_cast<char*> (&blksz64), sizeof(size_t) );
                if( blksz64 != nx*ny*sizeof(double) )
                    throw std::runtime_error("corrupt random number file");
            }
		}
	}
	
	double mean, var;
	mean = sum/count;
	var = sum2/count-mean*mean;
	
	LOGINFO("Random numbers in file have \n     mean = %f and var = %f", mean, var);
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
		LOGINFO("Generating a coarse white noise field by k-space degrading");
		//... initialize properties of container		
		res_		= rc.res_/2;
		cubesize_	= res_;
		ncubes_		= 1;
		baseseed_	= -2;
		
		if( sizeof(fftw_real)!=sizeof(T) )
		{
            LOGERR("type mismatch with fftw_real in k-space averaging");
            throw std::runtime_error("type mismatch with fftw_real in k-space averaging");
		}
            
		fftw_real 
		*rfine = new fftw_real[(size_t)rc.res_*(size_t)rc.res_*2*((size_t)rc.res_/2+1)],
		*rcoarse = new fftw_real[(size_t)res_*(size_t)res_*2*((size_t)res_/2+1)];
		
		fftw_complex
		*ccoarse = reinterpret_cast<fftw_complex*> (rcoarse),
		*cfine = reinterpret_cast<fftw_complex*> (rfine);
		
		int nx(rc.res_), ny(rc.res_), nz(rc.res_), nxc(res_), nyc(res_), nzc(res_);
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
		fftwf_plan
		pf = fftwf_plan_dft_r2c_3d(nx, ny, nz, rfine, cfine, FFTW_ESTIMATE),
		ipc= fftwf_plan_dft_c2r_3d(nxc, nyc, nzc, ccoarse, rcoarse, FFTW_ESTIMATE);
	#else
		fftw_plan
		pf = fftw_plan_dft_r2c_3d(nx, ny, nz, rfine, cfine, FFTW_ESTIMATE),
		ipc= fftw_plan_dft_c2r_3d(nxc, nyc, nzc, ccoarse, rcoarse, FFTW_ESTIMATE);
	#endif
		
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
	#ifdef SINGLE_PRECISION
		fftwf_execute( pf );
	#else
		fftw_execute( pf );
	#endif
#else
#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), pf, rfine, NULL );
#else
		rfftwnd_one_real_to_complex( pf, rfine, NULL );
#endif
#endif
		
		double fftnorm = 1.0/((double)nxc*(double)nyc*(double)nzc);
		
#pragma omp parallel for
		for( int i=0; i<nxc; i++ )
			for( int j=0; j<nyc; j++ )
				for( int k=0; k<nzc/2+1; k++ )
				{
					int ii(i),jj(j),kk(k);
					
					if( i > nxc/2 ) ii += nx/2;
					if( j > nyc/2 ) jj += ny/2;
					
					size_t qc,qf;

					double kx = (i <= (int)nxc/2)? (double)i : (double)(i-(int)nxc);
					double ky = (j <= (int)nyc/2)? (double)j : (double)(j-(int)nyc);
					double kz = (k <= (int)nzc/2)? (double)k : (double)(k-(int)nzc);
					
					
					qc = ((size_t)i*nyc+(size_t)j)*(nzc/2+1)+(size_t)k;
					qf = ((size_t)ii*ny+(size_t)jj)*(nz/2+1)+(size_t)kk;

					std::complex<double> val_fine(RE(cfine[qf]),IM(cfine[qf]));
					double phase =  (kx/nxc + ky/nyc + kz/nzc) * 0.5 * M_PI;
					std::complex<double> val_phas( cos(phase), sin(phase) );

					val_fine *= val_phas * fftnorm/sqrt(8.0);
					
					RE(ccoarse[qc]) = val_fine.real();
					IM(ccoarse[qc]) = val_fine.imag();
				}
		
		delete[] rfine;
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
		fftwf_execute( ipc );
	#else
		fftw_execute( ipc );
	#endif
#else
#ifndef SINGLETHREAD_FFTW		
		rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), ipc, ccoarse, NULL );
#else
		rfftwnd_one_complex_to_real( ipc, ccoarse, NULL );
#endif
#endif
		rnums_.push_back( new Meshvar<T>( res_, 0, 0, 0 ) );
        cubemap_[0] = 0; // map all to single array
		
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
	else
	{
		LOGINFO("Generating a coarse white noise field by averaging");
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
            cubemap_[0] = 0; // map all to single array
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
            cubemap_[0] = 0;
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
	
	LOGINFO("Restricted random numbers have\n       mean = %f, var = %f", rmean, rvar);
}


template< typename T >
random_numbers<T>::random_numbers( random_numbers<T>& rc, unsigned cubesize, long baseseed, 
				   bool kspace,bool isolated, int *x0_, int *lx_, bool zeromean )
: res_( 2*rc.res_ ), cubesize_( cubesize ), ncubes_( 1 ), baseseed_( baseseed )
{
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
		
	  LOGINFO("Generating a constrained random number set with seed %ld\n    using coarse mode replacement...",baseseed);
	  assert(lx[0]%2==0 && lx[1]%2==0 && lx[2]%2==0);
	  size_t nx=lx[0], ny=lx[1], nz=lx[2],
	    nxc=lx[0]/2, nyc=lx[1]/2, nzc=lx[2]/2;
	  
		
	  fftw_real *rfine = new fftw_real[nx*ny*(nz+2l)];
	  fftw_complex *cfine = reinterpret_cast<fftw_complex*> (rfine);
		
#ifdef FFTW3
#ifdef SINGLE_PRECISION
	  fftwf_plan
	    pf  = fftwf_plan_dft_r2c_3d( nx, ny, nz, rfine, cfine, FFTW_ESTIMATE),
	    ipf	= fftwf_plan_dft_c2r_3d( nx, ny, nz, cfine, rfine, FFTW_ESTIMATE);
#else
	  fftw_plan
	    pf  = fftw_plan_dft_r2c_3d( nx, ny, nz, rfine, cfine, FFTW_ESTIMATE),
	    ipf	= fftw_plan_dft_c2r_3d( nx, ny, nz, cfine, rfine, FFTW_ESTIMATE);
#endif
#else
	  rfftwnd_plan 
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
	  //this->free_all_mem();	// temporarily free memory, allocate again later
		
		
		
	  fftw_real *rcoarse = new fftw_real[nxc*nyc*(nzc+2)];
	  fftw_complex *ccoarse = reinterpret_cast<fftw_complex*> (rcoarse);
		
#ifdef FFTW3
#ifdef SINGLE_PRECISION
	  fftwf_plan pc  = fftwf_plan_dft_r2c_3d( nxc, nyc, nzc, rcoarse, ccoarse, FFTW_ESTIMATE);
#else
	  fftw_plan pc  = fftw_plan_dft_r2c_3d( nxc, nyc, nzc, rcoarse, ccoarse, FFTW_ESTIMATE);
#endif
#else
	  rfftwnd_plan pc	= rfftw3d_create_plan( nxc, nyc, nzc, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE);
#endif		
	  
#pragma omp parallel for
	  for( int i=0; i<(int)nxc; i++ )
	    for( int j=0; j<(int)nyc; j++ )
	      for( int k=0; k<(int)nzc; k++ )
		{
		  size_t q = ((size_t)i*(size_t)nyc+(size_t)j)*(size_t)(nzc+2)+(size_t)k;
		  rcoarse[q] = rc(x0[0]/2+i,x0[1]/2+j,x0[2]/2+k);
		}
#ifdef FFTW3
#ifdef SINGLE_PRECISION
	  fftwf_execute( pc );
	  fftwf_execute( pf );
#else
	  fftw_execute( pc );
	  fftw_execute( pf );	
#endif
#else
#ifndef SINGLETHREAD_FFTW		
	  rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), pc, rcoarse, NULL );
	  rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), pf, rfine, NULL );
#else
	  rfftwnd_one_real_to_complex( pc, rcoarse, NULL );
	  rfftwnd_one_real_to_complex( pf, rfine, NULL );
#endif
#endif
	  
	  double fftnorm = 1.0/((double)nx*(double)ny*(double)nz);
	  double sqrt8 = sqrt(8.0);
        double phasefac = -0.5;//-1.0;//-0.125;

	  //if( isolated ) phasefac *= 1.5;

        // embedding of coarse white noise by fourier interpolation
     
#if 1
#pragma omp parallel for
        for( int i=0; i<(int)nxc; i++ )
            for( int j=0; j<(int)nyc; j++ )
                for( int k=0; k<(int)nzc/2+1; k++ )
                {
                    int ii(i),jj(j),kk(k);
                    
                    //if( i==(int)nxc/2 ) continue;
                    //if( j==(int)nyc/2 ) continue;
                    
                    if( i > (int)nxc/2 ) ii += (int)nx/2;
                    if( j > (int)nyc/2 ) jj += (int)ny/2;
                    
                    size_t qc,qf;
                    
                    double kx = (i <= (int)nxc/2)? (double)i : (double)(i-(int)nxc);
                    double ky = (j <= (int)nyc/2)? (double)j : (double)(j-(int)nyc);
                    double kz = (k <= (int)nzc/2)? (double)k : (double)(k-(int)nzc);
                    
                    
                    qc = ((size_t)i*nyc+(size_t)j)*(nzc/2+1)+(size_t)k;
                    qf = ((size_t)ii*ny+(size_t)jj)*(nz/2+1)+(size_t)kk;
                    
                    std::complex<double> val(RE(ccoarse[qc]),IM(ccoarse[qc]));
                    double phase =  (kx/nxc + ky/nyc + kz/nzc) * phasefac * M_PI;
                    
                    std::complex<double> val_phas( cos(phase), sin(phase) );
                    
                    val *= val_phas * sqrt8;
                    
                    if( i!=(int)nxc/2 && j!=(int)nyc/2 && k!=(int)nzc/2 )
                    {
                        RE(cfine[qf]) = val.real();
                        IM(cfine[qf]) = val.imag();
                    }
                    else
                    {
                        //RE(cfine[qf]) = val.real();
                        //IM(cfine[qf]) = 0.0;
                    }
                }
        
#else
        
        // 0 0
#pragma omp parallel for
	for( int i=0; i<(int)nxc/2+1; i++ )
	  for( int j=0; j<(int)nyc/2+1; j++ )
	    for( int k=0; k<(int)nzc/2+1; k++ )
	      {
		int ii(i),jj(j),kk(k);
		size_t qc,qf;
		qc = ((size_t)i*(size_t)nyc+(size_t)j)*(nzc/2+1)+(size_t)k;
		qf = ((size_t)ii*(size_t)ny+(size_t)jj)*(nz/2+1)+(size_t)kk;

		double kx = (i <= (int)nxc/2)? (double)i : (double)(i-(int)nxc);
		double ky = (j <= (int)nyc/2)? (double)j : (double)(j-(int)nyc);
		double kz = (k <= (int)nzc/2)? (double)k : (double)(k-(int)nzc);
					
		double phase =  phasefac * (kx/nxc + ky/nyc + kz/nzc) * M_PI;
		std::complex<double> val_phas( cos(phase), sin(phase) );

		std::complex<double> val(RE(ccoarse[qc]),IM(ccoarse[qc]));
		val *= sqrt8 * val_phas;
                
		RE(cfine[qf]) = val.real();
		IM(cfine[qf]) = val.imag();

		//if( k==0 & (i==(int)nxc/2 || j==(int)nyc/2) )
		//  IM(cfine[qf]) *= -1.0;
	      }
        // 1 0
#pragma omp parallel for
	for( int i=nxc/2; i<(int)nxc; i++ )
	  for( int j=0; j<(int)nyc/2+1; j++ )
	    for( int k=0; k<(int)nzc/2+1; k++ )
	      {
		int ii(i+nx/2),jj(j),kk(k);
		size_t qc,qf;
		qc = ((size_t)i*(size_t)nyc+(size_t)j)*(nzc/2+1)+(size_t)k;
		qf = ((size_t)ii*(size_t)ny+(size_t)jj)*(nz/2+1)+(size_t)kk;
                
		double kx = (i <= (int)nxc/2)? (double)i : (double)(i-(int)nxc);
		double ky = (j <= (int)nyc/2)? (double)j : (double)(j-(int)nyc);
		double kz = (k <= (int)nzc/2)? (double)k : (double)(k-(int)nzc);
					
		double phase =  phasefac * (kx/nxc + ky/nyc + kz/nzc) * M_PI;
		std::complex<double> val_phas( cos(phase), sin(phase) );

		std::complex<double> val(RE(ccoarse[qc]),IM(ccoarse[qc]));
		val *= sqrt8 * val_phas;
                
		RE(cfine[qf]) = val.real();
		IM(cfine[qf]) = val.imag();
                
		//if( k==0 & (i==(int)nxc/2 || j==(int)nyc/2) )
		//IM(cfine[qf]) *= -1.0;
	      }
        // 0 1
#pragma omp parallel for
	for( int i=0; i<(int)nxc/2+1; i++ )
	  for( int j=nyc/2; j<(int)nyc; j++ )
	    for( int k=0; k<(int)nzc/2+1; k++ )
	      {
		int ii(i),jj(j+ny/2),kk(k);
		size_t qc,qf;
		qc = ((size_t)i*(size_t)nyc+(size_t)j)*(nzc/2+1)+(size_t)k;
		qf = ((size_t)ii*(size_t)ny+(size_t)jj)*(nz/2+1)+(size_t)kk;
                
		double kx = (i <= (int)nxc/2)? (double)i : (double)(i-(int)nxc);
		double ky = (j <= (int)nyc/2)? (double)j : (double)(j-(int)nyc);
		double kz = (k <= (int)nzc/2)? (double)k : (double)(k-(int)nzc);
					
		double phase =  phasefac * (kx/nxc + ky/nyc + kz/nzc)  * M_PI;
		std::complex<double> val_phas( cos(phase), sin(phase) );

		std::complex<double> val(RE(ccoarse[qc]),IM(ccoarse[qc]));
		val *= sqrt8 * val_phas;
                
		RE(cfine[qf]) = val.real();
		IM(cfine[qf]) = val.imag();
                
		//if( k==0 && (i==(int)nxc/2 || j==(int)nyc/2) )
		//  IM(cfine[qf]) *= -1.0;
	      }
        
        // 1 1
#pragma omp parallel for
	for( int i=nxc/2; i<(int)nxc; i++ )
	  for( int j=nyc/2; j<(int)nyc; j++ )
	    for( int k=0; k<(int)nzc/2+1; k++ )
	      {
		int ii(i+nx/2),jj(j+ny/2),kk(k);
		size_t qc,qf;
		qc = ((size_t)i*(size_t)nyc+(size_t)j)*(nzc/2+1)+(size_t)k;
		qf = ((size_t)ii*(size_t)ny+(size_t)jj)*(nz/2+1)+(size_t)kk;
                
		double kx = (i <= (int)nxc/2)? (double)i : (double)(i-(int)nxc);
		double ky = (j <= (int)nyc/2)? (double)j : (double)(j-(int)nyc);
		double kz = (k <= (int)nzc/2)? (double)k : (double)(k-(int)nzc);
					
		double phase =  phasefac * (kx/nxc + ky/nyc + kz/nzc) * M_PI;
		std::complex<double> val_phas( cos(phase), sin(phase) );

		std::complex<double> val(RE(ccoarse[qc]),IM(ccoarse[qc]));
		val *= sqrt8 * val_phas;
                
		RE(cfine[qf]) = val.real();
		IM(cfine[qf]) = val.imag();
	      }
#endif
        
	delete[] rcoarse;
	
#pragma omp parallel for
	for( int i=0; i<(int)nx; i++ )
	  for( int j=0; j<(int)ny; j++ )
	    for( int k=0; k<(int)nz/2+1; k++ )
	      {
		size_t q = ((size_t)i*ny+(size_t)j)*(nz/2+1)+(size_t)k;
		
		RE(cfine[q]) *= fftnorm;
		IM(cfine[q]) *= fftnorm;
	      }

		
		
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
		fftwf_execute( ipf );
	#else
		fftw_execute( ipf );
	#endif
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
					(*this)(x0[0]+i,x0[1]+j,x0[2]+k,false) = rfine[q];
				}
		
		delete[] rfine;
		
#ifdef FFTW3
	#ifdef SINGLE_PRECISION
		fftwf_destroy_plan(pf);
		fftwf_destroy_plan(pc);
		fftwf_destroy_plan(ipf);
	#else
		fftw_destroy_plan(pf);
		fftw_destroy_plan(pc);
		fftw_destroy_plan(ipf);
	#endif
#else
		fftwnd_destroy_plan(pf);
		fftwnd_destroy_plan(pc);
		fftwnd_destroy_plan(ipf);
#endif
		
	}
	else
	{
		LOGINFO("Generating a constrained random number set with seed %ld\n    using Hoffman-Ribak constraints...", baseseed);
		
		double fac = 1.0/sqrt(8.0);//1./sqrt(8.0);
		
		for( int i=x0[0],ii=x0[0]/2; i<x0[0]+lx[0]; i+=2,++ii )
			for( int j=x0[1],jj=x0[1]/2; j<x0[1]+lx[1]; j+=2,++jj )
				for( int k=x0[2],kk=x0[2]/2; k<x0[2]+lx[2]; k+=2,++kk )
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
void random_numbers<T>::register_cube( int i, int j, int k)
{
	i = (i+ncubes_)%ncubes_;
	j = (j+ncubes_)%ncubes_;
	k = (k+ncubes_)%ncubes_;
    size_t icube = ((size_t)i*ncubes_+(size_t)j)*ncubes_+(size_t)k;
    
    cubemap_iterator it = cubemap_.find( icube );
    
    if( it == cubemap_.end() )
    {
        rnums_.push_back( NULL );
        cubemap_[icube] = rnums_.size()-1;
#ifdef DEBUG
        LOGDEBUG("registering new cube %d,%d,%d . ID = %ld, memloc = %ld",i,j,k,icube,cubemap_[icube]);
#endif
    }
}


template< typename T >
double random_numbers<T>::fill_cube( int i, int j, int k)
{
	
	gsl_rng	*RNG = gsl_rng_alloc( gsl_rng_mt19937 );
	
	i = (i+ncubes_)%ncubes_;
	j = (j+ncubes_)%ncubes_;
	k = (k+ncubes_)%ncubes_;
	
	size_t icube = ((size_t)i*ncubes_+(size_t)j)*ncubes_+(size_t)k;
	long cubeseed = baseseed_+icube; //... each cube gets its unique seed
	
	gsl_rng_set( RNG, cubeseed );
    
    cubemap_iterator it = cubemap_.find( icube );
    
    if( it == cubemap_.end() )
    {
        LOGERR("Attempt to access non-registered random number cube!");
        throw std::runtime_error("Attempt to access non-registered random number cube!");
    }
    
    size_t cubeidx = it->second;
	
	if( rnums_[cubeidx] == NULL )
        rnums_[cubeidx] =  new Meshvar<T>( cubesize_, 0, 0, 0 );

	double mean = 0.0;
	
	for( int ii=0; ii<(int)cubesize_; ++ii )
		for( int jj=0; jj<(int)cubesize_; ++jj )
			for( int kk=0; kk<(int)cubesize_; ++kk )
			{	
				(*rnums_[cubeidx])(ii,jj,kk) = gsl_ran_ugaussian_ratio_method( RNG );
				mean += (*rnums_[cubeidx])(ii,jj,kk);
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
	
	size_t icube = ((size_t)i*ncubes_+(size_t)j)*ncubes_+(size_t)k;
    
    cubemap_iterator it = cubemap_.find( icube );
    
    if( it == cubemap_.end() )
    {
        LOGERR("Attempt to access unallocated RND cube %d,%d,%d in random_numbers::subtract_from_cube",i,j,k);
        throw std::runtime_error("Attempt to access unallocated RND cube in random_numbers::subtract_from_cube");
    }
    
    size_t cubeidx = it->second;
	
	for( int ii=0; ii<(int)cubesize_; ++ii )
		for( int jj=0; jj<(int)cubesize_; ++jj )
			for( int kk=0; kk<(int)cubesize_; ++kk )
				(*rnums_[cubeidx])(ii,jj,kk) -= val;
	
}

template< typename T >
void random_numbers<T>::free_cube( int i, int j, int k ) 
{
	
	i = (i+ncubes_)%ncubes_;
	j = (j+ncubes_)%ncubes_;
	k = (k+ncubes_)%ncubes_;
    
	size_t icube = ((size_t)i*(size_t)ncubes_+(size_t)j)*(size_t)ncubes_+(size_t)k;
    
    cubemap_iterator it = cubemap_.find( icube );
    
    if( it == cubemap_.end() )
    {
        LOGERR("Attempt to access unallocated RND cube %d,%d,%d in random_numbers::free_cube",i,j,k);
        throw std::runtime_error("Attempt to access unallocated RND cube in random_numbers::free_cube");
    }
    
    size_t cubeidx = it->second;
	
	if( rnums_[cubeidx] != NULL )
	{
		delete rnums_[cubeidx];
		rnums_[cubeidx] = NULL;
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
	
	LOGINFO("Generating random numbers w/ sample cube size of %d", cubesize_ );
}

template< typename T >
double random_numbers<T>::fill_subvolume( int *i0, int *n )
{
	int i0cube[3], ncube[3];
	
	i0cube[0] = (int)((double)(res_+i0[0])/cubesize_);
	i0cube[1] = (int)((double)(res_+i0[1])/cubesize_);
	i0cube[2] = (int)((double)(res_+i0[2])/cubesize_);
	
	ncube[0] = (int)(n[0]/cubesize_) + 2;
	ncube[1] = (int)(n[1]/cubesize_) + 2;
	ncube[2] = (int)(n[2]/cubesize_) + 2;
    
#ifdef DEBUG
    LOGDEBUG("random numbers needed for region %d,%d,%d ..+ %d,%d,%d",i0[0],i0[1],i0[2],n[0],n[1],n[2]);
    LOGDEBUG("filling cubes %d,%d,%d ..+ %d,%d,%d",i0cube[0],i0cube[1],i0cube[2],ncube[0],ncube[1],ncube[2]);
#endif

	double mean = 0.0;
    
    for( int i=i0cube[0]; i<i0cube[0]+ncube[0]; ++i )
		for( int j=i0cube[1]; j<i0cube[1]+ncube[1]; ++j )
			for( int k=i0cube[2]; k<i0cube[2]+ncube[2]; ++k )
			{
                int ii(i),jj(j),kk(k);
				
				ii = (ii+ncubes_)%ncubes_;
				jj = (jj+ncubes_)%ncubes_;
				kk = (kk+ncubes_)%ncubes_;
                
                register_cube( ii,jj,kk );
            }
	
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
    
    for( int i=0; i<(int)ncubes_; ++i )
		for( int j=0; j<(int)ncubes_; ++j )
			for( int k=0; k<(int)ncubes_; ++k )
			{
				int ii(i),jj(j),kk(k);

                                ii = (ii+ncubes_)%ncubes_;
                                jj = (jj+ncubes_)%ncubes_;
                                kk = (kk+ncubes_)%ncubes_;

                register_cube(ii,jj,kk);
            }
	
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
	
	return sum/(ncubes_*ncubes_*ncubes_);
}


template< typename T >
void random_numbers<T>:: print_allocated( void )
{
	unsigned ncount = 0, ntot = rnums_.size();
	for( size_t i=0; i<rnums_.size(); ++i )
		if( rnums_[i]!=NULL ) ncount++;
	
	LOGINFO(" -> %d of %d random number cubes currently allocated",ncount,ntot);
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
	disk_cached_	= pcf_->getValueSafe<bool>("random","disk_cached",true);
	restart_		= pcf_->getValueSafe<bool>("random","restart",false);
	
	mem_cache_.assign(levelmax_-levelmin_+1, (std::vector<T>*)NULL);
	
	if( restart_ && !disk_cached_ )
	{
        LOGERR("Cannot restart from mem cached random numbers.");
        throw std::runtime_error("Cannot restart from mem cached random numbers.");
    }
	////disk_cached_ = false;
	
	//... determine seed/white noise file data to be applied
	parse_rand_parameters();
	
	if( !restart_ )
	{
		//... compute the actual random numbers
		compute_random_numbers();		
	}
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
			//			unlink(fname);
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
				//rngseeds_.push_back( -abs((unsigned)(ltemp-i)%123+(unsigned)(ltemp+827342523521*i)%123456789) );
                rngseeds_.push_back( -abs((long)(ltemp-i)%123+(long)(ltemp+7342523521*i)%123456789) );
			else
				rngseeds_.push_back( ltemp );
		}else{
			rngfnames_.push_back( tempstr );
			rngseeds_.push_back(-1);
			LOGINFO("Random numbers for level %3d will be read from file.",i);
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
        {
            LOGERR("White noise file mismatch. This should not happen. Notify a developer!");
            throw std::runtime_error("White noise file mismatch. This should not happen. Notify a developer!");
        }
        int nxd(nxf/2),nyd(nyf/2),nzd(nzf/2);
        std::vector<T> deg_rand( (size_t)nxd*(size_t)nyd*(size_t)nzd, 0.0 );
        double fac = 1.0/sqrt(8.0);
        
        for( int i=0, ic=0; i<nxf; i+=2, ic++ )
        {
            std::vector<T> fine_rand( 2*nyf*nzf, 0.0 );
            iffine.read( reinterpret_cast<char*> (&fine_rand[0]), 2*nyf*nzf*sizeof(T) );
            
#pragma omp parallel for
            for( int j=0; j<nyf; j+=2 )
                for( int k=0; k<nzf; k+=2 )
                {
                    int jc = j/2, kc = k/2;
                    //size_t qc = (((size_t)i/2)*(size_t)nyd+((size_t)j/2))*(size_t)nzd+((size_t)k/2);
                    size_t qc = ((size_t)(ic*nyd+jc))*(size_t)nzd+(size_t)kc;
                    
                    size_t qf[8];
                    qf[0] = (0*(size_t)nyf+(size_t)j+0)*(size_t)nzf+(size_t)k+0;
                    qf[1] = (0*(size_t)nyf+(size_t)j+0)*(size_t)nzf+(size_t)k+1;
                    qf[2] = (0*(size_t)nyf+(size_t)j+1)*(size_t)nzf+(size_t)k+0;
                    qf[3] = (0*(size_t)nyf+(size_t)j+1)*(size_t)nzf+(size_t)k+1;
                    qf[4] = (1*(size_t)nyf+(size_t)j+0)*(size_t)nzf+(size_t)k+0;
                    qf[5] = (1*(size_t)nyf+(size_t)j+0)*(size_t)nzf+(size_t)k+1;
                    qf[6] = (1*(size_t)nyf+(size_t)j+1)*(size_t)nzf+(size_t)k+0;
                    qf[7] = (1*(size_t)nyf+(size_t)j+1)*(size_t)nzf+(size_t)k+1;
                    
                    double d = 0.0;
                    for( int q=0; q<8; ++q )
                        d += fac*fine_rand[qf[q]];
                    
                    //deg_rand[qc] += d;
                    deg_rand[qc] = d;
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
                    
                    size_t qc = (((size_t)i+(size_t)di)*(size_t)nyc+((size_t)j+(size_t)dj))*(size_t)nzc+(size_t)(k+dk);
                    size_t qcd = (size_t)(i*nyd+j)*(size_t)nzd+(size_t)k;
                    
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
                    
                    size_t qf[8];
                    qf[0] = (size_t)((2*i+0)*nyf+2*j+0)*(size_t)nzf+(size_t)(2*k+0);
                    qf[1] = (size_t)((2*i+0)*nyf+2*j+0)*(size_t)nzf+(size_t)(2*k+1);
                    qf[2] = (size_t)((2*i+0)*nyf+2*j+1)*(size_t)nzf+(size_t)(2*k+0);
                    qf[3] = (size_t)((2*i+0)*nyf+2*j+1)*(size_t)nzf+(size_t)(2*k+1);
                    qf[4] = (size_t)((2*i+1)*nyf+2*j+0)*(size_t)nzf+(size_t)(2*k+0);
                    qf[5] = (size_t)((2*i+1)*nyf+2*j+0)*(size_t)nzf+(size_t)(2*k+1);
                    qf[6] = (size_t)((2*i+1)*nyf+2*j+1)*(size_t)nzf+(size_t)(2*k+0);
                    qf[7] = (size_t)((2*i+1)*nyf+2*j+1)*(size_t)nzf+(size_t)(2*k+1);
                    
                    double finesum = 0.0;
                    for( int q=0; q<8; ++q )
                        finesum += fac*(*mem_cache_[ifine-levelmin_])[qf[q]];
                    
                    size_t qc = ((size_t)(i+di)*nyc+(size_t)(j+dj))*(size_t)nzc+(size_t)(k+dk);
                    
                    (*mem_cache_[icoarse-levelmin_])[qc] = finesum;
                }
    }
    
    
}

template< typename rng, typename T >
void random_number_generator<rng,T>::compute_random_numbers( void )
{
	bool kavg = pcf_->getValueSafe<bool>("random","kaveraging",true);
	bool rndsign = pcf_->getValueSafe<bool>("random","grafic_sign",false);
	bool brealspace_tf = !pcf_->getValue<bool>("setup","kspace_TF");
	
	std::vector< rng* > randc(std::max(levelmax_,levelmin_seed_)+1,(rng*)NULL);
	
	//--- FILL ALL WHITE NOISE ARRAYS FOR WHICH WE NEED THE FULL FIELD ---//
	
	//... seeds are given for a level coarser than levelmin
	if( levelmin_seed_ < levelmin_ )
	{
		if( rngfnames_[levelmin_seed_].size() > 0 )
			randc[levelmin_seed_]
			= new rng( 1<<levelmin_seed_, rngfnames_[levelmin_seed_], rndsign );
		else
			randc[levelmin_seed_]
			= new rng( 1<<levelmin_seed_, ran_cube_size_, rngseeds_[levelmin_seed_], true );
		
		for( int i=levelmin_seed_+1; i<=levelmin_; ++i )
		{
//#warning add possibility to read noise from file also here!
            
            if( rngfnames_[i].size() > 0 )
                LOGINFO("Warning: Cannot use filenames for higher levels currently! Ignoring!");
            
			randc[i] = new rng( *randc[i-1], ran_cube_size_, rngseeds_[i], kavg );
			delete randc[i-1];
			randc[i-1] = NULL;
		}
	}
	
	//... seeds are given for a level finer than levelmin, obtain by averaging
	if( levelmin_seed_ > levelmin_ )
	{
	  if( rngfnames_[levelmin_seed_].size() > 0 )
            randc[levelmin_seed_] = new rng( 1<<levelmin_seed_, rngfnames_[levelmin_seed_], rndsign );
	  else
            randc[levelmin_seed_] = new rng( 1<<levelmin_seed_, ran_cube_size_, rngseeds_[levelmin_seed_], true );//, x0, lx );
		
	  for( int ilevel = levelmin_seed_-1; ilevel >= (int)levelmin_; --ilevel ){
	    if( rngseeds_[ilevel-levelmin_] > 0 )
	      LOGINFO("Warning: random seed for level %d will be ignored.\n" \
		      "            consistency requires that it is obtained by restriction from level %d", ilevel, levelmin_seed_ );
	    
	    //if( brealspace_tf && ilevel < levelmax_ )
	    //  randc[ilevel] = new rng( *randc[ilevel+1], false );
	    //else // do k-space averaging
	    randc[ilevel] = new rng( *randc[ilevel+1], kavg );
	    
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
	{
        if( rngfnames_[levelmin_].size() > 0 )
			randc[levelmin_] = new rng( 1<<levelmin_, rngfnames_[levelmin_], rndsign );
		else
			randc[levelmin_] = new rng( 1<<levelmin_, ran_cube_size_, rngseeds_[levelmin_], true );
	}
    
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
		  randc[ilevel] = new rng( *randc[ilevel-1], ran_cube_size_, rngseeds_[ilevel], kavg, ilevel==levelmin_+1, x0, lx );
		delete randc[ilevel-1];
		randc[ilevel-1] = NULL;
		
		//... apply constraints to this level, if any
		//if( ilevel == levelmax_ )
		//constraints.apply( ilevel, x0, lx, randc[ilevel] );
		
		//... store numbers
		store_rnd( ilevel, randc[ilevel] );
	}
	
	delete randc[levelmax_];
	randc[levelmax_] = NULL;
    
	//... make sure that the coarse grid contains oct averages where it overlaps with a fine grid
	//... this also ensures that constraints enforced on fine grids are carried to the coarser grids
	if( brealspace_tf )
	  {
	    for( int ilevel=levelmax_; ilevel>levelmin_; --ilevel )
	      correct_avg( ilevel-1, ilevel );
	  }

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
	
	
	bool grafic_out = false;
	
	
	if( grafic_out )
	{
		std::vector<float> data;
		if( ilevel == levelmin_ )
		{
			int N = 1<<levelmin_;
			int i0,j0,k0;
			i0 = -lfac*shift[0];
			j0 = -lfac*shift[1];
			k0 = -lfac*shift[2];
			
			char fname[128];
			sprintf(fname,"grafic_wnoise_%04d.bin",ilevel);
			
			LOGUSER("Storing white noise field for grafic in file \'%s\'...", fname );
			
			std::ofstream ofs(fname,std::ios::binary|std::ios::trunc);
			data.assign( N*N, 0.0 );
			
			int blksize = 4*sizeof(int);
			int iseed = 0;
			
			ofs.write(reinterpret_cast<char*> (&blksize), sizeof(int) );
			ofs.write(reinterpret_cast<char*> (&N), sizeof(int) );
			ofs.write(reinterpret_cast<char*> (&N), sizeof(int) );
			ofs.write(reinterpret_cast<char*> (&N), sizeof(int) );
			ofs.write(reinterpret_cast<char*> (&iseed), sizeof(int) );
			ofs.write(reinterpret_cast<char*> (&blksize), sizeof(int) );
			
			for( int k=0; k<N; ++k )
			{	
				#pragma omp parallel for
				for( int j=0; j<N; ++j )
					for( int i=0; i<N; ++i )
						data[j*N+i] = -(*prng)(i+i0,j+j0,k+k0);
			
				blksize = N*N*sizeof(float);
				ofs.write(reinterpret_cast<char*> (&blksize), sizeof(int) );
				ofs.write(reinterpret_cast<char*> (&data[0]), N*N*sizeof(float) );
				ofs.write(reinterpret_cast<char*> (&blksize), sizeof(int) );
			}
			
			ofs.close();
			
		}
		else {
			
			int nx,ny,nz;
			int i0,j0,k0;
			
			nx = prefh_->size(ilevel, 0);
			ny = prefh_->size(ilevel, 1);
			nz = prefh_->size(ilevel, 2);
			i0 = prefh_->offset_abs(ilevel, 0) - lfac*shift[0];
			j0 = prefh_->offset_abs(ilevel, 1) - lfac*shift[1];
			k0 = prefh_->offset_abs(ilevel, 2) - lfac*shift[2];
			
			char fname[128];
			sprintf(fname,"grafic_wnoise_%04d.bin",ilevel);
			
			LOGUSER("Storing white noise field for grafic in file \'%s\'...", fname );
            LOGDEBUG("(%d,%d,%d) -- (%d,%d,%d) -- lfac = %d",nx,ny,nz,i0,j0,k0,lfac);
			
			std::ofstream ofs(fname,std::ios::binary|std::ios::trunc);
			data.assign( nx*ny, 0.0 );
			
			int blksize = 4*sizeof(int);
			int iseed = 0;
			
			ofs.write(reinterpret_cast<char*> (&blksize), sizeof(int) );
			ofs.write( reinterpret_cast<char*> (&nz), sizeof(unsigned) );
			ofs.write( reinterpret_cast<char*> (&ny), sizeof(unsigned) );
			ofs.write( reinterpret_cast<char*> (&nx), sizeof(unsigned) );
			ofs.write(reinterpret_cast<char*> (&iseed), sizeof(int) );
			ofs.write(reinterpret_cast<char*> (&blksize), sizeof(int) );
			
			for( int k=0; k<nz; ++k )
			{	
				#pragma omp parallel for
				for( int j=0; j<ny; ++j )
					for( int i=0; i<nx; ++i )
						data[j*nx+i] = -(*prng)(i+i0,j+j0,k+k0);
				
				blksize = nx*ny*sizeof(float);
				ofs.write(reinterpret_cast<char*> (&blksize), sizeof(int) );
				ofs.write(reinterpret_cast<char*> (&data[0]), nx*ny*sizeof(float) );
				ofs.write(reinterpret_cast<char*> (&blksize), sizeof(int) );
			}
			ofs.close();
			
		}

	}
	
	
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
			j0 = prefh_->offset_abs(ilevel, 1) - lfac*shift[1] - ny/4; // was nx/4
			k0 = prefh_->offset_abs(ilevel, 2) - lfac*shift[2] - nz/4; // was nx/4
			
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
			j0 = prefh_->offset_abs(ilevel, 1) - lfac*shift[1] - ny/4; // was nx/4
			k0 = prefh_->offset_abs(ilevel, 2) - lfac*shift[2] - nz/4; // was nx/4
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
