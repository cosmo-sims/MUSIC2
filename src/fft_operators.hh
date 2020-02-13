#ifndef __FFT_OPERATORS_HH
#define __FFT_OPERATORS_HH
struct fft_interp{

  template< typename m1, typename m2 >
  void interpolate( m1& V, m2& v, bool fourier_splice = false ) const
  {
    int oxc = V.offset(0), oyc = V.offset(1), ozc = V.offset(2);
    int oxf = v.offset(0), oyf = v.offset(1), ozf = v.offset(2);
    
    size_t nxf = v.size(0), nyf = v.size(1), nzf = v.size(2), nzfp = nzf+2;

    // cut out piece of coarse grid that overlaps the fine:
    assert( nxf%2==0 && nyf%2==0 && nzf%2==0 );

    size_t nxc = nxf/2, nyc = nyf/2, nzc = nzf/2, nzcp = nzf/2+2;

    fftw_real *rcoarse = new fftw_real[ nxc * nyc * nzcp ];
    fftw_complex *ccoarse = reinterpret_cast<fftw_complex*> (rcoarse);

    fftw_real *rfine = new fftw_real[ nxf * nyf * nzfp];
    fftw_complex *cfine = reinterpret_cast<fftw_complex*> (rfine);

    #pragma omp parallel for
    for( int i=0; i<(int)nxc; ++i )
      for( int j=0; j<(int)nyc; ++j )
	for( int k=0; k<(int)nzc; ++k ) 
	  {
	    size_t q = ((size_t)i*nyc+(size_t)j)*nzcp+(size_t)k;
	    rcoarse[q] = V( oxf+i, oyf+j, ozf+k );
	  }

    if( fourier_splice )
      {
	#pragma omp parallel for
	for( int i=0; i<(int)nxf; ++i )
	  for( int j=0; j<(int)nyf; ++j )
	    for( int k=0; k<(int)nzf; ++k ) 
	      {
		size_t q = ((size_t)i*nyf+(size_t)j)*nzfp+(size_t)k;
		rfine[q] = v(i,j,k);
	      }
      }
    else
      {
	#pragma omp parallel for
	for( size_t i=0; i<nxf*nyf*nzfp; ++i )
	  rfine[i] = 0.0;
      }

#ifdef FFTW3
#ifdef SINGLE_PRECISION
    fftwf_plan
      pc  = fftwf_plan_dft_r2c_3d( nxc, nyc, nzc, rcoarse, ccoarse, FFTW_ESTIMATE),
      pf  = fftwf_plan_dft_r2c_3d( nxf, nyf, nzf, rfine, cfine, FFTW_ESTIMATE),
      ipf = fftwf_plan_dft_c2r_3d( nxf, nyf, nzf, cfine, rfine, FFTW_ESTIMATE);
    fftwf_execute( pc );
    if( fourier_splice )
      fftwf_execute( pf );
#else
    fftw_plan
      pc  = fftw_plan_dft_r2c_3d( nxc, nyc, nzc, rcoarse, ccoarse, FFTW_ESTIMATE),
      pf  = fftw_plan_dft_r2c_3d( nxf, nyf, nzf, rfine, cfine, FFTW_ESTIMATE),
      ipf = fftw_plan_dft_c2r_3d( nxf, nyf, nzf, cfine, rfine, FFTW_ESTIMATE);
    fftw_execute( pc );
    if( fourier_splice )
      fftwf_execute( pf );
#endif
#else
    rfftwnd_plan 
      pc  = rfftw3d_create_plan( nxc, nyc, nzc, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE),
      pf  = rfftw3d_create_plan( nxf, nyf, nzf, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE|FFTW_IN_PLACE),
      ipf = rfftw3d_create_plan( nxf, nyf, nzf, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE|FFTW_IN_PLACE);
    
#ifndef SINGLETHREAD_FFTW		
    rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), pc, rcoarse, NULL );
    if( fourier_splice )
      rfftwnd_threads_one_real_to_complex( omp_get_max_threads(), pf, rfine, NULL );
#else
    rfftwnd_one_real_to_complex( pc, rcoarse, NULL );
    if( fourier_splice )
      rfftwnd_one_real_to_complex( pf, rfine, NULL );
#endif
#endif

    /*************************************************/
    //.. perform actual interpolation
    double fftnorm = 1.0/((double)nxf*(double)nyf*(double)nzf);
    double sqrt8 = sqrt(8.0);

    // 0 0
    #pragma omp parallel for
    for( int i=0; i<(int)nxc/2+1; i++ )
      for( int j=0; j<(int)nyc/2+1; j++ )
	for( int k=0; k<(int)nzc/2+1; k++ )
	  {
	    int ii(i),jj(j),kk(k);
	    size_t qc,qf;
	    qc = ((size_t)i*(size_t)nyc+(size_t)j)*(nzc/2+1)+(size_t)k;
	    qf = ((size_t)ii*(size_t)nyf+(size_t)jj)*(nzf/2+1)+(size_t)kk;
            
	    RE(cfine[qf]) = sqrt8*RE(ccoarse[qc]);
	    IM(cfine[qf]) = sqrt8*IM(ccoarse[qc]);
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
            
	    RE(cfine[qf]) = sqrt8*RE(ccoarse[qc]);
	    IM(cfine[qf]) = sqrt8*IM(ccoarse[qc]);
            
	    //if( k==0 & (i==(int)nxc/2 || j==(int)nyc/2) )
	    //  IM(cfine[qf]) *= -1.0;
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
            
	    RE(cfine[qf]) = sqrt8*RE(ccoarse[qc]);
	    IM(cfine[qf]) = sqrt8*IM(ccoarse[qc]);
            
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
	    qf = ((size_t)ii*(size_t)nyf+(size_t)jj)*(nzf/2+1)+(size_t)kk;
            
	    RE(cfine[qf]) = sqrt8*RE(ccoarse[qc]);
	    IM(cfine[qf]) = sqrt8*IM(ccoarse[qc]);
	  }
        
    delete[] rcoarse;

    /*************************************************/    

#ifdef FFTW3
  #ifdef SINGLE_PRECISION
    fftwf_execute( ipf );
    fftwf_destroy_plan(pf);
    fftwf_destroy_plan(pc);
    fftwf_destroy_plan(ipf);
  #else
    fftw_execute( ipf );
    fftw_destroy_plan(pf);
    fftw_destroy_plan(pc);
    fftw_destroy_plan(ipf);
  #endif
#else
  #ifndef SINGLETHREAD_FFTW		
    rfftwnd_threads_one_complex_to_real( omp_get_max_threads(), ipf, cfine, NULL );
  #else
    rfftwnd_one_complex_to_real( ipf, cfine, NULL );
  #endif
    fftwnd_destroy_plan(pf);
    fftwnd_destroy_plan(pc);
    fftwnd_destroy_plan(ipf);
#endif

    // copy back and normalize
    #pragma omp parallel for
    for( int i=0; i<(int)nxf; ++i )
      for( int j=0; j<(int)nyf; ++j )
	for( int k=0; k<(int)nzf; ++k ) 
	  {
	    size_t q = ((size_t)i*nyf+(size_t)j)*nzfp+(size_t)k;
	    v(i,j,k) = rfine[q] * fftnorm;
	  }

    delete[] rcoarse;
    delete[] rfine;

  }



};


#endif //__FFT_OPERATORS_HH
