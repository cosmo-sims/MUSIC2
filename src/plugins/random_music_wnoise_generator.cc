
#include <complex>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "random.hh"
#include "random_music_wnoise_generator.hh"

template <typename T>
music_wnoise_generator<T>::music_wnoise_generator(unsigned res, unsigned cubesize, long baseseed, int *x0, int *lx)
    : res_(res), cubesize_(cubesize), ncubes_(1), baseseed_(baseseed)
{
  LOGINFO("Generating random numbers (1) with seed %ld", baseseed);

  initialize();
  fill_subvolume(x0, lx);
}

template <typename T>
music_wnoise_generator<T>::music_wnoise_generator(unsigned res, unsigned cubesize, long baseseed, bool zeromean)
    : res_(res), cubesize_(cubesize), ncubes_(1), baseseed_(baseseed)
{
  LOGINFO("Generating random numbers (2) with seed %ld", baseseed);

  double mean = 0.0;
  size_t res_l = res;

  bool musicnoise = true;
  if (!musicnoise)
    cubesize_ = res_;

  if (!musicnoise)
    LOGERR("This currently breaks compatibility. Need to disable by hand! Make sure to not check into repo");

  initialize();

  if (musicnoise)
    mean = fill_all();
  else
  {
    rnums_.push_back(new Meshvar<T>(res, 0, 0, 0));
    cubemap_[0] = 0; // create dummy map index
    register_cube(0, 0, 0);
    // rapid_proto_ngenic_rng( res_, baseseed_, *this );
  }

  /*

    if( musicnoise )
    mean = fill_all();
    else
    {
    rnums_.push_back( new Meshvar<T>( res, 0, 0, 0 ) );
    cubemap_[0] = 0; // create dummy map index
    register_cube(0,0,0);
    rapid_proto_ngenic_rng( res_, baseseed_, *this );
    }

  */

  if (zeromean)
  {
    mean = 0.0;

#pragma omp parallel for reduction(+ \
                                   : mean)
    for (int i = 0; i < (int)res_; ++i)
      for (unsigned j = 0; j < res_; ++j)
        for (unsigned k = 0; k < res_; ++k)
          mean += (*this)(i, j, k);

    mean *= 1.0 / (double)(res_l * res_l * res_l);

#pragma omp parallel for
    for (int i = 0; i < (int)res_; ++i)
      for (unsigned j = 0; j < res_; ++j)
        for (unsigned k = 0; k < res_; ++k)
          (*this)(i, j, k) = (*this)(i, j, k) - mean;
  }
}

template <typename T>
music_wnoise_generator<T>::music_wnoise_generator(unsigned res, std::string randfname, bool randsign)
    : res_(res), cubesize_(res), ncubes_(1)
{
  rnums_.push_back(new Meshvar<T>(res, 0, 0, 0));
  cubemap_[0] = 0; // create dummy map index

  std::ifstream ifs(randfname.c_str(), std::ios::binary);
  if (!ifs)
  {
    LOGERR("Could not open random number file \'%s\'!", randfname.c_str());
    throw std::runtime_error(std::string("Could not open random number file \'") + randfname + std::string("\'!"));
  }

  unsigned vartype;
  unsigned nx, ny, nz, blksz32;
  size_t blksz64;
  int iseed;
  // long seed;

  float sign4 = -1.0f;
  double sign8 = -1.0;

  int addrtype = 32;

  if (randsign) // use grafic2 sign convention
  {
    sign4 = 1.0f;
    sign8 = 1.0;
  }

  //... read header and check if 32bit or 64bit block size .../
  ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
  ifs.read(reinterpret_cast<char *>(&nx), sizeof(unsigned));
  if (blksz32 != 4 * sizeof(int) || nx != res_)
  {
    addrtype = 64;

    ifs.seekg(0);
    ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));
    ifs.read(reinterpret_cast<char *>(&nx), sizeof(unsigned));

    if (blksz64 != 4 * sizeof(int) || nx != res_)
      addrtype = -1;
  }
  ifs.seekg(0);

  if (addrtype < 0)
    throw std::runtime_error("corrupt random number file");

  if (addrtype == 32)
    ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
  else
    ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));

  ifs.read(reinterpret_cast<char *>(&nx), sizeof(unsigned));
  ifs.read(reinterpret_cast<char *>(&ny), sizeof(unsigned));
  ifs.read(reinterpret_cast<char *>(&nz), sizeof(unsigned));
  ifs.read(reinterpret_cast<char *>(&iseed), sizeof(int));
  // seed = (long)iseed;

  if (nx != res_ || ny != res_ || nz != res_)
  {
    char errmsg[128];
    sprintf(errmsg, "White noise file dimensions do not match level dimensions: %ux%ux%u vs. %u**3", nx, ny, nz, res_);
    throw std::runtime_error(errmsg);
  }

  if (addrtype == 32)
    ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
  else
    ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));

  //... read data ...//
  // check whether random numbers are single or double precision numbers
  if (addrtype == 32)
  {
    ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
    if (blksz32 == nx * ny * sizeof(float))
      vartype = 4;
    else if (blksz32 == nx * ny * sizeof(double))
      vartype = 8;
    else
      throw std::runtime_error("corrupt random number file");
  }
  else
  {

    ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));
    if (blksz64 == nx * ny * sizeof(float))
      vartype = 4;
    else if (blksz64 == nx * ny * sizeof(double))
      vartype = 8;
    else
      throw std::runtime_error("corrupt random number file");
  }

  // rewind to beginning of block
  if (addrtype == 32)
    ifs.seekg(-sizeof(int), std::ios::cur);
  else
    ifs.seekg(-sizeof(size_t), std::ios::cur);

  std::vector<float> in_float;
  std::vector<double> in_double;

  LOGINFO("Random number file \'%s\'\n   contains %ld numbers. Reading...", randfname.c_str(), nx * ny * nz);

  long double sum = 0.0, sum2 = 0.0;
  size_t count = 0;

  // perform actual reading
  if (vartype == 4)
  {
    for (int ii = 0; ii < (int)nz; ++ii)
    {

      if (addrtype == 32)
      {
        ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
        if (blksz32 != nx * ny * sizeof(float))
          throw std::runtime_error("corrupt random number file");
      }
      else
      {
        ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));
        if (blksz64 != nx * ny * sizeof(float))
          throw std::runtime_error("corrupt random number file");
      }

      in_float.assign(nx * ny, 0.0f);
      ifs.read((char *)&in_float[0], nx * ny * sizeof(float));

      for (int jj = 0, q = 0; jj < (int)ny; ++jj)
        for (int kk = 0; kk < (int)nx; ++kk)
        {
          sum += in_float[q];
          sum2 += in_float[q] * in_float[q];
          ++count;

          (*rnums_[0])(kk, jj, ii) = sign4 * in_float[q++];
        }

      if (addrtype == 32)
      {
        ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
        if (blksz32 != nx * ny * sizeof(float))
          throw std::runtime_error("corrupt random number file");
      }
      else
      {
        ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));
        if (blksz64 != nx * ny * sizeof(float))
          throw std::runtime_error("corrupt random number file");
      }
    }
  }
  else if (vartype == 8)
  {
    for (int ii = 0; ii < (int)nz; ++ii)
    {
      if (addrtype == 32)
      {
        ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
        if (blksz32 != nx * ny * sizeof(double))
          throw std::runtime_error("corrupt random number file");
      }
      else
      {
        ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));
        if (blksz64 != nx * ny * sizeof(double))
          throw std::runtime_error("corrupt random number file");
      }

      in_double.assign(nx * ny, 0.0f);
      ifs.read((char *)&in_double[0], nx * ny * sizeof(double));

      for (int jj = 0, q = 0; jj < (int)ny; ++jj)
        for (int kk = 0; kk < (int)nx; ++kk)
        {
          sum += in_double[q];
          sum2 += in_double[q] * in_double[q];
          ++count;
          (*rnums_[0])(kk, jj, ii) = sign8 * in_double[q++];
        }

      if (addrtype == 32)
      {
        ifs.read(reinterpret_cast<char *>(&blksz32), sizeof(int));
        if (blksz32 != nx * ny * sizeof(double))
          throw std::runtime_error("corrupt random number file");
      }
      else
      {
        ifs.read(reinterpret_cast<char *>(&blksz64), sizeof(size_t));
        if (blksz64 != nx * ny * sizeof(double))
          throw std::runtime_error("corrupt random number file");
      }
    }
  }

  double mean, var;
  mean = sum / count;
  var = sum2 / count - mean * mean;

  LOGINFO("Random numbers in file have \n     mean = %f and var = %f", mean, var);
}

//... copy construct by averaging down
template <typename T>
music_wnoise_generator<T>::music_wnoise_generator(/*const*/ music_wnoise_generator<T> &rc)
{
  // if( res > rc.m_res || (res/rc.m_res)%2 != 0 )
  //			throw std::runtime_error("Invalid restriction in random number container copy constructor.");

  long double sum = 0.0, sum2 = 0.0;
  size_t count = 0;


  LOGINFO("Generating a coarse white noise field by k-space degrading");
  //... initialize properties of container
  res_ = rc.res_ / 2;
  cubesize_ = res_;
  ncubes_ = 1;
  baseseed_ = -2;

  if (sizeof(fftw_real) != sizeof(T))
  {
    LOGERR("type mismatch with fftw_real in k-space averaging");
    throw std::runtime_error("type mismatch with fftw_real in k-space averaging");
  }

  fftw_real
      *rfine = new fftw_real[(size_t)rc.res_ * (size_t)rc.res_ * 2 * ((size_t)rc.res_ / 2 + 1)],
      *rcoarse = new fftw_real[(size_t)res_ * (size_t)res_ * 2 * ((size_t)res_ / 2 + 1)];

  fftw_complex
      *ccoarse = reinterpret_cast<fftw_complex *>(rcoarse),
      *cfine = reinterpret_cast<fftw_complex *>(rfine);

  int nx(rc.res_), ny(rc.res_), nz(rc.res_), nxc(res_), nyc(res_), nzc(res_);
#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_plan
      pf = fftwf_plan_dft_r2c_3d(nx, ny, nz, rfine, cfine, FFTW_ESTIMATE),
      ipc = fftwf_plan_dft_c2r_3d(nxc, nyc, nzc, ccoarse, rcoarse, FFTW_ESTIMATE);
#else
  fftw_plan
      pf = fftw_plan_dft_r2c_3d(nx, ny, nz, rfine, cfine, FFTW_ESTIMATE),
      ipc = fftw_plan_dft_c2r_3d(nxc, nyc, nzc, ccoarse, rcoarse, FFTW_ESTIMATE);
#endif

#else
  rfftwnd_plan
      pf = rfftw3d_create_plan(nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE),
      ipc = rfftw3d_create_plan(nxc, nyc, nzc, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
#endif

#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
      {
        size_t q = ((size_t)i * ny + (size_t)j) * (nz + 2) + (size_t)k;
        rfine[q] = rc(i, j, k);
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

  double fftnorm = 1.0 / ((double)nxc * (double)nyc * (double)nzc);

#pragma omp parallel for
  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc / 2 + 1; k++)
      {
        int ii(i), jj(j), kk(k);

        if (i > nxc / 2)
          ii += nx / 2;
        if (j > nyc / 2)
          jj += ny / 2;

        size_t qc, qf;

        double kx = (i <= (int)nxc / 2) ? (double)i : (double)(i - (int)nxc);
        double ky = (j <= (int)nyc / 2) ? (double)j : (double)(j - (int)nyc);
        double kz = (k <= (int)nzc / 2) ? (double)k : (double)(k - (int)nzc);

        qc = ((size_t)i * nyc + (size_t)j) * (nzc / 2 + 1) + (size_t)k;
        qf = ((size_t)ii * ny + (size_t)jj) * (nz / 2 + 1) + (size_t)kk;

        std::complex<double> val_fine(RE(cfine[qf]), IM(cfine[qf]));
        double phase = (kx / nxc + ky / nyc + kz / nzc) * 0.5 * M_PI;
        std::complex<double> val_phas(cos(phase), sin(phase));

        val_fine *= val_phas * fftnorm / sqrt(8.0);

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
  rnums_.push_back(new Meshvar<T>(res_, 0, 0, 0));
  cubemap_[0] = 0; // map all to single array

#pragma omp parallel for reduction(+ \
                                  : sum, sum2, count)
  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++)
      {
        size_t q = ((size_t)i * nyc + (size_t)j) * (nzc + 2) + (size_t)k;
        (*rnums_[0])(i, j, k) = rcoarse[q];
        sum += (*rnums_[0])(i, j, k);
        sum2 += (*rnums_[0])(i, j, k) * (*rnums_[0])(i, j, k);
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
  
  double rmean, rvar;
  rmean = sum / count;
  rvar = sum2 / count - rmean * rmean;

  LOGINFO("Restricted random numbers have\n       mean = %f, var = %f", rmean, rvar);
}

template <typename T>
music_wnoise_generator<T>::music_wnoise_generator(music_wnoise_generator<T> &rc, unsigned cubesize, long baseseed, int *x0_, int *lx_, bool zeromean)
    : res_(2 * rc.res_), cubesize_(cubesize), ncubes_(1), baseseed_(baseseed)
{
  initialize();

  int x0[3], lx[3];
  if (x0_ == NULL || lx_ == NULL)
  {
    for (int i = 0; i < 3; ++i)
    {
      x0[i] = 0;
      lx[i] = res_;
    }
    fill_all();
  }
  else
  {
    for (int i = 0; i < 3; ++i)
    {
      x0[i] = x0_[i];
      lx[i] = lx_[i];
    }
    fill_subvolume(x0, lx);
  }


  LOGINFO("Generating a constrained random number set with seed %ld\n    using coarse mode replacement...", baseseed);
  assert(lx[0] % 2 == 0 && lx[1] % 2 == 0 && lx[2] % 2 == 0);
  size_t nx = lx[0], ny = lx[1], nz = lx[2],
          nxc = lx[0] / 2, nyc = lx[1] / 2, nzc = lx[2] / 2;

  fftw_real *rfine = new fftw_real[nx * ny * (nz + 2l)];
  fftw_complex *cfine = reinterpret_cast<fftw_complex *>(rfine);

#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_plan
      pf = fftwf_plan_dft_r2c_3d(nx, ny, nz, rfine, cfine, FFTW_ESTIMATE),
      ipf = fftwf_plan_dft_c2r_3d(nx, ny, nz, cfine, rfine, FFTW_ESTIMATE);
#else
  fftw_plan
      pf = fftw_plan_dft_r2c_3d(nx, ny, nz, rfine, cfine, FFTW_ESTIMATE),
      ipf = fftw_plan_dft_c2r_3d(nx, ny, nz, cfine, rfine, FFTW_ESTIMATE);
#endif
#else
  rfftwnd_plan
      pf = rfftw3d_create_plan(nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE),
      ipf = rfftw3d_create_plan(nx, ny, nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
#endif

#pragma omp parallel for
  for (int i = 0; i < (int)nx; i++)
    for (int j = 0; j < (int)ny; j++)
      for (int k = 0; k < (int)nz; k++)
      {
        size_t q = ((size_t)i * (size_t)ny + (size_t)j) * (size_t)(nz + 2) + (size_t)k;
        rfine[q] = (*this)(x0[0] + i, x0[1] + j, x0[2] + k);
      }
  // this->free_all_mem();	// temporarily free memory, allocate again later

  fftw_real *rcoarse = new fftw_real[nxc * nyc * (nzc + 2)];
  fftw_complex *ccoarse = reinterpret_cast<fftw_complex *>(rcoarse);

#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_plan pc = fftwf_plan_dft_r2c_3d(nxc, nyc, nzc, rcoarse, ccoarse, FFTW_ESTIMATE);
#else
  fftw_plan pc = fftw_plan_dft_r2c_3d(nxc, nyc, nzc, rcoarse, ccoarse, FFTW_ESTIMATE);
#endif
#else
  rfftwnd_plan pc = rfftw3d_create_plan(nxc, nyc, nzc, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
#endif

#pragma omp parallel for
  for (int i = 0; i < (int)nxc; i++)
    for (int j = 0; j < (int)nyc; j++)
      for (int k = 0; k < (int)nzc; k++)
      {
        size_t q = ((size_t)i * (size_t)nyc + (size_t)j) * (size_t)(nzc + 2) + (size_t)k;
        rcoarse[q] = rc(x0[0] / 2 + i, x0[1] / 2 + j, x0[2] / 2 + k);
      }
#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_execute(pc);
  fftwf_execute(pf);
#else
  fftw_execute(pc);
  fftw_execute(pf);
#endif
#else
#ifndef SINGLETHREAD_FFTW
  rfftwnd_threads_one_real_to_complex(omp_get_max_threads(), pc, rcoarse, NULL);
  rfftwnd_threads_one_real_to_complex(omp_get_max_threads(), pf, rfine, NULL);
#else
  rfftwnd_one_real_to_complex(pc, rcoarse, NULL);
  rfftwnd_one_real_to_complex(pf, rfine, NULL);
#endif
#endif

  double fftnorm = 1.0 / ((double)nx * (double)ny * (double)nz);
  double sqrt8 = sqrt(8.0);
  double phasefac = -0.5; 

  // embedding of coarse white noise by fourier interpolation

#pragma omp parallel for
  for (int i = 0; i < (int)nxc; i++)
    for (int j = 0; j < (int)nyc; j++)
      for (int k = 0; k < (int)nzc / 2 + 1; k++)
      {
        int ii(i), jj(j), kk(k);

        // if( i==(int)nxc/2 ) continue;
        // if( j==(int)nyc/2 ) continue;

        if (i > (int)nxc / 2)
          ii += (int)nx / 2;
        if (j > (int)nyc / 2)
          jj += (int)ny / 2;

        size_t qc, qf;

        double kx = (i <= (int)nxc / 2) ? (double)i : (double)(i - (int)nxc);
        double ky = (j <= (int)nyc / 2) ? (double)j : (double)(j - (int)nyc);
        double kz = (k <= (int)nzc / 2) ? (double)k : (double)(k - (int)nzc);

        qc = ((size_t)i * nyc + (size_t)j) * (nzc / 2 + 1) + (size_t)k;
        qf = ((size_t)ii * ny + (size_t)jj) * (nz / 2 + 1) + (size_t)kk;

        std::complex<double> val(RE(ccoarse[qc]), IM(ccoarse[qc]));
        double phase = (kx / nxc + ky / nyc + kz / nzc) * phasefac * M_PI;

        std::complex<double> val_phas(cos(phase), sin(phase));

        val *= val_phas * sqrt8;

        if (i != (int)nxc / 2 && j != (int)nyc / 2 && k != (int)nzc / 2)
        {
          RE(cfine[qf]) = val.real();
          IM(cfine[qf]) = val.imag();
        }
        else
        {
          // RE(cfine[qf]) = val.real();
          // IM(cfine[qf]) = 0.0;
        }
      }

    delete[] rcoarse;

#pragma omp parallel for
    for (int i = 0; i < (int)nx; i++)
      for (int j = 0; j < (int)ny; j++)
        for (int k = 0; k < (int)nz / 2 + 1; k++)
        {
          size_t q = ((size_t)i * ny + (size_t)j) * (nz / 2 + 1) + (size_t)k;

          RE(cfine[q]) *= fftnorm;
          IM(cfine[q]) *= fftnorm;
        }

#ifdef FFTW3
#ifdef SINGLE_PRECISION
    fftwf_execute(ipf);
#else
    fftw_execute(ipf);
#endif
#else
#ifndef SINGLETHREAD_FFTW
    rfftwnd_threads_one_complex_to_real(omp_get_max_threads(), ipf, cfine, NULL);
#else
    rfftwnd_one_complex_to_real(ipf, cfine, NULL);
#endif
#endif

#pragma omp parallel for
    for (int i = 0; i < (int)nx; i++)
      for (int j = 0; j < (int)ny; j++)
        for (int k = 0; k < (int)nz; k++)
        {
          size_t q = ((size_t)i * ny + (size_t)j) * (nz + 2) + (size_t)k;
          (*this)(x0[0] + i, x0[1] + j, x0[2] + k, false) = rfine[q];
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

template <typename T>
void music_wnoise_generator<T>::register_cube(int i, int j, int k)
{
  i = (i + ncubes_) % ncubes_;
  j = (j + ncubes_) % ncubes_;
  k = (k + ncubes_) % ncubes_;
  size_t icube = ((size_t)i * ncubes_ + (size_t)j) * ncubes_ + (size_t)k;

  cubemap_iterator it = cubemap_.find(icube);

  if (it == cubemap_.end())
  {
    rnums_.push_back(NULL);
    cubemap_[icube] = rnums_.size() - 1;
#ifdef DEBUG
    LOGDEBUG("registering new cube %d,%d,%d . ID = %ld, memloc = %ld", i, j, k, icube, cubemap_[icube]);
#endif
  }
}

template <typename T>
double music_wnoise_generator<T>::fill_cube(int i, int j, int k)
{

  gsl_rng *RNG = gsl_rng_alloc(gsl_rng_mt19937);

  i = (i + ncubes_) % ncubes_;
  j = (j + ncubes_) % ncubes_;
  k = (k + ncubes_) % ncubes_;

  size_t icube = ((size_t)i * ncubes_ + (size_t)j) * ncubes_ + (size_t)k;
  long cubeseed = baseseed_ + icube; //... each cube gets its unique seed

  gsl_rng_set(RNG, cubeseed);

  cubemap_iterator it = cubemap_.find(icube);

  if (it == cubemap_.end())
  {
    LOGERR("Attempt to access non-registered random number cube!");
    throw std::runtime_error("Attempt to access non-registered random number cube!");
  }

  size_t cubeidx = it->second;

  if (rnums_[cubeidx] == NULL)
    rnums_[cubeidx] = new Meshvar<T>(cubesize_, 0, 0, 0);

  double mean = 0.0;

  for (int ii = 0; ii < (int)cubesize_; ++ii)
    for (int jj = 0; jj < (int)cubesize_; ++jj)
      for (int kk = 0; kk < (int)cubesize_; ++kk)
      {
        (*rnums_[cubeidx])(ii, jj, kk) = gsl_ran_ugaussian_ratio_method(RNG);
        mean += (*rnums_[cubeidx])(ii, jj, kk);
      }

  gsl_rng_free(RNG);

  return mean / (cubesize_ * cubesize_ * cubesize_);
}

template <typename T>
void music_wnoise_generator<T>::subtract_from_cube(int i, int j, int k, double val)
{
  i = (i + ncubes_) % ncubes_;
  j = (j + ncubes_) % ncubes_;
  k = (k + ncubes_) % ncubes_;

  size_t icube = ((size_t)i * ncubes_ + (size_t)j) * ncubes_ + (size_t)k;

  cubemap_iterator it = cubemap_.find(icube);

  if (it == cubemap_.end())
  {
    LOGERR("Attempt to access unallocated RND cube %d,%d,%d in music_wnoise_generator::subtract_from_cube", i, j, k);
    throw std::runtime_error("Attempt to access unallocated RND cube in music_wnoise_generator::subtract_from_cube");
  }

  size_t cubeidx = it->second;

  for (int ii = 0; ii < (int)cubesize_; ++ii)
    for (int jj = 0; jj < (int)cubesize_; ++jj)
      for (int kk = 0; kk < (int)cubesize_; ++kk)
        (*rnums_[cubeidx])(ii, jj, kk) -= val;
}

template <typename T>
void music_wnoise_generator<T>::free_cube(int i, int j, int k)
{

  i = (i + ncubes_) % ncubes_;
  j = (j + ncubes_) % ncubes_;
  k = (k + ncubes_) % ncubes_;

  size_t icube = ((size_t)i * (size_t)ncubes_ + (size_t)j) * (size_t)ncubes_ + (size_t)k;

  cubemap_iterator it = cubemap_.find(icube);

  if (it == cubemap_.end())
  {
    LOGERR("Attempt to access unallocated RND cube %d,%d,%d in music_wnoise_generator::free_cube", i, j, k);
    throw std::runtime_error("Attempt to access unallocated RND cube in music_wnoise_generator::free_cube");
  }

  size_t cubeidx = it->second;

  if (rnums_[cubeidx] != NULL)
  {
    delete rnums_[cubeidx];
    rnums_[cubeidx] = NULL;
  }
}

template <typename T>
void music_wnoise_generator<T>::initialize(void)
{

  ncubes_ = std::max((int)((double)res_ / cubesize_), 1);
  if (res_ < cubesize_)
  {
    ncubes_ = 1;
    cubesize_ = res_;
  }

  LOGINFO("Generating random numbers w/ sample cube size of %d", cubesize_);
}

template <typename T>
double music_wnoise_generator<T>::fill_subvolume(int *i0, int *n)
{
  int i0cube[3], ncube[3];

  i0cube[0] = (int)((double)(res_ + i0[0]) / cubesize_);
  i0cube[1] = (int)((double)(res_ + i0[1]) / cubesize_);
  i0cube[2] = (int)((double)(res_ + i0[2]) / cubesize_);

  ncube[0] = (int)(n[0] / cubesize_) + 2;
  ncube[1] = (int)(n[1] / cubesize_) + 2;
  ncube[2] = (int)(n[2] / cubesize_) + 2;

#ifdef DEBUG
  LOGDEBUG("random numbers needed for region %d,%d,%d ..+ %d,%d,%d", i0[0], i0[1], i0[2], n[0], n[1], n[2]);
  LOGDEBUG("filling cubes %d,%d,%d ..+ %d,%d,%d", i0cube[0], i0cube[1], i0cube[2], ncube[0], ncube[1], ncube[2]);
#endif

  double mean = 0.0;

  for (int i = i0cube[0]; i < i0cube[0] + ncube[0]; ++i)
    for (int j = i0cube[1]; j < i0cube[1] + ncube[1]; ++j)
      for (int k = i0cube[2]; k < i0cube[2] + ncube[2]; ++k)
      {
        int ii(i), jj(j), kk(k);

        ii = (ii + ncubes_) % ncubes_;
        jj = (jj + ncubes_) % ncubes_;
        kk = (kk + ncubes_) % ncubes_;

        register_cube(ii, jj, kk);
      }

  #pragma omp parallel for reduction(+ : mean)
  for (int i = i0cube[0]; i < i0cube[0] + ncube[0]; ++i)
    for (int j = i0cube[1]; j < i0cube[1] + ncube[1]; ++j)
      for (int k = i0cube[2]; k < i0cube[2] + ncube[2]; ++k)
      {
        int ii(i), jj(j), kk(k);

        ii = (ii + ncubes_) % ncubes_;
        jj = (jj + ncubes_) % ncubes_;
        kk = (kk + ncubes_) % ncubes_;

        mean += fill_cube(ii, jj, kk);
      }
  return mean / (ncube[0] * ncube[1] * ncube[2]);
}

template <typename T>
double music_wnoise_generator<T>::fill_all(void)
{
  double sum = 0.0;

  for (int i = 0; i < (int)ncubes_; ++i)
    for (int j = 0; j < (int)ncubes_; ++j)
      for (int k = 0; k < (int)ncubes_; ++k)
      {
        int ii(i), jj(j), kk(k);

        ii = (ii + ncubes_) % ncubes_;
        jj = (jj + ncubes_) % ncubes_;
        kk = (kk + ncubes_) % ncubes_;

        register_cube(ii, jj, kk);
      }

#pragma omp parallel for reduction(+ \
                                   : sum)
  for (int i = 0; i < (int)ncubes_; ++i)
    for (int j = 0; j < (int)ncubes_; ++j)
      for (int k = 0; k < (int)ncubes_; ++k)
      {
        int ii(i), jj(j), kk(k);

        ii = (ii + ncubes_) % ncubes_;
        jj = (jj + ncubes_) % ncubes_;
        kk = (kk + ncubes_) % ncubes_;

        sum += fill_cube(ii, jj, kk);
      }

//... subtract mean
#pragma omp parallel for reduction(+ \
                                   : sum)
  for (int i = 0; i < (int)ncubes_; ++i)
    for (int j = 0; j < (int)ncubes_; ++j)
      for (int k = 0; k < (int)ncubes_; ++k)
      {
        int ii(i), jj(j), kk(k);

        ii = (ii + ncubes_) % ncubes_;
        jj = (jj + ncubes_) % ncubes_;
        kk = (kk + ncubes_) % ncubes_;
        subtract_from_cube(ii, jj, kk, sum / (ncubes_ * ncubes_ * ncubes_));
      }

  return sum / (ncubes_ * ncubes_ * ncubes_);
}

template <typename T>
void music_wnoise_generator<T>::print_allocated(void)
{
  unsigned ncount = 0, ntot = rnums_.size();
  for (size_t i = 0; i < rnums_.size(); ++i)
    if (rnums_[i] != NULL)
      ncount++;

  LOGINFO(" -> %d of %d random number cubes currently allocated", ncount, ntot);
}

template class music_wnoise_generator<float>;
template class music_wnoise_generator<double>;
