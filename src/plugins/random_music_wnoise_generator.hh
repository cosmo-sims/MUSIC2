#pragma once

#include <vector>
#include <map>
#include "general.hh"
#include "mesh.hh"

#define DEF_RAN_CUBE_SIZE 32

/*!
 * @brief encapsulates all things random number generator related
 */
template <typename T>
class music_wnoise_generator
{
public:
  unsigned
      res_,       //!< resolution of the full mesh
      cubesize_,  //!< size of one independent random number cube
      ncubes_;    //!< number of random number cubes to cover the full mesh
  long baseseed_; //!< base seed from which cube seeds are computed

protected:
  //! vector of 3D meshes (the random number cubes) with random numbers
  std::vector<Meshvar<T> *> rnums_;

  //! map of 3D indices to cube index
  std::map<size_t, size_t> cubemap_;

  typedef std::map<size_t, size_t>::iterator cubemap_iterator;

protected:
  //! register a cube with the hash map
  void register_cube(int i, int j, int k);

  //! fills a subcube with random numbers
  double fill_cube(int i, int j, int k);

  //! subtract a constant from an entire cube
  void subtract_from_cube(int i, int j, int k, double val);

  //! copy random numbers from a cube to a full grid array
  template <class C>
  void copy_cube(int i, int j, int k, C &dat)
  {
    int offi, offj, offk;

    offi = i * cubesize_;
    offj = j * cubesize_;
    offk = k * cubesize_;

    i = (i + ncubes_) % ncubes_;
    j = (j + ncubes_) % ncubes_;
    k = (k + ncubes_) % ncubes_;

    size_t icube = (i * ncubes_ + j) * ncubes_ + k;
    cubemap_iterator it = cubemap_.find(icube);

    if (it == cubemap_.end())
    {
      LOGERR("attempting to copy data from non-existing RND cube %d,%d,%d", i, j, k);
      throw std::runtime_error("attempting to copy data from non-existing RND cube");
    }

    size_t cubeidx = it->second;

    for (int ii = 0; ii < (int)cubesize_; ++ii)
      for (int jj = 0; jj < (int)cubesize_; ++jj)
        for (int kk = 0; kk < (int)cubesize_; ++kk)
          dat(offi + ii, offj + jj, offk + kk) = (*rnums_[cubeidx])(ii, jj, kk);
  }

  //! free the memory associated with a subcube
  void free_cube(int i, int j, int k);

  //! initialize member variables and allocate memory
  void initialize(void);

  //! fill a cubic subvolume of the full grid with random numbers
  double fill_subvolume(int *i0, int *n);

  //! fill an entire grid with random numbers
  double fill_all(void);

  //! fill an external array instead of the internal field
  template <class C>
  double fill_all(C &dat)
  {
    double sum = 0.0;

    for (int i = 0; i < (int)ncubes_; ++i)
      for (int j = 0; j < (int)ncubes_; ++j)
        for (int k = 0; k < (int)ncubes_; ++k)
        {
          int ii(i), jj(j), kk(k);
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
          copy_cube(ii, jj, kk, dat);
          free_cube(ii, jj, kk);
        }

    return sum / (ncubes_ * ncubes_ * ncubes_);
  }

  //! write the number of allocated random number cubes to stdout
  void print_allocated(void);

public:
  //! constructor
  music_wnoise_generator(unsigned res, unsigned cubesize, long baseseed, int *x0, int *lx);

  //! constructor for constrained fine field
  music_wnoise_generator(music_wnoise_generator<T> &rc, unsigned cubesize, long baseseed, int *x0_ = NULL, int *lx_ = NULL, bool zeromean = true);

  //! constructor
  music_wnoise_generator(unsigned res, unsigned cubesize, long baseseed, bool zeromean = true);

  //! constructor to read white noise from file
  music_wnoise_generator(unsigned res, std::string randfname, bool rndsign);

  //! copy constructor for averaged field (not copying) hence explicit!
  explicit music_wnoise_generator(/*const*/ music_wnoise_generator<T> &rc);

  //! destructor
  ~music_wnoise_generator()
  {
    for (unsigned i = 0; i < rnums_.size(); ++i)
      if (rnums_[i] != NULL)
        delete rnums_[i];
    rnums_.clear();
  }

  //! access a random number, this allocates a cube and fills it with consistent random numbers
  inline T &operator()(int i, int j, int k, bool fillrand = true)
  {
    int ic, jc, kc, is, js, ks;

    if (ncubes_ == 0)
      throw std::runtime_error("music_wnoise_generator: internal error, not properly initialized");

    //... determine cube
    ic = (int)((double)i / cubesize_ + ncubes_) % ncubes_;
    jc = (int)((double)j / cubesize_ + ncubes_) % ncubes_;
    kc = (int)((double)k / cubesize_ + ncubes_) % ncubes_;

    size_t icube = ((size_t)ic * ncubes_ + (size_t)jc) * ncubes_ + (size_t)kc;

    cubemap_iterator it = cubemap_.find(icube);

    if (it == cubemap_.end())
    {
      LOGERR("Attempting to copy data from non-existing RND cube %d,%d,%d @ %d,%d,%d", ic, jc, kc, i, j, k);
      throw std::runtime_error("attempting to copy data from non-existing RND cube");
    }

    size_t cubeidx = it->second;

    if (rnums_[cubeidx] == NULL)
    {
      LOGERR("Attempting to access data from non-allocated RND cube %d,%d,%d", ic, jc, kc);
      throw std::runtime_error("attempting to access data from non-allocated RND cube");
    }

    //... determine cell in cube
    is = (i - ic * cubesize_ + cubesize_) % cubesize_;
    js = (j - jc * cubesize_ + cubesize_) % cubesize_;
    ks = (k - kc * cubesize_ + cubesize_) % cubesize_;

    return (*rnums_[cubeidx])(is, js, ks);
  }

  //! free all cubes
  void free_all_mem(void)
  {
    for (unsigned i = 0; i < rnums_.size(); ++i)
      if (rnums_[i] != NULL)
      {
        delete rnums_[i];
        rnums_[i] = NULL;
      }
  }
};

