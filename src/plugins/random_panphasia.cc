#ifdef HAVE_PANPHASIA
#include "random.hh"
#include <cctype>
#include <cstring>
#include <stdint.h>

#include "densities.hh"
#include "HDF_IO.hh"

const int maxdim = 60, maxlev = 50, maxpow = 3 * maxdim;
typedef int rand_offset_[5];
typedef struct {
  int state[133]; // Nstore = Nstate (=5) + Nbatch (=128)
  int need_fill;
  int pos;
} rand_state_; 

/* pan_state_ struct -- corresponds to respective fortran module in panphasia_routines.f
 * data structure that contains all panphasia state variables
 * it needs to get passed between the fortran routines to enable
 * thread-safe execution.
 */
typedef struct {
  int base_state[5], base_lev_start[5][maxdim + 1];
  rand_offset_ poweroffset[maxpow + 1], superjump;
  rand_state_ current_state[maxpow + 2];

  int layer_min, layer_max, indep_field;

  long long xorigin_store[2][2][2], yorigin_store[2][2][2], zorigin_store[2][2][2];
  int lev_common, layer_min_store, layer_max_store;
  long long ix_abs_store, iy_abs_store, iz_abs_store, ix_per_store, iy_per_store, iz_per_store, ix_rel_store,
      iy_rel_store, iz_rel_store;
  double exp_coeffs[8][8][maxdim + 2];
  long long xcursor[maxdim + 1], ycursor[maxdim + 1], zcursor[maxdim + 1];
  int ixshift[2][2][2], iyshift[2][2][2], izshift[2][2][2];

  double cell_data[9][8];
  int ixh_last, iyh_last, izh_last;
  int init;

  int init_cell_props;
  int init_lecuyer_state;
  long long p_xcursor[62], p_ycursor[62], p_zcursor[62];

} pan_state_;

extern "C" {
void start_panphasia_(pan_state_ *lstate, const char *descriptor, int *ngrid, int *bverbose);

void parse_descriptor_(const char *descriptor, int16_t *l, int32_t *ix, int32_t *iy, int32_t *iz, int16_t *side1,
                       int16_t *side2, int16_t *side3, int32_t *check_int, char *name);

void panphasia_cell_properties_(pan_state_ *lstate, int *ixcell, int *iycell, int *izcell, double *cell_prop);

void adv_panphasia_cell_properties_(pan_state_ *lstate, int *ixcell, int *iycell, int *izcell, int *layer_min,
                                    int *layer_max, int *indep_field, double *cell_prop);

void set_phases_and_rel_origin_(pan_state_ *lstate, const char *descriptor, int *lev, long long *ix_rel,
                                long long *iy_rel, long long *iz_rel, int *VERBOSE);
/*void set_local_box_( pan_state_ *lstate, int lev, int8_t ix_abs, int8_t iy_abs, int8_t iz_abs,
                     int8_t ix_per, int8_t iy_per, int8_t iz_per, int8_t ix_rel, int8_t iy_rel,
                     int8_t iz_rel, int wn_level_base, int8_t check_rand, char *phase_name, int MYID);*/
/*extern struct {
  int layer_min, layer_max, hoswitch;
  }oct_range_;
*/
}

class RNG_panphasia : public RNG_plugin {
private:
  void forward_transform_field(real_t *field, int n0, int n1, int n2);
  void forward_transform_field(real_t *field, int n) { forward_transform_field(field, n, n, n); }

  void backward_transform_field(real_t *field, int n0, int n1, int n2);
  void backward_transform_field(real_t *field, int n) { backward_transform_field(field, n, n, n); }

protected:
  std::string descriptor_string_;
  int num_threads_;
  int levelmin_, levelmin_final_, levelmax_, ngrid_;
  bool incongruent_fields_;
  double inter_grid_phase_adjustment_;
  // double translation_phase_;
  pan_state_ *lstate;
  int grid_p_,grid_m_;
  double grid_rescale_fac_;
  int coordinate_system_shift_[3];
  int ix_abs_[3], ix_per_[3], ix_rel_[3], level_p_, lextra_;
  const refinement_hierarchy *prefh_;

  struct panphasia_descriptor {
    int16_t wn_level_base;
    int32_t i_xorigin_base, i_yorigin_base, i_zorigin_base;
    int16_t i_base, i_base_y, i_base_z;
    int32_t check_rand;
    std::string name;

    explicit panphasia_descriptor(std::string dstring) {
      char tmp[100];
      memset(tmp, ' ', 100);
      parse_descriptor_(dstring.c_str(), &wn_level_base, &i_xorigin_base, &i_yorigin_base, &i_zorigin_base, &i_base,
                        &i_base_y, &i_base_z, &check_rand, tmp);
      for (int i = 0; i < 100; i++)
        if (tmp[i] == ' ') {
          tmp[i] = '\0';
          break;
        }
      name = tmp;
      name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
    }
  };

  void clear_panphasia_thread_states(void) {
    for (int i = 0; i < num_threads_; ++i) {
      lstate[i].init = 0;
      lstate[i].init_cell_props = 0;
      lstate[i].init_lecuyer_state = 0;
    }
  }

  // greatest common divisor
  int gcd(int a, int b) {
    if (b == 0)
      return a;
    return gcd(b, a % b);
  }

  // least common multiple
  int lcm(int a, int b) { return abs(a * b) / gcd(a, b); }

// Two or largest power of 2 less than the argument
  int largest_power_two_lte(int b) {
    int a = 1;
    if (b<=a) return a;
    while (2*a < b) a = 2*a;
    return a;
   }


  panphasia_descriptor *pdescriptor_;

public:
  explicit RNG_panphasia(config_file &cf) : RNG_plugin(cf) {
    descriptor_string_ = pcf_->getValue<std::string>("random", "descriptor");

#ifdef _OPENMP
    num_threads_ = omp_get_max_threads();
#else
    num_threads_ = 1;
#endif

    // create independent state descriptions for each thread
    lstate = new pan_state_[num_threads_];

    // parse the descriptor for its properties
    pdescriptor_ = new panphasia_descriptor(descriptor_string_);
    LOGINFO("PANPHASIA: descriptor \'%s\' is base %d,", pdescriptor_->name.c_str(), pdescriptor_->i_base);

    // write panphasia base size into config file for the grid construction
    // as the gridding unit we use the least common multiple of 2 and i_base
    std::stringstream ss;
    //ARJ  ss << lcm(2, pdescriptor_->i_base);
    //ss <<  two_or_largest_power_two_less_than(pdescriptor_->i_base);//ARJ
    ss << 2; //ARJ - set gridding unit to two
    pcf_->insertValue("setup", "gridding_unit", ss.str());
    ss.str(std::string());
    ss <<  pdescriptor_->i_base ;
    pcf_->insertValue("random","base_unit", ss.str());
  }

  void initialize_for_grid_structure(const refinement_hierarchy &refh) {
    prefh_ = &refh;
    levelmin_ = prefh_->levelmin();
    levelmin_final_ = pcf_->getValue<unsigned>("setup", "levelmin");
    levelmax_ = prefh_->levelmax();

    clear_panphasia_thread_states();
    LOGINFO("PANPHASIA: running with %d threads", num_threads_);

    // if ngrid is not a multiple of i_base, then we need to enlarge and then sample down
    ngrid_ = 1 << levelmin_;
    
    grid_p_ = pdescriptor_->i_base;
    grid_m_ = largest_power_two_lte(grid_p_);

    lextra_ = (log10((double)ngrid_ / (double)pdescriptor_->i_base) + 0.001) / log10(2.0);
    int ratio  = 1 << lextra_;
    grid_rescale_fac_ = 1.0;
 
    coordinate_system_shift_[0] = -pcf_->getValue<int>("setup", "shift_x");
    coordinate_system_shift_[1] = -pcf_->getValue<int>("setup", "shift_y");
    coordinate_system_shift_[2] = -pcf_->getValue<int>("setup", "shift_z");

    incongruent_fields_ = false;
    if (ngrid_ != ratio * pdescriptor_->i_base) {
      incongruent_fields_ = true;
      ngrid_ = 2 * ratio * pdescriptor_->i_base;
      grid_rescale_fac_ = (double)ngrid_ / (1 << levelmin_);
      LOGINFO("PANPHASIA: will use a higher resolution:\n"
              "     (%d -> %d) * 2**ref compatible with PANPHASIA\n"
              "     will Fourier interpolate after.",
              grid_m_, grid_p_);
    } 
  }

  ~RNG_panphasia() { delete[] lstate; }

  void fill_grid(int level, DensityGrid<real_t> &R);

  bool is_multiscale() const { return true; }
};

void RNG_panphasia::forward_transform_field(real_t *field, int nx, int ny, int nz) {

  fftw_real *rfield = reinterpret_cast<fftw_real *>(field);
  fftw_complex *cfield = reinterpret_cast<fftw_complex *>(field);

#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_plan pf = fftwf_plan_dft_r2c_3d(nx, ny, nz, rfield, cfield, FFTW_ESTIMATE);
#else
  fftw_plan pf = fftw_plan_dft_r2c_3d(nx, ny, nz, rfield, cfield, FFTW_ESTIMATE);
#endif
#else
  rfftwnd_plan pf = rfftw3d_create_plan(nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
#endif

#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_execute(pf);
#else
  fftw_execute(pf);
#endif
#else
#ifndef SINGLETHREAD_FFTW
  rfftwnd_threads_one_real_to_complex(num_threads_, pf, rfield, NULL);
#else
  rfftwnd_one_real_to_complex(pf, rfield, NULL);
#endif
#endif
}

void RNG_panphasia::backward_transform_field(real_t *field, int nx, int ny, int nz) {

  fftw_real *rfield = reinterpret_cast<fftw_real *>(field);
  fftw_complex *cfield = reinterpret_cast<fftw_complex *>(field);

#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_plan ipf = fftwf_plan_dft_c2r_3d(nx, ny, nz, cfield, rfield, FFTW_ESTIMATE);
#else
  fftw_plan ipf = fftw_plan_dft_c2r_3d(nx, ny, nz, cfield, rfield, FFTW_ESTIMATE);
#endif
#else
  rfftwnd_plan ipf = rfftw3d_create_plan(nx, ny, nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
#endif

#ifdef FFTW3
#ifdef SINGLE_PRECISION
  fftwf_execute(ipf);
#else
  fftw_execute(ipf);
#endif
#else
#ifndef SINGLETHREAD_FFTW
  rfftwnd_threads_one_complex_to_real(num_threads_, ipf, cfield, NULL);
#else
  rfftwnd_one_complex_to_real(ipf, cfield, NULL);
#endif
#endif
}

#include <sys/time.h>
inline double get_wtime(void) {
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return (double)clock() / CLOCKS_PER_SEC;
#endif
}

void RNG_panphasia::fill_grid(int level, DensityGrid<real_t> &R) {
  fftw_real *pr0, *pr1, *pr2, *pr3, *pr4;
  fftw_complex *pc0, *pc1, *pc2, *pc3, *pc4;

  
  // determine resolution and offset so that we can do proper resampling
  int ileft[3], ileft_corner[3], nx[3], nxremap[3];
  int iexpand_left[3];

  for (int k = 0; k < 3; ++k) {
    ileft[k] = prefh_->offset_abs(level, k); 
    nx[k] = R.size(k);
    assert(nx[k] % 4 == 0); 
    if (level == levelmin_) {
        ileft_corner[k] = ileft[k];   // Top level - periodic
    }else{
        ileft_corner[k] = (ileft[k] - nx[k]/4 + (1 << level))%(1 << level); // Isolated
    }
    iexpand_left[k] =  (ileft_corner[k]%grid_m_ ==0) ? 0 : ileft_corner[k]%grid_m_;
    fprintf(stderr, "dim=%c : ileft = %d, ileft_corner %d, nx = %d\n", 'x' + k, ileft[k],ileft_corner[k],nx[k]);
  };

  int ileft_corner_m[3], ileft_corner_p[3],nx_m[3];
  int ileft_max_expand = std::max(iexpand_left[0],std::max(iexpand_left[1],iexpand_left[2]));

  for (int k = 0; k < 3; ++k) {
    ileft_corner_m[k] = ((ileft_corner[k] - iexpand_left[k]) + 
 		       coordinate_system_shift_[k] * (1 << (level - levelmin_final_)) + (1 << level)) % (1 << level);

    ileft_corner_p[k] =  grid_p_ * ileft_corner_m[k]/grid_m_;
    nx_m[k] = (ileft_max_expand!=0)? nx[k] + ileft_max_expand: nx[k];
    if (nx_m[k]%grid_m_ !=0) nx_m[k] = nx_m[k] + grid_m_ - nx_m[k]%grid_m_;
    nxremap[k] = grid_p_ * nx_m[k]/grid_m_;
    if (nxremap[k]%2==1){
      nx_m[k] = nx_m[k] + grid_m_;
      nxremap[k] = grid_p_ * nx_m[k]/grid_m_; 
     }
  }


  if ( (nx_m[0]!=nx_m[1]) || (nx_m[0]!=nx_m[2])) LOGERR("Fatal error: non-cubic refinement being requested");

  inter_grid_phase_adjustment_ = M_PI * (1.0 / (double)nx_m[0] - 1.0 / (double)nxremap[0]);
  LOGINFO("The value of the phase adjustement is %f\n", inter_grid_phase_adjustment_);


  LOGINFO("ileft[0],ileft[1],ileft[2] %d %d %d", ileft[0], ileft[1], ileft[2]);
  LOGINFO("ileft_corner[0,1,2] %d %d %d", ileft_corner[0], ileft_corner[1], ileft_corner[2]);

  LOGINFO("iexpand_left[1,2,3] = (%d, %d, %d) Max %d ",iexpand_left[0],iexpand_left[1],iexpand_left[2],
	  ileft_max_expand);

  LOGINFO("ileft_corner_m[0,1,2]  = (%d,%d,%d)",ileft_corner_m[0],ileft_corner_m[1],ileft_corner_m[2]); 
  LOGINFO("grid_m_ %d grid_p_ %d",grid_m_,grid_p_);
  LOGINFO("nx_m[0,1,2]  = (%d,%d,%d)",nx_m[0],nx_m[1],nx_m[2]); 
  LOGINFO("ileft_corner_p[0,1,2]  = (%d,%d,%d)",ileft_corner_p[0],ileft_corner_p[1],ileft_corner_p[2]);   
  LOGINFO("nxremap[0,1,2]  = (%d,%d,%d)",nxremap[0],nxremap[1],nxremap[2]);

  size_t ngp = nxremap[0] * nxremap[1] * (nxremap[2] + 2);

  pr0 = new fftw_real[ngp];
  pr1 = new fftw_real[ngp];
  pr2 = new fftw_real[ngp];
  pr3 = new fftw_real[ngp];
  pr4 = new fftw_real[ngp];

  pc0 = reinterpret_cast<fftw_complex *>(pr0);
  pc1 = reinterpret_cast<fftw_complex *>(pr1);
  pc2 = reinterpret_cast<fftw_complex *>(pr2);
  pc3 = reinterpret_cast<fftw_complex *>(pr3);
  pc4 = reinterpret_cast<fftw_complex *>(pr4);

  LOGINFO("calculating PANPHASIA random numbers for level %d...", level);
  clear_panphasia_thread_states();

  double t1 = get_wtime();
  double tp = t1;


#pragma omp parallel
  {
#ifdef _OPENMP
    const int mythread = omp_get_thread_num();
#else
    const int mythread = 0;
#endif
    int odd_x, odd_y, odd_z;
    int ng_level = ngrid_ * (1 << (level - levelmin_)); // full resolution of current level
  
    int verbosity = (mythread == 0);
    char descriptor[100];
    memset(descriptor, 0, 100);
    memcpy(descriptor, descriptor_string_.c_str(), descriptor_string_.size());

    if (level == levelmin_) {
      start_panphasia_(&lstate[mythread], descriptor, &ng_level, &verbosity);
    }

    {
      int level_p, lextra;
      long long ix_rel[3];
      panphasia_descriptor d(descriptor_string_);

      lextra = (log10((double)ng_level / (double)d.i_base) + 0.001) / log10(2.0);
      level_p = d.wn_level_base + lextra;
      int ratio = 1 << lextra;
      assert(ng_level == ratio * d.i_base);


 
        ix_rel[0] = ileft_corner_p[0];
        ix_rel[1] = ileft_corner_p[1];
        ix_rel[2] = ileft_corner_p[2];



// Code above ignores the coordinate_system_shift_ - but currently this is set to zero // 


      lstate[mythread].layer_min = 0;
      lstate[mythread].layer_max = level_p;
      lstate[mythread].indep_field = 1;

      set_phases_and_rel_origin_(&lstate[mythread], descriptor, &level_p, &ix_rel[0], &ix_rel[1], &ix_rel[2],
                                 &verbosity);

      LOGUSER(" called set_phases_and_rel_origin level %d ix_rel iy_rel iz_rel %d %d %d\n", level_p, ix_rel[0],
              ix_rel[1], ix_rel[2]);

      odd_x = ix_rel[0] % 2;
      odd_y = ix_rel[1] % 2;
      odd_z = ix_rel[2] % 2;
    }

    if (verbosity)
      t1 = get_wtime();

//***************************************************************
// Process Panphasia values: p000, p001, p010, p100 and indep field
//****************************************************************
//                START                                         //

#pragma omp for //nowait
    for (int i = 0; i < nxremap[0] / 2 + odd_x; ++i) {
      double cell_prop[9];
      pan_state_ *ps = &lstate[mythread];

      for (int j = 0; j < nxremap[1] / 2 + odd_y; ++j)
        for (int k = 0; k < nxremap[2] / 2 + odd_z; ++k) {

          // ARJ - added inner set of loops to speed up evaluation of Panphasia

          for (int ix = 0; ix < 2; ++ix)
            for (int iy = 0; iy < 2; ++iy)
              for (int iz = 0; iz < 2; ++iz) {
                int ii = 2 * i + ix - odd_x;
                int jj = 2 * j + iy - odd_y;
                int kk = 2 * k + iz - odd_z;

                if (((ii >= 0) && (ii < nxremap[0])) && ((jj >= 0) && (jj < nxremap[1])) &&
                    ((kk >= 0) && (kk < nxremap[2]))) {

                  size_t idx = ((size_t)ii * nxremap[1] + (size_t)jj) * (nxremap[2] + 2) + (size_t)kk;
                  adv_panphasia_cell_properties_(ps, &ii, &jj, &kk, &ps->layer_min, &ps->layer_max, &ps->indep_field,
                                                 cell_prop);

                  pr0[idx] = cell_prop[0];
                  pr1[idx] = cell_prop[4];
                  pr2[idx] = cell_prop[2];
                  pr3[idx] = cell_prop[1];
                  pr4[idx] = cell_prop[8];
                }
              }
        }
    }
  }
  LOGINFO("time for calculating PANPHASIA for level %d : %f s, %f µs/cell", level, get_wtime() - t1,
          1e6 * (get_wtime() - t1) / ((double)nxremap[2] * (double)nxremap[1] * (double)nxremap[0]));
  LOGINFO("time for calculating PANPHASIA for level %d : %f s, %f µs/cell", level, get_wtime() - t1,
          1e6 * (get_wtime() - t1) / ((double)nxremap[2] * (double)nxremap[1] * (double)nxremap[0]));

  //////////////////////////////////////////////////////////////////////////////////////////////

  LOGINFO("\033[31mtiming level %d [adv_panphasia_cell_properties]: %f s\033[0m", level, get_wtime() - tp);
  tp = get_wtime();

  /////////////////////////////////////////////////////////////////////////
  // transform and convolve with Legendres

  forward_transform_field(pr0, nxremap[0], nxremap[1], nxremap[2]);
  forward_transform_field(pr1, nxremap[0], nxremap[1], nxremap[2]);
  forward_transform_field(pr2, nxremap[0], nxremap[1], nxremap[2]);
  forward_transform_field(pr3, nxremap[0], nxremap[1], nxremap[2]);
  forward_transform_field(pr4, nxremap[0], nxremap[1], nxremap[2]);

#pragma omp parallel for
  for (int i = 0; i < nxremap[0]; i++)
    for (int j = 0; j < nxremap[1]; j++)
      for (int k = 0; k < nxremap[2] / 2 + 1; k++) {
        size_t idx = ((size_t)i * nxremap[1] + (size_t)j) * (nxremap[2] / 2 + 1) + (size_t)k;

        double fx(1.0), fy(1.0), fz(1.0), arg = 0.;
        complex gx(0., 0.), gy(0., 0.), gz(0., 0.);

        int ii(i), jj(j), kk(k);
        if (i > nxremap[0] / 2)
          ii -= nxremap[0];
        if (j > nxremap[1] / 2)
          jj -= nxremap[1];

        // int kkmax = std::max(abs(ii),std::max(abs(jj),abs(kk)));


        if (ii != 0) {
          arg = M_PI * (double)ii / (double)nxremap[0];
          fx = sin(arg) / arg;
          gx = complex(0.0, (arg * cos(arg) - sin(arg)) / (arg * arg));
        } else {
          fx = 1.0;
          gx = 0.0;
        }

        if (jj != 0) {
          arg = M_PI * (double)jj / (double)nxremap[1];
          fy = sin(arg) / arg;
          gy = complex(0.0, (arg * cos(arg) - sin(arg)) / (arg * arg));
        } else {
          fy = 1.0;
          gy = 0.0;
        }

        if (kk != 0) {
          arg = M_PI * (double)kk / (double)nxremap[2];
          fz = sin(arg) / arg;
          gz = complex(0.0, (arg * cos(arg) - sin(arg)) / (arg * arg));
        } else {
          fz = 1.0;
          gz = 0.0;
        }

        complex temp_comp = (fx + sqrt(3.0) * gx) * (fy + sqrt(3.0) * gy) * (fz + sqrt(3.0) * gz);
        double magnitude = sqrt(1.0 - std::abs(temp_comp * temp_comp));

        if (abs(ii) != nxremap[0] / 2 && abs(jj) != nxremap[1] / 2 &&
            abs(kk) != nxremap[2] / 2) { // kkmax != nxremap[2]/2 ){
          complex x, y0(RE(pc0[idx]), IM(pc0[idx])), y1(RE(pc1[idx]), IM(pc1[idx])), y2(RE(pc2[idx]), IM(pc2[idx])),
              y3(RE(pc3[idx]), IM(pc3[idx])), y4(RE(pc4[idx]), IM(pc4[idx]));

          x = y0 * fx * fy * fz + sqrt(3.0) * (y1 * gx * fy * fz + y2 * fx * gy * fz + y3 * fx * fy * gz) +
              y4 * magnitude;

          RE(pc0[idx]) = x.real();
          IM(pc0[idx]) = x.imag();
        }
      }

  //                END

  LOGINFO("\033[31mtiming level %d [build panphasia field]: %f s\033[0m", level, get_wtime() - tp);
  tp = get_wtime();

//***************************************************************
// Process Panphasia values: p000, p001, p010, p100 and indep field
//****************************************************************

#pragma omp parallel
  {
#ifdef _OPENMP
    const int mythread = omp_get_thread_num();
#else
    const int mythread = 0;
#endif
    int odd_x, odd_y, odd_z;
    int ng_level = ngrid_ * (1 << (level - levelmin_)); // full resolution of current level
    int verbosity = (mythread == 0);
    char descriptor[100];
    memset(descriptor, 0, 100);
    memcpy(descriptor, descriptor_string_.c_str(), descriptor_string_.size());

    if (level == levelmin_) {
      start_panphasia_(&lstate[mythread], descriptor, &ng_level, &verbosity);
    }

    {
      int level_p, lextra;
      long long ix_rel[3];
      panphasia_descriptor d(descriptor_string_);

      lextra = (log10((double)ng_level / (double)d.i_base) + 0.001) / log10(2.0);
      level_p = d.wn_level_base + lextra;
      int ratio = 1 << lextra;
      assert(ng_level == ratio * d.i_base);

      ix_rel[0] = ileft_corner_p[0];
      ix_rel[1] = ileft_corner_p[1];
      ix_rel[2] = ileft_corner_p[2];

// Code above ignores the coordinate_system_shift_ - but currently this is set to zero // 

      lstate[mythread].layer_min = 0;
      lstate[mythread].layer_max = level_p;
      lstate[mythread].indep_field = 1;

      set_phases_and_rel_origin_(&lstate[mythread], descriptor, &level_p, &ix_rel[0], &ix_rel[1], &ix_rel[2],
                                 &verbosity);

      LOGUSER(" called set_phases_and_rel_origin level %d ix_rel iy_rel iz_rel %d %d %d\n", level_p, ix_rel[0],
              ix_rel[1], ix_rel[2]);

      odd_x = ix_rel[0] % 2;
      odd_y = ix_rel[1] % 2;
      odd_z = ix_rel[2] % 2;
    }

    if (verbosity)
      t1 = get_wtime();

//                START                                         //
//***************************************************************
// Process Panphasia values: p110, p011, p101, p111
//****************************************************************
#pragma omp for //nowait
    for (int i = 0; i < nxremap[0] / 2 + odd_x; ++i) {
      double cell_prop[9];
      pan_state_ *ps = &lstate[mythread];

      for (int j = 0; j < nxremap[1] / 2 + odd_y; ++j)
        for (int k = 0; k < nxremap[2] / 2 + odd_z; ++k) {

          // ARJ - added inner set of loops to speed up evaluation of Panphasia

          for (int ix = 0; ix < 2; ++ix)
            for (int iy = 0; iy < 2; ++iy)
              for (int iz = 0; iz < 2; ++iz) {
                int ii = 2 * i + ix - odd_x;
                int jj = 2 * j + iy - odd_y;
                int kk = 2 * k + iz - odd_z;

                if (((ii >= 0) && (ii < nxremap[0])) && ((jj >= 0) && (jj < nxremap[1])) &&
                    ((kk >= 0) && (kk < nxremap[2]))) {

                  size_t idx = ((size_t)ii * nxremap[1] + (size_t)jj) * (nxremap[2] + 2) + (size_t)kk;
                  adv_panphasia_cell_properties_(ps, &ii, &jj, &kk, &ps->layer_min, &ps->layer_max, &ps->indep_field,
                                                 cell_prop);

                  pr1[idx] = cell_prop[6];
                  pr2[idx] = cell_prop[3];
                  pr3[idx] = cell_prop[5];
                  pr4[idx] = cell_prop[7];
                }
              }
        }
    }
  }
  LOGINFO("time for calculating PANPHASIA for level %d : %f s, %f µs/cell", level, get_wtime() - t1,
          1e6 * (get_wtime() - t1) / ((double)nxremap[2] * (double)nxremap[1] * (double)nxremap[0]));

  LOGINFO("\033[31mtiming level %d [adv_panphasia_cell_properties2]: %f s \033[0m", level, get_wtime() - tp);
  tp = get_wtime();

  /////////////////////////////////////////////////////////////////////////
  // transform and convolve with Legendres

  forward_transform_field(pr1, nxremap[0], nxremap[1], nxremap[2]);
  forward_transform_field(pr2, nxremap[0], nxremap[1], nxremap[2]);
  forward_transform_field(pr3, nxremap[0], nxremap[1], nxremap[2]);
  forward_transform_field(pr4, nxremap[0], nxremap[1], nxremap[2]);

#pragma omp parallel for
  for (int i = 0; i < nxremap[0]; i++)
    for (int j = 0; j < nxremap[1]; j++)
      for (int k = 0; k < nxremap[2] / 2 + 1; k++) {
        size_t idx = ((size_t)i * nxremap[1] + (size_t)j) * (nxremap[2] / 2 + 1) + (size_t)k;

        double fx(1.0), fy(1.0), fz(1.0), arg = 0.;
        complex gx(0., 0.), gy(0., 0.), gz(0., 0.);

        int ii(i), jj(j), kk(k);
        if (i > nxremap[0] / 2)
          ii -= nxremap[0];
        if (j > nxremap[1] / 2)
          jj -= nxremap[1];

        // int kkmax = std::max(abs(ii),std::max(abs(jj),abs(kk)));

        if (ii != 0) {
          arg = M_PI * (double)ii / (double)nxremap[0];
          fx = sin(arg) / arg;
          gx = complex(0.0, (arg * cos(arg) - sin(arg)) / (arg * arg));
        } else {
          fx = 1.0;
          gx = 0.0;
        }

        if (jj != 0) {
          arg = M_PI * (double)jj / (double)nxremap[1];
          fy = sin(arg) / arg;
          gy = complex(0.0, (arg * cos(arg) - sin(arg)) / (arg * arg));
        } else {
          fy = 1.0;
          gy = 0.0;
        }

        if (kk != 0) {
          arg = M_PI * (double)kk / (double)nxremap[2];
          fz = sin(arg) / arg;
          gz = complex(0.0, (arg * cos(arg) - sin(arg)) / (arg * arg));
        } else {
          fz = 1.0;
          gz = 0.0;
        }

        if (abs(ii) != nxremap[0] / 2 && abs(jj) != nxremap[1] / 2 &&
            abs(kk) != nxremap[2] / 2) { // kkmax != nxremap[2]/2 ){
          complex x, y1(RE(pc1[idx]), IM(pc1[idx])), y2(RE(pc2[idx]), IM(pc2[idx])), y3(RE(pc3[idx]), IM(pc3[idx])),
              y4(RE(pc4[idx]), IM(pc4[idx]));

          x = 3.0 * (y1 * gx * gy * fz + y2 * fx * gy * gz + y3 * gx * fy * gz) + sqrt(27.0) * y4 * gx * gy * gz;

          RE(pc0[idx]) = RE(pc0[idx]) + x.real();
          IM(pc0[idx]) = IM(pc0[idx]) + x.imag();
        }
      }

  LOGINFO("\033[31mtiming level %d [build panphasia field2]: %f s\033[0m", level, get_wtime() - tp);
  tp = get_wtime();

  //                END
  //***************************************************************
  // Compute Panphasia values of p011, p101, p110, p111 coefficients
  // and combine with p000, p001, p010, p100 and indep field.
  //****************************************************************

  /////////////////////////////////////////////////////////////////////////
  // do we need to cut off the small scales?
  //   int nn = 1<<level;

  if (incongruent_fields_) {

    LOGINFO("Remapping fields from dimension %d -> %d",nxremap[0],nx_m[0]);
    memset(pr1, 0, ngp * sizeof(fftw_real));

#pragma omp parallel for
    for (int i = 0; i < nxremap[0]; i++)
      for (int j = 0; j < nxremap[1]; j++)
        for (int k = 0; k < nxremap[2] / 2 + 1; k++) {

          int ii = (i > nxremap[0] / 2) ? i - nxremap[0] : i, jj = (j > nxremap[1] / 2) ? j - nxremap[1] : j, kk = k;

          int ia(abs(ii)), ja(abs(jj)), ka(abs(kk));

          if (ia < nx_m[0] / 2 && ja < nx_m[1] / 2 && ka < nx_m[2] / 2) {

            size_t idx = ((size_t)(i)*nxremap[1] + (size_t)(j)) * (nxremap[2] / 2 + 1) + (size_t)(k);

            int ir = (ii < 0) ? ii + nx_m[0] : ii, jr = (jj < 0) ? jj + nx_m[1] : jj, kr = kk; // never negative

            size_t idx2 = ((size_t)ir * nx_m[1] + (size_t)jr) * ((size_t)nx_m[2] / 2 + 1) + (size_t)kr;


           complex x(RE(pc0[idx]),IM(pc0[idx]));
           double total_phase_shift;
           total_phase_shift =  inter_grid_phase_adjustment_ * (double)(ii+jj+kk);
           x = x * exp(complex(0.0, total_phase_shift)); 
           RE(pc1[idx2]) = x.real(); 
           IM(pc1[idx2]) = x.imag();
           }
        }

    memcpy(pr0, pr1, ngp * sizeof(fftw_real));
  }


  if( level == 9 ){
    LOGINFO("DC mode of level is %g",RE(pc0[0]));
    //RE(pc0[0]) = 1e8;
    //IM(pc0[0]) = 0.0;
  }

  LOGINFO("\033[31mtiming level %d [remap noncongruent]: %f s\033[0m", level, get_wtime() - tp);
  tp = get_wtime();
  /////////////////////////////////////////////////////////////////////////
  // transform back

  backward_transform_field(pr0, nx_m[0], nx_m[1], nx_m[2]);

  /////////////////////////////////////////////////////////////////////////
  // copy to random data structure
  delete[] pr1;
  delete[] pr2;
  delete[] pr3;
  delete[] pr4;

  LOGINFO("Copying random field data %d,%d,%d -> %d,%d,%d", nxremap[0], nxremap[1], nxremap[2], nx[0], nx[1], nx[2]);
  
  //    n = 1<<level;
  //    ng = n;
  //    ngp = ng*ng*2*(ng/2+1);

  double sum = 0.0, sum2 = 0.0;
  size_t count = 0;

  /*double norm = 1.0 / sqrt((double)nxremap[0] * (double)nxremap[1] * (double)nxremap[2] * (double)nx[0] *
                           (double)nx[1] * (double)nx[2]);*/

 double norm = 1.0 / sqrt((double)nxremap[0] * (double)nxremap[1] * (double)nxremap[2] * (double)nx_m[0] *
                          (double)nx_m[1] * (double)nx_m[2]);

#pragma omp parallel for reduction(+ : sum, sum2, count)
  for (int k = 0; k < nx[2]; ++k) // ARJ - swapped roles of i,k, and reverse ordered loops
    for (int j = 0; j < nx[1]; ++j)
      for (int i = 0; i < nx[0]; ++i) {
        size_t idx = ((size_t)(i+iexpand_left[0])*nx_m[1] + (size_t)(j+iexpand_left[1])) * (nx_m[2] + 2) 
     + (size_t)(k+iexpand_left[2]);
        R(i, j, k) = pr0[idx] * norm;

        sum += R(i, j, k);
        sum2 += R(i, j, k) * R(i, j, k);
        ++count;
      }

  delete[] pr0;

  sum /= (double)count;
  sum2 /= (double)count;

  sum2 = (sum2 - sum * sum);

  LOGUSER("done with PANPHASIA for level %d:\n       mean=%g, var=%g", level, sum, sum2);
  LOGUSER("Copying into R array: nx[0],nx[1],nx[2] %d %d %d \n", nx[0], nx[1], nx[2]);

  LOGINFO("PANPHASIA level %d mean and variance are\n       <p> = %g | var(p) = %g", level, sum, sum2);
}

namespace {
RNG_plugin_creator_concrete<RNG_panphasia> creator("PANPHASIA");
}

#endif
