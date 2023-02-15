#include "random.hh"
#include "random_music_wnoise_generator.hh"

typedef music_wnoise_generator<real_t> rng;

class RNG_music : public RNG_plugin
{
protected:
  std::vector<long> rngseeds_;
  std::vector<std::string> rngfnames_;
  unsigned ran_cube_size_;

  int levelmin_, levelmax_, levelmin_seed_;

  bool disk_cached_;
  bool restart_;
  bool initialized_;

  std::vector<std::vector<real_t> *> mem_cache_;

  //! checks if the specified string is numeric
  bool is_number(const std::string &s);

  //! parses the random number parameters in the conf file
  void parse_random_parameters(void);

  //! computes the white noise fields and keeps them either in memory or on disk
  void compute_random_numbers(void);

  //! adjusts averages
  void correct_avg(int icoarse, int ifine);

  //! store the white noise fields in memory or on disk
  void store_rnd(int ilevel, rng *prng);

public:
  explicit RNG_music(config_file &cf) : RNG_plugin(cf), initialized_(false) {}

  ~RNG_music() {}

  bool is_multiscale() const { return true; }

  void initialize_for_grid_structure(const refinement_hierarchy &refh)
  {
    prefh_ = &refh;
    levelmin_ = prefh_->levelmin();
    levelmax_ = prefh_->levelmax();

    ran_cube_size_ = pcf_->get_value_safe<unsigned>("random", "cubesize", DEF_RAN_CUBE_SIZE);
    disk_cached_ = pcf_->get_value_safe<bool>("random", "disk_cached", true);
    restart_ = pcf_->get_value_safe<bool>("random", "restart", false);

    pcf_->insert_value("setup","fourier_splicing","true");

    mem_cache_.assign(levelmax_ - levelmin_ + 1, (std::vector<real_t> *)NULL);

    if (restart_ && !disk_cached_)
    {
      music::elog.Print("Cannot restart from mem cached random numbers.");
      throw std::runtime_error("Cannot restart from mem cached random numbers.");
    }

    //... determine seed/white noise file data to be applied
    parse_random_parameters();

    if (!restart_)
    {
      //... compute the actual random numbers
      compute_random_numbers();
    }

    initialized_ = true;
  }

  void fill_grid(int level, DensityGrid<real_t> &R);
};

bool RNG_music::is_number(const std::string &s)
{
  for (size_t i = 0; i < s.length(); i++)
    if (!std::isdigit(s[i]) && s[i] != '-')
      return false;

  return true;
}

void RNG_music::parse_random_parameters(void)
{
  //... parse random number options
  for (int i = 0; i <= 100; ++i)
  {
    char seedstr[128];
    std::string tempstr;
    bool noseed = false;
    snprintf(seedstr,128, "seed[%d]", i);
    if (pcf_->contains_key("random", seedstr))
      tempstr = pcf_->get_value<std::string>("random", seedstr);
    else
    {
      // "-2" means that no seed entry was found for that level
      tempstr = std::string("-2");
      noseed = true;
    }

    if (is_number(tempstr))
    {
      long ltemp;
      pcf_->convert(tempstr, ltemp);
      rngfnames_.push_back("");
      if (noseed) // ltemp < 0 )
        //... generate some dummy seed which only depends on the level, negative so we know it's not
        //... an actual seed and thus should not be used as a constraint for coarse levels
        // rngseeds_.push_back( -abs((unsigned)(ltemp-i)%123+(unsigned)(ltemp+827342523521*i)%123456789) );
        rngseeds_.push_back(-abs((long)(ltemp - i) % 123 + (long)(ltemp + 7342523521 * i) % 123456789));
      else
      {
        if (ltemp <= 0)
        {
          music::elog.Print("Specified seed [random]/%s needs to be a number >0!", seedstr);
          throw std::runtime_error("Seed values need to be >0");
        }
        rngseeds_.push_back(ltemp);
      }
    }
    else
    {
      rngfnames_.push_back(tempstr);
      rngseeds_.push_back(-1);
      music::ilog.Print("Random numbers for level %3d will be read from file.", i);
    }
  }

  //.. determine for which levels random seeds/random number files are given
  levelmin_seed_ = -1;
  for (unsigned ilevel = 0; ilevel < rngseeds_.size(); ++ilevel)
  {
    if (levelmin_seed_ < 0 && (rngfnames_[ilevel].size() > 0 || rngseeds_[ilevel] >= 0))
      levelmin_seed_ = ilevel;
  }
}

void RNG_music::compute_random_numbers(void)
{
  bool rndsign = pcf_->get_value_safe<bool>("random", "grafic_sign", false);

  std::vector<rng *> randc(std::max(levelmax_, levelmin_seed_) + 1, (rng *)NULL);

  //--- FILL ALL WHITE NOISE ARRAYS FOR WHICH WE NEED THE FULL FIELD ---//

  //... seeds are given for a level coarser than levelmin
  if (levelmin_seed_ < levelmin_)
  {
    if (rngfnames_[levelmin_seed_].size() > 0)
      randc[levelmin_seed_] = new rng(1 << levelmin_seed_, rngfnames_[levelmin_seed_], rndsign);
    else
      randc[levelmin_seed_] = new rng(1 << levelmin_seed_, ran_cube_size_, rngseeds_[levelmin_seed_]);

    for (int i = levelmin_seed_ + 1; i <= levelmin_; ++i)
    {
      // #warning add possibility to read noise from file also here!

      if (rngfnames_[i].size() > 0)
        music::ilog.Print("Warning: Cannot use filenames for higher levels currently! Ignoring!");

      randc[i] = new rng(*randc[i - 1], ran_cube_size_, rngseeds_[i]);
      delete randc[i - 1];
      randc[i - 1] = NULL;
    }
  }

  //... seeds are given for a level finer than levelmin, obtain by averaging
  if (levelmin_seed_ > levelmin_)
  {
    if (rngfnames_[levelmin_seed_].size() > 0)
      randc[levelmin_seed_] = new rng(1 << levelmin_seed_, rngfnames_[levelmin_seed_], rndsign);
    else
      randc[levelmin_seed_] =
          new rng(1 << levelmin_seed_, ran_cube_size_, rngseeds_[levelmin_seed_]); //, x0, lx );

    for (int ilevel = levelmin_seed_ - 1; ilevel >= (int)levelmin_; --ilevel)
    {
      if (rngseeds_[ilevel - levelmin_] > 0)
        music::ilog.Print("Warning: random seed for level %d will be ignored.\n"
                "            consistency requires that it is obtained by restriction from level %d",
                ilevel, levelmin_seed_);

      randc[ilevel] = new rng(*randc[ilevel + 1]);

      if (ilevel + 1 > levelmax_)
      {
        delete randc[ilevel + 1];
        randc[ilevel + 1] = NULL;
      }
    }
  }

  //--- GENERATE AND STORE ALL LEVELS, INCLUDING REFINEMENTS ---//

  //... levelmin
  if (randc[levelmin_] == NULL)
  {
    if (rngfnames_[levelmin_].size() > 0)
      randc[levelmin_] = new rng(1 << levelmin_, rngfnames_[levelmin_], rndsign);
    else
      randc[levelmin_] = new rng(1 << levelmin_, ran_cube_size_, rngseeds_[levelmin_]);
  }

  store_rnd(levelmin_, randc[levelmin_]);

  //... refinement levels
  for (int ilevel = levelmin_ + 1; ilevel <= levelmax_; ++ilevel)
  {
    int lx[3], x0[3];
    int shift[3], levelmin_poisson;
    shift[0] = pcf_->get_value<int>("setup", "shift_x");
    shift[1] = pcf_->get_value<int>("setup", "shift_y");
    shift[2] = pcf_->get_value<int>("setup", "shift_z");

    levelmin_poisson = pcf_->get_value<unsigned>("setup", "levelmin");

    int lfac = 1 << (ilevel - levelmin_poisson);

    std::array<int,3> margin;
    if( prefh_->get_margin()>0 ){
      margin[0] = prefh_->get_margin();
      margin[1] = prefh_->get_margin();
      margin[2] = prefh_->get_margin();
    }else{
      margin[0] = prefh_->size(ilevel, 0)/2;
      margin[1] = prefh_->size(ilevel, 1)/2;
      margin[2] = prefh_->size(ilevel, 2)/2;
    }

    lx[0] = prefh_->size(ilevel, 0) + 2 * margin[0];
    lx[1] = prefh_->size(ilevel, 1) + 2 * margin[1];
    lx[2] = prefh_->size(ilevel, 2) + 2 * margin[2];
    x0[0] = prefh_->offset_abs(ilevel, 0) - lfac * shift[0] - margin[0];
    x0[1] = prefh_->offset_abs(ilevel, 1) - lfac * shift[1] - margin[1];
    x0[2] = prefh_->offset_abs(ilevel, 2) - lfac * shift[2] - margin[2];

    if (randc[ilevel] == NULL)
      randc[ilevel] =
          new rng(*randc[ilevel - 1], ran_cube_size_, rngseeds_[ilevel], x0, lx);
    delete randc[ilevel - 1];
    randc[ilevel - 1] = NULL;

    //... store numbers
    store_rnd(ilevel, randc[ilevel]);
  }

  delete randc[levelmax_];
  randc[levelmax_] = NULL;
}

void RNG_music::store_rnd(int ilevel, rng *prng)
{
  int shift[3], levelmin_poisson;
  shift[0] = pcf_->get_value<int>("setup", "shift_x");
  shift[1] = pcf_->get_value<int>("setup", "shift_y");
  shift[2] = pcf_->get_value<int>("setup", "shift_z");

  levelmin_poisson = pcf_->get_value<unsigned>("setup", "levelmin");

  int lfac = 1 << (ilevel - levelmin_poisson);

  bool grafic_out = false;

  if (grafic_out)
  {
    std::vector<float> data;
    if (ilevel == levelmin_)
    {
      int N = 1 << levelmin_;
      int i0, j0, k0;
      i0 = -lfac * shift[0];
      j0 = -lfac * shift[1];
      k0 = -lfac * shift[2];

      char fname[128];
      snprintf(fname, 128, "grafic_wnoise_%04d.bin", ilevel);

      music::ulog.Print("Storing white noise field for grafic in file \'%s\'...", fname);

      std::ofstream ofs(fname, std::ios::binary | std::ios::trunc);
      data.assign(N * N, 0.0);

      int blksize = 4 * sizeof(int);
      int iseed = 0;

      ofs.write(reinterpret_cast<char *>(&blksize), sizeof(int));
      ofs.write(reinterpret_cast<char *>(&N), sizeof(int));
      ofs.write(reinterpret_cast<char *>(&N), sizeof(int));
      ofs.write(reinterpret_cast<char *>(&N), sizeof(int));
      ofs.write(reinterpret_cast<char *>(&iseed), sizeof(int));
      ofs.write(reinterpret_cast<char *>(&blksize), sizeof(int));

      for (int k = 0; k < N; ++k)
      {
        #pragma omp parallel for
        for (int j = 0; j < N; ++j)
          for (int i = 0; i < N; ++i)
            data[j * N + i] = -(*prng)(i + i0, j + j0, k + k0);

        blksize = N * N * sizeof(float);
        ofs.write(reinterpret_cast<char *>(&blksize), sizeof(int));
        ofs.write(reinterpret_cast<char *>(&data[0]), N * N * sizeof(float));
        ofs.write(reinterpret_cast<char *>(&blksize), sizeof(int));
      }

      ofs.close();
    }
    else
    {

      int nx, ny, nz;
      int i0, j0, k0;

      nx = prefh_->size(ilevel, 0);
      ny = prefh_->size(ilevel, 1);
      nz = prefh_->size(ilevel, 2);
      i0 = prefh_->offset_abs(ilevel, 0) - lfac * shift[0];
      j0 = prefh_->offset_abs(ilevel, 1) - lfac * shift[1];
      k0 = prefh_->offset_abs(ilevel, 2) - lfac * shift[2];

      char fname[128];
      snprintf(fname, 128,"grafic_wnoise_%04d.bin", ilevel);

      music::ulog.Print("Storing white noise field for grafic in file \'%s\'...", fname);
      music::dlog.Print("(%d,%d,%d) -- (%d,%d,%d) -- lfac = %d", nx, ny, nz, i0, j0, k0, lfac);

      std::ofstream ofs(fname, std::ios::binary | std::ios::trunc);
      data.assign(nx * ny, 0.0);

      int blksize = 4 * sizeof(int);
      int iseed = 0;

      ofs.write(reinterpret_cast<char *>(&blksize), sizeof(int));
      ofs.write(reinterpret_cast<char *>(&nz), sizeof(unsigned));
      ofs.write(reinterpret_cast<char *>(&ny), sizeof(unsigned));
      ofs.write(reinterpret_cast<char *>(&nx), sizeof(unsigned));
      ofs.write(reinterpret_cast<char *>(&iseed), sizeof(int));
      ofs.write(reinterpret_cast<char *>(&blksize), sizeof(int));

      for (int k = 0; k < nz; ++k)
      {
        #pragma omp parallel for
        for (int j = 0; j < ny; ++j)
          for (int i = 0; i < nx; ++i)
            data[j * nx + i] = -(*prng)(i + i0, j + j0, k + k0);

        blksize = nx * ny * sizeof(float);
        ofs.write(reinterpret_cast<char *>(&blksize), sizeof(int));
        ofs.write(reinterpret_cast<char *>(&data[0]), nx * ny * sizeof(float));
        ofs.write(reinterpret_cast<char *>(&blksize), sizeof(int));
      }
      ofs.close();
    }
  }

  if (disk_cached_)
  {
    std::vector<real_t> data;
    if (ilevel == levelmin_)
    {
      int N = 1 << levelmin_;
      int i0, j0, k0;

      i0 = -lfac * shift[0];
      j0 = -lfac * shift[1];
      k0 = -lfac * shift[2];

      char fname[128];
      snprintf(fname, 128,"wnoise_%04d.bin", ilevel);

      music::ulog.Print("Storing white noise field in file \'%s\'...", fname);

      std::ofstream ofs(fname, std::ios::binary | std::ios::trunc);

      ofs.write(reinterpret_cast<char *>(&N), sizeof(unsigned));
      ofs.write(reinterpret_cast<char *>(&N), sizeof(unsigned));
      ofs.write(reinterpret_cast<char *>(&N), sizeof(unsigned));

      data.assign(N * N, 0.0);
      for (int i = 0; i < N; ++i)
      {
        #pragma omp parallel for
        for (int j = 0; j < N; ++j)
          for (int k = 0; k < N; ++k)
            data[j * N + k] = (*prng)(i + i0, j + j0, k + k0);

        ofs.write(reinterpret_cast<char *>(&data[0]), N * N * sizeof(real_t));
      }
      ofs.close();
    }
    else
    {
      int nx, ny, nz;
      int i0, j0, k0;

      std::array<int,3> margin;
      if( prefh_->get_margin()>0 ){
        margin[0] = prefh_->get_margin();
        margin[1] = prefh_->get_margin();
        margin[2] = prefh_->get_margin();
      }else{
        margin[0] = prefh_->size(ilevel, 0)/2;
        margin[1] = prefh_->size(ilevel, 1)/2;
        margin[2] = prefh_->size(ilevel, 2)/2;
      }

      nx = prefh_->size(ilevel, 0) + 2 * margin[0];
      ny = prefh_->size(ilevel, 1) + 2 * margin[1];
      nz = prefh_->size(ilevel, 2) + 2 * margin[2];
      i0 = prefh_->offset_abs(ilevel, 0) - lfac * shift[0] - margin[0];
      j0 = prefh_->offset_abs(ilevel, 1) - lfac * shift[1] - margin[1];
      k0 = prefh_->offset_abs(ilevel, 2) - lfac * shift[2] - margin[2];

      char fname[128];
      snprintf(fname, 128,"wnoise_%04d.bin", ilevel);

      music::ulog.Print("Storing white noise field in file \'%s\'...", fname);

      std::ofstream ofs(fname, std::ios::binary | std::ios::trunc);

      ofs.write(reinterpret_cast<char *>(&nx), sizeof(unsigned));
      ofs.write(reinterpret_cast<char *>(&ny), sizeof(unsigned));
      ofs.write(reinterpret_cast<char *>(&nz), sizeof(unsigned));

      data.assign(ny * nz, 0.0);
      for (int i = 0; i < nx; ++i)
      {
        #pragma omp parallel for
        for (int j = 0; j < ny; ++j)
          for (int k = 0; k < nz; ++k)
            data[j * nz + k] = (*prng)(i + i0, j + j0, k + k0);

        ofs.write(reinterpret_cast<char *>(&data[0]), ny * nz * sizeof(real_t));
      }
      ofs.close();
    }
  }
  else
  {
    int nx, ny, nz;
    int i0, j0, k0;

    if (ilevel == levelmin_)
    {
      i0 = -lfac * shift[0];
      j0 = -lfac * shift[1];
      k0 = -lfac * shift[2];

      nx = ny = nz = 1 << levelmin_;
    }
    else
    {
      std::array<int,3> margin;
      if( prefh_->get_margin()>0 ){
        margin[0] = prefh_->get_margin();
        margin[1] = prefh_->get_margin();
        margin[2] = prefh_->get_margin();
      }else{
        margin[0] = prefh_->size(ilevel, 0)/2;
        margin[1] = prefh_->size(ilevel, 1)/2;
        margin[2] = prefh_->size(ilevel, 2)/2;
      }
      nx = prefh_->size(ilevel, 0) + 2 * margin[0];
      ny = prefh_->size(ilevel, 1) + 2 * margin[1];
      nz = prefh_->size(ilevel, 2) + 2 * margin[2];
      i0 = prefh_->offset_abs(ilevel, 0) - lfac * shift[0] - margin[0];
      j0 = prefh_->offset_abs(ilevel, 1) - lfac * shift[1] - margin[1]; 
      k0 = prefh_->offset_abs(ilevel, 2) - lfac * shift[2] - margin[2];
    }

    mem_cache_[ilevel - levelmin_] = new std::vector<real_t>(nx * ny * nz, 0.0);

    music::ulog.Print("Copying white noise to mem cache...");

    #pragma omp parallel for
    for (int i = 0; i < nx; ++i)
      for (int j = 0; j < ny; ++j)
        for (int k = 0; k < nz; ++k)
          (*mem_cache_[ilevel - levelmin_])[((size_t)i * ny + (size_t)j) * nz + (size_t)k] =
              (*prng)(i + i0, j + j0, k + k0);
  }
}

void RNG_music::fill_grid(int ilevel, DensityGrid<real_t> &A)
{
  if (!initialized_)
  {
    music::elog.Print("Call to RNG_music::fill_grid before call to RNG_music::initialize_for_grid_structure");
    throw std::runtime_error("invalid call order for random number generator");
  }

  if (restart_)
    music::ilog.Print("Attempting to restart using random numbers for level %d\n     from file \'wnoise_%04d.bin\'.", ilevel,
            ilevel);

  if (disk_cached_)
  {
    char fname[128];
    snprintf(fname, 128,"wnoise_%04d.bin", ilevel);

    music::ulog.Print("Loading white noise from file \'%s\'...", fname);

    std::ifstream ifs(fname, std::ios::binary);
    if (!ifs.good())
    {
      music::elog.Print("White noise file \'%s\'was not found.", fname);
      throw std::runtime_error("A white noise file was not found. This is an internal inconsistency and bad.");
    }

    int nx, ny, nz;
    ifs.read(reinterpret_cast<char *>(&nx), sizeof(int));
    ifs.read(reinterpret_cast<char *>(&ny), sizeof(int));
    ifs.read(reinterpret_cast<char *>(&nz), sizeof(int));

    if (nx != (int)A.size(0) || ny != (int)A.size(1) || nz != (int)A.size(2))
    {
      std::array<int,3> margin;
      if( prefh_->get_margin()>0 ){
        margin[0] = prefh_->get_margin();
        margin[1] = prefh_->get_margin();
        margin[2] = prefh_->get_margin();
      }else{
        margin[0] = prefh_->size(ilevel, 0)/2;
        margin[1] = prefh_->size(ilevel, 1)/2;
        margin[2] = prefh_->size(ilevel, 2)/2;
      }

      if (nx == (int)A.size(0) + 2 * margin[0] && ny == (int)A.size(1) + 2 * margin[1] && nz == (int)A.size(2) + 2 * margin[2])
      {
        int ox = margin[0], oy = margin[1], oz = margin[2];
        std::vector<real_t> slice(ny * nz, 0.0);

        for (int i = 0; i < nx; ++i)
        {
          ifs.read(reinterpret_cast<char *>(&slice[0]), ny * nz * sizeof(real_t));

          if (i < ox)
            continue;
          if (i >= 3 * ox)
            break;

          #pragma omp parallel for
          for (int j = oy; j < 3 * oy; ++j)
            for (int k = oz; k < 3 * oz; ++k)
              A(i - ox, j - oy, k - oz) = slice[j * nz + k];
        }

        ifs.close();
      }
      else
      {
        music::elog.Print("White noise file is not aligned with array. File: [%d,%d,%d]. Mem: [%d,%d,%d].", nx, ny, nz, A.size(0),
               A.size(1), A.size(2));
        throw std::runtime_error(
            "White noise file is not aligned with array. This is an internal inconsistency and bad.");
      }
    }
    else
    {

      for (int i = 0; i < nx; ++i)
      {
        std::vector<real_t> slice(ny * nz, 0.0);
        ifs.read(reinterpret_cast<char *>(&slice[0]), ny * nz * sizeof(real_t));

        #pragma omp parallel for
        for (int j = 0; j < ny; ++j)
          for (int k = 0; k < nz; ++k)
            A(i, j, k) = slice[j * nz + k];
      }
      ifs.close();
    }
  }
  else
  {
    music::ulog.Print("Copying white noise from memory cache...");

    if (mem_cache_[ilevel - levelmin_] == NULL)
      music::elog.Print("Tried to access mem-cached random numbers for level %d. But these are not available!\n", ilevel);

    int nx(A.size(0)), ny(A.size(1)), nz(A.size(2));

    if ((size_t)nx * (size_t)ny * (size_t)nz != mem_cache_[ilevel - levelmin_]->size())
    {
      music::elog.Print("White noise file is not aligned with array. File: [%d,%d,%d]. Mem: [%d,%d,%d].", nx, ny, nz, A.size(0),
             A.size(1), A.size(2));
      throw std::runtime_error("White noise file is not aligned with array. This is an internal inconsistency and bad");
    }

    #pragma omp parallel for
    for (int i = 0; i < nx; ++i)
      for (int j = 0; j < ny; ++j)
        for (int k = 0; k < nz; ++k)
          A(i, j, k) = (*mem_cache_[ilevel - levelmin_])[((size_t)i * ny + (size_t)j) * nz + (size_t)k];

    std::vector<real_t>().swap(*mem_cache_[ilevel - levelmin_]);
    delete mem_cache_[ilevel - levelmin_];
    mem_cache_[ilevel - levelmin_] = NULL;
  }
}

namespace
{
  RNG_plugin_creator_concrete<RNG_music> creator("MUSIC");
}
