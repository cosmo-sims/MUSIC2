/*

 output_gamer2.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions
 for cosmological simulations

 Copyright (C) 2010-13  Oliver Hahn

 Plugin: Yuri Oku (yuri.oku.astro@gmail.com)

 */

#include "output.hh"

class gamer2_output_plugin : public output_plugin
{
  protected:
    enum iofields
    {
        id_dm_mass,
        id_dm_pos,
        id_dm_vel,
        id_gas_rho,
        id_gas_vel,
    };

    int PATCH_SIZE, PS2;
    std::vector<int> refine_offset; // offset from left or right edge of simulation volume (the number of patch in parent level)
    std::vector<int> refine_size;   // the size of refinement region (the number of cell on the level)

    // round up to the nearest multiple of b
    int updiv(int a, int b)
    {
        if (a % b == 0)
            return a / b;
        else
            return a / b + 1;
    }

    std::string par_ic_;
    std::string um_ic_;

    double boxlength_;
    double zstart_;
    double omegam_;
    double omegab_;
    double hubble_;
    bool bbaryons_;

    double gamer_unit_length_;
    double gamer_unit_time_;
    double gamer_unit_density_;
    double gamer_unit_mass_;
    double gamer_unit_velocity_;

    double music_unit_length_;
    double music_unit_mass_;
    double music_unit_velocity_;

    size_t write2tempfile_par(std::string fname, const grid_hierarchy &gh, unsigned ilevel, real_t fac = 1.0, real_t shift = 0.0)
    {
        const MeshvarBnd<real_t> *data = gh.get_grid(ilevel);

        int n0 = data->size(0), n1 = data->size(1), n2 = data->size(2);
        std::vector<real_t> vdata;
        vdata.reserve((unsigned)(n0) * (n1) * (n2));

        size_t count = 0;
        for (int i = 0; i < n0; ++i)
            for (int j = 0; j < n1; ++j)
                for (int k = 0; k < n2; ++k)
                    if (gh.is_in_mask(ilevel, i, j, k) && !gh.is_refined(ilevel, i, j, k))
                    {
                        vdata.push_back((*data)(i, j, k) * fac + shift);

                        count++;
                    }

        std::ofstream ofs_temp(fname, std::ios::binary | std::ios::trunc);
        ofs_temp.write((char *)&vdata[0], vdata.size() * sizeof(real_t));

        ofs_temp.flush();
        ofs_temp.close();

        return count;
    }

    size_t write2tempfile_grid(std::string fname, const grid_hierarchy &gh, unsigned ilevel, real_t fac = 1.0, real_t shift = 0.0)
    {
        const MeshvarBnd<real_t> *data = gh.get_grid(ilevel);

        int n0 = data->size(0), n1 = data->size(1), n2 = data->size(2);
        std::vector<real_t> vdata;
        vdata.reserve((unsigned)(n0) * (n1) * (n2));

        size_t count = 0;
        for (int k = 0; k < n2; ++k)
            for (int j = 0; j < n1; ++j)
                for (int i = 0; i < n0; ++i)
                {
                    vdata.push_back((*data)(i, j, k) * fac + shift);

                    count++;
                }

        std::ofstream ofs_temp(fname, std::ios::binary | std::ios::trunc);
        ofs_temp.write((char *)&vdata[0], vdata.size() * sizeof(real_t));

        ofs_temp.flush();
        ofs_temp.close();

        return count;
    }

    size_t copy2outputfile(std::ofstream &ofs, std::string temp_fname)
    {
        std::ifstream ifs(temp_fname, std::ios::binary);
        if (!ifs)
        {
            std::cerr << "Error: Could not open file " << temp_fname << std::endl;
            exit(1);
        }

        ofs << ifs.rdbuf();

        /* count elements */
        ifs.seekg(0, std::ios::end);
        size_t count = ifs.tellg() / sizeof(real_t);

        ifs.close();
        return count;
    }

    void assemble_gamer2_file()
    {
        /* ----- write particle file ----- */
        std::ofstream ofs_par(par_ic_.c_str(), std::ios::binary | std::ios::trunc);
        std::cout << "=============================================================\n";
        music::ilog.Print("GAMER-2 plugin: writing %s", par_ic_.c_str());

        /* write mass */
        size_t npart = 0;
        for (unsigned ilevel = levelmin_; ilevel <= levelmax_; ++ilevel)
        {
            std::string temp_fname = "___ic_temp_" + std::to_string(1000 * id_dm_mass + ilevel) + ".bin";

            size_t npart_ilevel = copy2outputfile(ofs_par, temp_fname);

            music::ilog.Print("Level %2d: %12llu particles", ilevel, npart_ilevel);
            npart += npart_ilevel;

            remove(temp_fname.c_str());
        }

        music::ilog.Print("Total: %12llu particles", npart);

        /* write position */
        for (unsigned coord = 0; coord < 3; ++coord)
        {
            for (unsigned ilevel = levelmin_; ilevel <= levelmax_; ++ilevel)
            {
                std::string temp_fname = "___ic_temp_" + std::to_string(1000 * id_dm_pos + 100 * coord + ilevel) + ".bin";
                copy2outputfile(ofs_par, temp_fname);
                remove(temp_fname.c_str());
            }
        }

        /* write velocity */
        for (unsigned coord = 0; coord < 3; ++coord)
        {
            for (unsigned ilevel = levelmin_; ilevel <= levelmax_; ++ilevel)
            {
                std::string temp_fname = "___ic_temp_" + std::to_string(1000 * id_dm_vel + 100 * coord + ilevel) + ".bin";
                copy2outputfile(ofs_par, temp_fname);
                remove(temp_fname.c_str());
            }
        }

        ofs_par.flush();
        ofs_par.close();

        /* ----- write mesh file ----- */
        std::ofstream ofs_um(um_ic_.c_str(), std::ios::binary | std::ios::trunc);
        music::ilog.Print("GAMER-2 plugin: writing %s", um_ic_.c_str());

        // write out mesh data compatible with GAMER format UM_IC_FORMAT_VZYX
        // [field][z][y][x] in a row-major order, level by level from the coarsest to the finest

        size_t nmesh = 0;
        for (unsigned ilevel = levelmin_; ilevel <= levelmax_; ++ilevel)
        {
            /* output density */
            std::string temp_fname = "___ic_temp_" + std::to_string(1000 * id_gas_rho + ilevel) + ".bin";
            size_t nmesh_level     = copy2outputfile(ofs_um, temp_fname);
            remove(temp_fname.c_str());

            music::ilog.Print("Level %2d: %12llu cells", ilevel, nmesh_level);
            nmesh += nmesh_level;

            /* output velocity */
            for (unsigned coord = 0; coord < 3; ++coord)
            {
                temp_fname = "___ic_temp_" + std::to_string(1000 * id_gas_vel + 100 * coord + ilevel) + ".bin";
                copy2outputfile(ofs_um, temp_fname);
                remove(temp_fname.c_str());
            }
        }

        music::ilog.Print("Total: %12llu cells", nmesh);

        ofs_um.flush();
        ofs_um.close();

        music::ilog.Print("GAMER-2 plugin: writing recommended parameters for this initial condition in Input__Parameter_MUSIC_IC");
        std::string fname = "Input__Parameter_MUSIC_IC";
        std::ofstream ofs(fname, std::ios::trunc);
        ofs << "# Recommended parameters\n";
        ofs << "\n# simulation scale\n";
        ofs << "BOX_SIZE                " << std::setw(12) << boxlength_ << "   # box size along the longest side (in Mpc/h)\n";
        ofs << "NX0_TOT_X               " << std::setw(12) << (1 << levelmin_) << "   # number of base-level cells along x\n";
        ofs << "NX0_TOT_Y               " << std::setw(12) << (1 << levelmin_) << "   # number of base-level cells along y\n";
        ofs << "NX0_TOT_Z               " << std::setw(12) << (1 << levelmin_) << "   # number of base-level cells along z\n";
        ofs << "\n# cosmology (COMOVING only)\n";
        ofs << "A_INIT                  " << std::setw(12) << (1.0 / (zstart_ + 1)) << "   # initial scale factor\n";
        ofs << "OMEGA_M0                " << std::setw(12) << omegam_ << "   # omega matter at the present time\n";
        ofs << "HUBBLE0                 " << std::setw(12) << hubble_
            << "   # dimensionless Hubble parameter (currently only for converting ELBDM_MASS to code units)\n";
        ofs << "\n# particle\n";
        ofs << "PAR_NPAR                " << std::setw(12) << npart
            << "   # total number of particles (must be set for PAR_INIT==1/3; must be an integer)\n";
        ofs << "PAR_INIT                " << std::setw(12) << 3
            << "   # initialization option for particles: (1=FUNCTION, 2=RESTART, 3=FILE->'PAR_IC')\n";
        ofs << "PAR_IC_FORMAT           " << std::setw(12) << 1
            << "   # data format of PAR_IC: (1=[attribute][id], 2=[id][attribute]; row-major) [1]\n";
        ofs << "PAR_IC_FLOAT8           " << std::setw(12) << ((sizeof(real_t) == 8) ? 1 : 0)
            << "   # floating-point precision for PAR_IC (<0: default, 0: single, 1: double) [default: same as FLOAT8_PAR]\n";
        ofs << "PAR_IC_MASS             " << std::setw(12) << -1.0 << "   # mass of all particles for PAR_INIT==3 (<0=off) [-1.0]\n";
        ofs << "PAR_IC_TYPE             " << std::setw(12) << 2
            << "   # type of all particles for PAR_INIT==3 (<0=off, 2=dark matter) [-1]\n";
        ofs << "\n# initialization\n";
        ofs << "OPT__INIT               " << std::setw(12) << 3 << "   # initialization option: (1=FUNCTION, 2=RESTART, 3=FILE->'UM_IC')\n";
        ofs << "OPT__UM_IC_LEVEL        " << std::setw(12) << 0 << "   # AMR level corresponding to UM_IC (must >= 0) [0]\n";
        ofs << "OPT__UM_IC_NLEVEL       " << std::setw(12) << (levelmax_ - levelmin_ + 1)
            << "   # number of AMR refinement levels in UM_IC \n";
        ofs << "OPT__UM_IC_NVAR         " << std::setw(12) << 4 << "   # number of variables in UM_IC: (0:Dens, 1:VelX, 2:VelY, 3:VelZ)\n";
        ofs << "OPT__UM_IC_FORMAT       " << std::setw(12) << 1
            << "   # data format of UM_IC: (1=vzyx, 2=zyxv; row-major and v=field) [1]\n";
        ofs << "OPT__UM_IC_FLOAT8       " << std::setw(12) << ((sizeof(real_t) == 8) ? 1 : 0)
            << "   # floating-point precision for PAR_IC (<0: default, 0: single, 1: double) [default: same as FLOAT8_PAR]\n";
    }

  public:
    gamer2_output_plugin(config_file &cf) // std::string fname, Cosmology cosm, Parameters param )
        : output_plugin(cf)               // fname, cosm, param )
    {
        par_ic_ = cf.get_value<std::string>("output", "parfilename");
        um_ic_  = cf.get_value<std::string>("output", "filename");

        boxlength_ = cf.get_value<double>("setup", "boxlength");
        zstart_    = cf.get_value<double>("setup", "zstart");
        bbaryons_  = cf.get_value<bool>("setup", "baryons");
        omegam_    = cf.get_value<double>("cosmology", "Omega_m");
        omegab_    = cf.get_value<double>("cosmology", "Omega_b");
        hubble_    = cf.get_value<double>("cosmology", "H0") / 100.0;

        // necessary to cover refinement region by patches (blocks)
        int blocking_factor = cf.get_value_safe<int>("setup", "blocking_factor", -1);
        if (blocking_factor == -1)
        {
            music::elog.Print("GAMER-2 plugin: require [setup]/blocking_factor to set the patch size (default = 8)");
            throw std::runtime_error("GAMER-2 plugin: require blocking_factor");
        }
        // necessary for the mesh resolution of neighbouring patchs to be within factor of two
        int padding = cf.get_value_safe<int>("setup", "padding", 0);
        if (padding <= 0)
        {
            music::elog.Print("GAMER-2 plugin: require [setup]/padding > 0 for the mesh resolution of neighbouring patchs to be within factor of two");
            throw std::runtime_error("GAMER-2 plugin: require padding > 0");
        }

        PATCH_SIZE = blocking_factor;
        PS2        = 2 * PATCH_SIZE;

        // "GAMER" physical constants (/include/PhysicalConstant.h), input parameters (Input__Parameter), and COMOVING units
        // (src/Init/Init_Unit.cpp)
        double const_cm   = 1.0;
        double const_km   = 1.0e5 * const_cm;
        double const_pc   = 3.08567758149e18 * const_cm; // parsec
        double const_Mpc  = 1.0e6 * const_pc;
        double const_s    = 1.0;                                                // second
        double const_Msun = 1.9885e33;                                          // solar mass
        double const_G    = 6.6738e-8;                                          // gravitational constant in cgs
        double H0         = 100.0 * hubble_ * const_km / (const_s * const_Mpc); // H0 = 100*h*km/(s*Mpc)
        // see https://github.com/gamer-project/gamer/wiki/Runtime-Parameters%3A-Units#units-in-cosmological-simulations
        gamer_unit_length_   = const_Mpc / hubble_;
        gamer_unit_time_     = 1.0 / H0;
        gamer_unit_density_  = 3.0 * omegam_ * H0 * H0 / (8.0 * M_PI * const_G);
        gamer_unit_mass_     = gamer_unit_density_ * gamer_unit_length_ * gamer_unit_length_ * gamer_unit_length_;
        gamer_unit_velocity_ = gamer_unit_length_ / gamer_unit_time_;

        // MUSIC units
        music_unit_length_   = const_Mpc / hubble_;
        music_unit_mass_     = const_Msun / hubble_;
        music_unit_velocity_ = 1.0e5 * boxlength_;
    }

    ~gamer2_output_plugin()
    {
    }

    void write_dm_density(const grid_hierarchy &gh)
    {
        // write refinement mask. the AMR grid should satisfy the following requirements.
        // REQUIREMENT: 1. the size of refinement region must be multiple of patch size
        //              2. the refinement region must be aligned with the patch in parent level
        //              3. the mesh resolution of neighbouring patchs must be within factor of two
        // the MUSIC hierarchical grid should satisfy the above requirements (when blocking_factor and padding are set), but here we
        // check it again
        refine_offset = {0, 0, 0, 0, 0, 0};
        refine_size   = {(1 << levelmin_), (1 << levelmin_), (1 << levelmin_)};

        for (unsigned ilevel = levelmin_ + 1; ilevel <= levelmax_; ++ilevel)
        {
            unsigned dlv = ilevel - levelmin_;
            for (int dim = 0; dim < 3; dim++)
            {
                // left edge
                int left = updiv(gh.offset_abs(ilevel, dim), PS2);

                // right edge
                int right = updiv((1 << ilevel) - gh.offset_abs(ilevel, dim) - gh.size(ilevel, dim), PS2);

                // mesh resolution of neighbouring patchs must be within factor of two
                if (left != 0 || right != 0)
                {
                    int left_coarse = refine_offset[6 * (dlv - 1) + 2 * dim];
                    if (left <= 2 * left_coarse)
                        left = 2 * left_coarse + 1;

                    int right_course = refine_offset[6 * (dlv - 1) + 2 * dim + 1];
                    if (right <= 2 * right_course)
                        right = 2 * right_course + 1;
                }

                int size = (1 << ilevel) - (left + right) * PS2;

                if (size <= 0)
                {
                    music::elog.Print("zero size refinement region on level %d", ilevel);
                    throw std::runtime_error("Fatal: zero size refinement region");
                }

                int left_music  = gh.offset_abs(ilevel, dim);
                int right_music = (1 << ilevel) - gh.offset_abs(ilevel, dim) - gh.size(ilevel, dim);
                if (left_music != left * PS2 || right_music != right * PS2)
                {
                    music::elog.Print("MUSIC hierarchical grid does not satisfy GAMER-2 AMR requirement on level %d", ilevel);
                    throw std::runtime_error("Fatal: MUSIC hierarchical grid does not satisfy GAMER-2 AMR requirement");
                }

                refine_offset.push_back(left);
                refine_offset.push_back(right);

                refine_size.push_back(size);
            }
        }

        /* write refinement mask */
        music::ilog.Print("GAMER-2 plugin: writing refinement mask Input__UM_IC_RefineRegion");
        std::string fname = "Input__UM_IC_RefineRegion";
        std::ofstream ofs(fname, std::ios::trunc);
        ofs << "# RefineRegion\n";
        ofs << "# Offset describes the shift of a refinement grid with respect to its coarser parent grid in number of patches\n";
        ofs << "# lv, offset (xl, xr, yl, yr, zl, zr)\n";
        for (unsigned ilevel = levelmin_ + 1; ilevel <= levelmax_; ++ilevel)
        {
            unsigned dlv = ilevel - levelmin_;
            ofs << std::setw(3) << ilevel << " \t";
            ofs << std::setw(5) << refine_offset[6 * dlv + 0] - 2 * refine_offset[6 * (dlv - 1) + 0] << " \t";
            ofs << std::setw(5) << refine_offset[6 * dlv + 1] - 2 * refine_offset[6 * (dlv - 1) + 1] << " \t";
            ofs << std::setw(5) << refine_offset[6 * dlv + 2] - 2 * refine_offset[6 * (dlv - 1) + 2] << " \t";
            ofs << std::setw(5) << refine_offset[6 * dlv + 3] - 2 * refine_offset[6 * (dlv - 1) + 3] << " \t";
            ofs << std::setw(5) << refine_offset[6 * dlv + 4] - 2 * refine_offset[6 * (dlv - 1) + 4] << " \t";
            ofs << std::setw(5) << refine_offset[6 * dlv + 5] - 2 * refine_offset[6 * (dlv - 1) + 5] << " \n";
        }

        ofs.flush();
        ofs.close();

        if (!bbaryons_)
        {
            // music::ilog.Print("GAMER-2 plugin: writing vacuum mesh to temporary files");
            // write vacuum mesh (floor will be applied in GAMER)
            for (unsigned ilevel = levelmin_; ilevel <= levelmax_; ++ilevel)
            {
                std::string temp_fname = "___ic_temp_" + std::to_string(1000 * id_gas_rho + ilevel) + ".bin";
                write2tempfile_grid(temp_fname, gh, ilevel, 0, 0);

                for (int coord = 0; coord < 3; ++coord)
                {
                    temp_fname = "___ic_temp_" + std::to_string(1000 * id_gas_vel + 100 * coord + ilevel) + ".bin";
                    write2tempfile_grid(temp_fname, gh, ilevel, 0, 0);
                }
            }
        }
    }

    void write_dm_mass(const grid_hierarchy &gh)
    {
        // music::ilog.Print("GAMER-2 plugin: writing particle mass to temporary files");

        // write particle mass
        for (unsigned ilevel = levelmin_; ilevel <= levelmax_; ++ilevel)
        {
            /* uniform particle mass for each level */
            double dx   = boxlength_ / (double)(1 << ilevel);
            double dx3  = dx * dx * dx;
            double rhom = 2.77519737e11; // h-1 M_o / (h-1 Mpc)**3
            real_t cmass;

            if (bbaryons_)
                cmass = (omegam_ - omegab_) * rhom * dx3 * music_unit_mass_ / gamer_unit_mass_;
            else
                cmass = omegam_ * rhom * dx3 * music_unit_mass_ / gamer_unit_mass_;

            std::string temp_fname = "___ic_temp_" + std::to_string(1000 * id_dm_mass + ilevel) + ".bin";

            write2tempfile_par(temp_fname, gh, ilevel, 0, cmass);
        }
    }

    void write_dm_position(int coord, const grid_hierarchy &gh)
    {
        // music::ilog.Print("GAMER-2 plugin: writing particle position to temporary files");

        for (unsigned ilevel = levelmin_; ilevel <= levelmax_; ++ilevel)
        {
            const MeshvarBnd<real_t> *data = gh.get_grid(ilevel);

            int n0 = data->size(0), n1 = data->size(1), n2 = data->size(2);
            std::vector<real_t> vdata;
            vdata.reserve((unsigned)(n0) * (n1) * (n2));

            size_t count = 0;
            real_t fac   = music_unit_length_ / gamer_unit_length_;
            for (int i = 0; i < n0; ++i)
                for (int j = 0; j < n1; ++j)
                    for (int k = 0; k < n2; ++k)
                        if (gh.is_in_mask(ilevel, i, j, k) && !gh.is_refined(ilevel, i, j, k))
                        {
                            /* convert displacement to particle position */
                            double xx[3];
                            gh.cell_pos(ilevel, i, j, k, xx);
                            real_t pos = (xx[coord] + (*data)(i, j, k)) * boxlength_;

                            /* apply periodic boundary condition */
                            while (pos < 0)
                                pos += boxlength_;
                            while (pos >= boxlength_)
                                pos -= boxlength_;

                            vdata.push_back(pos * fac);

                            count++;
                        }

            std::string temp_fname = "___ic_temp_" + std::to_string(1000 * id_dm_pos + 100 * coord + ilevel) + ".bin";

            std::ofstream ofs_temp(temp_fname, std::ios::binary | std::ios::trunc);
            ofs_temp.write((char *)&vdata[0], vdata.size() * sizeof(real_t));

            ofs_temp.flush();
            ofs_temp.close();
        }
    }

    void write_dm_velocity(int coord, const grid_hierarchy &gh)
    {
        // music::ilog.Print("GAMER-2 plugin: writing particle velocity to temporary files");

        for (unsigned ilevel = levelmin_; ilevel <= levelmax_; ++ilevel)
        {
            std::string temp_fname = "___ic_temp_" + std::to_string(1000 * id_dm_vel + 100 * coord + ilevel) + ".bin";

            // Code velocity is v_code = a * v_peculiar. See Eq. (16) of SCHIVE, TSAI, & CHIUEH (2010)
            real_t a_start = 1.0 / (1.0 + zstart_);
            real_t fac     = a_start * music_unit_velocity_ / gamer_unit_velocity_;
            write2tempfile_par(temp_fname, gh, ilevel, fac, 0);
        }
    }

    void write_dm_potential(const grid_hierarchy &gh)
    { /* skip */
    }

    void write_gas_potential(const grid_hierarchy &gh)
    { /* skip */
    }

    void write_gas_position(int coord, const grid_hierarchy &gh)
    { /* not used */
    }

    void write_gas_density(const grid_hierarchy &gh)
    {
        // music::ilog.Print("GAMER-2 plugin: writing gas density to temporary files");

        for (unsigned ilevel = levelmin_; ilevel <= levelmax_; ++ilevel)
        {
            std::string temp_fname = "___ic_temp_" + std::to_string(1000 * id_gas_rho + ilevel) + ".bin";

            real_t fac = omegab_ / omegam_;
            write2tempfile_grid(temp_fname, gh, ilevel, fac, fac);
        }
    }

    void write_gas_velocity(int coord, const grid_hierarchy &gh)
    {
        // music::ilog.Print("GAMER-2 plugin: writing gas velocity to temporary files");

        for (unsigned ilevel = levelmin_; ilevel <= levelmax_; ++ilevel)
        {
            std::string temp_fname = "___ic_temp_" + std::to_string(1000 * id_gas_vel + 100 * coord + ilevel) + ".bin";

            // Code velocity is v_code = a * v_peculiar. See Eq. (16) of SCHIVE, TSAI, & CHIUEH (2010)
            real_t a_start = 1.0 / (1.0 + zstart_);
            real_t fac     = a_start * music_unit_velocity_ / gamer_unit_velocity_;
            write2tempfile_grid(temp_fname, gh, ilevel, fac, 0);
        }
    }

    void finalize(void)
    {
        assemble_gamer2_file();
    }
};

namespace
{
output_plugin_creator_concrete<gamer2_output_plugin> creator("gamer2");
}
