#include <vector>
/*
 
 region_multibox.cc - A plugin for MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2015 Ben Keller
 
 */

#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>

#include "region_generator.hh"
#include "point_file_reader.hh"


typedef std::vector<unsigned> col;
typedef std::vector<col> slice;
typedef std::vector<slice> grid;

typedef struct {
    unsigned x;
    unsigned y;
    unsigned z;
} coord;

typedef std::vector<coord> region;

class region_multibox_plugin : public region_generator_plugin{
private:
    double cen[3];
    unsigned res;
    double vfac_;
    grid refgrid;
    
    region where(unsigned level)
    {
        slice curslice;
        col curcol;
        region equal_region;
        coord cp;
        for(unsigned i=0;i<refgrid.size();i++)
        {
            cp.x = i;
            curslice = refgrid[i];
            for(unsigned j=0;j<curslice.size();j++)
            {
                cp.y = j;
                curcol = curslice[j];
                for(unsigned k=0;k<curcol.size();k++)
                {
                    cp.z = k;
                    if(curcol[k] == level)
                    {
                        equal_region.push_back(cp);
                    }
                }
            }
        }
        return equal_region;
    }
    //This function takes the grid, a vector of doubles containing the particle
    //coordinates, and sets each particle-containing gridcell to levelmax_.
    void deposit_particles(std::vector<double> pp, unsigned res)
    {
        unsigned i,j,k;
        std::vector<double>::iterator cp = pp.begin();
        while(cp != pp.end())
        {
            i = (*(cp))*res;
            cp++;
            j = (*(cp))*res;
            cp++;
            k = (*(cp))*res;
            cp++;
            refgrid[i][j][k] = levelmax_;
        }
    }

    //This function takes the grid, which has been already had particles
    //deposited onto it, and set to the maximum refinement level.  It then
    //fills the remaining refinement levels
    void build_refgrid()
    {
        region curregion;
        //Build an extra layer about the maxlevel layer.
        curregion = where(levelmax_);
        for(region::iterator cp= curregion.begin(); cp != curregion.end(); ++cp)
        {
            for(int i=-1; i<2; i++)
            {
                for(int j=-1; j<2; j++)
                {
                    for(int k=-1; k<2; k++)
                    {
                        if(refgrid[cp->x+i][cp->y+j][cp->z+k] == levelmin_)
                        {
                            refgrid[cp->x+i][cp->y+j][cp->z+k] = levelmax_;
                        }
                    }
                }
            }
        }
        for(unsigned curlevel=levelmax_; curlevel>(levelmin_+1); curlevel--)
        {
            curregion = where(curlevel);
            for(region::iterator cp= curregion.begin(); cp != curregion.end(); ++cp)
            {
                for(int i=-1; i<2; i++)
                {
                    for(int j=-1; j<2; j++)
                    {
                        for(int k=-1; k<2; k++)
                        {
                            if(refgrid[cp->x+i][cp->y+j][cp->z+k] == levelmin_)
                            {
                                refgrid[cp->x+i][cp->y+j][cp->z+k] = curlevel-1;
                            }
                        }
                    }
                }
            }
        }
    }
    
public:
    explicit region_multibox_plugin( config_file& cf )
    : region_generator_plugin( cf )
    {
        res = 1<<(levelmin_-1);
        //check parameters
        if ( !cf.containsKey("setup", "region_point_file"))
        {
            LOGERR("Missing parameter \'region_point_file\' needed for region=multibox");
            throw std::runtime_error("Missing parameter for region=multibox");
        }
        if ( cf.containsKey("setup", "region_point_level"))
        {
        res = 1<<(levelmin_-1);
        res = 1<<(cf.getValue<int>("setup","region_point_level")-1);
        }
        vfac_ = cf.getValue<double>("cosmology","vfact");
        if (levelmin_ >= levelmax_)
        {
            LOGERR("Why are you specifying a region if you aren't refining?");
            throw std::runtime_error("region=multibox needs levelmax>levelmin");
        }
        std::vector<double> pp;
        point_reader pfr;
        pfr.read_points_from_file(cf.getValue<std::string>("setup", "region_point_file"), vfac_, pp);
        LOGINFO("Multibox Resolution: %d", res);
        //Initialize the grid with zeros, the base level
        refgrid = grid(res,slice(res,col(res,levelmin_)));
        //Build the highest refinement region
        deposit_particles(pp, res);
        LOGINFO("Deposited Multigrid Particles");
        //Build the intermediate refinement regions
        build_refgrid();
        LOGINFO("Built Multigrid");
        get_center(cen);
        dump_grid();
    }
    
    
    void get_AABB( double *left, double *right, unsigned level )
    {
        left[0] = left[1] = left[2] = 1.0;
        right[0] = right[1] = right[2] = 0.0;
        if( level <= levelmin_ )
        {
            left[0] = left[1] = left[2] = 0.0;
            right[0] = right[1] = right[2] = 1.0;
            return;
        }
        region myres = where(level);
        std::vector<coord>::iterator cp = myres.begin();
        while (cp != myres.end())
        {
            //Check left side
            if (float(cp->x)/res < left[0])
                left[0] = float(cp->x)/res;
            if (float(cp->y)/res < left[1])
                left[1] = float(cp->y)/res;
            if (float(cp->z)/res < left[2])
                left[2] = float(cp->z)/res;
            //Check right side
            if (float(cp->x)/res > right[0])
                right[0] = float(cp->x)/res;
            if (float(cp->y)/res > right[1])
                right[1] = float(cp->y)/res;
            if (float(cp->z)/res > right[2])
                right[2] = float(cp->z)/res;
            cp++;
        }
    }
    void dump_grid()
    {
        FILE * gridfile = fopen("griddump.dat", "w");
        for(unsigned i=0;i<res;i++)
        {
            for(unsigned j=0; j<res; j++)
            {
                for(unsigned k=0;k<res;k++)
                {
                    fprintf(gridfile, "%03d %03d %03d %1d\n", i-int(res*cen[0]),j-int(res*cen[1]),k-int(res*cen[2]), refgrid[i][j][k]);
                }
            }
        }
        fclose(gridfile);
    }
  
    void update_AABB( double *left, double *right )
    {
        //no need for this I think?
    }
  
    bool query_point( double *x, int level )
    {
        for(int i=0; i<3; ++i)
        {
            if(x[i] >= 1.0 || x[i] <= 0.0) return false;
        }
        return (level <= int(refgrid[(x[0])*res][(x[1])*res][(x[2])*res]));
    }
    
    bool is_grid_dim_forced( size_t* ndims )
    {   
        return false; //is this true?
    }
    
    void get_center( double *xc )
    {
        double xc2[3];
        get_AABB(xc, xc2, levelmax_);
        for(int i=0; i<3; ++i)
        {
            xc[i] += xc2[i];
            xc[i] *= 0.5;
        }
    }

  void get_center_unshifted( double *xc )
  {
      get_center(xc);
  }

};

namespace{
    region_generator_plugin_creator_concrete< region_multibox_plugin > creator("multibox");
}
