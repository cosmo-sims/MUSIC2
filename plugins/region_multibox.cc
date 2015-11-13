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


typedef std::vector<int> col;
typedef std::vector<col> slice;
typedef std::vector<slice> grid;

typedef struct {
    int x;
    int y;
    int z;
} coord;

typedef std::vector<coord> region;

class region_multibox_plugin : public region_generator_plugin{
private:
    int res, reflevel;
    grid refgrid;
    region find_equal(grid ingrid, int level)
    {
        slice curslice;
        col curcol;
        region equal_region;
        coord cp;
        for(int i=0;i<ingrid.size();i++)
        {
            cp.x = i;
            curslice = ingrid[i];
            for(int j=0;j<curslice.size();j++)
            {
                cp.y = j;
                curcol = curslice[j];
                for(int k=0;k<curcol.size();k++)
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
    //This function takes the grid, which has been already had particles
    //deposited onto it, and set to the maximum refinement level.  It then
    //fills the remaining refinement levels
    void build_refgrid(grid ingrid, int maxlevel)
    {
        region curregion;
        for(int curlevel=maxlevel; curlevel>2; curlevel--)
        {
            curregion = find_equal(ingrid, curlevel);
            for(region::iterator cp= curregion.begin(); cp != curregion.end(); ++cp)
            {
                for(int i=-1; i<2; i++)
                {
                    for(int j=-1; j<2; j++)
                    {
                        for(int k=-1; k<2; k++)
                        {
                            if(ingrid[cp->x+i][cp->y+j][cp->z+k] == 0)
                            {
                                ingrid[cp->x+i][cp->y+j][cp->z+k] = curlevel-1;
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
        res = 128;
        //Initialize the grid with zeros, the base level
        refgrid = grid(res,slice(res,col(res,0)));
    }
    
    
    void get_AABB( double *left, double *right, unsigned level )
    {
    }
  
    void update_AABB( double *left, double *right )
    {
    }
  
    bool query_point( double *x, int level )
    {
    }
    
    bool is_grid_dim_forced( size_t* ndims )
    {   
    }
    
    void get_center( double *xc )
    {
    }

  void get_center_unshifted( double *xc )
  {
  }

};

namespace{
    region_generator_plugin_creator_concrete< region_multibox_plugin > creator("multibox");
}
