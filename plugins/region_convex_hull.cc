/*
 
 region_convex_hull.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010-13  Oliver Hahn
 
 */

#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>

#include "region_generator.hh"
#include "convex_hull.hh"
#include "point_file_reader.hh"


//! Convex hull region plugin
class region_convex_hull_plugin : public region_generator_plugin{
private:
    
    convex_hull<double> *phull_;
    int shift[3], shift_level, padding_;
    double vfac_;
    bool do_extra_padding_;
    
    std::vector<float> level_dist_;
    
    void apply_shift( size_t Np, double *p, int *shift, int levelmin )
    {
        double dx = 1.0/(1<<levelmin);
        LOGINFO("unapplying shift of previous zoom region to region particles :\n" \
                "\t [%d,%d,%d] = (%f,%f,%f)",shift[0],shift[1],shift[2],shift[0]*dx,shift[1]*dx,shift[2]*dx);
        
        for( size_t i=0,i3=0; i<Np; i++,i3+=3 )
            for( size_t j=0; j<3; ++j )
                p[i3+j] = p[i3+j]-shift[j]*dx;
    }
    
public:
    explicit region_convex_hull_plugin( config_file& cf )
    : region_generator_plugin( cf )
    {
        std::vector<double> pp;
        
        vfac_ = cf.getValue<double>("cosmology","vfact");
        padding_ = cf.getValue<int>("setup","padding");
        
        
        std::string point_file = cf.getValue<std::string>("setup","region_point_file");
        
        point_reader pfr;
        pfr.read_points_from_file( point_file, vfac_, pp );
        
        if( cf.containsKey("setup","region_point_shift") )
        {
            std::string point_shift = cf.getValue<std::string>("setup","region_point_shift");
            if(sscanf( point_shift.c_str(), "%d,%d,%d", &shift[0],&shift[1],&shift[2] )!=3){
	      LOGERR("Error parsing triple for region_point_shift");
	      throw std::runtime_error("Error parsing triple for region_point_shift");
	    }
            unsigned point_levelmin = cf.getValue<unsigned>("setup","region_point_levelmin");
            
            apply_shift( pp.size()/3, &pp[0], shift, point_levelmin );
            shift_level = point_levelmin;
        }
        
        // compute the convex hull
        phull_ =  new convex_hull<double>(  &pp[0], pp.size()/3 );
        
        //expand the ellipsoid by one grid cell
        unsigned levelmax = cf.getValue<unsigned>("setup","levelmax");
        double dx = 1.0/(1ul<<levelmax);
        phull_->expand( sqrt(3.)*dx );
        
        // output the center
        double c[3] = { phull_->centroid_[0], phull_->centroid_[1], phull_->centroid_[2] };
        LOGINFO("Region center from convex hull centroid determined at\n\t (%f,%f,%f)",c[0],c[1],c[2]);
        
        //-----------------------------------------------------------------
        // when querying the bounding box, do we need extra padding?
        do_extra_padding_ = false;
        
        // conditions should be added here
        {
            std::string output_plugin = cf.getValue<std::string>("output","format");
            if( output_plugin == std::string("grafic2") )
                do_extra_padding_ = true;
        }
        
        level_dist_.assign( levelmax+1, 0.0 );
        // generate the higher level ellipsoids
        for( int ilevel = levelmax_-1; ilevel > 0; --ilevel )
        {
            dx = 1.0/(1ul<<(ilevel));
            level_dist_[ilevel] = level_dist_[ilevel+1] + padding_ * dx;
        }
    }
    
    ~region_convex_hull_plugin()
    {
        delete phull_;
    }
    
    void get_AABB( double *left, double *right, unsigned level )
    {
        for( int i=0; i<3; ++i )
        {
            left[i] = phull_->left_[i];
            right[i] = phull_->right_[i];
        }
        double dx = 1.0/(1ul<<level);
        double pad = (double)(padding_+1) * dx;
        
        if( ! do_extra_padding_ ) pad = 0.0;
        
        double ext = sqrt(3)*dx + pad;
        
        for( int i=0;i<3;++i )
        {
            left[i]  -= ext;
            right[i] += ext;
        }
        
    }
    
    void update_AABB( double *left, double *right )
    {
        // we ignore this, the grid generator must have generated a grid that contains the ellipsoid
        // it might have enlarged it, but who cares...
    }

    bool query_point( double *x, int ilevel )
    {   return phull_->check_point( x, level_dist_[ilevel] );   }
    
    bool is_grid_dim_forced( size_t* ndims )
    {   return false;   }
    
    void get_center( double *xc )
    {
        xc[0] = phull_->centroid_[0];
        xc[1] = phull_->centroid_[1];
        xc[2] = phull_->centroid_[2];
    }
    
    void get_center_unshifted( double *xc )
    {
        double dx = 1.0/(1<<shift_level);
        double c[3] = { phull_->centroid_[0], phull_->centroid_[1], phull_->centroid_[2] };
        xc[0] = c[0]+shift[0]*dx;
        xc[1] = c[1]+shift[1]*dx;
        xc[2] = c[2]+shift[2]*dx;
        
    }
};

namespace{
    region_generator_plugin_creator_concrete< region_convex_hull_plugin > creator("convex_hull");
}

