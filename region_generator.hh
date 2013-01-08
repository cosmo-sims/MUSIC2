#ifndef __REGION_GENERATOR_HH
#define __REGION_GENERATOR_HH

#include <vector>
#include "config_file.hh"

//! Abstract base class for region generators
/*!
 This class implements a purely virtual interface that can be
 used to derive instances implementing various region generators.
 */
class region_generator_plugin{
public:
    config_file *pcf_;
public:
    region_generator_plugin( config_file& cf )
    : pcf_( &cf )
    {
    }
    
    //! destructor
    virtual ~region_generator_plugin() { };
    
    //! compute the bounding box of the region
    virtual void get_AABB( double *left, double *right, unsigned level) = 0;
    
    //! query whether a point intersects the region
    virtual bool query_point( double *x ) = 0;
    
    //! query whether the region generator explicitly forces the grid dimensions
    virtual bool is_grid_dim_forced( size_t *ndims ) = 0;
    
    //! get the center of the region
    virtual void get_center( double *xc ) = 0;
    
};

//! Implements abstract factory design pattern for region generator plug-ins
struct region_generator_plugin_creator
{
	//! create an instance of a transfer function plug-in
	virtual region_generator_plugin * create( config_file& cf ) const = 0;
	
	//! destroy an instance of a plug-in
	virtual ~region_generator_plugin_creator() { }
};

//! Write names of registered region generator plug-ins to stdout
std::map< std::string, region_generator_plugin_creator *>& get_region_generator_plugin_map();
void print_region_generator_plugins( void );

//! Concrete factory pattern for region generator plug-ins
template< class Derived >
struct region_generator_plugin_creator_concrete : public region_generator_plugin_creator
{
	//! register the plug-in by its name
	region_generator_plugin_creator_concrete( const std::string& plugin_name )
	{
		get_region_generator_plugin_map()[ plugin_name ] = this;
	}
	
	//! create an instance of the plug-in
	region_generator_plugin * create( config_file& cf ) const
	{
		return new Derived( cf );
	}
};

typedef region_generator_plugin region_generator;

region_generator_plugin *select_region_generator_plugin( config_file& cf );

extern region_generator_plugin *the_region_generator;

#endif
