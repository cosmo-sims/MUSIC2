/*
 
 poisson.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 */

#ifndef __POISSON_HH
#define __POISSON_HH

#include <string>
#include <map>

#include "general.hh"

class poisson_plugin
{
protected:
	
	//! reference to the config_file object that holds all configuration options
	config_file& cf_;
	
public:
	
	//! constructor
	explicit poisson_plugin( config_file& cf )
	: cf_(cf)
	{ }
	
	//! destructor
	virtual ~poisson_plugin()
	{ }
	
	virtual double solve( grid_hierarchy& f, grid_hierarchy& u ) = 0;
	virtual double gradient( int dir, grid_hierarchy& u, grid_hierarchy& Du ) = 0;
	virtual double gradient_add( int dir, grid_hierarchy& u, grid_hierarchy& Du ) = 0;
	
};

#pragma mark -

/*!
 * @brief implements abstract factory design pattern for poisson solver plug-ins
 */
struct poisson_plugin_creator
{
	//! create an instance of a plug-in
	virtual poisson_plugin * create( config_file& cf ) const = 0;
	
	//! destroy an instance of a plug-in
	virtual ~poisson_plugin_creator() { }
};

//! maps the name of a plug-in to a pointer of the factory pattern 
std::map< std::string, poisson_plugin_creator *>& get_poisson_plugin_map();

//! print a list of all registered output plug-ins
void print_poisson_plugins();


/*!
 * @brief concrete factory pattern for output plug-ins
 */
template< class Derived >
struct poisson_plugin_creator_concrete : public poisson_plugin_creator
{
	//! register the plug-in by its name
	poisson_plugin_creator_concrete( const std::string& plugin_name )
	{
		get_poisson_plugin_map()[ plugin_name ] = this;
	}
	
	//! create an instance of the plug-in
	poisson_plugin * create( config_file& cf ) const
	{
		return new Derived( cf );
	}
};

/**************************************************************************************/
/**************************************************************************************/
#pragma mark -

class multigrid_poisson_plugin : public poisson_plugin
{
public:
	
	explicit multigrid_poisson_plugin( config_file& cf )
	: poisson_plugin( cf )
	{ }
	
	double solve( grid_hierarchy& f, grid_hierarchy& u );
	double gradient( int dir, grid_hierarchy& u, grid_hierarchy& Du );
	double gradient_add( int dir, grid_hierarchy& u, grid_hierarchy& Du );
	
protected:
	
	struct implementation
	{
		double solve_O2( grid_hierarchy& f, grid_hierarchy& u );
		double solve_O4( grid_hierarchy& f, grid_hierarchy& u );
		double solve_O6( grid_hierarchy& f, grid_hierarchy& u );
		void gradient_O2( int dir, grid_hierarchy& u, grid_hierarchy& Du );
		void gradient_add_O2( int dir, grid_hierarchy& u, grid_hierarchy& Du );
		void gradient_O4( int dir, grid_hierarchy& u, grid_hierarchy& Du );
		void gradient_add_O4( int dir, grid_hierarchy& u, grid_hierarchy& Du );
		void gradient_O6( int dir, grid_hierarchy& u, grid_hierarchy& Du );
		void gradient_add_O6( int dir, grid_hierarchy& u, grid_hierarchy& Du );
	};
};

/**************************************************************************************/
/**************************************************************************************/
#pragma mark -

class fft_poisson_plugin : public poisson_plugin
{
public:
	
	explicit fft_poisson_plugin( config_file& cf )
	: poisson_plugin( cf )
	{ }
	
	double solve( grid_hierarchy& f, grid_hierarchy& u );
	double gradient( int dir, grid_hierarchy& u, grid_hierarchy& Du );
	double gradient_add( int dir, grid_hierarchy& u, grid_hierarchy& Du ){ return 0.0; }
	
	
};












#endif // __POISSON_HH

