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
#include "mesh.hh"

//! abstract base class for Poisson solvers and gradient calculations
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
	
	//! solve Poisson's equation Du=f
	virtual double solve( grid_hierarchy& f, grid_hierarchy& u ) = 0;
	
	//! compute the gradient of u
	virtual double gradient( int dir, grid_hierarchy& u, grid_hierarchy& Du ) = 0;
	
	//! compute the gradient and add
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

//! adaptive FAS multigrid implementation of abstract base class poisson_plugin
class multigrid_poisson_plugin : public poisson_plugin
{
public:
	
	//! constructor
	explicit multigrid_poisson_plugin( config_file& cf )
	: poisson_plugin( cf )
	{ }
	
	//! solve Poisson's equation Du=f
	double solve( grid_hierarchy& f, grid_hierarchy& u );
	
	//! compute the gradient of u
	double gradient( int dir, grid_hierarchy& u, grid_hierarchy& Du );
	
	//! compute the gradient and add
	double gradient_add( int dir, grid_hierarchy& u, grid_hierarchy& Du );
	
protected:
	
	//! various FD approximation implementations
	struct implementation
	{
		//! solve 2nd order FD approximation to Poisson's equation
		double solve_O2( grid_hierarchy& f, grid_hierarchy& u );
		
		//! solve 4th order FD approximation to Poisson's equation
		double solve_O4( grid_hierarchy& f, grid_hierarchy& u );
		
		//! solve 6th order FD approximation to Poisson's equation		
		double solve_O6( grid_hierarchy& f, grid_hierarchy& u );
		
		//! compute 2nd order FD gradient
		void gradient_O2( int dir, grid_hierarchy& u, grid_hierarchy& Du );

		//! compute and add 2nd order FD gradient
		void gradient_add_O2( int dir, grid_hierarchy& u, grid_hierarchy& Du );
		
		//! compute 4th order FD gradient
		void gradient_O4( int dir, grid_hierarchy& u, grid_hierarchy& Du );
		
		//! compute and add 4th order FD gradient
		void gradient_add_O4( int dir, grid_hierarchy& u, grid_hierarchy& Du );
		
		//! compute 6th order FD gradient
		void gradient_O6( int dir, grid_hierarchy& u, grid_hierarchy& Du );
		
		//! compute and add 6th order FD gradient
		void gradient_add_O6( int dir, grid_hierarchy& u, grid_hierarchy& Du );
	};
};

/**************************************************************************************/
/**************************************************************************************/
#pragma mark -

//! FFT based implementation of abstract base class poisson_plugin
class fft_poisson_plugin : public poisson_plugin
{
public:
	
	//! constructor
	explicit fft_poisson_plugin( config_file& cf )
	: poisson_plugin( cf )
	{ }
	
	//! solve Poisson's equation Du=f
	double solve( grid_hierarchy& f, grid_hierarchy& u );
	
	//! compute the gradient of u
	double gradient( int dir, grid_hierarchy& u, grid_hierarchy& Du );
	
	//! compute the gradient and add
	double gradient_add( int dir, grid_hierarchy& u, grid_hierarchy& Du ){ return 0.0; }
	
	
};

/**************************************************************************************/
/**************************************************************************************/
#pragma mark -

void poisson_hybrid( MeshvarBnd<double>& f, int idir, int order, bool periodic );









#endif // __POISSON_HH

