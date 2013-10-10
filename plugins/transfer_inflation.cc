/*
 
 tranfer_inflation.cc - This file is part of MUSIC -
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

#include "transfer_function.hh"

class transfer_inflation_plugin : public transfer_function_plugin
{
protected:
	
	double ns2_;
	
public:
	
	transfer_inflation_plugin( config_file& cf ) 
	: transfer_function_plugin( cf )
	{ 
		ns2_ = 0.5*cf.getValue<double>("cosmology","nspec");
		tf_distinct_ = true;
	}
	
	~transfer_inflation_plugin(){ };
	
	double compute( double k, tf_type type=baryon)
	{
		return pow(k,ns2_);
	}
	
	double get_kmax( void )
	{
		return 1e10;
	}
	
	double get_kmin( void )
	{
		return 1e-30;
	}
	
};


namespace{
	transfer_function_plugin_creator_concrete< transfer_inflation_plugin > creator("inflation");
}

