/*
 
 output.cc - This file is part of MUSIC -
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

#include "output.hh"


std::map< std::string, output_plugin_creator *>& 
get_output_plugin_map()
{
	static std::map< std::string, output_plugin_creator* > output_plugin_map;
	return output_plugin_map;
}

void print_output_plugins()
{
	std::map< std::string, output_plugin_creator *>& m = get_output_plugin_map();
	
	std::map< std::string, output_plugin_creator *>::iterator it;
	it = m.begin();
	std::cout << " - Available output plug-ins:\n";
	while( it!=m.end() )
	{
		if( (*it).second )
			std::cout << "\t\'" << (*it).first << "\'\n";
		++it;
	}
		
}

output_plugin *select_output_plugin( config_file& cf )
{
	std::string formatname = cf.getValue<std::string>( "output", "format" );
	
	output_plugin_creator *the_output_plugin_creator 
	= get_output_plugin_map()[ formatname ];
	
	if( !the_output_plugin_creator )
	{	
		std::cerr << " - Error: output plug-in \'" << formatname << "\' not found." << std::endl;
		print_output_plugins();
		throw std::runtime_error("Unknown output plug-in");
		
	}else
		std::cout << " - Selecting output plug-in \'" << formatname << "\'..." << std::endl;
	
	output_plugin *the_output_plugin 
	= the_output_plugin_creator->create( cf );
	
	return the_output_plugin;
}



