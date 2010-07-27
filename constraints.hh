/*
 
 constraints.hh - This file is part of MUSIC -
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

#ifndef __CONSTRAINTS_HH
#define __CONSTRAINTS_HH

#include <vector>
#include "config_file.hh"

class constraint_set
{
	
protected:
	
	struct constraint{
		double x,y,z;
		unsigned level;
		double sigma;
	};
	
	
	std::vector<constraint> cset_;
	
public:
	
	constraint_set( config_file& cf )
	{
		float tf0, tf1,tf2,tf3;
		unsigned ti;
		unsigned i=0;
		
		while(true)
		{
			char temp1[128];
			std::string temp2;
			sprintf(temp1,"constraint[%d]",i);
			if( cf.containsKey( "constraints", temp1 ) )
			{
				temp2				= cf.getValue<std::string>( "constraints", temp1 );
				sscanf( temp2.c_str(), "%f,%f,%f,%d,%f", &tf0, &tf1, &tf2, &ti, &tf3 ); 
				constraint new_c;
				new_c.x = tf0;
				new_c.y = tf1;
				new_c.z = tf2;
				new_c.level = ti;
				new_c.sigma = tf3;
				
				cset_.push_back( new_c );
				
			}		
			else
				break;
			++i;
		}
		
		
		std::cout << " - Found " << cset_.size() << " density constraint(s) to be obeyed.\n";
		
	}
	
};


#endif // __CONSTRAINTS_HH
