/*
 
 defaults.hh - This file is part of MUSIC -
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


#ifndef __DEFAULTS_HH
#define __DEFAULTS_HH

#include <string>
#include <map>

struct default_conf{
	std::string sec;
	std::string tag;
	std::string val;
	default_conf( std::string sec_, std::string tag_, std::string val_ )
	: sec(sec_), tag(tag_), val(val_)
	{ }
};


class default_options{
protected:
	std::map<std::string,default_conf> def;
public:
	default_options();
	
	template<typename T> 
	void query( std::string tag )
	{}
	
};

extern default_options defaults;


#endif //__DEFAULTS_HH

