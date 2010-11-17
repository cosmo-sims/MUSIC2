/*
 
 defaults.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
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

