/*

 random.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions for cosmological simulations

 Copyright (C) 2010-23  Oliver Hahn

 */

#include "random.hh"

std::map<std::string, RNG_plugin_creator *> &
get_RNG_plugin_map()
{
	static std::map<std::string, RNG_plugin_creator *> RNG_plugin_map;
	return RNG_plugin_map;
}

void print_RNG_plugins()
{
	std::map<std::string, RNG_plugin_creator *> &m = get_RNG_plugin_map();
	std::map<std::string, RNG_plugin_creator *>::iterator it;
	it = m.begin();
	std::cout << " - Available random number generator plug-ins:\n";
	while (it != m.end())
	{
		if ((*it).second)
			std::cout << "\t\'" << (*it).first << "\'\n";
		++it;
	}
}

RNG_plugin *select_RNG_plugin(config_file &cf)
{
	std::string rngname = cf.get_value_safe<std::string>("random", "generator", "MUSIC");

	RNG_plugin_creator *the_RNG_plugin_creator = get_RNG_plugin_map()[rngname];

	if (!the_RNG_plugin_creator)
	{
		std::cerr << " - Error: random number generator plug-in \'" << rngname << "\' not found." << std::endl;
		music::elog.Print("Invalid/Unregistered random number generator plug-in encountered : %s", rngname.c_str());
		print_RNG_plugins();
		throw std::runtime_error("Unknown random number generator plug-in");
	}
	else
	{
		std::cout << " - Selecting random number generator plug-in \'" << rngname << "\'..." << std::endl;
		music::ulog.Print("Selecting random number generator plug-in  : %s", rngname.c_str());
	}

	RNG_plugin *the_RNG_plugin = the_RNG_plugin_creator->create(cf);

	return the_RNG_plugin;
}
