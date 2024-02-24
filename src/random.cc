// This file is part of MUSIC
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2010-2024 by Oliver Hahn
// 
// monofonIC is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// monofonIC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
