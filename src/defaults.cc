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

#include "defaults.hh"

#define ADD_DEF(x,a,b,c) def.insert(std::pair<std::string,default_conf>(x,default_conf(a,b,c)));

default_options defaults;

default_options::default_options()
{
	//... [setup] ...
	ADD_DEF("align_top",		"setup",	"align_top",		"yes");
	ADD_DEF("baryons",			"setup",	"baryons",			"no");
	ADD_DEF("center_v",			"setup",	"center_velocities","no");
	ADD_DEF("deconvolve",		"setup",	"deconvolve",		"yes");
	ADD_DEF("exact_shotnoise",	"setup",	"exact_shotnoise",	"yes");
	ADD_DEF("overlap",			"setup",	"overlap",			"8");
	ADD_DEF("padding",			"setup",	"padding",			"16");
	ADD_DEF("periodic_TF",		"setup",	"periodic_TF",		"yes");
	ADD_DEF("use_2LPT",			"setup",	"use_2LPT",			"yes");
	ADD_DEF("use_LLA",			"setup",	"use_LLA",			"no");

	//... [poisson] ...
	ADD_DEF("mgacc",			"poisson",	"accuracy",			"1e-4");
	ADD_DEF("mggrad",			"poisson",	"grad_order",		"6");
	ADD_DEF("mglapl",			"poisson",	"laplce_order",		"6");
	ADD_DEF("fft_fine",			"poisson",	"fft_fine",			"yes");
	ADD_DEF("kspace_poisson",	"poisson",	"kspace",			"no");
	
	
	//... deprecated
	ADD_DEF("avg_fine",			"setup",	"avg_fine",			"no");
}




#undef ADD_DEF







