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

#ifndef __DENSITIES_HH
#define __DENSITIES_HH

#include <assert.h>

#include "general.hh"
#include "config_file.hh"
#include "random.hh"
#include "transfer_function.hh"
#include "general.hh"

void GenerateDensityHierarchy(config_file &cf, const cosmology::calculator* cc, tf_type type,
															refinement_hierarchy &refh, noise_generator &rand, grid_hierarchy &delta, bool smooth, bool shift);

void GenerateDensityUnigrid(config_file &cf, const cosmology::calculator*, tf_type type,
														refinement_hierarchy &refh, noise_generator &rand, grid_hierarchy &delta, bool smooth, bool shift);

void normalize_density(grid_hierarchy &delta);

void coarsen_density(const refinement_hierarchy &rh, GridHierarchy<real_t> &u, bool kspace);

#endif
