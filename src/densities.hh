/*

 densities.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions
 for cosmological simulations

 Copyright (C) 2010  Oliver Hahn

 */

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
