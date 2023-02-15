/*
 
 cosmology.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#pragma once

#include "mesh.hh"
#include "general.hh"


//! compute the jeans sound speed
/*! given a density in g/cm^-3 and a mass in g it gives back the sound
 *  speed in cm/s for which the input mass is equal to the jeans mass
 *  @param rho density 
 *  @param mass mass scale
 *  @returns jeans sound speed
 */
inline double jeans_sound_speed( double rho, double mass )
{
	const double G = 6.67e-8;
	return pow( 6.0*mass/M_PI*sqrt(rho)*pow(G,1.5), 1.0/3.0 );
}

//! computes the density from the potential using the Laplacian
void compute_Lu_density( const grid_hierarchy& u, grid_hierarchy& fnew, unsigned order=4 );

//! computes the 2nd order density perturbations using also off-diagonal terms in the potential Hessian 
void compute_LLA_density( const grid_hierarchy& u, grid_hierarchy& fnew, unsigned order=4 );

//! computes the source term for the 2nd order perturbations in the displacements
void compute_2LPT_source( const grid_hierarchy& u, grid_hierarchy& fnew, unsigned order=4 );

void compute_2LPT_source_FFT( config_file& cf_, const grid_hierarchy& u, grid_hierarchy& fnew );

