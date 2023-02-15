// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2020 by Oliver Hahn
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

#include <cosmology_parameters.hh>

/**
 * @brief namespace encapsulating all things cosmology
 * 
 */
namespace cosmology{

//! we store here the preset cosmological paramters
parameters::defaultmmap_t parameters::default_pmaps_
{
  //=============================================================================
  // Planck 2018 baseline cosmologies
  // cf. https://wiki.cosmos.esa.int/planck-legacy-archive/images/b/be/Baseline_params_table_2018_68pc.pdf
  //=============================================================================

  // baseline 2.17 base_plikHM_TTTEEE_lowl_lowE_lensing
  {"Planck2018EE", {
    {"h",           0.67321},
    {"Omega_m",     0.3158},
    {"Omega_b",     0.04938898},
    {"Omega_DE",    0.6842},
    {"w_0",         -1.0},
    {"w_a",         0.0},
    {"n_s",         0.96605},
    {"A_s",         2.1005e-9},
    {"k_p",         0.05},
    {"YHe",         0.245401},
    {"N_ur",        2.046},
    {"m_nu1",       0.06},
    {"m_nu2",       0.0},
    {"m_nu3",       0.0},
    {"Tcmb",        2.7255}}},

  // baseline 2.18 base_plikHM_TTTEEE_lowl_lowE_lensing_post_BAO
  {"Planck2018EE+BAO", {
    {"h",           0.67702},
    {"Omega_m",     0.3106},
    {"Omega_b",     0.04897284},
    {"Omega_DE",    0.6894},
    {"w_0",         -1.0},
    {"w_a",         0.0},
    {"n_s",         0.96824},
    {"A_s",         2.1073e-9},
    {"k_p",         0.05},
    {"YHe",         0.245425},
    {"N_ur",        2.046},
    {"m_nu1",       0.06},
    {"m_nu2",       0.0},
    {"m_nu3",       0.0},
    {"Tcmb",        2.7255}}},

  // baseline 2.19 base_plikHM_TTTEEE_lowl_lowE_lensing_post_Pantheon
  {"Planck2018EE+SN", {
    {"h",           0.6749},
    {"Omega_m",     0.3134},
    {"Omega_b",     0.04919537},
    {"Omega_DE",    0.6866},
    {"w_0",         -1.0},
    {"w_a",         0.0},
    {"n_s",         0.96654},
    {"A_s",         2.1020e-9},
    {"k_p",         0.05},
    {"YHe",         0.245411},
    {"N_ur",        2.046},
    {"m_nu1",       0.06},
    {"m_nu2",       0.0},
    {"m_nu3",       0.0},
    {"Tcmb",        2.7255}}},

  // baseline 2.20 base_plikHM_TTTEEE_lowl_lowE_lensing_post_BAO_Pantheon
  {"Planck2018EE+BAO+SN", {
    {"h",           0.67742},
    {"Omega_m",     0.3099},
    {"Omega_b",     0.048891054},
    {"Omega_DE",    0.6901},
    {"w_0",         -1.0},
    {"w_a",         0.0},
    {"n_s",         0.96822},
    {"A_s",         2.1064e-9},
    {"k_p",         0.05},
    {"YHe",         0.245421},
    {"N_ur",        2.046},
    {"m_nu1",       0.06},
    {"m_nu2",       0.0},
    {"m_nu3",       0.0},
    {"Tcmb",        2.7255}}}
};

/**
 * @brief Output all sets of cosmological parameters that we store internally
 * 
 */
void print_ParameterSets( void ){
  music::ilog << "Available cosmology parameter sets:" << std::endl;
  parameters p;
  p.print_available_sets();
  music::ilog << std::endl;
}

}// end namespace cosmology