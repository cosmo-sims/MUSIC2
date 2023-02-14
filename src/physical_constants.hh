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
#pragma once

//------------------------------------------------------------------------------------
// physical constants for convenience, all values have been taken from
// the 2018 edition of the Particle Data Group Booklet,
// http://pdg.lbl.gov/2019/mobile/reviews/pdf/rpp2018-rev-phys-constants-m.pdf

namespace phys_const
{
// helper value of pi so that we don't need to include any other header just for this
static constexpr double pi_ = 3.141592653589793115997963468544185161590576171875;

//--- unit conversions ---------------------------------------------------

// 1 Mpc in m
static constexpr double Mpc_SI = 3.0857e22;

// 1 Gyr in s
static constexpr double Gyr_SI = 3.1536e16;

// 1 eV in J
static constexpr double eV_SI = 1.602176487e-19;

// 1 erg in J
static constexpr double erg_SI = 1e-7;

//--- physical constants ------------------------------------------------

// speed of light c in m/s
static constexpr double c_SI = 2.99792458e8;

// gravitational constant G in m^3/s^2/kg
static constexpr double G_SI = 6.6740800e-11;

// Boltzmann constant k_B in kg m^2/s^2/K
static constexpr double kB_SI = 1.38064852e-23;

// reduced Planck's quantum \hbar in kg m^2/s
static constexpr double hbar_SI = 1.054571800e-34;

// Stefan-Boltzmann constant sigma in J/m^2/s/K^-4
static constexpr double sigma_SI = (pi_ * pi_) * (kB_SI * kB_SI * kB_SI * kB_SI) / 60. / (hbar_SI * hbar_SI * hbar_SI) / (c_SI * c_SI);

// electron mass in kg
static constexpr double me_SI = 9.10938356e-31;

// proton mass in kg
static constexpr double mp_SI = 1.672621898e-27;

// unified atomic mass unit (u) in kg
static constexpr double u_SI = 1.660539040e-27;

// critical density of the Universe in h^2 kg/m^3
static constexpr double rhocrit_h2_SI = 3 * 1e10 / (8 * pi_ * G_SI) / Mpc_SI / Mpc_SI;

// solar mass in kg
static constexpr double Msun_SI = 1.98847e30;

} // namespace phys_const