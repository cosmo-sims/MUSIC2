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

#include <cmath>

inline double Meyer_scaling_function( double k, double kmax )
{
	constexpr double twopithirds{2.0*M_PI/3.0};
	constexpr double fourpithirds{4.0*M_PI/3.0};
	auto nu = []( double x ){ return x<0.0?0.0:(x<1.0?x:1.0); };

    k = std::abs(k)/kmax * 2 * M_PI;

	if( k < twopithirds ) return 1.0;
	else if( k< fourpithirds ){
		return std::cos( 0.5*M_PI * nu(3*k/(2*M_PI)-1.0) );
	}
	return 0.0;
}

inline double Shannon_scaling_function( double k, double kmax )
{
	if( std::abs(k) < kmax ) return 1.0;
	return 0.0;
}