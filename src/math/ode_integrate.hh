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

//! ODE integration namespace
namespace ode_integrate
{

// simple Runge-Kutta 4th order step without error estimate
template <typename vector_t, typename function_t>
inline void rk4_step(double h, double &t, vector_t &y, function_t f)
{
    vector_t k1(h * f(t, y));
    vector_t k2(h * f(t + h / 2, y + k1 / 2));
    vector_t k3(h * f(t + h / 2, y + k2 / 2));
    vector_t k4(h * f(t + h, y + k3));
    y += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    t += h;
}

// Cash-Karp modified Runge-Kutta scheme, 5th order with 4th order error estimate
// see Press & Teukolsky (1992): "Adaptive Stepsize Runge-Kutta Integration"
// in Computers in Physics 6, 188 (1992); doi: 10.1063/1.4823060
template <typename vector_t, typename function_t>
inline vector_t ckrk5_step(double h, double &t, vector_t &y, function_t f)
{
  static constexpr double
      a2 = 0.20,
      a3 = 0.30, a4 = 0.60, a5 = 1.0, a6 = 0.8750,
      b21 = 0.20,
      b31 = 3.0 / 40.0, b32 = 9.0 / 40.0,
      b41 = 0.30, b42 = -0.90, b43 = 1.20,
      b51 = -11.0 / 54.0, b52 = 2.50, b53 = -70.0 / 27.0, b54 = 35.0 / 27.0,
      b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0, b63 = 575.0 / 13824.0, b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0,
      c1 = 37.0 / 378.0, c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0,
      dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
      dc4 = c4 - 13525.0 / 55296.0, dc5 = -277.0 / 14336.0, dc6 = c6 - 0.250;

  vector_t k1(h * f(t, y));
  vector_t k2(h * f(t + a2 * h, y + b21 * k1));
  vector_t k3(h * f(t + a3 * h, y + b31 * k1 + b32 * k2));
  vector_t k4(h * f(t + a4 * h, y + b41 * k1 + b42 * k2 + b43 * k3));
  vector_t k5(h * f(t + a5 * h, y + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4));
  vector_t k6(h * f(t + a6 * h, y + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5));

  y += c1 * k1 + c3 * k3 + c4 * k4 + c6 * k6;

  return dc1 * k1 + dc3 * k3 + dc4 * k4 + dc5 * k5 + dc6 * k6;
}

// Adaptive step-size quality-controlled routine for ckrk5_step, see
// Press & Teukolsky (1992): "Adaptive Stepsize Runge-Kutta Integration"
// in Computers in Physics 6, 188 (1992); doi: 10.1063/1.4823060
template <typename vector_t, typename function_t>
inline void rk_step_qs(double htry, double &t, vector_t &y, vector_t &yscale, function_t f, double eps, double &hdid, double &hnext)
{
  static constexpr double SAFETY{0.9};
  static constexpr double PSHRNK{-0.25};
  static constexpr double PGROW{-0.2};
  static constexpr double ERRCON{1.89e-4};

  auto h(htry);
  vector_t ytemp(y);
  vector_t yerr;
  double errmax;

do_ckrk5trialstep:
  yerr = ckrk5_step(h, t, ytemp, f);
  errmax = 0.0;
  for (size_t i = 0; i < yerr.size(); ++i)
  {
    errmax = std::max(errmax, std::fabs(yerr[i] / yscale[i]));
  }
  errmax = errmax / eps;
  if (errmax > 1.0)
  {
    h *= std::max(0.1, SAFETY*std::pow(errmax, PSHRNK));
    if (t + h == t)
    {
      std::cerr << "stepsize underflow in rkqs" << std::endl;
      abort();
    }
    goto do_ckrk5trialstep;
  }
  else
  {
    if( errmax > ERRCON ){
      hnext = h * SAFETY * std::pow(errmax, PGROW);
    }else{
      hnext = 5*h;
    }
    hdid = h;
    t += h;
    y = ytemp;
  }
}


} // namespace ode_integrate