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

#include <array>

//! implements general N-dim vectors of arbitrary primtive type with some arithmetic ops
template <int N, typename T = double>
struct vec_t
{
  std::array<T, N> data_;

  vec_t() {}

  vec_t(const vec_t<N, T> &v)
      : data_(v.data_) {}

  vec_t(vec_t<N, T> &&v)
      : data_(std::move(v.data_)) {}

  template <typename... E>
  vec_t(E... e)
      : data_{{std::forward<E>(e)...}}
  {
    static_assert(sizeof...(E) == N, "Brace-enclosed initialiser list doesn't match vec_t length!");
  }

  //! bracket index access to vector components
  T &operator[](size_t i) noexcept { return data_[i]; }

  //! const bracket index access to vector components
  const T &operator[](size_t i) const noexcept { return data_[i]; }

  // assignment operator
  vec_t<N, T> &operator=(const vec_t<N, T> &v) noexcept
  {
    data_ = v.data_;
    return *this;
  }

  //! implementation of summation of vec_t
  vec_t<N, T> operator+(const vec_t<N, T> &v) const noexcept
  {
    vec_t<N, T> res;
    for (int i = 0; i < N; ++i)
      res[i] = data_[i] + v[i];
    return res;
  }

  //! implementation of difference of vec_t
  vec_t<N, T> operator-(const vec_t<N, T> &v) const noexcept
  {
    vec_t<N, T> res;
    for (int i = 0; i < N; ++i)
      res[i] = data_[i] - v[i];
    return res;
  }

  //! implementation of unary negative
  vec_t<N, T> operator-() const noexcept
  {
    vec_t<N, T> res;
    for (int i = 0; i < N; ++i)
      res[i] = -data_[i];
    return res;
  }

  //! implementation of scalar multiplication
  template <typename T2>
  vec_t<N, T> operator*(T2 s) const noexcept
  {
    vec_t<N, T> res;
    for (int i = 0; i < N; ++i)
      res[i] = data_[i] * s;
    return res;
  }

  //! implementation of scalar division
  vec_t<N, T> operator/(T s) const noexcept
  {
    vec_t<N, T> res;
    for (int i = 0; i < N; ++i)
      res[i] = data_[i] / s;
    return res;
  }

  //! takes the absolute value of each element
  vec_t<N, T> abs(void) const noexcept
  {
    vec_t<N, T> res;
    for (int i = 0; i < N; ++i)
      res[i] = std::fabs(data_[i]);
    return res;
  }

  //! implementation of implicit summation of vec_t
  vec_t<N, T> &operator+=(const vec_t<N, T> &v) noexcept
  {
    for (int i = 0; i < N; ++i)
      data_[i] += v[i];
    return *this;
  }

  //! implementation of implicit subtraction of vec_t
  vec_t<N, T> &operator-=(const vec_t<N, T> &v) noexcept
  {
    for (int i = 0; i < N; ++i)
      data_[i] -= v[i];
    return *this;
  }

  //! implementation of implicit scalar multiplication of vec_t
  vec_t<N, T> &operator*=(T s) noexcept
  {
    for (int i = 0; i < N; ++i)
      data_[i] *= s;
    return *this;
  }

  //! implementation of implicit scalar division of vec_t
  vec_t<N, T> &operator/=(T s) noexcept
  {
    for (int i = 0; i < N; ++i)
      data_[i] /= s;
    return *this;
  }

  size_t size(void) const noexcept { return N; }
};

//! multiplication with scalar
template <typename T2, int N, typename T = double>
inline vec_t<N, T> operator*(T2 s, const vec_t<N, T> &v)
{
  vec_t<N, T> res;
  for (int i = 0; i < N; ++i)
    res[i] = v[i] * s;
  return res;
}
