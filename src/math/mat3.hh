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

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include <math/vec3.hh>

//! class for 3x3 matrix calculations
template<typename T>
class mat3_t{
protected:
    std::array<T,9> data_;
    std::array<double,9> data_double_;
    gsl_matrix_view m_;
    gsl_vector *eval_;
    gsl_matrix *evec_;
	gsl_eigen_symmv_workspace * wsp_;
    bool bdid_alloc_gsl_;
						
    void init_gsl(){
        // allocate memory for GSL operations if we haven't done so yet
        if( !bdid_alloc_gsl_ )
        {
            if( typeid(T)!=typeid(double) ){
                m_ = gsl_matrix_view_array ( &data_double_[0], 3, 3);
            }else{
                m_ = gsl_matrix_view_array ( (double*)(&data_[0]), 3, 3); // this should only ever be called for T==double so cast is to avoid constexpr ifs from C++17
            }
            eval_ = gsl_vector_alloc (3);
            evec_ = gsl_matrix_alloc (3, 3);
            wsp_ = gsl_eigen_symmv_alloc (3);
            bdid_alloc_gsl_ = true;
        }

        if( typeid(T)!=typeid(double) ){
            for( int i=0; i<9; ++i ) data_double_[i] = double(data_[i]);
        }
    }

    void free_gsl(){
        // free memory for GSL operations if it was allocated
        if( bdid_alloc_gsl_ )
        {
            gsl_eigen_symmv_free (wsp_);
            gsl_vector_free (eval_);
            gsl_matrix_free (evec_);
        }
    }

public:

    mat3_t()
    : bdid_alloc_gsl_(false) 
    {}

    //! copy constructor
    mat3_t( const mat3_t<T> &m)
    : data_(m.data_), bdid_alloc_gsl_(false) 
    {}
    
    //! move constructor
    mat3_t( mat3_t<T> &&m)
    : data_(std::move(m.data_)), bdid_alloc_gsl_(false) 
    {}

    //! construct mat3_t from initializer list
    template<typename ...E>
    mat3_t(E&&...e) 
    : data_{{std::forward<E>(e)...}}, bdid_alloc_gsl_(false)
    {}

    mat3_t<T>& operator=(const mat3_t<T>& m) noexcept{
        data_ = m.data_;
        return *this;
    }

    mat3_t<T>& operator=(const mat3_t<T>&& m) noexcept{
        data_ = std::move(m.data_);
        return *this;
    }

    //! destructor
    ~mat3_t(){
        this->free_gsl();
    }
    
    //! bracket index access to vector components
    T &operator[](size_t i) noexcept { return data_[i];}
    
    //! const bracket index access to vector components
    const T &operator[](size_t i) const noexcept { return data_[i]; }

    //! matrix 2d index access
    T &operator()(size_t i, size_t j) noexcept { return data_[3*i+j]; }

    //! const matrix 2d index access
    const T &operator()(size_t i, size_t j) const noexcept { return data_[3*i+j]; }

    //! in-place addition
    mat3_t<T>& operator+=( const mat3_t<T>& rhs ) noexcept{
        for (size_t i = 0; i < 9; ++i) {
           (*this)[i] += rhs[i];
        }
        return *this;
    }

    //! in-place subtraction
    mat3_t<T>& operator-=( const mat3_t<T>& rhs ) noexcept{
        for (size_t i = 0; i < 9; ++i) {
           (*this)[i] -= rhs[i];
        }
        return *this;
    }

    void zero() noexcept{
        for (size_t i = 0; i < 9; ++i) data_[i]=0;
    }

    void eigen( vec3_t<T>& evals, vec3_t<T>& evec1, vec3_t<T>& evec2, vec3_t<T>& evec3_t )
    {
        this->init_gsl();

        gsl_eigen_symmv (&m_.matrix, eval_, evec_, wsp_);
        gsl_eigen_symmv_sort (eval_, evec_, GSL_EIGEN_SORT_VAL_ASC);

        for( int i=0; i<3; ++i ){
            evals[i] = gsl_vector_get( eval_, i );
            evec1[i] = gsl_matrix_get( evec_, i, 0 );
            evec2[i] = gsl_matrix_get( evec_, i, 1 );
            evec3_t[i] = gsl_matrix_get( evec_, i, 2 );
        }
    }
};

template<typename T>
constexpr const mat3_t<T> operator+(const mat3_t<T> &lhs, const mat3_t<T> &rhs) noexcept
{
    mat3_t<T> result;
    for (size_t i = 0; i < 9; ++i) {
        result[i] = lhs[i] + rhs[i];
    }
    return result;
}

// matrix - vector multiplication
template<typename T>
inline vec3_t<T> operator*( const mat3_t<T> &A, const vec3_t<T> &v ) noexcept
{
    vec3_t<T> result;
    for( int mu=0; mu<3; ++mu ){
        result[mu] = 0.0;
        for( int nu=0; nu<3; ++nu ){
            result[mu] += A(mu,nu)*v[nu];
        }
    }
    return result;
}

