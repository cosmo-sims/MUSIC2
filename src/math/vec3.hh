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

//! implements a simple class of 3-vectors of arbitrary scalar type
template< typename T >
class vec3_t{
private:
    //! holds the data
    std::array<T,3> data_;
    
public: 
    //! expose access to elements via references
    T &x,&y,&z;

    //! empty constructor
    vec3_t()
    : data_{{T(0),T(0),T(0)}},x(data_[0]),y(data_[1]),z(data_[2]){}

    //! copy constructor
    vec3_t( const vec3_t<T> &v)
    : data_(v.data_), x(data_[0]),y(data_[1]),z(data_[2]){}

    //! copy constructor for non-const reference, needed to avoid variadic template being called for non-const reference
    vec3_t( vec3_t<T>& v)
    : data_(v.data_), x(data_[0]),y(data_[1]),z(data_[2]){}

    //! move constructor
    vec3_t( vec3_t<T> &&v)
    : data_(std::move(v.data_)), x(data_[0]), y(data_[1]), z(data_[2]){}

    //! construct vec3_t from initializer list
//    template<typename ...E>
//    vec3_t(E&&...e) 
//    : data_{{std::forward<E>(e)...}}, x{data_[0]}, y{data_[1]}, z{data_[2]}
//    {}
    vec3_t( T a, T b, T c )
        : data_{{a,b,c}}, x(data_[0]), y(data_[1]), z(data_[2]){}
  
    explicit vec3_t( T a)
        : data_{{a,a,a}}, x(data_[0]), y(data_[1]), z(data_[2]){}
    
    //! bracket index access to vector components
    T &operator[](size_t i) noexcept{ return data_[i];}
    
    //! const bracket index access to vector components
    const T &operator[](size_t i) const noexcept { return data_[i]; }

    // assignment operator
    vec3_t<T>& operator=( const vec3_t<T>& v ) noexcept { data_=v.data_; return *this; }

    //! implementation of summation of vec3_t
    vec3_t<T> operator+( const vec3_t<T>& v ) const noexcept{ return vec3_t<T>({x+v.x,y+v.y,z+v.z}); }

    //! implementation of difference of vec3_t
    vec3_t<T> operator-( const vec3_t<T>& v ) const noexcept{ return vec3_t<T>({x-v.x,y-v.y,z-v.z}); }

    //! implementation of unary negative
    vec3_t<T> operator-() const noexcept{ return vec3_t<T>({-x,-y,-z}); }

    //! implementation of scalar multiplication
    vec3_t<T> operator*( T s ) const noexcept{ return vec3_t<T>({x*s,y*s,z*s}); }

    //! implementation of scalar division
    vec3_t<T> operator/( T s ) const noexcept{ return vec3_t<T>({x/s,y/s,z/s}); }

    //! implementation of += operator
    vec3_t<T>& operator+=( const vec3_t<T>& v ) noexcept{ x+=v.x; y+=v.y; z+=v.z; return *this; }

    //! implementation of -= operator
    vec3_t<T>& operator-=( const vec3_t<T>& v ) noexcept{ x-=v.x; y-=v.y; z-=v.z; return *this; }

    //! multiply with scalar
    vec3_t<T>& operator*=( T s ) noexcept{ x*=s; y*=s; z*=s; return *this; }
    
    //! divide by scalar
    vec3_t<T>& operator/=( T s ) noexcept{ x/=s; y/=s; z/=s; return *this; }

    //! compute dot product with another vector
    T dot(const vec3_t<T> &a) const noexcept
    {
        return data_[0] * a.data_[0] + data_[1] * a.data_[1] + data_[2] * a.data_[2];
    }
    
    //! returns 2-norm squared of vector
    T norm_squared(void) const noexcept { return this->dot(*this); }

    //! returns 2-norm of vector
    T norm(void) const noexcept { return std::sqrt( this->norm_squared() ); }

    //! wrap absolute vector to box of size p
    vec3_t<T>& wrap_abs( T p = 1.0 ) noexcept{
        for( auto& x : data_ ) x = std::fmod( 2*p + x, p );
        return *this;
    }

    //! wrap relative vector to box of size p
    vec3_t<T>& wrap_rel( T p = 1.0 ) noexcept{
        for( auto& x : data_ ) x = (x<-p/2)? x+p : (x>=p/2)? x-p : x;
        return *this;
    }

    //! takes the absolute value of each element
    vec3_t<T> abs(void) const noexcept
    {
        vec3_t<T> res;
        for (int i = 0; i < 3; ++i)
            res[i] = std::fabs(data_[i]);
        return res;
    }

    //! ordering, allows 3d sorting of vec3_ts
    bool operator<( const vec3_t<T>& o ) const noexcept{
        if( x!=o.x ) return x<o.x?true:false;
        if( y!=o.y ) return y<o.y?true:false;
        if( z!=o.z ) return z<o.z?true:false;
        return false;
    }
};

//! multiplication with scalar
template<typename T>
vec3_t<T> operator*( T s, const vec3_t<T>& v ){
    return vec3_t<T>({v.x*s,v.y*s,v.z*s});
}
