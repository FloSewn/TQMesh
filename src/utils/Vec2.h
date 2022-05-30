/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once
#include <algorithm>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <float.h>

namespace CppUtils {

namespace Vec2Def
{
constexpr int coutprec  = 3;
constexpr int coutwidth = 6;
};

/***********************************************************
* Vec2d Class
***********************************************************/
template <class T>
class Vec2
{
  T e[2] { T{}, T{} };

public:
  T& x = e[0];
  T& y = e[1];

  // Construct
  Vec2() {}
  Vec2(T e0, T e1) : e{e0,e1} {}

  // Copy 
  Vec2(const Vec2<T>& v) : e{ v.x,v.y } {}
  Vec2<T>& operator=(const Vec2<T>& v) 
  { x=v.x; y=v.y; return *this; }

  // Move
  Vec2(Vec2<T>&& v) : e{v.x,v.y} {}
  Vec2<T>& operator=(Vec2<T>&& v) 
  { x=v.x; y=v.y; return *this; }

  // Vector access
  T operator[](int i) const 
  { 
    if (i == 0) return x;
    if (i == 1) return y;
    else throw std::runtime_error("Invalid access to Vec2");
  } 
  T& operator[](int i) 
  { 
    if (i == 0) return x;
    if (i == 1) return y;
    else throw std::runtime_error("Invalid access to Vec2");
  } 

  // Vector negation
  Vec2<T> operator-() const { return Vec2(-x,-y); }

  // Vector-vector addition
  Vec2<T>& operator+=(const Vec2<T> &v)
  { x += v.x; y += v.y; return *this; }

  // Vector-vector substraction
  Vec2<T>& operator-=(const Vec2<T> &v)
  { y -= v.x; y -= v.y; return *this; }

  // Vector-vector multiplication
  Vec2<T>& operator*=(const Vec2<T> &v)
  { x *= v.x; y *= v.y; return *this; }

  // Vector-vector division
  Vec2<T>& operator/=(const Vec2<T> &v)
  { x /= v.x; y /= v.y; return *this; }

  // Vector-scalar addition
  Vec2<T>& operator+=(const T t)
  { x += t; y += t; return *this; }

  // Vector-scalar substraction
  Vec2<T>& operator-=(const T t)
  { x -= t; y -= t; return *this; }

  // Vector-scalar multiplication
  Vec2<T>& operator*=(const T t)
  { x *= t; y *= t; return *this; }

  // Vector-scalar division
  Vec2<T>& operator/=(const T t)
  { return *this *= 1./t; }

  // Estimate vector length
  double length() const 
  { return std::sqrt(length_squared()); }

  // Estimate squared vector length
  double length_squared() const 
  { return x*x + y*y; }

  // Return true if each vector-value is close to zero 
  bool near_zero_values(const T s=DBL_EPSILON) const
  { return (fabs(x) < s) && (fabs(y) < s); }

  // Return true if vector-length is close to zero 
  bool near_zero_length(const T s=DBL_EPSILON) const
  { return (length() < s); }

};

/***********************************************************
* Vec2 utility functions
***********************************************************/

// Output 
template <typename T>
inline std::ostream& operator<<(std::ostream &out, 
                                const Vec2<T> &v)
{ return out 
  << std::fixed << std::setprecision(Vec2Def::coutprec)  
  << '(' << std::setw(Vec2Def::coutwidth) << v.x 
  << ',' << std::setw(Vec2Def::coutwidth) << v.y << ')'; 
}

// Equality
template <typename T>
inline bool operator==(const Vec2<T>& a, const Vec2<T>& b)
{ return (a-b).near_zero_values(); }

// Un-Equality
template <typename T>
inline bool operator!=(const Vec2<T>& a, const Vec2<T>& b)
{ return !(a==b); }

// Vector-Vector addition
template <typename T>
inline Vec2<T> operator+(const Vec2<T> &u,const Vec2<T> &v)
{ return Vec2<T>(u.x+v.x, u.y+v.y); }

// Vector-Vector substraction
template <typename T>
inline Vec2<T> operator-(const Vec2<T> &u,const Vec2<T> &v)
{ return Vec2<T>(u.x-v.x, u.y-v.y); }

// Vector-Vector multiplication
template <typename T>
inline Vec2<T> operator*(const Vec2<T> &u,const Vec2<T> &v)
{ return Vec2<T>(u.x*v.x, u.y*v.y); }

// Vector-Vector division
template <typename T>
inline Vec2<T> operator/(const Vec2<T> &u,const Vec2<T> &v)
{ return Vec2<T>(u.x/v.x, u.y/v.y); }


// Vector-scalar addition
template <typename T>
inline Vec2<T> operator+(const Vec2<T> &u, const T v)
{ return Vec2<T>(u.x+v, u.y+v); }
template <typename T>
inline Vec2<T> operator+(const T v, const Vec2<T> &u)
{ return Vec2<T>(u.x+v, u.y+v); }

// Vector-scalar substraction
template <typename T>
inline Vec2<T> operator-(const Vec2<T> &u, const T v)
{ return Vec2<T>(u.x-v, u.y-v); }
template <typename T>
inline Vec2<T> operator-(const T v, const Vec2<T> &u)
{ return Vec2<T>(v-u.x, v-u.y); }

// Vector-scalar multiplication
template <typename T>
inline Vec2<T> operator*(const Vec2<T> &u, const T v)
{ return Vec2<T>(u.x*v, u.y*v); }
template <typename T>
inline Vec2<T> operator*(const T v, const Vec2<T> &u)
{ return Vec2<T>(u.x*v, u.y*v); }

// Vector-scalar division
template <typename T>
inline Vec2<T> operator/(const Vec2<T> &u, const T v)
{ return Vec2<T>(u.x/v, u.y/v); }
template <typename T>
inline Vec2<T> operator/(const T v, const Vec2<T> &u)
{ return Vec2<T>(v/u.x, v/u.y); }


// Scalar product
template <typename T>
inline T dot(const Vec2<T> &u, const Vec2<T> &v)
{ return u.x * v.x + u.y * v.y; }
// Cross product
template <typename T>
inline T cross(const Vec2<T> &u, const Vec2<T> &v)
{ return u.x * v.y - u.y * v.x; }

// Angle
template <typename T>
inline T angle(const Vec2<T> &u, const Vec2<T> &v)
{ 
  const T cos_a = dot(u,v) / ( u.length() * v.length() );
  return acos( std::clamp( cos_a ,-1.0, 1.0 ) ); 
}

// Unit-vector
template <typename T>
inline Vec2<T> unit_vector(const Vec2<T>& u)
{ return u / u.length(); }


/***********************************************************
* integer / double vector aliases
***********************************************************/
using Vec2i = Vec2<int>;
using Vec2d = Vec2<double>;
using Vec2f = Vec2<float>;


} // namespace CppUtils
