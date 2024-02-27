/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once
#include <algorithm>
#include <functional>
#include <array>
#include <initializer_list>
#include <numeric>
#include <limits>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <float.h>

#include "MathUtility.h"

namespace CppUtils {

/**********************************************************************
* VecND class:
* This class describes an n-dimensional vector
**********************************************************************/
template <class T, std::size_t N>
class VecND
{
  std::array<T,N> entries_ {};  

public:

  /*------------------------------------------------------------------
  | Variables for ostream formatting / precision 
  ------------------------------------------------------------------*/
  static inline std::size_t os_precision = 3;
  static inline std::size_t os_width = 0; 
  static inline std::size_t cmp_ulp = 2;

  /*------------------------------------------------------------------
  | Iterators
  ------------------------------------------------------------------*/
  using iterator       = typename std::array<T,N>::iterator;
  using const_iterator = typename std::array<T,N>::const_iterator;

  iterator begin() { return entries_.begin(); }
  iterator end() { return entries_.end(); }

  const_iterator begin() const { return entries_.begin(); }
  const_iterator end() const { return entries_.end(); }

  const_iterator cbegin() const { return entries_.cbegin(); }
  const_iterator cend() const { return entries_.cend(); }

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  VecND() {}
  VecND(const std::array<T,N>& e) : entries_ {e} {}
  VecND(std::array<T,N>&& e) : entries_ {std::move(e)} {}

  /*------------------------------------------------------------------
  | Constructor with varyadic argument list
  | Requires C++-17
  | 
  | Source: https://stackoverflow.com/questions/37351061/how-can-i-\
  |         create-a-c-constructor-that-accepts-a-variable-number-\
  |         of-ints
  ------------------------------------------------------------------*/
  template<
    class... TT,
    class E = std::enable_if_t<(std::is_same_v<TT, T> && ...)>
  >
  VecND(TT... tt) : entries_ { tt... } {}

  /*------------------------------------------------------------------
  | Copy
  ------------------------------------------------------------------*/
  VecND(const VecND<T,N>& v) : entries_ { v.entries_ } {}
  VecND<T,N>& operator=(const VecND<T,N>& v)
  { entries_ = v.entries_; return *this; }

  /*------------------------------------------------------------------
  | Move
  ------------------------------------------------------------------*/
  VecND(const VecND<T,N>&& v) : entries_ { std::move(v.entries_) } {}
  VecND<T,N>& operator=(VecND<T,N>&& v)
  { entries_ = std::move(v.entries_); return *this; }

  /*------------------------------------------------------------------
  | Access
  ------------------------------------------------------------------*/
  T operator[](std::size_t i) const
  { return entries_[i]; }

  T& operator[](std::size_t i)
  { return entries_[i]; }

  T& get_x() { return entries_[0]; }
  const T& get_x() const { return entries_[0]; }

  T& get_y() 
  { 
    if constexpr (N < 2)
      return entries_[0];
    else
      return entries_[1]; 
  }
  const T& get_y() const 
  { 
    if constexpr (N < 2)
      return entries_[0];
    else
      return entries_[1]; 
  }

  T& get_z() 
  { 
    if constexpr (N < 3)
      return entries_[0];
    else
      return entries_[2]; 
  }
  const T& get_z() const 
  { 
    if constexpr (N < 3)
      return entries_[0];
    else
      return entries_[2]; 
  }

  T& x = get_x(); 
  T& y = get_y(); 
  T& z = get_z(); 

  /*------------------------------------------------------------------
  | Negation
  ------------------------------------------------------------------*/
  VecND<T,N> operator-() const 
  {
    std::array<T,N> out {};
    std::transform(cbegin(), cend(), out.begin(),
                   std::negate<T>());
    return VecND<T,N> { std::move(out) };
  }

  /*------------------------------------------------------------------
  | VecND-VecND elementwise addition
  ------------------------------------------------------------------*/
  VecND<T,N> operator+=(const VecND<T,N> &v) 
  {
    std::transform(cbegin(), cend(), v.cbegin(), begin(),
                   std::plus<T>());
    return *this;
  }

  /*------------------------------------------------------------------
  | VecND-VecND elementwise subtraction
  ------------------------------------------------------------------*/
  VecND<T,N> operator-=(const VecND<T,N> &v)  
  {
    std::transform(cbegin(), cend(), v.cbegin(), begin(),
                   std::minus<T>());
    return *this;
  }

  /*------------------------------------------------------------------
  | VecND-VecND elementwise multiplication
  ------------------------------------------------------------------*/
  VecND<T,N> operator*=(const VecND<T,N> &v)  
  {
    std::transform(cbegin(), cend(), v.cbegin(), begin(),
                   std::multiplies<T>());
    return *this;
  }

  /*------------------------------------------------------------------
  | VecND-VecND elementwise division
  ------------------------------------------------------------------*/
  VecND<T,N> operator/=(const VecND<T,N> &v)  
  {
    std::transform(cbegin(), cend(), v.cbegin(), begin(),
                   std::divides<T>());
    return *this;
  }

  /*------------------------------------------------------------------
  | VecND-Scalar addition 
  ------------------------------------------------------------------*/
  VecND<T,N> operator+=(const T &v) 
  {
    std::transform(cbegin(), cend(), begin(),
                   [&](auto const& elem) { return elem + v; });
    return *this;
  }

  /*------------------------------------------------------------------
  | VecND-Scalar subtraction 
  ------------------------------------------------------------------*/
  VecND<T,N> operator-=(const T &v) 
  {
    std::transform(cbegin(), cend(), begin(),
                   [&](auto const& elem) { return elem - v; });
    return *this;
  }

  /*------------------------------------------------------------------
  | VecND-Scalar multiplication 
  ------------------------------------------------------------------*/
  VecND<T,N> operator*=(const T &v) 
  {
    std::transform(cbegin(), cend(), begin(),
                   [&](auto const& elem) { return elem * v; });
    return *this;
  }

  /*------------------------------------------------------------------
  | VecND-Scalar division 
  ------------------------------------------------------------------*/
  VecND<T,N> operator/=(const T &v) 
  {
    std::transform(cbegin(), cend(), begin(),
                   [&](auto const& elem) { return elem / v; });
    return *this;
  }

  /*------------------------------------------------------------------
  | L2-Norm of VecND
  ------------------------------------------------------------------*/
  T norm_sqr() const
  { return std::inner_product(cbegin(), cend(), cbegin(), T {}); }

  T norm() const
  { return std::sqrt(norm_sqr()); }

  /*------------------------------------------------------------------
  | Dot product 
  ------------------------------------------------------------------*/
  T dot(const VecND<T,N>& v) const 
  { return std::inner_product(cbegin(), cend(), v.cbegin(), T {}); }

  /*------------------------------------------------------------------
  | Cross product 
  | SFINAE -> Returns T {} for N != 2 and N != 3
  | Sources: https://stackoverflow.com/questions/57449491/correct-way\
  |          -to-use-stdenable-if
  |          https://stackoverflow.com/questions/13786479/using-c11-\
  |          stdenable-if-to-enable-member-function-if-vector-is-\
  |          specific-lengt
  ------------------------------------------------------------------*/
  template<std::size_t NN=N>
  typename std::enable_if_t<NN==2, T>
  cross(const VecND<T,N>& v) const
  { return (x * v.y) - (y * v.x); }

  template<std::size_t NN=N>
  typename std::enable_if_t<NN==3, VecND<T,N>>
  cross(const VecND<T,N>& v) const
  { return {
      y * v.z - z * v.y,
      z * v.x - x * v.z,
      x * v.y - y * v.x,
  };}

  template<std::size_t NN=N>
  typename std::enable_if_t<(NN!=2 && NN!=3), T>
  cross(const VecND<T,N>& v) const
  { return {}; }

  /*------------------------------------------------------------------
  | Angle  
  ------------------------------------------------------------------*/
  T angle(const VecND<T,N>& v) const 
  {
    const T cos_ang = this->dot(v) / ( this->norm() * v.norm() );
    return acos( CLAMP( cos_ang ,-1.0, 1.0 ) ); 
  }

  /*------------------------------------------------------------------
  | Check if vector entries are near zero
  | ulp == units in the last place 
  ------------------------------------------------------------------*/
  bool is_zero(std::size_t ulp = VecND<T,N>::cmp_ulp) const
  {
    T v = norm_sqr();
    return (v <= std::numeric_limits<T>::epsilon() * v * ulp);
  }

  /*------------------------------------------------------------------
  | Returns the sum of all vector entries
  ------------------------------------------------------------------*/
  T sum() const
  { 
    return std::accumulate(cbegin(), cend(), 0);
  }

  /*------------------------------------------------------------------
  | Returns the product of all vector entries
  ------------------------------------------------------------------*/
  T product() const 
  {
    return std::accumulate(cbegin(), cend(), 1.0, 
                           std::multiplies<T>());
  }

  /*------------------------------------------------------------------
  | Extrema  
  ------------------------------------------------------------------*/
  T min() const { return *std::min_element(cbegin(), cend()); }
  T max() const { return *std::max_element(cbegin(), cend()); }

}; // VecND

/*********************************************************************
* VecND utility functions
*********************************************************************/

/*--------------------------------------------------------------------
| Output to ostream
--------------------------------------------------------------------*/
template <typename T, std::size_t N>
inline std::ostream& operator<<(std::ostream &os, const VecND<T,N> &v)
{
  os << std::fixed << std::setprecision(VecND<T,N>::os_precision) << '(';

  for (std::size_t i = 0; i < N; ++i)
    os << std::setw(VecND<T,N>::os_width) << v[i]
       << (i < N-1 ? ',' : ')');

  return os;
}

/*--------------------------------------------------------------------
| Equality / ineuality
--------------------------------------------------------------------*/
template <typename T, std::size_t N>
inline bool operator==(const VecND<T,N>& a, const VecND<T,N>& b)
{ return (a-b).is_zero(); }

template <typename T, std::size_t N>
inline bool operator!=(const VecND<T,N>& a, const VecND<T,N>& b)
{ return !(a==b); }

/*--------------------------------------------------------------------
| VecND-VecND elementwise addition
--------------------------------------------------------------------*/
template <typename T, std::size_t N>
inline VecND<T,N> operator+(const VecND<T,N>& a, const VecND<T,N>& b)
{ 
  std::array<T,N> out {};
  std::transform(a.cbegin(), a.cend(), b.cbegin(), out.begin(),
                 std::plus<T>());
  return VecND<T,N> { std::move(out) };
}

/*--------------------------------------------------------------------
| VecND-VecND elementwise subtraction
--------------------------------------------------------------------*/
template <typename T, std::size_t N>
inline VecND<T,N> operator-(const VecND<T,N>& a, const VecND<T,N>& b)
{ 
  std::array<T,N> out {};
  std::transform(a.cbegin(), a.cend(), b.cbegin(), out.begin(),
                 std::minus<T>());
  return VecND<T,N> { std::move(out) };
}

/*--------------------------------------------------------------------
| VecND-VecND elementwise multiplication
--------------------------------------------------------------------*/
template <typename T, std::size_t N>
inline VecND<T,N> operator*(const VecND<T,N>& a, const VecND<T,N>& b)
{ 
  std::array<T,N> out {};
  std::transform(a.cbegin(), a.cend(), b.cbegin(), out.begin(),
                 std::multiplies<T>());
  return VecND<T,N> { std::move(out) };
}

/*--------------------------------------------------------------------
| VecND-VecND elementwise division
--------------------------------------------------------------------*/
template <typename T, std::size_t N>
inline VecND<T,N> operator/(const VecND<T,N>& a, const VecND<T,N>& b)
{ 
  std::array<T,N> out {};
  std::transform(a.cbegin(), a.cend(), b.cbegin(), out.begin(),
                 std::divides<T>());
  return VecND<T,N> { std::move(out) };
}


/*--------------------------------------------------------------------
| VecND-scalar addition 
--------------------------------------------------------------------*/
template <typename T, std::size_t N>
inline VecND<T,N> operator+(const VecND<T,N>& a, const T& b)
{
  std::array<T,N> out {};
  std::transform(a.cbegin(), a.cend(), out.begin(),
                 [&](auto const& elem) { return elem + b; });
  return VecND<T,N> { std::move(out) };
}
template <typename T, std::size_t N>
inline VecND<T,N> operator+(const T& b, const VecND<T,N>& a)
{
  std::array<T,N> out {};
  std::transform(a.cbegin(), a.cend(), out.begin(),
                 [&](auto const& elem) { return elem + b; });
  return VecND<T,N> { std::move(out) };
}

/*--------------------------------------------------------------------
| VecND-scalar subtraction 
--------------------------------------------------------------------*/
template <typename T, std::size_t N>
inline VecND<T,N> operator-(const VecND<T,N>& a, const T& b)
{
  std::array<T,N> out {};
  std::transform(a.cbegin(), a.cend(), out.begin(),
                 [&](auto const& elem) { return elem - b; });
  return VecND<T,N> { std::move(out) };
}
template <typename T, std::size_t N>
inline VecND<T,N> operator-(const T& b, const VecND<T,N>& a)
{
  std::array<T,N> out {};
  std::transform(a.cbegin(), a.cend(), out.begin(),
                [&](T elem) { return b - elem; });
  return VecND<T,N> { std::move(out) };
}

/*--------------------------------------------------------------------
| VecND-scalar multiplication 
--------------------------------------------------------------------*/
template <typename T, std::size_t N>
inline VecND<T,N> operator*(const VecND<T,N>& a, const T& b)
{
  std::array<T,N> out {};
  std::transform(a.cbegin(), a.cend(), out.begin(),
                 [&](auto const& elem) { return elem * b; });
  return VecND<T,N> { std::move(out) };
}
template <typename T, std::size_t N>
inline VecND<T,N> operator*(const T& b, const VecND<T,N>& a)
{
  std::array<T,N> out {};
  std::transform(a.cbegin(), a.cend(), out.begin(),
                 [&](auto const& elem) { return elem * b; });
  return VecND<T,N> { std::move(out) };
}

/*--------------------------------------------------------------------
| VecND-scalar division 
--------------------------------------------------------------------*/
template <typename T, std::size_t N>
inline VecND<T,N> operator/(const VecND<T,N>& a, const T& b)
{
  std::array<T,N> out {};
  std::transform(a.cbegin(), a.cend(), out.begin(),
                 [&](auto const& elem) { return elem / b; });
  return VecND<T,N> { std::move(out) };
}

template <typename T, std::size_t N>
inline VecND<T,N> operator/(const T& b, const VecND<T,N>& a)
{
  std::array<T,N> out {};
  std::transform(a.cbegin(), a.cend(), out.begin(),
                [&](T elem) { return b / elem; } );
  return VecND<T,N> { std::move(out) };
}

/*--------------------------------------------------------------------
| Dot product
--------------------------------------------------------------------*/
template <typename T, std::size_t N>
inline T dot(const VecND<T,N>& a, const VecND<T,N>& b)
{ return std::inner_product(a.cbegin(), a.cend(), b.cbegin(), T {}); }


/*--------------------------------------------------------------------
| Cross product
--------------------------------------------------------------------*/
template <typename T, std::size_t N, std::size_t NN=N>
typename std::enable_if_t<NN==2, T>
cross(const VecND<T,N>& a, const VecND<T,N>& b)
{ return (a.x * b.y) - (a.y * b.x); }

template <typename T, std::size_t N, std::size_t NN=N>
typename std::enable_if_t<NN==3, VecND<T,N>>
cross(const VecND<T,N>& a, const VecND<T,N>& b)
{ return {
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x,
};}

template <typename T, std::size_t N, std::size_t NN=N>
typename std::enable_if_t<(NN!=2 && NN!=3), T>
cross(const VecND<T,N>& a, const VecND<T,N>& b)
{ return {}; }

/*--------------------------------------------------------------------
| Angle
--------------------------------------------------------------------*/
template <typename T, std::size_t N>
inline T angle(const VecND<T,N>& a, const VecND<T,N>& b)
{
  const T cos_ang = dot(a,b) / ( a.norm() * b.norm() );
  return acos( CLAMP( cos_ang ,-1.0, 1.0 ) ); 
}

/**********************************************************************
* Some commont type definitions
**********************************************************************/
using Vec2i = VecND<int,2>;
using Vec2d = VecND<double,2>;
using Vec2f = VecND<float,2>;

using Vec3i = VecND<int,3>;
using Vec3d = VecND<double,3>;
using Vec3f = VecND<float,3>;

using Vec4i = VecND<int,4>;
using Vec4d = VecND<double,4>;
using Vec4f = VecND<float,4>;

template <typename T>
using Vec2 = VecND<T,2>;
template <typename T>
using Vec3 = VecND<T,3>;
template <typename T>
using Vec4 = VecND<T,4>;


} // namepsace CppUtils
