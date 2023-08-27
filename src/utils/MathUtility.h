/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <functional>
#include <float.h>
#include <limits>
#include <vector>
#include <numeric>
#include <algorithm>

namespace CppUtils {

/*********************************************************************
* GLOBAL CONSTANTS
*********************************************************************/
constexpr double CPPUTILS_SMALL   = std::numeric_limits<double>::epsilon();
constexpr double CPPUTILS_MAX     = std::numeric_limits<double>::max();
constexpr double CPPUTILS_MIN     = std::numeric_limits<double>::min();
constexpr double CPPUTILS_PI      = 3.14159265358979323846;
constexpr double CPPUTILS_PI_HALF = 0.5 * 3.14159265358979323846;

/*********************************************************************
* USEFUL FUNCTIONS
*
* Sources:
* --------
*  - https://randomascii.wordpress.com/2012/02/25/comparing-\
*    floating-point-numbers-2012-edition/
*********************************************************************/
template <typename T> 
static inline T ABS(T a)
{ return ( (a > T{}) ? a : -a); }

template <typename T>
static inline T MIN(T a, T b) 
{ return a < b ? a : b; }

template <typename T>
static inline T MAX(T a, T b) 
{ return a > b ? a : b; }

template <typename T>
static inline T MOD(T n, T M)
{ return ((n % M) + M) % M;  }

template <typename T>
static inline T CLAMP(T x, T lower, T upper)
{ return MAX(lower, MIN(upper, x)); }

template <typename T>
static inline bool EQ(T a, T b, 
                      T max_diff=std::numeric_limits<T>::min(),
                      T max_rel_diff=std::numeric_limits<T>::epsilon())
{ 
  // Check if the numbers are really close
  // -> for comparing numbers near zero
  const T diff = ABS(a-b);

  if ( diff <= max_diff )
    return true;

  const T a_abs = ABS(a);
  const T b_abs = ABS(b);
  const T largest = MAX(a_abs, b_abs);
  
  if ( diff <= largest * max_rel_diff )
    return true;

  return false; 
}

template <typename T>
static inline bool EQ0(T a,
                       T max_diff=std::numeric_limits<T>::min(),
                       T max_rel_diff=std::numeric_limits<T>::epsilon())
{ return EQ(a, {0}, max_diff, max_rel_diff); }


/*********************************************************************
* Argsort
*
* Reference: 
* ---------
*   https://stackoverflow.com/questions/1577475/c-sorting-\
*   and-keeping-track-of-indexes
*********************************************************************/
template <typename T>
std::vector<size_t> argsort(const std::vector<T> &v)
{
  // Initial index locations
  std::vector<size_t> index( v.size() );
  std::iota(index.begin(), index.end(), 0);

  // Use std::stable_sort() instead of sort to avoid unnecessary
  // index re-orderings when v contains elements of equal values
  std::stable_sort(index.begin(), index.end(), 
    [&v](size_t i1, size_t i2) 
  { 
    return v[i1] < v[i2]; 
  });

  return index;
}


} // namespace CppUtils
