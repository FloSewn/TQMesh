/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <functional>
#include <cassert>
#include <float.h>
#include <limits.h>

namespace TQMesh {
namespace TQUtils {

/*********************************************************************
* Asserts messages
*********************************************************************/
#ifndef NDEBUG
#define ASSERT(cond, str)                               \
  do { if (!cond) {                                     \
        std::cerr << "# [ERROR] " << str << std::endl;  \
        assert(cond);                                   \
        } else {  } } while(false)
#else
#define ASSERT(cond, str) \
  do { } while (false)
#endif

/*********************************************************************
* Debug messages
*********************************************************************/
#ifndef NDEBUG
#define DBG_MSG(str) \
  do { std::clog << "# " << str << std::endl; } while(false)
#else
#define DBG_MSG(str) \
  do { } while (false)
#endif

/*********************************************************************
* Output messages
*********************************************************************/
#define MSG(str) \
  do { std::clog << "# " << str << std::endl; } while(false)

/*********************************************************************
* DEFAULT VALUES
*********************************************************************/
// *** Size function parameters ***
#define TQ_DOMAIN_MIN_SCALING  (0.001)
#define TQ_DOMAIN_MIN_ELEMSIZE (0.001)
// *** QuadTree parameters ***
#define TQ_QTREE_SCALE         (10000.0)
#define TQ_QTREE_ITEMS         (100)
#define TQ_QTREE_DEPTH         (25)
// *** Marker for advancing front edges ***
#define TQ_INTR_EDGE_MARKER    (-1)


#define TQ_RANGE_FACTOR        (1.0)
#define TQ_WIDE_SEARCH_FACTOR  (10.0)
#define TQ_QUAD_RANGE_FACTOR   (0.50)
#define TQ_QUAD_WEDGE_ANGLE    (2.35619449) // = 3/4 pi
#define TQ_QUAD_LAYER_ANGLE    (1.57079633) // = 1/2 pi
#define TQ_QUAD_LAYER_RANGE    (0.75)

/*********************************************************************
* GLOBAL CONSTANTS
*********************************************************************/
constexpr double TQ_SMALL  = 1.0E-13; //DBL_EPSILON;
constexpr double TQ_MAX    = DBL_MAX;
constexpr double TQ_MIN    = DBL_MIN;

/*********************************************************************
* USEFUL FUNCTIONS
*********************************************************************/
template <typename T> 
static inline T ABS(T a)
{ return ( (a > 0) ? a : -a); }

template <typename T>
static inline bool EQ(T a, T b)
{ return ( ABS(a-b) < TQ_SMALL ); }

template <typename T>
static inline bool EQ0(T a)
{ return ( ABS(a) < TQ_SMALL ); }

template <typename T>
static inline T MIN(T a, T b) 
{ return a < b ? a : b; }

template <typename T>
static inline T MAX(T a, T b) 
{ return a > b ? a : b; }

template <typename T>
static inline T MOD(T n, T M)
{ return ((n % M) + M) % M;  }


} // namespace TQUtils
} // namespace TQMesh 
