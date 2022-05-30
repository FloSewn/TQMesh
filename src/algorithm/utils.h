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

#include "Helpers.h"

namespace TQMesh {
namespace TQAlgorithm {


/*********************************************************************
* Output messages
*********************************************************************/
static CppUtils::SimpleLogger TQLogger { std::clog, "# " };
#define MSG(str) \
  do { TQMesh::TQAlgorithm::TQLogger << str << std::endl; } while(false)

#ifndef NDEBUG
#define DBG_MSG(str) \
  MSG(str) 
#else
#define DBG_MSG(str) \
  do { } while (false)
#endif


/*********************************************************************
* GLOBALS
*********************************************************************/
// *** Marker for interior mesh edges ***
//define TQ_INTR_EDGE_MARKER    (-1)

/*********************************************************************
* CONSTANTS
*********************************************************************/
constexpr double TQ_SMALL  = 1.0E-13; //DBL_EPSILON;
constexpr double TQ_MAX    = DBL_MAX;
constexpr double TQ_MIN    = DBL_MIN;

constexpr double TQMeshMinimumElementSize = 0.001;
constexpr double TQMeshMinimumElementScaling = 0.001;

constexpr int TQMeshInteriorEdgeMarker = -1;

constexpr double TQMeshQuadLayerAngle = 1.57079633; // = 1/2 pi
constexpr double TQMeshQuadLayerRange = 0.75;

constexpr double TQMeshQuadRangeFactor = 0.50;
constexpr double TQMeshRangeFactor = 1.0;

constexpr double TQMeshWideSearchFactor = 10.0;

} // namespace TQAlgorithm
} // namespace TQMesh 
