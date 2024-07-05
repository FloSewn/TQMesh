/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <TQMeshConfig.h>
#include "STLHeaders.h"
#include "CppUtils.h"

namespace TQMesh {

/*********************************************************************
* CONSTANTS
*********************************************************************/
constexpr int    INTERIOR_EDGE_COLOR   = -1;
constexpr int    DEFAULT_EDGE_COLOR    =  0;
constexpr int    DEFAULT_ELEMENT_COLOR =  0;
constexpr int    DEFAULT_MESH_ID       =  0;
constexpr double TQ_SMALL              = 1.0E-13;
constexpr double TQ_MAX                = DBL_MAX;
constexpr double TQ_MIN                = DBL_MIN;

} // namespace TQMesh 
