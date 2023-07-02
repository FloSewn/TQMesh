/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <algorithm>
#include <limits.h>

#include "VecND.h"
#include "utils.h"
#include "VtkIO.h"
#include "ProgressBar.h"

#include "Vertex.h"
#include "Triangle.h"
#include "Quad.h"
#include "Front.h"
#include "Boundary.h"
#include "Domain.h"
#include "QuadLayer.h"
#include "Mesh.h"


namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* The actual interface to generate meshes
*********************************************************************/
class MeshGenerator
{

public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  MeshGenerator() {}

private:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/

};

} // namespace TQAlgorithm
} // namespace TQMesh
