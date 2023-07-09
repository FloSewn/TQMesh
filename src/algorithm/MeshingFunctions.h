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
#include "Edge.h"
#include "Facet.h"
#include "Boundary.h"
#include "Front.h"
#include "Domain.h"
#include "Mesh.h"
#include "MeshValidator.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* 
*********************************************************************/
class MeshingFunctions
{
public:
  using VertexVector   = std::vector<Vertex*>;
  using TriVector      = std::vector<Triangle*>;


}; // MeshingFunctions

} // namespace TQAlgorithm
} // namespace TQMesh
