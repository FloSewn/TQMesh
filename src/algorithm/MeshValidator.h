/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "Vertex.h"
#include "Edge.h"
#include "Facet.h"
#include "Boundary.h"
#include "Front.h"
#include "Domain.h"
#include "Mesh.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* This class manages the removal of mesh entities
*********************************************************************/
class MeshValidator 
{
public:
  /*------------------------------------------------------------------
  | Constructor / Destructor 
  ------------------------------------------------------------------*/
  MeshValidator(Mesh& mesh, const Domain& domain, Front& front)
  : mesh_ { mesh }, domain_ { domain }, front_ { front } {}

  ~MeshValidator() {}


private:


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Mesh&           mesh_;
  const Domain&   domain_;
  Front&          front_;


}; // MeshValidator

} // namespace TQAlgorithm
} // namespace TQMesh
