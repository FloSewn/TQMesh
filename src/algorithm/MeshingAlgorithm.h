/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <vector>

#include "Domain.h"
#include "Mesh.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* 
*********************************************************************/
class MeshingAlgorithm
{
public:
  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  MeshingAlgorithm() = default;
  virtual ~MeshingAlgorithm() {}

  virtual bool generate_elements(Mesh& mesh,  const Domain& domain,
                                 int n_elements=0) = 0;

private:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/

}; // MeshingAlgorithm
 

} // namespace TQAlgorithm
} // namespace TQMesh
