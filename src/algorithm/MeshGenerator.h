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

#include "Domain.h"
#include "Mesh.h"
#include "MeshBuilder.h"


namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* The actual interface to generate meshes
*********************************************************************/
class MeshGenerator
{
  using MeshVector = std::vector<std::unique_ptr<Mesh>>;

public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  MeshGenerator() {}

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  Mesh& mesh(std::size_t i_mesh)
  {
    if ( i_mesh >= meshes_.size() )
      TERMINATE("MeshGenerator::mesh(): Mesh index out of range");
    return *meshes_[i_mesh];
  }

  /*------------------------------------------------------------------
  | Initialize a new mesh for a given domain
  ------------------------------------------------------------------*/
  void define_mesh(Domain& domain,
                   int     mesh_id = DEFAULT_MESH_ID,
                   int     element_color = DEFAULT_ELEMENT_COLOR)
  {
    meshes_.push_back(
      std::make_unique<Mesh>(mesh_id, element_color,
                             domain.vertices().quad_tree().scale(),
                             domain.vertices().quad_tree().max_items(),
                             domain.vertices().quad_tree().max_depth() )  
    );

    Mesh& new_mesh = *( meshes_.back() );
    mesh_builder_.prepare_mesh( new_mesh, domain );
    mesh_builder_.add_mesh_and_domain( new_mesh, domain );
  }

  /*------------------------------------------------------------------
  | Generate the actual mesh elements for all defined meshes
  ------------------------------------------------------------------*/
  bool generate_mesh_elements()
  {
    return false;
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  bool generate_quad_layers()
  {
    return false;
  }


private:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  MeshVector  meshes_ {};
  MeshBuilder mesh_builder_ {};

}; // MeshGenerator

} // namespace TQAlgorithm
} // namespace TQMesh
