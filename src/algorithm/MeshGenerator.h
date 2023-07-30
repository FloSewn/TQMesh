/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <algorithm>
#include <memory>
#include <limits.h>

#include "VecND.h"

#include "Domain.h"
#include "Mesh.h"
#include "MeshBuilder.h"
#include "MeshingStrategy.h"
#include "TriangulationStrategy.h"
#include "QuadLayerStrategy.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* Different element generation algorithms that utlize an advancing
* front approach
*********************************************************************/
enum class Algorithm { 
  None,
  Triangulation,
  QuadLayer,
};


/*********************************************************************
* The actual interface to generate meshes
*********************************************************************/
class MeshGenerator
{
  using MeshVector         = std::vector<std::unique_ptr<Mesh>>;
  using MeshingStrategyPtr = std::unique_ptr<MeshingStrategy>;

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
  Mesh& new_mesh(Domain& domain,
                 int     mesh_id = DEFAULT_MESH_ID,
                 int     element_color = DEFAULT_ELEMENT_COLOR)
  {
    meshes_.push_back( 
      mesh_builder_.create_empty_mesh_ptr(domain, mesh_id, element_color)
    );

    Mesh& new_mesh = *( meshes_.back() );
    mesh_builder_.prepare_mesh( new_mesh, domain );
    mesh_builder_.add_mesh_and_domain( new_mesh, domain );
    
    return new_mesh;
  }

  /*------------------------------------------------------------------
  | Set an advancing front algorithm for a specified mesh.
  | Returns false if the mesh is not connected to this MeshGenerator
  | or if the algorithm was not found.
  ------------------------------------------------------------------*/
  bool set_meshing_algorithm(Mesh& mesh, Algorithm algorithm_type)
  {
    Domain* domain = mesh_builder_.get_domain( mesh );

    if ( !domain ) 
      return false;

    switch (algorithm_type)
    {
      case Algorithm::Triangulation:
        front_algorithm_ 
          = std::make_unique<TriangulationStrategy>(mesh, *domain);
        algorithm_ = Algorithm::Triangulation;
        break;

      case Algorithm::QuadLayer:
        front_algorithm_ 
          = std::make_unique<QuadLayerStrategy>(mesh, *domain);
        algorithm_ = Algorithm::QuadLayer;
        break;

      default:
        algorithm_ = Algorithm::None;
    }

    return true;
  }

  /*------------------------------------------------------------------
  | Use the current front algorithm to generate mesh elements in the
  | associated mesh
  ------------------------------------------------------------------*/
  bool generate_mesh_elements()
  {
    if ( algorithm_ == Algorithm::None )
      return false;
    return front_algorithm_->generate_elements();
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  QuadLayerStrategy& quad_layer()
  {
    if ( algorithm_ != Algorithm::QuadLayer )
      TERMINATE("MeshGenerator::quad_layer(): "
        "QuadLayer is not chosen as current algorithm.");
    return *dynamic_cast<QuadLayerStrategy*>(front_algorithm_.get());
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  TriangulationStrategy& triangulation()
  {
    if ( algorithm_ != Algorithm::Triangulation )
      TERMINATE("MeshGenerator::triangulation(): "
        "Triangulation is not chosen as current algorithm.");
    return *dynamic_cast<TriangulationStrategy*>(front_algorithm_.get());
  }


private:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  MeshVector         meshes_ {};
  MeshBuilder        mesh_builder_ {};

  MeshingStrategyPtr front_algorithm_;
  Algorithm          algorithm_ { Algorithm::None };

}; // MeshGenerator

} // namespace TQAlgorithm
} // namespace TQMesh
