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
#include "SmoothingStrategy.h"
#include "TriangulationStrategy.h"
#include "QuadLayerStrategy.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* Different element generation / manipulation algorithms 
*********************************************************************/
enum class MeshingAlgorithm { 
  None,
  Triangulation,
  QuadLayer,
};

enum class SmoothingAlgorithm {
  None,
  Laplace,
  Torsion,
  Mixed,
};


/*********************************************************************
* The actual interface to generate meshes
*********************************************************************/
class MeshGenerator
{
  using MeshVector           = std::vector<std::unique_ptr<Mesh>>;
  using MeshingStrategyPtr   = std::unique_ptr<MeshingStrategy>;
  using SmoothingStrategyPtr = std::unique_ptr<SmoothingStrategy>;

public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  MeshGenerator() {}

  /*------------------------------------------------------------------
  | Getters
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
  | Set a mesh generation algorithm for a specified mesh.
  | Returns false if the mesh is not connected to this MeshGenerator
  | or if the algorithm was not found.
  ------------------------------------------------------------------*/
  bool set_algorithm(Mesh& mesh, MeshingAlgorithm algorithm_type)
  {
    Domain* domain = mesh_builder_.get_domain( mesh );

    if ( !domain ) 
      return false;

    switch (algorithm_type)
    {
      case MeshingAlgorithm::Triangulation:
        meshing_algorithm_ 
          = std::make_unique<TriangulationStrategy>(mesh, *domain);
        meshing_algorithm_type_ = MeshingAlgorithm::Triangulation;
        break;

      case MeshingAlgorithm::QuadLayer:
        meshing_algorithm_ 
          = std::make_unique<QuadLayerStrategy>(mesh, *domain);
        meshing_algorithm_type_ = MeshingAlgorithm::QuadLayer;
        break;

      default:
        meshing_algorithm_type_ = MeshingAlgorithm::None;
    }

    return true;
  }

  /*------------------------------------------------------------------
  | Set a mesh smoothing algorithm for a specified mesh.
  | Returns false if the mesh is not connected to this MeshGenerator
  | or if the algorithm was not found.
  ------------------------------------------------------------------*/
  bool set_algorithm(Mesh& mesh, SmoothingAlgorithm algorithm_type)
  {
    Domain* domain = mesh_builder_.get_domain( mesh );

    if ( !domain ) 
      return false;

    switch (algorithm_type)
    {
      case SmoothingAlgorithm::Laplace:
        smoothing_algorithm_ 
          = std::make_unique<LaplaceSmoothingStrategy>(mesh, *domain);
        smoothing_algorithm_type_ = SmoothingAlgorithm::Laplace;
        break;

      case SmoothingAlgorithm::Torsion:
        smoothing_algorithm_ 
          = std::make_unique<TorsionSmoothingStrategy>(mesh, *domain);
        smoothing_algorithm_type_ = SmoothingAlgorithm::Torsion;
        break;

      case SmoothingAlgorithm::Mixed:
        smoothing_algorithm_ 
          = std::make_unique<MixedSmoothingStrategy>(mesh, *domain);
        smoothing_algorithm_type_ = SmoothingAlgorithm::Mixed;
        break;

      default:
        smoothing_algorithm_type_ = SmoothingAlgorithm::None;
    }

    return true;
  }

  /*------------------------------------------------------------------
  | Use the current front algorithm to generate mesh elements in the
  | associated mesh
  ------------------------------------------------------------------*/
  bool generate_mesh_elements()
  {
    if ( meshing_algorithm_type_ == MeshingAlgorithm::None )
      return false;
    return meshing_algorithm_->generate_elements();
  }

  /*------------------------------------------------------------------
  | Use the current smoothing algorithm to smooth the grid elements
  ------------------------------------------------------------------*/
  bool smooth_mesh_elements(int iterations)
  {
    if ( smoothing_algorithm_type_ == SmoothingAlgorithm::None )
      return false;
    return smoothing_algorithm_->smooth(iterations);
  }



  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  QuadLayerStrategy& quad_layer()
  {
    if ( meshing_algorithm_type_ != MeshingAlgorithm::QuadLayer )
      TERMINATE("MeshGenerator::quad_layer(): "
        "QuadLayer is not chosen as current algorithm.");
    return *dynamic_cast<QuadLayerStrategy*>(meshing_algorithm_.get());
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  TriangulationStrategy& triangulation()
  {
    if ( meshing_algorithm_type_ != MeshingAlgorithm::Triangulation )
      TERMINATE("MeshGenerator::triangulation(): "
        "Triangulation is not chosen as current algorithm.");
    return *dynamic_cast<TriangulationStrategy*>(meshing_algorithm_.get());
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  LaplaceSmoothingStrategy& laplace_smoothing()
  {
    if ( smoothing_algorithm_type_ != SmoothingAlgorithm::Laplace )
      TERMINATE("MeshGenerator::laplace_smoothing(): "
        "Laplace is not chosen as current smoothing algorithm.");
    return *dynamic_cast<LaplaceSmoothingStrategy*>(smoothing_algorithm_.get());
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  TorsionSmoothingStrategy& torsion_smoothing()
  {
    if ( smoothing_algorithm_type_ != SmoothingAlgorithm::Torsion )
      TERMINATE("MeshGenerator::torsion_smoothing(): "
        "Torsion is not chosen as current smoothing algorithm.");
    return *dynamic_cast<TorsionSmoothingStrategy*>(smoothing_algorithm_.get());
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  MixedSmoothingStrategy& mixed_smoothing()
  {
    if ( smoothing_algorithm_type_ != SmoothingAlgorithm::Mixed )
      TERMINATE("MeshGenerator::mixed_smoothing(): "
        "Mixed is not chosen as current smoothing algorithm.");
    return *dynamic_cast<MixedSmoothingStrategy*>(smoothing_algorithm_.get());
  }



private:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  MeshVector           meshes_ {};
  MeshBuilder          mesh_builder_ {};

  MeshingStrategyPtr   meshing_algorithm_;
  MeshingAlgorithm     meshing_algorithm_type_ { MeshingAlgorithm::None };

  SmoothingStrategyPtr smoothing_algorithm_;
  SmoothingAlgorithm   smoothing_algorithm_type_ { SmoothingAlgorithm::None };

}; // MeshGenerator

} // namespace TQAlgorithm
} // namespace TQMesh
