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
#include "MeshWriter.h"
#include "MeshMerger.h"
#include "MeshingStrategy.h"
#include "SmoothingStrategy.h"
#include "RefinementStrategy.h"
#include "ModificationStrategy.h"
#include "TriangulationStrategy.h"
#include "QuadLayerStrategy.h"

namespace TQMesh {

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

enum class RefinementAlgorithm {
  None,
  Quad,
};

enum class ModificationAlgorithm {
  None,
  Tri2Quad,
  QMorph,
};


/*********************************************************************
* The actual interface to generate meshes
*********************************************************************/
class MeshGenerator
{
  using MeshVector              = std::vector<std::unique_ptr<Mesh>>;
  using MeshingStrategyPtr      = std::unique_ptr<MeshingStrategy>;
  using SmoothingStrategyPtr    = std::unique_ptr<SmoothingStrategy>;
  using RefinementStrategyPtr   = std::unique_ptr<RefinementStrategy>;
  using ModificationStrategyPtr = std::unique_ptr<ModificationStrategy>;

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
  
  bool is_valid(Mesh& mesh)
  { return ( mesh_builder_.get_domain(mesh) != nullptr ); }

  std::size_t size() const { return meshes_.size(); }


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
  | Merge all defined meshes 
  ------------------------------------------------------------------*/
  bool merge_meshes(Mesh& receiver, Mesh& donor)
  {
    if ( !mesh_builder_.get_domain( receiver ) || 
         !mesh_builder_.get_domain( donor ) )
      return false;

    MeshMerger merger ( receiver, donor );
    bool success = merger.merge();

    if ( success ) 
    {
      auto it = std::find_if(meshes_.begin(), meshes_.end(), 
        [&donor](const std::unique_ptr<Mesh>& ptr)
        { return ptr.get() ==&donor; }
      );

      if (it == meshes_.end())
        return false;

      if ( !mesh_builder_.remove_mesh_and_domain(donor) )
        return false;

      meshes_.erase( it );
    }

    return true;

  } // MeshGenerator::merge_meshes()

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  bool write_mesh(Mesh& mesh, const std::string& filename,
                  MeshExportType export_type)
  {
    Domain* domain = mesh_builder_.get_domain( mesh );

    if ( !domain )
      return false;

    MeshWriter writer { mesh, *domain };

    return writer.write(filename, export_type);
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  QuadLayerStrategy& quad_layer_generation(Mesh& mesh)
  {
    if ( meshing_algorithm_type_ != MeshingAlgorithm::QuadLayer ||
         &meshing_algorithm_->mesh() != &mesh )
      if ( !set_algorithm(mesh, MeshingAlgorithm::QuadLayer) )
        TERMINATE("MeshGenerator::quad_layer_generation(): Invalid mesh provided.");

    return *dynamic_cast<QuadLayerStrategy*>(meshing_algorithm_.get());
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  TriangulationStrategy& triangulation(Mesh& mesh)
  {
    if ( meshing_algorithm_type_ != MeshingAlgorithm::Triangulation ||
         &meshing_algorithm_->mesh() != &mesh )
      if ( !set_algorithm(mesh, MeshingAlgorithm::Triangulation) )
        TERMINATE("MeshGenerator::triangulation(): Invalid mesh provided.");

    return *dynamic_cast<TriangulationStrategy*>(meshing_algorithm_.get());
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  LaplaceSmoothingStrategy& laplace_smoothing(Mesh& mesh)
  {
    if ( smoothing_algorithm_type_ != SmoothingAlgorithm::Laplace ||
         &smoothing_algorithm_->mesh() != &mesh )
      if ( !set_algorithm(mesh, SmoothingAlgorithm::Laplace) )
        TERMINATE("MeshGenerator::laplace_smoothing(): Invalid mesh provided.");

    return *dynamic_cast<LaplaceSmoothingStrategy*>(smoothing_algorithm_.get());
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  TorsionSmoothingStrategy& torsion_smoothing(Mesh& mesh)
  {
    if ( smoothing_algorithm_type_ != SmoothingAlgorithm::Torsion ||
         &smoothing_algorithm_->mesh() != &mesh )
      if ( !set_algorithm(mesh, SmoothingAlgorithm::Torsion) )
        TERMINATE("MeshGenerator::torsion_smoothing(): Invalid mesh provided.");

    return *dynamic_cast<TorsionSmoothingStrategy*>(smoothing_algorithm_.get());
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  MixedSmoothingStrategy& mixed_smoothing(Mesh& mesh)
  {
    if ( smoothing_algorithm_type_ != SmoothingAlgorithm::Mixed ||
         &smoothing_algorithm_->mesh() != &mesh )
      if ( !set_algorithm(mesh, SmoothingAlgorithm::Mixed) )
        TERMINATE("MeshGenerator::mixed_smoothing(): Invalid mesh provided.");

    return *dynamic_cast<MixedSmoothingStrategy*>(smoothing_algorithm_.get());
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  QuadRefinementStrategy& quad_refinement(Mesh& mesh)
  {
    if ( refinement_algorithm_type_ != RefinementAlgorithm::Quad ||
         &refinement_algorithm_->mesh() != &mesh )
      if ( !set_algorithm(mesh, RefinementAlgorithm::Quad) )
        TERMINATE("MeshGenerator::quad_refinement(): Invalid mesh provided.");

    return *dynamic_cast<QuadRefinementStrategy*>(refinement_algorithm_.get());
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  Tri2QuadStrategy& tri2quad_modification(Mesh& mesh)
  {
    if ( modification_algorithm_type_ != ModificationAlgorithm::Tri2Quad ||
         &modification_algorithm_->mesh() != &mesh )
      if ( !set_algorithm(mesh, ModificationAlgorithm::Tri2Quad) )
        TERMINATE("MeshGenerator::tri2quad_modification(): Invalid mesh provided.");

    return *dynamic_cast<Tri2QuadStrategy*>(modification_algorithm_.get());
  }


private:

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
  | Set a mesh refinement algorithm for a specified mesh.
  | Returns false if the mesh is not connected to this MeshGenerator
  | or if the algorithm was not found.
  ------------------------------------------------------------------*/
  bool set_algorithm(Mesh& mesh, RefinementAlgorithm algorithm_type)
  {
    Domain* domain = mesh_builder_.get_domain( mesh );

    if ( !domain ) 
      return false;

    switch (algorithm_type)
    {
      case RefinementAlgorithm::Quad:
        refinement_algorithm_ 
          = std::make_unique<QuadRefinementStrategy>(mesh, *domain);
        refinement_algorithm_type_ = RefinementAlgorithm::Quad;
        break;

      default:
        refinement_algorithm_type_ = RefinementAlgorithm::None;
    }

    return true;
  }

  /*------------------------------------------------------------------
  | Set a mesh modification algorithm for a specified mesh.
  | Returns false if the mesh is not connected to this MeshGenerator
  | or if the algorithm was not found.
  ------------------------------------------------------------------*/
  bool set_algorithm(Mesh& mesh, ModificationAlgorithm algorithm_type)
  {
    Domain* domain = mesh_builder_.get_domain( mesh );

    if ( !domain ) 
      return false;

    switch (algorithm_type)
    {
      case ModificationAlgorithm::Tri2Quad:
        modification_algorithm_ 
          = std::make_unique<Tri2QuadStrategy>(mesh, *domain);
        modification_algorithm_type_ = ModificationAlgorithm::Tri2Quad;
        break;

      default:
        modification_algorithm_type_ = ModificationAlgorithm::None;
    }

    return true;
  }



  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  MeshVector              meshes_ {};
  MeshBuilder             mesh_builder_ {};

  MeshingStrategyPtr      meshing_algorithm_;
  MeshingAlgorithm        meshing_algorithm_type_ { MeshingAlgorithm::None };

  SmoothingStrategyPtr    smoothing_algorithm_;
  SmoothingAlgorithm      smoothing_algorithm_type_ { SmoothingAlgorithm::None };

  RefinementStrategyPtr   refinement_algorithm_;
  RefinementAlgorithm     refinement_algorithm_type_ { RefinementAlgorithm::None };

  ModificationStrategyPtr modification_algorithm_;
  ModificationAlgorithm   modification_algorithm_type_ { ModificationAlgorithm::None };

}; // MeshGenerator

} // namespace TQMesh
