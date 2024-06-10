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

#include "Domain.h"
#include "Mesh.h"
#include "MeshCleanup.h"
#include "EntityChecks.h"

namespace TQMesh {

using namespace CppUtils;

/*********************************************************************
* This class can be used to check a mesh for its completeness and 
* for certain quality measures
*********************************************************************/
class MeshChecker
{
public:

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  MeshChecker(Mesh& mesh, const Domain& domain)
  : mesh_   { &mesh }
  , domain_ { &domain }
  {}

  virtual ~MeshChecker() {}

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  Mesh& mesh() { return *mesh_; }

  /*------------------------------------------------------------------
  | This function checks for the completeness of the given mesh
  ------------------------------------------------------------------*/
  bool check_completeness(bool mesh_cleanup=true) const
  {
    // Make sure that the mesh structure is up to date
    if ( mesh_cleanup )
      perform_mesh_cleanup();

    bool complete = true;

    // The mesh may not contain any advancing front edges
    complete &= ( mesh_->get_front_edges().size() == 0 );

    // Check the mesh connectivity
    complete &= EntityChecks::check_mesh_validity( *mesh_ );

    // Check meshed area 
    const double mesh_area = mesh_->area();
    const double domain_area = domain_->area();
    auto extent = domain_->extent();
    const double scale_x = extent.first[1] - extent.first[0];
    const double scale_y = extent.second[1] - extent.second[0];
    const double area_difference = ABS(mesh_area - domain_area);
    complete &= (area_difference / (scale_x*scale_y)) < area_threshold_;

    return complete;
  }

private:
  /*------------------------------------------------------------------
  | This function calls some MeshCleanup functions which might be
  | required to update the internal mesh structure
  ------------------------------------------------------------------*/
  void perform_mesh_cleanup() const
  {
    MeshCleanup::assign_mesh_indices( *mesh_ );
    MeshCleanup::setup_facet_connectivity( *mesh_ );
  }

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Mesh*         mesh_;
  const Domain* domain_;

  const double  area_threshold_ { 1.0E-10 };

}; // MeshChecker

} // namespace TQMesh
