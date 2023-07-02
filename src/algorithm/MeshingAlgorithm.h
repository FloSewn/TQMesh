/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <vector>

#include "utils.h"

#include "Domain.h"
#include "Mesh.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* Base class for different meshing algorithm implementations 
*********************************************************************/
class MeshingAlgorithm
{
public:
  using MeshVector     = std::vector<Mesh*>;
  using DomainVector   = std::vector<Domain*>;

  /*------------------------------------------------------------------
  | Constructor / Destructor
  ------------------------------------------------------------------*/
  MeshingAlgorithm() = default;
  virtual ~MeshingAlgorithm() {}

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  size_t n_meshes() const { return meshes_.size(); }
  size_t n_domains() const { return domains_.size(); }

  /*------------------------------------------------------------------
  | Create a new empty mesh entity, based on the extent of a given
  | domain
  ------------------------------------------------------------------*/
  static inline Mesh
  create_empty_mesh(Domain& domain, 
                    int mesh_id=CONSTANTS.default_mesh_id(),
                    int element_color=CONSTANTS.default_element_color())
  { return { mesh_id, element_color, 
             domain.vertices().quad_tree().scale(),
             domain.vertices().quad_tree().max_items(),
             domain.vertices().quad_tree().max_depth() }; }

  /*------------------------------------------------------------------
  | Add meshes and corresponding domains
  ------------------------------------------------------------------*/
  void add_mesh_and_domain(Mesh& mesh, Domain& domain)
  {
    ASSERT( n_meshes() == n_domains(), 
      "MeshingAlgorithm: Invalid mesh-domain structure.");

    meshes_.push_back( &mesh );
    domains_.push_back( &domain );

  } // add_mesh_and_domain()

  /*------------------------------------------------------------------
  | Remove a mesh and its corresponding domain 
  ------------------------------------------------------------------*/
  bool remove_mesh_and_domain(Mesh& mesh)
  {
    ASSERT( n_meshes() == n_domains(), 
      "MeshingAlgorithm: Invalid mesh-domain structure.");

    auto it = std::find(meshes_.begin(), meshes_.end(), &mesh);

    if (it == meshes_.end())
      return false;

    size_t index = std::distance(meshes_.begin(), it);

    domains_.erase(domains_.begin() + index);
    meshes_.erase(meshes_.begin() + index);

  } // remove_mesh_and_domain()

protected:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  MeshVector   meshes_  {};
  DomainVector domains_ {};


}; // MeshingAlgorithm
 

} // namespace TQAlgorithm
} // namespace TQMesh
