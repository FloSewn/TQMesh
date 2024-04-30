/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "TQMesh.h"

namespace TQMesh {

using namespace CppUtils;

/*********************************************************************
* This class contains the functionality to initialize a set of meshes
* from given set of corresponding domains.
* The initialization comprises the generation of boundary edges, 
* whereas no internal edges, triangles or quads are constructed.
*********************************************************************/
class MeshBuilder
{
public:
  using MeshVector     = std::vector<Mesh*>;
  using DomainVector   = std::vector<Domain*>;
  using VertexVector   = std::vector<Vertex*>;
  using EdgeVector     = std::vector<Edge*>;
  using BoolVector     = std::vector<bool>;
  using IntVector      = std::vector<int>;

  /*------------------------------------------------------------------
  | Constructor / Destructor
  ------------------------------------------------------------------*/
  MeshBuilder() = default;
  virtual ~MeshBuilder() {}

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  size_t n_meshes() const { return meshes_.size(); }
  size_t n_domains() const { return domains_.size(); }

  const MeshVector& meshes() const { return meshes_; }
  MeshVector& meshes() { return meshes_; }

  const DomainVector& domains() const { return domains_; }
  DomainVector& domains() { return domains_; }

  /*------------------------------------------------------------------
  | Get the domain that corresponds to a given mesh
  ------------------------------------------------------------------*/
  Domain* get_domain(Mesh& mesh)
  {
    for ( std::size_t i = 0; i < meshes_.size(); ++i )
      if ( meshes_[i] == &mesh ) 
        return domains_[i];
    return nullptr;
  }

  /*------------------------------------------------------------------
  | Create a new empty mesh entity, based on the extent of a given
  | domain
  ------------------------------------------------------------------*/
  static inline Mesh
  create_empty_mesh(Domain& domain, 
                    int mesh_id=DEFAULT_MESH_ID,
                    int element_color=DEFAULT_ELEMENT_COLOR,
                    double domain_extent=0.0)
  { 
    // Obtain domain extents 
    if ( EQ0(domain_extent, 0.0) )
    {
      Vec2d extrema { DBL_MAX, -DBL_MAX };

      for ( const auto& v_ptr : domain.vertices() )
      {
        const Vec2d& v_xy = v_ptr->xy();
        domain_extent = MAX(domain_extent, ABS(v_xy.x));
        domain_extent = MAX(domain_extent, ABS(v_xy.y));
      }

      // Double extent and enlarge slightly
      domain_extent *= domain_enlargement_;
    }

    return { mesh_id, element_color, domain_extent,
             domain.vertices().quad_tree().max_items(),
             domain.vertices().quad_tree().max_depth() }; 

  } // create_empty_mesh()

  /*------------------------------------------------------------------
  | Create a new unique_ptr to an empty mesh entity, 
  | based on the extent of a given domain
  ------------------------------------------------------------------*/
  static inline std::unique_ptr<Mesh>
  create_empty_mesh_ptr(Domain& domain, 
                        int mesh_id=DEFAULT_MESH_ID,
                        int element_color=DEFAULT_ELEMENT_COLOR)
  {
    return std::move(
      std::make_unique<Mesh>(mesh_id, element_color,
                             domain.vertices().quad_tree().scale(),
                             domain.vertices().quad_tree().max_items(),
                             domain.vertices().quad_tree().max_depth() )  
    );

  } // create_empty_mesh_ptr()

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

    return true;

  } // remove_mesh_and_domain()

  /*------------------------------------------------------------------
  | Prepare a given mesh entity for the advancing front triangulation.
  | The given mesh entity must be empty.
  | 
  | After the function call, it will consists of boundary edges
  | and the initial advancing front. The latter is stored indirectly
  | in terms of interior mesh edges, which do not feature any
  | neighboring facets. 
  | The mesh does not contain any elements.
  ------------------------------------------------------------------*/
  bool prepare_mesh(Mesh& mesh, Domain& domain)
  {
    // We require an empty mesh
    if ( !mesh.is_empty() )
    {
      LOG(ERROR) << "Failed mesh preparation: Mesh is not empty.";
      return false;
    }

    // Check if the domain is valid for the mesh generation
    if ( !EntityChecks::check_domain_validity(domain) )
    {
      LOG(ERROR) << "Failed mesh preparation: Invalid domain.";
      return false;
    }

    // Get all edges that will define the advancing front
    // Here, we also consider other meshes that have been defined with 
    // this mesh builder, in order to maintain vertex-adjacency
    FrontInitData front_data { domain, meshes_ };

    // Initialize edges in a temporary advancing front
    Front tmp_front {};
    tmp_front.init_front(domain, front_data, mesh.vertices());

    // Check if the initial advancing front is valid
    if ( !EntityChecks::check_front_validity( tmp_front ) )
    {
      LOG(ERROR) << "Failed mesh preparation: Invalid advancing front.";

      // Remove all previously generated vertices
      for ( auto& v_ptr : mesh.vertices() )
        mesh.remove_vertex( *v_ptr );
      mesh.clear_waste();
      
      ASSERT( mesh.is_empty(), 
        "Failed to get mesh back into initial state.");

      return false;
    }

    // Setup the mesh's boundary edges from the initial advancing front
    for ( auto& e : tmp_front )
    {
      Vertex& v1 = e->v1();
      Vertex& v2 = e->v2();
      int marker = e->marker();
      Edge& e_new = mesh.boundary_edges().add_edge( v1, v2, marker );

      ASSERT( v1.has_property( VertexProperty::on_boundary ),
        "MeshBuilder::prepare_mesh(): Missing "
        "vertex property \"on_boundary\".");
      ASSERT( v2.has_property( VertexProperty::on_boundary ),
        "MeshBuilder::prepare_mesh(): Missing "
        "vertex property \"on_boundary\".");

      // Connect boundary edges of this mesh and its parner mesh
      Edge* e_twin = e->twin_edge();

      if ( e_twin )
      {
        e_twin->twin_edge( &e_new );
        e_new.twin_edge( e_twin );
        e->twin_edge( nullptr );
      }
    }

    return true;

  } // prepare_mesh()

private:

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  MeshVector   meshes_  {};
  DomainVector domains_ {};

  // Constant factor to double the domain extent and to slightly 
  // enlarge it
  static constexpr double domain_enlargement_ { 2.1 };

}; // MeshBuilder
 

} // namespace TQMesh
