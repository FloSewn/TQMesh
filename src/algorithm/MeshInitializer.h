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
* This class contains the functionality to initialize a set of meshes
* from given set of corresponding domains.
* The initialization comprises the generation of boundary edges, 
* whereas no internal edges, triangles or quads are constructed.
*********************************************************************/
class MeshInitializer
{
public:
  using MeshVector     = std::vector<Mesh*>;
  using DomainVector   = std::vector<Domain*>;
  using EdgeVector     = std::vector<Edge*>;
  using BoolVector     = std::vector<bool>;
  using IntVector      = std::vector<int>;

  /*------------------------------------------------------------------
  | Constructor / Destructor
  ------------------------------------------------------------------*/
  MeshInitializer() = default;
  virtual ~MeshInitializer() {}

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
  { 
    return { mesh_id, element_color, 
             domain.vertices().quad_tree().scale(),
             domain.vertices().quad_tree().max_items(),
             domain.vertices().quad_tree().max_depth() }; 

  } // create_empty_mesh()

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


  /*------------------------------------------------------------------
  | Check if the provided domain is valid for the mesh generation
  ------------------------------------------------------------------*/
  static inline bool check_domain_validity(const Domain& domain) 
  {
    // Check if boundaries are defined 
    if ( domain.size() < 1 )
    {
      LOG(ERROR) << 
      "Can not initialize the mesh's advancing front, "
      "since no domain boundaries are defined.";
      return false;
    }

    // Check if an exterior boundary exists
    bool extr_bdry_found = false;

    for ( const auto& boundary : domain )
    {
      if ( boundary->is_exterior() )
      {
        extr_bdry_found = true;
        break;
      }
    }

    if ( !extr_bdry_found )
    {
      LOG(ERROR) << 
      "Can not initialize advancing front. "
      "No exterior mesh boundary found.";
      return false;
    }

    // Check if domain boundaries are traversable
    for ( const auto& boundary : domain )
    {
      Edge& e_start = boundary->edges()[0];

      if ( !boundary->is_traversable(e_start, e_start) )
      {
        LOG(ERROR) << 
        "Can not initialize advancing front. "
        "Mesh boundaries are not traversable.";
        return false;
      }
    }

    return true;

  } // check_domain_validity()


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
    if ( !mesh.is_empty() )
    {
      LOG(ERROR) << "Failed mesh preparation: Mesh not empty.";
      return false;
    }

    if ( !MeshInitializer::check_domain_validity(domain) )
    {
      LOG(ERROR) << "Failed mesh preparation: Invalid domain.";
      return false;
    }

    // Count the number of boundary edge overlaps between the current
    // domain and all domains that are already treated by the meshing 
    // algorithm structure
    size_t n_overlaps = 0;
    for (size_t i_domain = 0; i_domain < n_domains(); ++i_domain)
      n_overlaps += domain.get_overlaps( *domains_[i_domain] );

    // Get all edges that will define the advancing front
    FrontInitData front_data = collect_front_edges(domain);

    // Initialize edges in a temporary advancing front
    Front tmp_front {};
    tmp_front.init_front_edges(domain, front_data, mesh.vertices());

    // Setup the mesh's boundary edges from the initial advancing front
    for ( auto& e : tmp_front )
    {
      Vertex& v1 = e->v1();
      Vertex& v2 = e->v2();
      int marker = e->marker();
      Edge& e_new = mesh.boundary_edges().add_edge( v1, v2, marker );

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
  | For a given domain boundary edge <e>, 
  | search all boundary edges of a 
  | neighboring mesh, that are contained within <e>.
  | This function is used upon the initialization of the advancing
  | front.
  ------------------------------------------------------------------*/
  EdgeVector get_neighbor_mesh_edges(const Edge& e) const
  {
    EdgeVector   nbr_edges {};

    const Vec2d& c = e.xy();
    double       r = CONSTANTS.edge_search_factor() * e.length();

    const Vec2d& v = e.v1().xy();
    const Vec2d& w = e.v2().xy();

    for ( Mesh* m : meshes_ )
    {
      ASSERT( m, "Neighbor mesh data structure is invalid." );

      // Get edges in the vicinity of the given edge
      EdgeVector edges_in_vicinity = m->get_bdry_edges(c, r);

      // Get all edges, that are actually located on the segment
      // of the current domain boundary edge 
      for ( Edge* e_vicinity : edges_in_vicinity )
      {
        const Vec2d& p = e_vicinity->v1().xy();
        const Vec2d& q = e_vicinity->v2().xy();

        if ( in_on_segment(v,w,p) && in_on_segment(v,w,q) )
          nbr_edges.push_back( e_vicinity );
      }

      // In case edges have been found, sort them in 
      // ascending direction to the stating vertex
      std::sort(nbr_edges.begin(), nbr_edges.end(), 
      [v](Edge* e1, Edge* e2)
      {
        double delta_1 = (e1->xy() - v).norm_sqr();
        double delta_2 = (e2->xy() - v).norm_sqr();
        return (delta_1 < delta_2);
      });

      // If edges have been located, terminate the search
      // --> First come, first serve
      if ( nbr_edges.size() > 0 )
        break;
    }

    return std::move(nbr_edges);
  }


  /*------------------------------------------------------------------
  | Collect all boundary / neighbor-mesh edges that will be used 
  | for the creation of the advancing front
  | Additionally, return also an array of booleans, which indicate
  | the orientation of the given edges.
  | Edges which result from neighbor mesh are oriented in the 
  | opposite direction. This can be used to identify these edges
  | from normal boundary edges (e.g. for the advancing front 
  | refinement).
  ------------------------------------------------------------------*/
  FrontInitData collect_front_edges(Domain& domain) const
  {
    FrontInitData front_data {};

    for ( const auto& boundary : domain )
    {
      EdgeVector edges {};
      IntVector  markers {};
      BoolVector is_oriented {};

      for ( const auto& e : boundary->edges() )
      {
        EdgeVector nbr_edges = get_neighbor_mesh_edges(*e);

        if ( nbr_edges.size() > 0 )
        {
          for ( Edge* nbr_e : nbr_edges )
          {
            edges.push_back( nbr_e );
            is_oriented.push_back( false );
            markers.push_back( nbr_e->marker() );
          }
        }
        else
        {
          edges.push_back( e.get() ) ;
          is_oriented.push_back( true );
          markers.push_back( e->marker() );
        }
      }

      front_data.edges.push_back( edges );
      front_data.is_oriented.push_back( is_oriented );
      front_data.markers.push_back( markers );
    }

    return std::move(front_data); 

  } // Mesh::collect_front_edges()


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  MeshVector   meshes_  {};
  DomainVector domains_ {};


}; // MeshInitializer
 

} // namespace TQAlgorithm
} // namespace TQMesh