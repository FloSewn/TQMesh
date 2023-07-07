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
#include "utils.h"
#include "VtkIO.h"
#include "ProgressBar.h"

#include "Vertex.h"
#include "Triangle.h"
#include "Quad.h"
#include "NullFacet.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* This class contains all functions that are required to cleanup and
* prepare the mesh after the generation process
*********************************************************************/
class Cleanup
{
public:

  /*------------------------------------------------------------------
  | Check the facet-vertex-edge connectivtiy of a given mesh
  ------------------------------------------------------------------*/
  template <typename Mesh>
  static inline bool check_mesh_validity(Mesh& mesh)
  { 
    // Check connectivity for interior edges
    for ( const auto& e_ptr : mesh.interior_edges() )
    {
      const Vertex& v1 = e_ptr->v1();
      const Vertex& v2 = e_ptr->v2();

      bool check_1 = false; bool check_2 = false;

      // Traverse adjacents facets of both vertices 
      // and search for the current interior edge
      for ( auto f : v1.facets() )
      {
        check_1 = ( f->get_edge_index(v1, v2) >= 0 );
        if (check_1) break;
      }

      for ( auto f : v2.facets() )
      {
        check_2 = ( f->get_edge_index(v1, v2) >= 0 );
        if (check_2) break;
      }

      // Edge was not found 
      if (!check_1 || !check_2) 
        return false;
    }

    // Check connectivity for boundary edges
    for ( const auto& e_ptr : mesh.boundary_edges() )
    {
      const Vertex& v1 = e_ptr->v1();
      const Vertex& v2 = e_ptr->v2();

      bool check_1 = false; bool check_2 = false;

      // Traverse adjacents facets of both vertices 
      // and search for the current interior edge
      for ( auto f : v1.facets() )
      {
        check_1 = ( f->get_edge_index(v1, v2) >= 0 );
        if (check_1) break;
      }

      for ( auto f : v2.facets() )
      {
        check_2 = ( f->get_edge_index(v1, v2) >= 0 );
        if (check_2) break;
      }

      // Edge was not found 
      if (!check_1 || !check_2) 
        return false;
    }

    return true;

  } // Cleanup::check_mesh_validity()


  /*------------------------------------------------------------------
  | This function assigns a corresponding global index to the 
  | following entities of a given mesh:
  | - Vertices
  | - Quads, Triangles (-> Quads are stored prior to triangles)
  ------------------------------------------------------------------*/
  template <typename Mesh>
  static inline void assign_mesh_indices(Mesh& mesh)
  {
    unsigned int vertex_index = 0; 
    unsigned int element_index = 0;

    // Vertices 
    for ( auto& v_ptr : mesh.vertices() )
      v_ptr->index( vertex_index++ );

    // Quads
    for ( auto& q_ptr : mesh.quads() )
      q_ptr->index( element_index++ );

    // Triangles 
    for ( auto& t_ptr : mesh.triangles() )
      t_ptr->index( element_index++ );

  } // Cleanup::assign_mesh_indices()


  /*------------------------------------------------------------------
  | Each mesh vertex gets the value of the domain's size function 
  | at its position. This is required for a proper mesh output.
  ------------------------------------------------------------------*/
  template <typename Mesh, typename Domain>
  static inline void assign_size_function_to_vertices(Mesh& mesh,
                                                      const Domain& domain)
  {
    for ( auto& v_ptr : mesh.vertices() )
      v_ptr->mesh_size( domain.size_function(v_ptr->xy()) );

  } // Cleanup:assign_size_function_to_vertices()


  /*------------------------------------------------------------------
  | Every vertex gets assigned its neighboring vertices and these 
  | are then sorted by means of ascending angles
  | This function requires, that all the interior edges and 
  | boundary edges of the mesh have been generated 
  | -> element adjacency of the edges is not required here
  ------------------------------------------------------------------*/
  template <typename Mesh>
  static inline void setup_vertex_connectivity(Mesh& mesh)
  {
    // Remove all current vertex-to-vertex connectivities
    for ( auto& v_ptr : mesh.vertices() )
      v_ptr->vertices().clear();

    // Get vertex-to-vertex connectivity from interior edges
    for ( const auto& e : mesh.interior_edges() )
    {
      Vertex* v1 = &(e->v1());
      Vertex* v2 = &(e->v2());

      v1->vertices().push_back( v2 );
      v2->vertices().push_back( v1 );
    }

    // Get vertex-to-vertex connectivity from boundary edges
    for ( const auto& e : mesh.boundary_edges() )
    {
      Vertex* v1 = &(e->v1());
      Vertex* v2 = &(e->v2());

      v1->vertices().push_back( v2 );
      v2->vertices().push_back( v1 );
    }

    // Sort local vertex connectivities by ascending angle
    for ( auto& v_ptr : mesh.vertices() )
    {
      const Vec2d xy = v_ptr->xy();

      std::sort( v_ptr->vertices().begin(), v_ptr->vertices().end(),
      [xy] ( Vertex* v1, Vertex* v2 )
      {
        const Vec2d dxy1 = v1->xy() - xy;
        const Vec2d dxy2 = v2->xy() - xy;
        const double a1 = std::atan2(dxy1.y, dxy1.x);
        const double a2 = std::atan2(dxy2.y, dxy2.x);

        return ( a1 < a2 );
      });
    }

  } // Cleanup::setup_vertex_connectivity()


  /*------------------------------------------------------------------
  | Initialize the connectivity between facets and facets, as well  
  | as between edges and facets
  ------------------------------------------------------------------*/
  template <typename Mesh>
  static inline void setup_facet_connectivity(Mesh& mesh)
  {
    // Setup connectivity for interior edges
    for ( const auto& e_ptr : mesh.interior_edges() )
    {
      Facet* f1 = &NullFacet::get_instance();
      Facet* f2 = &NullFacet::get_instance();

      const Vertex& v1 = e_ptr->v1();
      const Vertex& v2 = e_ptr->v2();

      int idx1 = -1; int idx2 = -1;

      // Search the two facets (f1,f2), that share the current edge
      for ( auto f : v1.facets() )
      {
        int idx = f->get_edge_index(v1, v2);

        if ( idx < 0 ) continue;

        if ( idx1 < 0 )
        {
          idx1 = idx;
          f1 = f;
        }
        else if ( idx2 < 0 )
        {
          idx2 = idx;
          f2 = f;
          break;
        }
      }

      ASSERT( NullFacet::is_not_null(f1), 
        "Cleanup::setup_facet_connectivity(): "
        "First facet is NullFacet.");

      ASSERT( NullFacet::is_not_null(f2), 
        "Cleanup::setup_facet_connectivity(): "
        "Second facet is NullFacet.");

      // Setup connectivity between facets (f1,f2)
      f1->neighbor( idx1, f2 );
      f2->neighbor( idx2, f1 );

      // Setup connectivity between internal edge and facets (f1,f2)
      if ( is_left( v1.xy(), v2.xy(), f1->xy() ) )
      {
        e_ptr->facet_l( f1 );
        e_ptr->facet_r( f2 );
      }
      else
      {
        e_ptr->facet_l( f2 );
        e_ptr->facet_r( f1 );
      }
    }

    // Setup connectivity for boundary edges
    for ( const auto& e_ptr : mesh.boundary_edges() )
    {
      Facet* f1 = &NullFacet::get_instance();

      const Vertex& v1 = e_ptr->v1();
      const Vertex& v2 = e_ptr->v2();

      // Search adjacent facets to the current boundary edge
      for ( auto f : v1.facets() )
      {
        int idx = f->get_edge_index(v1, v2);

        if ( idx < 0 ) continue;

        f1 = f;
        break;
      }

      // Setup connectivity between boundary edge and facet
      if ( is_left( v1.xy(), v2.xy(), f1->xy() ) )
        e_ptr->facet_l( f1 );
      else
        ASSERT(false, "Invalid mesh boundary edge");
    }
    
  } // Mesh::setup_facet_connectivity


  /*------------------------------------------------------------------
  | Clean up quad elements:
  |
  | It may be, that some adjacent quads share two internal edges
  | Such elements are merged by this function.
  | 
  |                     v3                       vp   
  |                       x---------------------x
  |                       | \                   |
  |                       |   \ e2     q_nbr    |
  |                       |     \               |
  |                       |       \  v2         |
  |                       |         x           |
  |                       |           \         |
  |                       |   q_cur     \  e1   |
  |                       |               \     |
  |                       |                 \   |
  |                       |                   \ |
  |                       x---------------------x
  |                     v4                       v1
  |
  | -> This function requires a finalized mesh and that the functions
  |    assign_mesh_indices(), setup_vertex_connectivity() 
  |    and setup_facet_connectivity() have been called before.
  |
  ------------------------------------------------------------------*/
  template <typename Mesh>
  static inline void clear_double_quad_edges(Mesh& mesh, bool init=true)
  {
    if (init)
    {
      Cleanup::assign_mesh_indices(mesh);
      Cleanup::setup_vertex_connectivity(mesh);
      Cleanup::setup_facet_connectivity(mesh);
    }

    // Initialize arrays to store entities temporally
    std::vector<std::pair<Quad*,Quad*>>     quads_to_remove;
    std::vector<std::pair<Edge*,Edge*>>     edges_to_remove;
    std::vector<Vertex*>                    verts_to_remove;
    std::vector<std::pair<Vertex*,Vertex*>> opposing_vertices;

    // Array used to mark elements that will be merged
    // --> Quads are assumed to be located before triangles
    std::vector<bool> element_marker(mesh.n_elements(), false);
    size_t i_quad = 0;

    for ( auto& q_cur : mesh.quads() )
    {
      if ( element_marker[i_quad] ) 
        continue;

      ++i_quad;

      // Loop over all edges of the current quad
      for ( int i_vert = 0; i_vert < 4; ++i_vert )
      {
        // Vertices of the current quad
        Vertex& v1 = q_cur->vertex( i_vert );
        Vertex& v2 = q_cur->vertex( MOD(i_vert + 1, 4) );
        Vertex& v3 = q_cur->vertex( MOD(i_vert + 2, 4) );
        Vertex& v4 = q_cur->vertex( MOD(i_vert + 3, 4) );

        // Get the neighboring facets of the current successive
        // quad edges (v1,v2) and (v2,v3)
        int idx_1 = q_cur->get_edge_index(v1, v2);
        int idx_2 = q_cur->get_edge_index(v2, v3);

        Facet* nbr_1 = q_cur->neighbor(idx_1);
        Facet* nbr_2 = q_cur->neighbor(idx_2);

        // Proceed, if no neighbors are found (NullFacet)
        if ( NullFacet::is_null(nbr_1) || NullFacet::is_null(nbr_2) )
          continue;

        // Proceed if neighbors of the adjacent edges differ
        if ( nbr_1 != nbr_2 )
          continue;

        ASSERT( nbr_1 != nullptr, 
          "Cleanup::clear_double_quad_edges(): "
          "First edge neighbor is nullptr.");

        ASSERT( nbr_2 != nullptr, 
          "Cleanup::clear_double_quad_edges(): "
          "Second edge neighbor is nullptr.");

        // Proceed if neighbor has already been added
        if ( element_marker[nbr_1->index()] )
          continue;

        // In this stage, we address only quad / quad connections
        if ( nbr_1->n_vertices() < 4 )
          continue;

        // Now we can cast the facet to a quad
        Quad* q_nbr = static_cast<Quad*>(nbr_1);

        // Get the internal edges adjacent to both current quads
        Edge* e1 = mesh.interior_edges().get_edge(v1, v2);
        Edge* e2 = mesh.interior_edges().get_edge(v2, v3);

        ASSERT( (e1 != e2), 
          "Cleanup::clear_double_quad_edges(): "
          "Invalid edge data structure.");

        ASSERT( ( e1->v1() == e2->v1() || e1->v1() == e2->v2() 
               || e1->v2() == e2->v1() || e1->v2() == e2->v2() ),
          "Cleanup::clear_double_quad_edges(): "
          "Wrong edges found for quad merging.");

        // Mark the current quads, such that these won't get chosen
        // in upcoming loops
        element_marker[ q_cur->index() ] = true;
        element_marker[ q_nbr->index() ] = true;

        // Add elements to the removal vectors
        quads_to_remove.push_back( {q_cur.get(), q_nbr} ); 
        edges_to_remove.push_back( {e1, e2} );
        verts_to_remove.push_back( &v2 );
        
        // We still need the vertex of the neighboring quad, that 
        // is located on the opposite of the current edge segments
        // --> Use internal edge definition of quads
        int idx_op = q_nbr->get_edge_index(v2,v3);
        Vertex& v_op = q_nbr->vertex( idx_op );

        opposing_vertices.push_back( {&v4, &v_op} );

        ASSERT( (v_op != v1), 
          "Cleanup::clear_double_quad_edges(): "
          "Messed up quad data structure (1).");
        ASSERT( (v_op != v2), 
          "Cleanup::clear_double_quad_edges(): "
          "Messed up quad data structure (2).");
        ASSERT( (v_op != v3), 
          "Cleanup::clear_double_quad_edges(): "
          "Messed up quad data structure (3).");

        // At this point we can break the inner loop over 
        // the quad edges
        break;
      }
    }

    // Merge the marked edges and created triangles from the 
    // adjacent marked quads
    for (size_t i = 0; i < quads_to_remove.size(); ++i)
    {
      Edge* e1 = edges_to_remove[i].first;
      Edge* e2 = edges_to_remove[i].second;

      Vertex* v2 = verts_to_remove[i];

      // Get correct vertex order of edge segments
      Vertex& v1 = (e1->v2() == *v2) 
                 ? e1->v1() : e1->v2();
      Vertex& v3 = (e2->v1() == *v2) 
                 ? e2->v2() : e2->v1();

      Vertex* o1 = opposing_vertices[i].first;
      Vertex* o2 = opposing_vertices[i].second;

      // Create new quad 
      Quad& q_new = mesh.add_quad( *o1, v1, *o2, v3 );
      q_new.is_active(true);
    }

    // Removal of old quads
    for (size_t i = 0; i < quads_to_remove.size(); ++i)
    {
      Quad* q1 = quads_to_remove[i].first;
      Quad* q2 = quads_to_remove[i].second;

      mesh.remove_quad( *q1 );
      mesh.remove_quad( *q2 );
    }

    // Removal of old interior edges
    for (size_t i = 0; i < edges_to_remove.size(); ++i)
    {
      Edge* e1 = edges_to_remove[i].first;
      Edge* e2 = edges_to_remove[i].second;

      mesh.remove_interior_edge( *e1 );
      mesh.remove_interior_edge( *e2 );
    }

    // Removal of old vertices
    for (size_t i = 0; i < verts_to_remove.size(); ++i)
    {
      Vertex* v = verts_to_remove[i];
      mesh.remove_vertex( *v );
    }

    // Re-initialize facet-to-facet connectivity
    Cleanup::assign_mesh_indices(mesh);
    Cleanup::setup_vertex_connectivity(mesh);
    Cleanup::setup_facet_connectivity(mesh);

  } // Cleanup::clear_double_quad_edges() 


  /*------------------------------------------------------------------
  | Clean up triangle elements:
  | 
  | It may be, that some triangles share two internal edges with a 
  | single quad element. 
  | These triangles and quads will be removed in this function and then
  | replaced by a single triangle.
  |
  |               v3                                  v1
  |                x---------------------------------x
  |                 \                               /
  |                  \\                           //
  |                   \ \         t_nbr         / /
  |                    \  \                   /  /
  |                     \   \               /   /
  |                      \    \           /    /
  |                       \     \       /     /
  |                        \      \   /      /
  |                         \       x       /
  |                          \     v2      /
  |                           \           /
  |                            \  q_cur  /
  |                             \       /
  |                              \     /
  |                               \   /
  |                                \ /
  |                                 x
  |                                  v4
  |
  | -> This function requires a finalized mesh and that the functions
  |    assign_mesh_indices(), setup_vertex_connectivity() 
  |    and setup_facet_connectivity() have been called before.
  |
  | -> This function should be called after the function 
  |    clean_double_quad_edges(), because such triangles might be 
  |    genereated during the latter function.
  |
  ------------------------------------------------------------------*/
  template <typename Mesh>
  static inline void clear_double_triangle_edges(Mesh& mesh, bool init=true)
  {
    if (init)
    {
      Cleanup::assign_mesh_indices(mesh);
      Cleanup::setup_vertex_connectivity(mesh);
      Cleanup::setup_facet_connectivity(mesh);
    }

    // Initialize arrays to store entities temporally
    std::vector<std::pair<Quad*,Triangle*>> elements_to_remove;
    std::vector<std::pair<Edge*,Edge*>>     edges_to_remove;
    std::vector<Vertex*>                    verts_to_remove;

    // Array used to mark quads that will be merged
    // --> Quads are assumed to be located before triangles
    std::vector<bool> element_marker(mesh.n_elements(), false);
    size_t i_quad = 0;

    for ( auto& q_cur : mesh.quads() )
    {
      if ( element_marker[i_quad] ) 
        continue;

      ++i_quad;

      // Loop over all edges of the current quad
      for ( int i_vert = 0; i_vert < 4; ++i_vert )
      {
        // Vertices of the current quad
        Vertex& v1 = q_cur->vertex( i_vert );
        Vertex& v2 = q_cur->vertex( MOD(i_vert + 1, 4) );
        Vertex& v3 = q_cur->vertex( MOD(i_vert + 2, 4) );

        // Get the neighboring facets of the current successive
        // quad edges (v1,v2) and (v2,v3)
        int idx_1 = q_cur->get_edge_index(v1, v2);
        int idx_2 = q_cur->get_edge_index(v2, v3);

        Facet* nbr_1 = q_cur->neighbor(idx_1);
        Facet* nbr_2 = q_cur->neighbor(idx_2);

        // Proceed, if no neighbors are found (NullFacet)
        if ( NullFacet::is_null(nbr_1) || NullFacet::is_null(nbr_2) )
          continue;

        // Proceed if neighbors of the adjacent edges differ
        if ( nbr_1 != nbr_2 )
          continue;

        ASSERT( nbr_1 != nullptr, 
          "Cleanup::clear_double_quad_edges(): "
          "First edge neighbor is nullptr.");

        ASSERT( nbr_2 != nullptr, 
          "Cleanup::clear_double_quad_edges(): "
          "Second edge neighbor is nullptr.");

        // Proceed if neighbor has already been added
        if ( element_marker[nbr_1->index()] )
          continue;

        // In this stage, we address only quad / triangle connections
        if ( nbr_1->n_vertices() > 3 )
          continue;

        // Now we can cast the facet to a triangle
        Triangle* t_nbr = static_cast<Triangle*>(nbr_1);

        // Get the internal edges adjacent to both current elements
        Edge* e1 = mesh.interior_edges().get_edge(v1, v2);
        Edge* e2 = mesh.interior_edges().get_edge(v2, v3);

        ASSERT( (e1 != e2), 
          "Cleanup::clear_double_triangle_edges(): "
          "Invalid edge data structure.");

        ASSERT( ( e1->v1() == e2->v1() || e1->v1() == e2->v2() 
               || e1->v2() == e2->v1() || e1->v2() == e2->v2() ),
          "Cleanup::clear_double_triangle_edges(): "
          "Wrong edges found for element merging.");

        // Mark the current quads, such that the won't get chosen
        // in upcoming loops
        element_marker[ q_cur->index() ] = true;
        element_marker[ t_nbr->index() ] = true;

        // Add elements to the removal vectors
        elements_to_remove.push_back( {q_cur.get(), t_nbr} ); 
        edges_to_remove.push_back( {e1, e2} );
        verts_to_remove.push_back( &v2 );

        // At this point we can break the inner loop over 
        // the quad edges
        break;
      }
    }

    // Create a new triangle from the remaining vertices
    for ( size_t i = 0; i < elements_to_remove.size(); ++i )
    {
      Quad*   q = elements_to_remove[i].first;
      Vertex* v = verts_to_remove[i];

      int id_v = q->get_vertex_index( *v );

      ASSERT( (id_v > -1), 
        "Cleanup::clear_double_triangle_edges(): "
        "Invalid quad data structure.");

      // Get the remaining three vertices of the quad 
      Vertex& v1 = q->vertex( MOD(id_v+1, 4) );
      Vertex& v2 = q->vertex( MOD(id_v+2, 4) );
      Vertex& v3 = q->vertex( MOD(id_v+3, 4) );

      // Create new triangle 
      Triangle& t_new = mesh.add_triangle( v1, v2, v3 );
      t_new.is_active(true);
    }

    // Removal of old elements
    for ( size_t i = 0; i < elements_to_remove.size(); ++i )
    {
      Quad*     q = elements_to_remove[i].first;
      Triangle* t = elements_to_remove[i].second;
      
      mesh.remove_quad( *q );
      mesh.remove_triangle( *t );
    }

    // Removal of old interior edges
    for (size_t i = 0; i < edges_to_remove.size(); ++i)
    {
      Edge* e1 = edges_to_remove[i].first;
      Edge* e2 = edges_to_remove[i].second;

      mesh.remove_interior_edge( *e1 );
      mesh.remove_interior_edge( *e2 );
    }

    // Removal of old vertices
    for (size_t i = 0; i < verts_to_remove.size(); ++i)
    {
      Vertex* v = verts_to_remove[i];
      mesh.remove_vertex( *v );
    }

    // Re-initialize facet-to-facet connectivity
    Cleanup::assign_mesh_indices(mesh);
    Cleanup::setup_vertex_connectivity(mesh);
    Cleanup::setup_facet_connectivity(mesh);

  } // Cleanup::clear_double_triangle_edges()




private:
  /*------------------------------------------------------------------
  | We hide the constructor, since this class acts only as container
  | for static inline functions
  ------------------------------------------------------------------*/
  Cleanup() = default;
  ~Cleanup() {};



}; // Cleanup

} // namespace TQAlgorithm
} // namespace TQMesh
