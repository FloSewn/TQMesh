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

#include "Vertex.h"
#include "Edge.h"
#include "Triangle.h"
#include "Quad.h"
#include "Facet.h"

namespace TQMesh {

using namespace CppUtils;

/*********************************************************************
* This class contains all functions that are required to cleanup and
* prepare the mesh after the generation process
*********************************************************************/
class MeshCleanup
{
public:

  using VertexList = std::list<Vertex*>;
  using EdgeList   = std::list<Edge*>;

  /*------------------------------------------------------------------
  | Adjust the coordinate of a given vertex, while accounting for
  | the entire underlying mesh structure
  ------------------------------------------------------------------*/
  static inline void set_vertex_coordinates(Vertex& v, const Vec2d& xy)
  {
    v.adjust_xy( xy );

    for ( auto& e : v.edges() )
      e->update_metrics();

    for ( auto& f : v.facets() )
      f->update_metrics();
  } 


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

  } // MeshCleanup::assign_mesh_indices()


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

  } // MeshCleanup:assign_size_function_to_vertices()


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

      if ( NullFacet::is_not_null(f1) )
      {
        f1->neighbor( idx1, f2 );

        if ( is_left( v1.xy(), v2.xy(), f1->xy() ) )
          e_ptr->facet_l( f1 );
        else
          e_ptr->facet_r( f1 );
      }

      if ( NullFacet::is_not_null(f2) )
      {
        f2->neighbor( idx2, f1 );

        if ( is_left( v1.xy(), v2.xy(), f2->xy() ) )
          e_ptr->facet_l( f2 );
        else
          e_ptr->facet_r( f2 );
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

      if ( NullFacet::is_null(f1) )
        continue;

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
  |    assign_mesh_indices() and setup_facet_connectivity() 
  |    have been called before.
  |
  ------------------------------------------------------------------*/
  template <typename Mesh>
  static inline void clear_double_quad_edges(Mesh& mesh, bool init=true)
  {
    if (init)
    {
      MeshCleanup::assign_mesh_indices(mesh);
      MeshCleanup::setup_facet_connectivity(mesh);
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
          "MeshCleanup::clear_double_quad_edges(): "
          "First edge neighbor is nullptr.");

        ASSERT( nbr_2 != nullptr, 
          "MeshCleanup::clear_double_quad_edges(): "
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
          "MeshCleanup::clear_double_quad_edges(): "
          "Invalid edge data structure.");

        ASSERT( ( e1->v1() == e2->v1() || e1->v1() == e2->v2() 
               || e1->v2() == e2->v1() || e1->v2() == e2->v2() ),
          "MeshCleanup::clear_double_quad_edges(): "
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
          "MeshCleanup::clear_double_quad_edges(): "
          "Messed up quad data structure (1).");
        ASSERT( (v_op != v2), 
          "MeshCleanup::clear_double_quad_edges(): "
          "Messed up quad data structure (2).");
        ASSERT( (v_op != v3), 
          "MeshCleanup::clear_double_quad_edges(): "
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
    MeshCleanup::assign_mesh_indices(mesh);
    MeshCleanup::setup_facet_connectivity(mesh);

  } // MeshCleanup::clear_double_quad_edges() 


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
  |    assign_mesh_indices() and setup_facet_connectivity() have 
  |    been called before.
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
      MeshCleanup::assign_mesh_indices(mesh);
      MeshCleanup::setup_facet_connectivity(mesh);
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
          "MeshCleanup::clear_double_quad_edges(): "
          "First edge neighbor is nullptr.");

        ASSERT( nbr_2 != nullptr, 
          "MeshCleanup::clear_double_quad_edges(): "
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
          "MeshCleanup::clear_double_triangle_edges(): "
          "Invalid edge data structure.");

        ASSERT( ( e1->v1() == e2->v1() || e1->v1() == e2->v2() 
               || e1->v2() == e2->v1() || e1->v2() == e2->v2() ),
          "MeshCleanup::clear_double_triangle_edges(): "
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
        "MeshCleanup::clear_double_triangle_edges(): "
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
    MeshCleanup::assign_mesh_indices(mesh);
    MeshCleanup::setup_facet_connectivity(mesh);

  } // MeshCleanup::clear_double_triangle_edges()

  /*------------------------------------------------------------------
  | This function locates vertices that are adjacent to exactly 
  | three triangles of same color and merges them to a single one
  |           x                       x
  |          /|\                     / \
  |         / | \                   /   \
  |        /  x  \     ----->      /     \
  |       /  / \  \               /       \
  |      / /     \ \             /         \
  |     //         \\           /           \
  |    x-------------x         x-------------x
  |
  ------------------------------------------------------------------*/
  template <typename Mesh>
  static inline void merge_degenerate_triangles(Mesh& mesh, bool init=true)
  {
    // Collect all vertices, that are adjacent to exactly three 
    // triangles
    VertexList bad_vertices {};

    for ( const auto& vertex_ptr : mesh.vertices() )
    {
      // Skip fixed vertices
      if ( vertex_ptr->is_fixed() )
        continue;

      if ( vertex_ptr->facets().size() != 3 )
        continue;
      if ( vertex_ptr->edges().size() != 3 )
        continue;
      if ( vertex_ptr->facets(0).n_vertices() != 3 )
        continue;
      if ( vertex_ptr->facets(1).n_vertices() != 3 )
        continue;
      if ( vertex_ptr->facets(2).n_vertices() != 3 )
        continue;

      auto color = vertex_ptr->facets(0).color();

      if ( vertex_ptr->facets(1).color() != color )
        continue;
      if ( vertex_ptr->facets(2).color() != color )
        continue;

      bad_vertices.push_back( vertex_ptr.get() );
    }



    for ( auto& bad_vertex : bad_vertices )
    {
      const Triangle* t0 = static_cast<const Triangle*>(
          &bad_vertex->facets(0));

      const Triangle* t1 = static_cast<const Triangle*>(
          &bad_vertex->facets(1));

      const Triangle* t2 = static_cast<const Triangle*>(
          &bad_vertex->facets(2));

      auto color = t0->color();

      // Get surrounding vertices
      int i0 = t0->get_vertex_index( *bad_vertex );
      int i1 = t1->get_vertex_index( *bad_vertex );
      int i2 = t2->get_vertex_index( *bad_vertex );

      Vertex& v0 = const_cast<Vertex&>( t0->vertex(MOD(i0+1,3)) );
      Vertex& v1 = const_cast<Vertex&>( t1->vertex(MOD(i1+1,3)) );
      Vertex& v2 = const_cast<Vertex&>( t2->vertex(MOD(i2+1,3)) );

      // Remove interior edges 
      Edge* e0 = mesh.interior_edges().get_edge(*bad_vertex, v0);
      Edge* e1 = mesh.interior_edges().get_edge(*bad_vertex, v1);
      Edge* e2 = mesh.interior_edges().get_edge(*bad_vertex, v2);

      ASSERT( e0, "Interior edge not defined.\n" );
      mesh.remove_interior_edge( *e0 );

      ASSERT( e1, "Interior edge not defined.\n" );
      mesh.remove_interior_edge( *e1 );

      ASSERT( e2, "Interior edge not defined.\n" );
      mesh.remove_interior_edge( *e2 );

      // Remove triangles
      mesh.remove_triangle(const_cast<Triangle&>(*t0));
      mesh.remove_triangle(const_cast<Triangle&>(*t1));
      mesh.remove_triangle(const_cast<Triangle&>(*t2));

      // Remove vertex
      mesh.remove_vertex( *bad_vertex );

      // Add new triangle
      auto o = orientation(v0.xy(), v1.xy(), v2.xy());

      if ( o == Orientation::CCW ) 
        mesh.add_triangle( v0, v1, v2, color );
      else
        mesh.add_triangle( v0, v2, v1, color );
    }

    // Remove deleted entities
    mesh.clear_waste();

    // Re-initialize facet-to-facet connectivity
    MeshCleanup::assign_mesh_indices(mesh);
    MeshCleanup::setup_facet_connectivity(mesh);
   
  } // merge_degenerate_triangles()

  /*------------------------------------------------------------------
  | This function loops over all internal edges and, if possible,
  | merges two adjacent triangles to one quad element.
  ------------------------------------------------------------------*/
  template <typename Mesh>
  static inline void merge_triangles_to_quads(Mesh& mesh, bool init=true)
  {
    if (init)
      MeshCleanup::setup_facet_connectivity(mesh);

    // Pick all internal edges, which are adjacent to two triangles
    EdgeList tri_edges {};

    for ( const auto& e_ptr : mesh.interior_edges() )
    {
      Facet* f_l = e_ptr->facet_l();
      Facet* f_r = e_ptr->facet_r();

      if (  ( NullFacet::is_not_null(f_l) && f_l->n_vertices() == 3 )
         && ( NullFacet::is_not_null(f_r) && f_r->n_vertices() == 3 ) )
        tri_edges.push_back( e_ptr.get() );
    }

    // Sort edge list with increasing minimum edge length of their
    // adjacent triangles
    tri_edges.sort(
    []( Edge* a, Edge* b )
    {
      const double a_l = a->facet_l()->min_edge_length();
      const double a_r = a->facet_r()->min_edge_length();
      const double a_ang = MIN(a_l, a_r);

      const double b_l = b->facet_l()->min_edge_length();
      const double b_r = b->facet_r()->min_edge_length();
      const double b_ang = MIN(b_l, b_r);

      return a_ang < b_ang;
    });

    // Loop over all sorted edges and merge their adjacent triangles
    // to quadrilaterals. It may be, that the triangle of an 
    // upcoming edge has already been merged in this process. Thus
    // we have to make sure, that the triangles to merge still
    // exist.
    //
    //   v2               q_r
    //     *-------------*
    //     | \           |
    //     |   \   f_r   |
    //     |     \       |
    //     |       \     |
    //     |  f_l    \   |
    //     |           \ |
    //     *-------------*
    //    q_l             v1
    //
    for ( auto& e : tri_edges )
    {
      Facet* f_l = e->facet_l();
      Facet* f_r = e->facet_r();

      // Triangles have already been merged with prior edge
      if ( f_l->n_vertices() > 3 || f_r->n_vertices() > 3 ) 
        continue;

      Vertex& v1 = e->v1();
      Vertex& v2 = e->v2();

      int i_l = f_l->get_edge_index(v1, v2);
      int i_r = f_r->get_edge_index(v1, v2);

      Vertex& q_l = f_l->vertex(i_l);
      Vertex& q_r = f_r->vertex(i_r);

      // Remove internal edge
      mesh.remove_interior_edge( *e );

      // Create new quadrilateral element
      Quad& q_new = mesh.add_quad( q_l, v1, q_r, v2 );
      q_new.is_active( true );

      // Update internal edge connectivity to prevent bad memory 
      // access for upcoming internal edges that are no longer 
      // adjacent to two triangles due to the merge
      for (int i = 0; i < 4; ++i)
      {
        int i1 = i;
        int i2 = MOD(i+1, 4);

        Vertex& q1 = q_new.vertex(i1);
        Vertex& q2 = q_new.vertex(i2);

        Edge* e_share = mesh.interior_edges().get_edge(q1,q2);

        // Boundary edge found
        if ( e_share == nullptr ) 
          continue;

        Facet* t_l = e_share->facet_l();
        Facet* t_r = e_share->facet_r();

        if ( t_l && (t_l == f_l || t_l == f_r) )
          e_share->facet_l( &q_new );

        if ( t_r && (t_r == f_l || t_r == f_r) )
          e_share->facet_r( &q_new );
      }

      // Remove triangles
      mesh.remove_triangle( *(static_cast<Triangle*>(f_l)) );
      mesh.remove_triangle( *(static_cast<Triangle*>(f_r)) );
    }

    // Remove deleted entities
    mesh.clear_waste();

    // Bad elements may have been created up to this point
    // due to the merging of triangles to quads
    // The two upcoming function fix these bad elements
    MeshCleanup::clear_double_quad_edges(mesh, true);
    MeshCleanup::clear_double_triangle_edges(mesh, false);

  } // merge_triangles_to_quads()

private:
  /*------------------------------------------------------------------
  | We hide the constructor, since this class acts only as container
  | for static inline functions
  ------------------------------------------------------------------*/
  MeshCleanup() = default;
  ~MeshCleanup() {};

}; // MeshCleanup

} // namespace TQMesh
