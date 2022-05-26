/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <algorithm>

#include "Vec2.h"
#include "utils.h"
#include "ProgressBar.h"

#include "Vertex.h"
#include "Triangle.h"
#include "Quad.h"
#include "Front.h"
#include "Boundary.h"
#include "Domain.h"
#include "QuadLayer.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace TQUtils;

/*********************************************************************
* The mesh 
*********************************************************************/
class Mesh
{
  using EdgePair       = std::pair<const Edge*,const Edge*>;
  using EdgePairVector = std::vector<EdgePair>;
  using EdgeVector     = std::vector<Edge*>;
  using VertexVector   = std::vector<Vertex*>;
  using TriVector      = std::vector<Triangle*>;
  using QuadVector     = std::vector<Quad*>;
  using DoubleVector   = std::vector<double>;
  using Vec2dVector    = std::vector<Vec2d>;

public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Mesh(Domain&   domain,
       double    qtree_scale=TQ_QTREE_SCALE,
       size_t    qtree_items=TQ_QTREE_ITEMS, 
       size_t    qtree_depth=TQ_QTREE_DEPTH,
       bool      init_structure=true)
  : domain_ { &domain }
  , verts_  { qtree_scale, qtree_items, qtree_depth }
  , tris_   { qtree_scale, qtree_items, qtree_depth }
  , quads_  { qtree_scale, qtree_items, qtree_depth }
  , front_  { }
  {
    // Return if vertices & front should not be initialized
    if (!init_structure)
      return;

    // Copy vertices from domain
    for ( const auto& v_ptr : domain.vertices() )
    {
      Vertex& v = verts_.push_back( v_ptr->xy(), 
                                    v_ptr->sizing(), 
                                    v_ptr->range() );
      v.on_front( v_ptr->on_front() );
      v.on_boundary( v_ptr->on_boundary() );
      v.is_fixed( v_ptr->is_fixed() );
    }

    // Initialize edges in the advancing front
    front_.init_front_edges(domain, verts_);

    // Setup boundary edges from the initial advancing front
    for ( auto& e : front_ )
      bdry_edges_.add_edge( e->v1(), e->v2(), e->marker() );

  } // Mesh::Constructor()

  /*------------------------------------------------------------------
  | This function takes care of the garbage collection 
  ------------------------------------------------------------------*/
  void clear_waste()
  {
    verts_.clear_waste();
    tris_.clear_waste();
    quads_.clear_waste();
    front_.clear_waste();
    intr_edges_.clear_waste();
    bdry_edges_.clear_waste();
  }

  /*------------------------------------------------------------------
  | Getter
  ------------------------------------------------------------------*/
  const Front& front() const { return front_; }
  Front& front() { return front_; }

  const Triangles& triangles() const { return tris_; }
  Triangles& triangles() { return tris_; }

  const Quads& quads() const { return quads_; }
  Quads& quads() { return quads_; }

  const Vertices& vertices() const { return verts_; }
  Vertices& vertices() { return verts_; }

  const EdgeList& interior_edges() const { return intr_edges_; }
  EdgeList& interior_edges() { return intr_edges_; }

  const EdgeList& boundary_edges() const { return bdry_edges_; }
  EdgeList& boundary_edges() { return bdry_edges_; }

  double area() const { return mesh_area_; }


  /*------------------------------------------------------------------
  | This function assigns a corresponding global index to each entity
  | of the current mesh
  | --> Quads are stored prior to triangles
  ------------------------------------------------------------------*/
  void assign_mesh_indices()
  {
    unsigned int v_index = 0; 
    unsigned int e_index = 0;

    // Vertices 
    for ( const auto& v_ptr : verts_ )
      v_ptr->index( v_index++ );

    // Quads
    for ( const auto& q_ptr : quads_ )
      q_ptr->index( e_index++ );

    // Triangles 
    for ( const auto& t_ptr : tris_ )
      t_ptr->index( e_index++ );

  } // assign_mesh_indices()

  /*------------------------------------------------------------------
  | Every triangle is refined into three quads and every 
  | quads is refined into four quads.
  | This results in an all-quad mesh.
  |
  |                        
  |          (v7)         (v6)      (v5)      qi... sub-quads 
  |             o<--------o---------o         vi... vertices
  |             |         |         ^
  |             |  [q4]   |   [q3]  |
  |             |         |         |
  |             |         |(v9)     |  
  |         (v8)o---------o---------o(v4) 
  |             |         |         |
  |             |         |         |
  |             |  [q1]   |   [q2]  |
  |             v         |         |
  |             o---------o-------->o
  |          (v1)        (v2)       (v3)
  |                         
  |                         
  |                        (v5)                 
  |                         o                   
  |                       /   \
  |                     / [q3]  \
  |                   /           \
  |           (v6)  /      (v7)     \  (v4)
  |               o---------o---------o
  |             /           |           \
  |           /             |             \
  |         /     [q1]      |      [q2]     \
  |       /                 |                 \
  |     o-------------------o-------------------o
  |    (v1)                (v2)                 (v3)
  |
  ------------------------------------------------------------------*/
  bool refine_to_quads()
  {
    MSG("START MESH REFINEMENT");
    clear_waste();

    // Gather all coarse edges, quads and tris
    EdgeVector coarse_intr_edges {};
    EdgeVector coarse_bdry_edges {};
    TriVector  coarse_tris {};
    QuadVector coarse_quads {};

    for ( auto& e : intr_edges_ )
      coarse_intr_edges.push_back( e.get() );

    for ( auto& e : bdry_edges_ )
      coarse_bdry_edges.push_back( e.get() );

    for ( auto& t_ptr : tris_ )
      coarse_tris.push_back( t_ptr.get() );

    for ( auto& q_ptr : quads_ )
      coarse_quads.push_back( q_ptr.get() );

    // Refine interior edges
    for ( auto e : coarse_intr_edges )
    {
      Vertex& v = verts_.push_back( e->xy() );
      intr_edges_.add_edge( e->v1(), v );
      intr_edges_.add_edge( v, e->v2() );
      e->sub_vertex( &v );

      if ( e->v1().is_fixed() && e->v2().is_fixed() )
        v.is_fixed( true );
    }

    // Refine boundary edges
    for ( auto e : coarse_bdry_edges )
    {
      Vertex& v = verts_.push_back( e->xy() );
      bdry_edges_.add_edge( e->v1(), v, e->marker() );
      bdry_edges_.add_edge( v, e->v2(), e->marker() );
      e->sub_vertex( &v );
      v.on_boundary( true );

      if ( e->v1().is_fixed() && e->v2().is_fixed() )
        v.is_fixed( true );
    }

    // Refine all triangle elements
    for ( auto t : coarse_tris )
    {
      Vertex& v1 = t->v1();
      Vertex& v3 = t->v2();
      Vertex& v5 = t->v3();

      Edge* e13 = get_edge( v1, v3 );
      Edge* e35 = get_edge( v3, v5 );
      Edge* e51 = get_edge( v5, v1 );

      Vertex* v2 = e13->sub_vertex();
      Vertex* v4 = e35->sub_vertex();
      Vertex* v6 = e51->sub_vertex();

      // Create new vertex at center of quad
      Vertex& v7 = verts_.push_back( t->xy() );

      // Create new sub-quads 
      quads_.push_back( v1, *v2, v7, *v6 );
      quads_.push_back( v3, *v4, v7, *v2 );
      quads_.push_back( v5, *v6, v7, *v4 );

      // Create new interior edges
      intr_edges_.add_edge( *v2, v7 );
      intr_edges_.add_edge( *v4, v7 );
      intr_edges_.add_edge( *v6, v7 );

    }

    // Refine all quad elements
    for ( auto q : coarse_quads )
    {
      Vertex& v1 = q->v1();
      Vertex& v3 = q->v2();
      Vertex& v5 = q->v3();
      Vertex& v7 = q->v4();

      Edge* e13 = get_edge( v1, v3 );
      Edge* e35 = get_edge( v3, v5 );
      Edge* e57 = get_edge( v5, v7 );
      Edge* e71 = get_edge( v7, v1 );

      Vertex* v2 = e13->sub_vertex();
      Vertex* v4 = e35->sub_vertex();
      Vertex* v6 = e57->sub_vertex();
      Vertex* v8 = e71->sub_vertex();

      // Create new vertex at center of quad
      Vertex& v9 = verts_.push_back( q->xy() );

      // Create new sub-quads 
      quads_.push_back( v1, *v2, v9, *v8 );
      quads_.push_back( v3, *v4, v9, *v2 );
      quads_.push_back( v5, *v6, v9, *v4 );
      quads_.push_back( v7, *v8, v9, *v6 );

      // Create new interior edges
      intr_edges_.add_edge( *v2, v9 );
      intr_edges_.add_edge( *v4, v9 );
      intr_edges_.add_edge( *v6, v9 );
      intr_edges_.add_edge( *v8, v9 );

    }

    // Remove old entitires
    for ( auto e : coarse_intr_edges )
      check_remove_interior_edge( *e );
      
    for ( auto e : coarse_bdry_edges )
      check_remove_boundary_edge( *e );

    for ( auto t : coarse_tris )
      check_remove_triangle( *t );

    for ( auto q : coarse_quads )
      check_remove_quad( *q );


    // Re-initialize vertex-to-vertex connectivity
    setup_vertex_connectivity();

    // Update connectivity between elements and edges
    setup_facet_connectivity();

    MSG("DONE!");
      
    return true;
      
  } // Mesh::refine_to_quads()

  /*------------------------------------------------------------------
  | Algorithm to create a single quad layer for a connected list of
  | advancing front edges that start with v_start and end with v_end
  ------------------------------------------------------------------*/
  bool create_quad_layers(Vertex& v_start, Vertex& v_end, 
                          size_t n_layers, double first_height,
                          double growth_ratio)
  {
    bool quads_generated;

    double height = first_height;

    Vertex* v1 = &v_start;
    Vertex* v2 = &v_end;

    for ( size_t i = 0; i < n_layers; ++i )
    {
      quads_generated = add_quad_layer( v1, v2, height );

      if ( !quads_generated )
        return false;

      height *= growth_ratio;
    }

    setup_facet_connectivity();
    merge_triangles_to_quads();

    return true;

  } // Mesh::create_quad_layers()

  /*------------------------------------------------------------------
  | Create a triangular mesh using the advancing front algorithm
  ------------------------------------------------------------------*/
  bool triangulate()
  {
    unsigned int iter = 0;
    bool wide_search = false;

    ProgressBar progress_bar {};

    // Initialize base edge
    front_.set_base_first();
    Edge* base = &( front_.base() );

    // Invalid base definition
    if ( !base ) 
      return false;

    MSG("START MESHING");

    // Start advancing front loop
    while ( true )
    {
      // Try to advance the current base edge
      bool success = advance_front_triangle(*base, wide_search);

      // If it worked, reset iteration counter and wide search
      // and go to the next base edge
      if ( success )
      {
        iter = 0;
        wide_search = false;

        front_.set_base_first();
        base = &( front_.base() );

        clear_waste();
      }
      // If it failed, go to the next base edge
      else
      {
        front_.set_base_next();
        base = &( front_.base() );
        ++iter;
      }

      // All front edges failed to create new elements
      // --> Activate wide search for neighboring vertices
      //     and re-run the algorithm
      if ( iter == front_.size() && !wide_search )
      {
        wide_search = true;
        iter = 0;
      }

      // Update progress bar
      double state = std::ceil(100.0 * area() / domain_->area());
      progress_bar.update( static_cast<int>(state) );
      progress_bar.show( std::clog );

      // No more edges in the advancing front
      // --> Meshing algorithm succeeded
      if ( front_.size() == 0 )
      {
        MSG("[");
        MSG("MESHING SUCCEEDED!");
        break;
      }

      // All front edges faild to create new elements, even
      // when using the wide search 
      // --> Meshing algorithm failed
      if ( iter == front_.size() && wide_search )
      {
        MSG("[");
        MSG("MESHING FAILED.");
        return false;
      }
    }

    // Initialize vertex-to-vertex connectivity
    setup_vertex_connectivity();

    // Initialize facet-to-facet connectivity
    // as well as connectivity between edges and facets
    setup_facet_connectivity();

    return true;

  } // Mesh::triangulate()

  /*------------------------------------------------------------------
  | Create a quadrilateral mesh using the advancing front algorithm
  ------------------------------------------------------------------*/
  bool pave()
  {
    unsigned int iter = 0;
    bool wide_search = false;

    ProgressBar progress_bar {};

    // Initialize base edge
    front_.set_base_first();
    Edge* base = &( front_.base() );

    // Invalid base definition
    if ( !base ) 
      return false;

    MSG("START MESHING");

    // Start advancing front loop
    while ( true )
    {
      // Try to advance the current base edge
      bool success = advance_front_quad(*base, wide_search, -1.0);

      // If it worked, reset iteration counter and wide search
      // and go to the next base edge
      if ( success )
      {
        iter = 0;
        wide_search = false;

        front_.set_base_first();
        base = &( front_.base() );

        clear_waste();
      }
      // If it failed, go to the next base edge
      else
      {
        front_.set_base_next();
        base = &( front_.base() );
        ++iter;
      }

      // All front edges failed to create new elements
      // --> Activate wide search for neighboring vertices
      //     and re-run the algorithm
      if ( iter == front_.size() && !wide_search )
      {
        wide_search = true;
        iter = 0;
      }

      // Update progress bar
      double state = std::ceil(100.0 * area() / domain_->area());
      progress_bar.update( static_cast<int>(state) );
      progress_bar.show( std::clog );

      // No more edges in the advancing front
      // --> Meshing algorithm succeeded
      if ( front_.size() == 0 )
      {
        MSG("");
        MSG("MESHING SUCCEEDED!");
        break;
      }

      // All front edges faild to create new elements, even
      // when using the wide search 
      // --> Meshing algorithm failed
      if ( iter == front_.size() && wide_search )
      {
        MSG("");
        MSG("MESHING FAILED.");
        return false;
      }
    }

    // Initialize vertex-to-vertex connectivity
    setup_vertex_connectivity();

    // Initialize facet-to-facet connectivity
    // as well as connectivity between edges and facets
    setup_facet_connectivity();

    // Merge remaining triangles to quads
    merge_triangles_to_quads();

    return true;

  } // Mesh::pave()

  /*------------------------------------------------------------------
  | Let the advancing front create a new quad
  ------------------------------------------------------------------*/
  bool advance_front_quad(Edge& base, bool wide_search=false,
                          double height = -1.0)
  {
    Vertex& b1 = base.v1();
    Vertex& b2 = base.v2();

    double h_b1 = height;
    double h_b2 = height;

    Vec2d p1 = b1.xy() + h_b1 * base.normal();
    Vec2d p2 = b2.xy() + h_b2 * base.normal();

    if ( height <= 0 )
    {
      h_b1 = domain_->size_function( b1.xy() );
      h_b2 = domain_->size_function( b2.xy() );
      
      p1 = b1.xy() + h_b1 * base.normal();
      p2 = b2.xy() + h_b2 * base.normal();

      const double rho_q1 = domain_->size_function( p1 );
      const double rho_q2 = domain_->size_function( p2 );

      const double theta_1 = MAX(-0.25,MIN(0.25,1.0-(h_b1-rho_q1)/h_b1));
      const double theta_2 = MAX(-0.25,MIN(0.25,1.0-(h_b2-rho_q2)/h_b2));

      p1 += theta_1 * base.length() * base.tangent();
      p2 -= theta_2 * base.length() * base.tangent();
    }

    // Search radii
    double r1 = domain_->size_function( p1 );
    double r2 = domain_->size_function( p2 );

    if ( height > 0 )
    {
      r1 = TQ_QUAD_RANGE_FACTOR * height;
      r2 = TQ_QUAD_RANGE_FACTOR * height;
    }

    // ****** Create first triangle *******
    TriVector new_tris_p1 {};

    // Find all vertices in vicinity of p1 
    VertexVector vertex_candidates_p1 
      = find_local_vertices(p1, r1, wide_search );

    // Create potential triangles with all found vertices
    check_vertex_candidates(vertex_candidates_p1, base, new_tris_p1);

    // If potential triangles have been found, choose the best one
    // --> this function also updates the advancing front
    Triangle* t1 = choose_best_triangle(new_tris_p1, base);

    // No triangle has been found yet -> create new one
    if ( !t1 )
    {
      Vertex&   v_new = verts_.push_back( p1 );
      Triangle& t_new = tris_.push_back(b1, b2, v_new);

      // Algorithm fails if new vertex or new triangle is invalid
      if ( !check_vertex( v_new ) || !check_triangle( t_new) )
      {
        check_remove_triangle( t_new );
        check_remove_vertex( v_new );

        return false;
      }

      // Update the advancing front with new vertex
      update_front( base, v_new, t_new );

      // Keep track of the new triangle
      t1 = &t_new;
    }

    // ****** Create second triangle *******
    Vertex& d1 = t1->v3();
    Vertex& d2 = t1->v2();

    Edge* diag = front_.get_edge(d1, d2);

    // --> created triangle t1 is not located on the front
    if ( !diag )
      return true;

    TriVector new_tris_p2 {};

    // Find all vertices in vicinity of p2 
    VertexVector vertex_candidates_p2 
      = find_local_vertices(p2, r2, wide_search );

    // Create potential triangles with all found vertices
    check_vertex_candidates(vertex_candidates_p2, *diag, new_tris_p2);

    // If potential triangles have been found, choose the best one
    // --> this function also updates the advancing front
    Triangle* t2 = choose_best_triangle(new_tris_p2, *diag);

    // No triangle has been found yet -> create new one
    if ( !t2 )
    {
      Vertex&   v_new = verts_.push_back( p2 );
      Triangle& t_new = tris_.push_back(d1, d2, v_new);

      // Algorithm fails if new vertex or new triangle is invalid
      if ( !check_vertex( v_new ) || !check_triangle( t_new) )
      {
        check_remove_triangle( t_new );
        check_remove_vertex( v_new );

        return true;
      }

      // Update the advancing front with new vertex
      update_front( *diag, v_new, t_new );

      // Keep track of the new triangle
      t2 = &t_new;
    }

    // ****** Merge the newly created triangles  *******
    // Gather vertices
    Vertex& q1 = t1->v1();
    Vertex& q2 = t1->v2();
    Vertex& q3 = t2->v3();
    Vertex& q4 = t2->v1();

    // Remove internal edge
    Edge* e_rem = intr_edges_.get_edge(q2, q4);
    check_remove_interior_edge( *e_rem );

    // Remove old triangular elements
    check_remove_triangle( *t1 );
    check_remove_triangle( *t2 );

    // Create new quadrilateral element
    Quad& q_new = quads_.push_back( q1, q2, q3, q4 );
    q_new.is_active( true );

    return true;

  } // Mesh::advance_front_quad()

  /*------------------------------------------------------------------
  | Let the advancing front create a new triangle 
  ------------------------------------------------------------------*/
  bool advance_front_triangle(Edge& base, bool wide_search=false)
  {
    // Obtain the position of a potential new vertex and the radius 
    // to search around it for potential existing vertices
    const double height = domain_->size_function( base.xy() );
    const Vec2d  v_xy   = base.xy() + height * base.normal();
    const double r      = domain_->size_function( v_xy );

    // Find all vertices in the vicinity of the current base edge
    VertexVector vertex_candidates 
      = find_local_vertices(v_xy, TQ_RANGE_FACTOR * r, wide_search);

    // Create potential triangles with all found vertices
    TriVector new_triangles {};
    check_vertex_candidates(vertex_candidates, base, new_triangles);

    // If potential triangles have been found, choose the best one
    if ( choose_best_triangle(new_triangles, base) )
      return true;

    // Check if a potential triangle can be created with the base edge
    // and a newly created vertex
    Vertex& b1 = base.v1();
    Vertex& b2 = base.v2();

    Vertex&   v_new = get_base_vertex( base );
    Triangle& t_new = tris_.push_back(b1, b2, v_new);
    
    // Algorithm fails if new vertex or new triangle is invalid
    if ( !check_vertex( v_new ) || !check_triangle( t_new) )
    {
      check_remove_triangle( t_new );
      check_remove_vertex( v_new );

      return false;
    }
    
    // Update the advancing front with new vertex
    update_front( base, v_new, t_new );

    return true;

  } // Mesh::advance_front_triangle()

  /*------------------------------------------------------------------
  | This function loops over all internal edges and, if possible,
  | merges two adjacent triangles to one quad element.
  |
  | -> This function requires a finalized mesh and that the functions
  |    setup_vertex_connectivity() and setup_facet_connectivity()
  |    have been called before.
  ------------------------------------------------------------------*/
  void merge_triangles_to_quads()
  {
    // Pick all internal edges, which are adjacent to two triangles
    std::list<Edge*> tri_edges {};

    for ( const auto& e_ptr : intr_edges_ )
    {
      Facet* f_l = e_ptr->facet_l();
      Facet* f_r = e_ptr->facet_r();

      if (  (f_l && f_l->n_vertices() == 3) 
         && (f_r && f_r->n_vertices() == 3) )
        tri_edges.push_back( e_ptr.get() );
    }

    // Sort edge list with increasing minimum angles of the 
    // adjacent triangles
    tri_edges.sort(
    []( Edge* a, Edge* b )
    {
      const double a_l = ( a->facet_l() )
                       ? a->facet_l()->min_angle() 
                       : TQ_MAX;
      const double a_r = ( a->facet_r() )
                       ? a->facet_r()->min_angle() 
                       : TQ_MAX;
      const double a_ang = MIN(a_l, a_r);

      const double b_l = ( b->facet_l() )
                       ? b->facet_l()->min_angle() 
                       : TQ_MAX;
      const double b_r = ( b->facet_r() )
                       ? b->facet_r()->min_angle() 
                       : TQ_MAX;
      const double b_ang = MIN(b_l, b_r);

      return a_ang < b_ang;

    });

    // Loop over all chosen edges and merge the adjacent triangles
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

      if ( f_l->n_vertices() != 3 || f_r->n_vertices() != 3 )
        continue;

      Vertex& v1 = e->v1();
      Vertex& v2 = e->v2();


      int i_l = f_l->get_edge_index(v1, v2);
      int i_r = f_r->get_edge_index(v1, v2);

      Vertex& q_l = f_l->vertex(i_l);
      Vertex& q_r = f_r->vertex(i_r);

      // Remove internal edge
      check_remove_interior_edge( *e );
      e = nullptr;

      // Create new quadrilateral element
      Quad& q_new = quads_.push_back( q_l, v1, q_r, v2 );
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

        Edge* e_share = intr_edges_.get_edge(q1,q2);

        if ( !e_share ) continue;

        Facet* t_l = e_share->facet_l();
        Facet* t_r = e_share->facet_r();

        if ( t_l && (t_l == f_l || t_l == f_r) )
          e_share->facet_l( &q_new );

        if ( t_r && (t_r == f_l || t_r == f_r) )
          e_share->facet_r( &q_new );
      }

      // Remove triangles
      check_remove_triangle( *(static_cast<Triangle*>(f_l)) );
      check_remove_triangle( *(static_cast<Triangle*>(f_r)) );

    }

    // Re-initialize facet-to-facet connectivity
    setup_facet_connectivity();

    // Bad elements may have been created up to this point
    // due to the merging of triangles to quads
    // The two upcoming function fix these bad elements
    clean_double_quad_edges();
    clean_double_triangle_edges();
    
    // Remove deleted entities
    clear_waste();

    // Re-initialize vertex-to-vertex connectivity
    setup_vertex_connectivity();

    // Re-initialize facet-to-facet connectivity
    setup_facet_connectivity();

  } // Mesh::merge_triangles_to_quads()


private:

  /*------------------------------------------------------------------
  | This function is part of the advancing front loop.
  | For a given location <xy> and a respective search range <dist>,
  | all vertices that are located in this vicinity are put into a 
  | vector and the sorted ascendingly towards <xy>.
  | If the flag <wide_search> is activated, the search range is 
  | enlarged by a preset factor, in order to finde more vertex 
  | candidates.
  ------------------------------------------------------------------*/
  VertexVector find_local_vertices(const Vec2d& xy, double dist, 
                                   bool wide_search=false )
  {
    if (wide_search)
      dist *= TQ_WIDE_SEARCH_FACTOR;

    // Get vertices in vicinity of xy  
    VertexVector vertex_candidates = verts_.get_items(xy, dist);

    // Sort vertices in ascending order towards xy
    std::sort( vertex_candidates.begin(), vertex_candidates.end(), 
    [xy] ( const Vertex* a, const Vertex* b )
    {
      return ( (a->xy()-xy).length_squared() 
             < (b->xy()-xy).length_squared() );
    });

    return std::move( vertex_candidates ); 

  } // find_local_vertices() 

  /*------------------------------------------------------------------
  | This function is part of the advancing front loop.
  | We loop over a given set of <vertex_candidates> and check
  | if we can create a possible triangle with the current base edge
  | (<b1>,<b2>) and the given vertices.
  | If it is possible, new triangles are created and pushed back to 
  | the vector <new_triangles>.
  ------------------------------------------------------------------*/
  void check_vertex_candidates(const VertexVector& vertex_candidates,
                               Edge& base, TriVector& new_triangles)
  {
    for ( Vertex* v : vertex_candidates )
    {
      // Skip vertices that are not located on the advancing front
      if ( !v->on_front() )
        continue;

      // Skip vertices that are colinear to the current base edge
      if ( TQGeom::orientation( base.v1().xy(), base.v2().xy(), v->xy() )
          == TQGeom::Orientation::CL )
        continue;

      // Create new potential triangle 
      Triangle& t_new = tris_.push_back( base.v1(), base.v2(), *v );

      // Check if new potential triangle is valid
      if ( check_triangle( t_new ) )
        new_triangles.push_back( &t_new );
      else
        check_remove_triangle( t_new );

    }

  } // check_vertex_candidates()

  /*------------------------------------------------------------------
  | This function is part of the advancing front loop.
  | We sort a given vector of <new_triangles> in descending order
  | according to the triangle quality.
  | Finally, the advancing front is updated with the triangle of best
  | quality and all other triangles are removed.
  ------------------------------------------------------------------*/
  Triangle* choose_best_triangle(TriVector& new_triangles,
                                 Edge&      base)
  {
    if ( new_triangles.size() < 1 )
      return nullptr;

    DBG_MSG("VALID TRIANGLES IN NEIGHBORHOOD: " 
            << new_triangles.size());

    std::sort( new_triangles.begin(), new_triangles.end(),
    [this] ( Triangle* t1, Triangle* t2 )
    {
      const double h1 = domain_->size_function( t1->xy() );
      const double h2 = domain_->size_function( t2->xy() );
      const double q1 = t1->quality(h1);
      const double q2 = t2->quality(h2);

      return ( q1 > q2 );
    });

    Triangle* new_tri = new_triangles[0];
    Vertex&   v_adj   = new_tri->v3();

    update_front( base, v_adj, *new_tri );

    for (int i = 1; i < new_triangles.size(); i++)
      tris_.remove( *new_triangles[i] );

    return new_tri;

  } // choose_best_triangle()

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  bool check_triangle(const Triangle& tri)
  {
    const double rho   = domain_->size_function( tri.xy() );
    const double range = 2.0 * rho;

    DBG_MSG("CHECK NEW TRIANGLE: " << tri);

    if ( !tri.is_valid() )
      return false;

    if ( tri.intersects_front( front_, range ) )
    {
      DBG_MSG("  > FRONT INTERSECTION");
      return false;
    }

    if ( tri.intersects_vertex( verts_, range ) )
    {
      DBG_MSG("  > VERTEX INTERSECTION");
      return false;
    }

    if ( tri.intersects_triangle( tris_, range ) )
    {
      DBG_MSG("  > TRIANGLE INTERSECTION");
      return false;
    }

    if ( tri.intersects_quad( quads_, range ) )
    {
      DBG_MSG("  > QUAD INTERSECTION");
      return false;
    }

    DBG_MSG("  > VALID");

    return true;

  } // check_triangle()

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  bool check_vertex(const Vertex& v)
  {
    const double rho   = domain_->size_function( v.xy() );
    const double range = 2.0 * rho;

    DBG_MSG("CHECK NEW VERTEX: " << v);

    if ( !domain_->is_inside( v ) )
    {
      DBG_MSG("  > OUTSIDE DOMAIN");
      return false;
    }

    if ( v.intersects_facet(tris_, range) )
    {
      DBG_MSG("  > TRIANGLE INTERSECTION");
      return false;
    }

    if ( v.intersects_facet(quads_, range) )
    {
      DBG_MSG("  > QUAD INTERSECTION");
      return false;
    }

    DBG_MSG("  > VALID");

    return true;

  } // check_vertex()

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  Vertex& get_base_vertex(const Edge& base)
  {
    // Half of the factor h for height of equlateral triangle
    // h := sqrt(3) / 2  -   h_fac := h / 2
    constexpr double h_fac = 0.4330127019; 

    // Obtain size function value at the centroid of an equlateral
    // triangle, created from the current base edge
    Vec2d c = base.xy() + base.normal() * base.length() * h_fac;
    const double rho = domain_->size_function(c);

    // Coordinate of new vertex 
    Vec2d xy = base.xy() + base.normal() * rho;

    return verts_.push_back( xy.x, xy.y );

  } // get_base_vertex() 

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  void update_front( Edge& base, Vertex& v_new, Triangle& t_new )
  {
    // Get advancing front edges adjacent to vertex
    // -> First two vertices of new triangle tri are always
    //    the base edge vertices
    Edge* e1 = front_.get_edge(v_new, base.v1());
    Edge* e2 = front_.get_edge(v_new, base.v2());

    // *** Both edges are connected to vertex ***
    //     -> No new edge must be created
    //     -> Base vertex v1 no longer on the advancing front
    //     -> Base vertex v2 no longer on the advancing front
    //     -> vertex no longer on the advancing front
    if ( e1 && e2 )
    {
      ASSERT( v_new.on_front(), "Grid structure corrupted." );

      base.v1().on_front( false );
      base.v2().on_front( false );
      v_new.on_front( false );

      if ( e1->is_interior() )
        intr_edges_.add_edge(e1->v1(), e1->v2());
      if ( e2->is_interior() )
        intr_edges_.add_edge(e2->v1(), e2->v2());

      front_.remove( *e1 );
      front_.remove( *e2 );
    }
    // *** First edge is connected to vertex ***
    //     -> New edge between second base vertex and vertex
    //     -> Base vertex v1 no longer on the advancing front
    //     -> vertex is part of the advancing front
    else if ( e1 && !e2 )
    {
      ASSERT( v_new.on_front(), "Grid structure corrupted." );

      base.v1().on_front( false );

      if ( e1->is_interior() )
        intr_edges_.add_edge(e1->v1(), e1->v2());

      front_.remove( *e1 );
      front_.add_edge(v_new, base.v2());
    }
    // *** Second edge is connected to vertex ***
    //     -> New edge between first base vertex and vertex
    //     -> Base vertex v2 no longer on the advancing front
    //     -> vertex is part of the advancing front
    else if ( !e1 && e2 )
    {
      ASSERT( v_new.on_front(), "Grid structure corrupted." );

      base.v2().on_front( false );

      if ( e2->is_interior() )
        intr_edges_.add_edge(e2->v1(), e2->v2());

      front_.remove( *e2 );
      front_.add_edge(base.v1(), v_new);
    }
    // *** Both edges are not connected to vertex ***
    //     -> Create two new edges
    //     -> Vertex now part of the advancing front
    else
    {
      v_new.on_front( true );
      front_.add_edge(base.v1(), v_new);
      front_.add_edge(v_new, base.v2());
    }

    // If current base is not at the boundary, add it to the 
    // interior edge list
    if ( base.is_interior() )
      intr_edges_.add_edge(base.v1(), base.v2());

    // Remove base edge
    front_.remove( base );

    // Mark new triangle as active
    t_new.is_active( true );

    // Add element area to the total mesh area
    mesh_area_ += t_new.area();

  } // update_front() 

  /*------------------------------------------------------------------
  | Every vertex gets assigned its neighboring vertices and these 
  | are then sorted by means of ascending angles
  | This function requires, that all the interior edges and 
  | boundary edges of the mesh have been generated 
  | -> element adjacency of the edges is not required here
  ------------------------------------------------------------------*/
  void setup_vertex_connectivity()
  {
    // Remove all current vertex-to-vertex connectivities
    for ( auto& v_ptr : verts_ )
      v_ptr->vertices().clear();

    // Get vertex-to-vertex connectivity from interior edges
    for ( const auto& e : intr_edges_ )
    {
      Vertex* v1 = &(e->v1());
      Vertex* v2 = &(e->v2());

      v1->vertices().push_back( v2 );
      v2->vertices().push_back( v1 );
    }

    // Get vertex-to-vertex connectivity from boundary edges
    for ( const auto& e : bdry_edges_ )
    {
      Vertex* v1 = &(e->v1());
      Vertex* v2 = &(e->v2());

      v1->vertices().push_back( v2 );
      v2->vertices().push_back( v1 );
    }

    // Sort local vertex connectivities by ascending angle
    for ( auto& v_ptr : verts_ )
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

  } // Mesh::setup_vertex_connectivity()

  /*------------------------------------------------------------------
  | Initialize the connectivity between facets and facets, as well  
  | as between edges and facets
  ------------------------------------------------------------------*/
  void setup_facet_connectivity()
  {
    Facet* f1 = nullptr;
    Facet* f2 = nullptr;
    
    // Setup connectivity for interor edges
    for ( const auto& e_ptr : intr_edges_ )
    {
      const Vertex& v1 = e_ptr->v1();
      const Vertex& v2 = e_ptr->v2();

      int idx1 {-1};
      int idx2 {-1};

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

      // Setup connectivity between facets
      f1->neighbor( idx1, f2 );
      f2->neighbor( idx2, f1 );

      // Setup connectivity between internal edge and facets
      if ( TQGeom::is_left( v1.xy(), v2.xy(), f1->xy() ) )
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
    for ( const auto& e_ptr : bdry_edges_ )
    {
      const Vertex& v1 = e_ptr->v1();
      const Vertex& v2 = e_ptr->v2();

      for ( auto f : v1.facets() )
      {
        int idx = f->get_edge_index(v1, v2);

        if ( idx < 0 ) continue;

        f1 = f;
        break;
      }

      // Setup connectivity between internal edge and facets
      if ( TQGeom::is_left( v1.xy(), v2.xy(), f1->xy() ) )
        e_ptr->facet_l( f1 );

    }
    
  } // Mesh::setup_facet_connectivity

  /*------------------------------------------------------------------
  | Clean up quad elements:
  | It may be, that some adjacent quads share two internal edges
  | These elements will be merged in the next step.
  |
  | -> This function requires a finalized mesh and that the functions
  |    setup_vertex_connectivity() and setup_facet_connectivity()
  |    have been called before.
  | 
  |     v3                       vp   
  |       x---------------------x
  |       | \                   |
  |       |   \ e2     q_nbr    |
  |       |     \               |
  |       |       \  v2         |
  |       |         x           |
  |       |           \         |
  |       |   q_cur     \  e1   |
  |       |               \     |
  |       |                 \   |
  |       |                   \ |
  |       x---------------------x
  |     v4                       v1
  ------------------------------------------------------------------*/
  void clean_double_quad_edges()
  {
    // Initialize all quad colors
    for ( auto& q_cur : quads_ )
      q_cur->color(0);

    std::vector<std::pair<Quad*,Quad*>>     quads_to_remove;
    std::vector<std::pair<Edge*,Edge*>>     edges_to_remove;
    std::vector<Vertex*>                    verts_to_remove;
    std::vector<std::pair<Vertex*,Vertex*>> opposing_vertices;

    for ( auto& q_cur : quads_ )
    {
      if ( q_cur->color() > 0 )
        continue;

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

        // Proceed, if no neighbors are found (nullptr) or if 
        // neighbors of the adjacent edges differ, or if the neighbor
        // has already been added
        if ( !nbr_1 || !nbr_2 || nbr_1 != nbr_2 || nbr_1->color() > 0)
          continue;

        // In this stage, we address only quad / quad connections
        if ( nbr_1->n_vertices() < 4 )
          continue;

        // Now we can cast the facet to a quad
        Quad* q_nbr = static_cast<Quad*>(nbr_1);

        // Get the internal edges adjacent to both current quads
        Edge* e1 = intr_edges_.get_edge(v1, v2);
        Edge* e2 = intr_edges_.get_edge(v2, v3);

        ASSERT( (e1 != e2), "INVALID DATA STRUCTURE");
        ASSERT( ( e1->v1() == e2->v1() || e1->v1() == e2->v2() 
               || e1->v2() == e2->v1() || e1->v2() == e2->v2() ),
            "WRONG EDGES FOUND.");

        // Color the current quads, such that the won't get chosen
        // in upcoming loops
        q_cur->color(1);
        q_nbr->color(1);

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

        ASSERT( (v_op != v1), "BAD DATA STRUCTURE");
        ASSERT( (v_op != v2), "BAD DATA STRUCTURE");
        ASSERT( (v_op != v3), "BAD DATA STRUCTURE");

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
      Quad& q_new = quads_.push_back( *o1, v1, *o2, v3 );
      q_new.is_active(true);

    }

    // Removal of old quads
    for (size_t i = 0; i < quads_to_remove.size(); ++i)
    {
      Quad* q1 = quads_to_remove[i].first;
      Quad* q2 = quads_to_remove[i].second;

      check_remove_quad( *q1 );
      check_remove_quad( *q2 );
    }

    // Removal of old interior edges
    for (size_t i = 0; i < edges_to_remove.size(); ++i)
    {
      Edge* e1 = edges_to_remove[i].first;
      Edge* e2 = edges_to_remove[i].second;

      check_remove_interior_edge( *e1 );
      check_remove_interior_edge( *e2 );
    }

    // Removal of old vertices
    for (size_t i = 0; i < verts_to_remove.size(); ++i)
    {
      Vertex* v = verts_to_remove[i];
      check_remove_vertex( *v );
    }

    // Re-initialize facet-to-facet connectivity
    setup_facet_connectivity();

  } // Mesh::clean_double_quad_edges()

  /*------------------------------------------------------------------
  | Clean up triangle elements:
  | It may be, that some triangles share two internal edges with a 
  | single quad element. 
  | These triangles and quads will be removed in this step and then
  | replaced by a single triangle.
  |
  | -> This function requires a finalized mesh and that the functions
  |    setup_vertex_connectivity() and setup_facet_connectivity()
  |    have been called before.
  | -> This function should be called after the function 
  |    clean_double_quad_edges(), because such triangles might be 
  |    genereated during the latter function.
  |
  |        v3
  |       x
  |       | \
  |       |\  \
  |       | \   \
  |       |  \    \
  |       |   \     \
  |       |    \      \
  |       |     \ t_nbr \
  |       |      \        \
  |       |       x--.      \
  |       |      v2   --.     \
  |       |              --.    \
  |       |   q_cur         --.   \
  |       |                    ---  \
  |       x--------------------------x
  |     v4                          v1
  |
  ------------------------------------------------------------------*/
  void clean_double_triangle_edges()
  {
    // Initialize all quad colors
    for ( auto& q_cur : quads_ )
      q_cur->color(0);

    std::vector<std::pair<Quad*,Triangle*>> elements_to_remove;
    std::vector<std::pair<Edge*,Edge*>>     edges_to_remove;
    std::vector<Vertex*>                    verts_to_remove;

    for ( auto& q_cur : quads_ )
    {
      if ( q_cur->color() > 0 )
        continue;

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

        // Proceed, if no neighbors are found (nullptr) or if 
        // neighbors of the adjacent edges differ, or if the neighbor
        // has already been added
        if ( !nbr_1 || !nbr_2 || nbr_1 != nbr_2 || nbr_1->color() > 0)
          continue;

        // In this stage, we address only quad / triangle connections
        if ( nbr_1->n_vertices() > 3 )
          continue;

        // Now we can cast the facet to a triangle
        Triangle* t_nbr = static_cast<Triangle*>(nbr_1);

        // Get the internal edges adjacent to both current elements
        Edge* e1 = intr_edges_.get_edge(v1, v2);
        Edge* e2 = intr_edges_.get_edge(v2, v3);

        ASSERT( (e1 != e2), "INVALID DATA STRUCTURE");
        ASSERT( ( e1->v1() == e2->v1() || e1->v1() == e2->v2() 
               || e1->v2() == e2->v1() || e1->v2() == e2->v2() ),
            "WRONG EDGES FOUND.");

        // Color the current quads, such that the won't get chosen
        // in upcoming loops
        q_cur->color(1);
        t_nbr->color(1);

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

      ASSERT( (id_v > -1), "BAD DATA STRUCTURE" );

      // Get the remaining three vertices of the quad 
      Vertex& v1 = q->vertex( MOD(id_v+1, 4) );
      Vertex& v2 = q->vertex( MOD(id_v+2, 4) );
      Vertex& v3 = q->vertex( MOD(id_v+3, 4) );

      // Create new triangle 
      Triangle& t_new = tris_.push_back( v1, v2, v3 );
      t_new.is_active(true);
    }

    // Removal of old elements
    for ( size_t i = 0; i < elements_to_remove.size(); ++i )
    {
      Quad*     q = elements_to_remove[i].first;
      Triangle* t = elements_to_remove[i].second;
      
      check_remove_quad( *q );
      check_remove_triangle( *t );
    }

    // Removal of old interiord edges
    for (size_t i = 0; i < edges_to_remove.size(); ++i)
    {
      Edge* e1 = edges_to_remove[i].first;
      Edge* e2 = edges_to_remove[i].second;

      check_remove_interior_edge( *e1 );
      check_remove_interior_edge( *e2 );
    }

    // Removal of old vertices
    for (size_t i = 0; i < verts_to_remove.size(); ++i)
    {
      Vertex* v = verts_to_remove[i];
      check_remove_vertex( *v );
    }

    // Re-initialize facet-to-facet connectivity
    setup_facet_connectivity();

  } // Mesh::clean_double_triangle_edges()

  /*------------------------------------------------------------------
  | Algorithm to create a single quad layer for a connected list of
  | advancing front edges that start with v_start and end with v_end
  ------------------------------------------------------------------*/
  bool add_quad_layer(Vertex*& v_start_in, Vertex*& v_end_in, 
                      double height)
  {
    // Find closest vertices to given input vertices
    Vertex* v_start = nullptr;
    Vertex* v_end = nullptr;
    double d2_start_min = 1.0E+10;
    double d2_end_min = 1.0E+10;

    for ( const auto& v : verts_ )
    {
      double d2_start = (v_start_in->xy() - v->xy()).length_squared();
      double d2_end   = (v_end_in->xy() - v->xy()).length_squared();

      if (d2_start < d2_start_min)
      {
        v_start = v.get();
        d2_start_min = d2_start;
      }

      if (d2_end < d2_end_min)
      {
        v_end = v.get();
        d2_end_min = d2_end;
      }
    }

    if (!v_start || !v_end)
    {
      MSG("[ERROR]: Failed to create quad layer. "
        "Provided starting or ending vertex pointer is NULL.");
      return false;
    }

    // Get advancing front edges adjacent to input vertices
    Edge* e_start = front_.get_edge(*v_start, 1); 
    Edge* e_end   = front_.get_edge(*v_end, 2); 

    if ( !e_start || !e_end )
    {
      MSG("Failed to find an advancing front edge that is adjacent"
          " to the give input vertex for the quad layer creation.");
      return false;
    }

    // Check if given front segments can be traversed and if closed
    bool is_closed = (v_start == v_end);

    if ( !front_.is_traversable(*e_start, *e_end) )
    {
      MSG("Failed to traverse the advancing front "
          "for the given input vertices.");
      return false;
    }

    // For closed quad layers, try not to start at sharp angle edges
    if ( is_closed )
    {
      const Vec2d& v1 = e_end->v1().xy();
      const Vec2d& v2 = e_end->v2().xy();
      const Vec2d& v3 = e_start->v2().xy();

      const double ang = angle(v1-v2, v3-v2);

      Edge *e_next = e_start->get_next_edge();

      if ( e_next && ang <= TQ_QUAD_LAYER_ANGLE )
      {
        e_end = e_start;
        e_start = e_next; 
      }
    }

    // Triangulate front edges with critical angles
    //prepare_quad_layer_front(e_start, e_end, height);

    // Create the quad layer structure, which keeps track of the target
    // vertex coordinates, that are projected from the base vertex 
    // coordinates
    QuadLayer quad_layer { e_start, e_end, is_closed, height };
    quad_layer.smooth_heights( *domain_ );
    quad_layer.setup_vertex_projection( verts_, front_, bdry_edges_ );

    // For each base edge in the quad layer, try to create a quad
    // element with its given projected coordinates
    create_quad_layer_elements( quad_layer );

    // Triangulate the quad layer based edges, where the generation
    // of quads did not succeed
    finish_quad_layer( quad_layer );

    // Remove deleted entities
    clear_waste();

    // Set pointers to new start and ending vertices
    int i = 0;
    int n = quad_layer.n_bases();

    do 
    {
      v_start_in = quad_layer.p1()[i];

      if ( is_closed )
        v_end_in = v_start_in;
      else
        v_end_in = quad_layer.p2()[MOD(i-1,n)]; 

      ++i;

    } while ( !(v_start_in->on_front()) && !(v_end_in->on_front()) );

    return true;

  } // add_quad_layer()

  /*------------------------------------------------------------------
  | For each QuadProjection, create a triangle with its base 
  | vertices (b1,b2) and a vertex p1, which is either located in 
  | the vicinity of the base edge or which is otherwise generated at 
  | the projected coordinate of the base vertex b1
  |   
  |           p1            p2
  |          x-------------x-------------
  |          | \           | \          |
  |          |   \         |   \        |
  |          |     \       |     \      |
  |          |       \     |       \    |
  |          |         \   |         \  |
  |          |    base   \ |           \|
  | ---------x-------------x------------x-------
  |           b1            b2
  |   
  ------------------------------------------------------------------*/
  void create_quad_layer_elements(QuadLayer& quad_layer)
  {
    auto& b1        = quad_layer.b1();
    auto& b2        = quad_layer.b2();

    auto& p1        = quad_layer.p1();
    auto& p2        = quad_layer.p2();

    auto& p1_xy     = quad_layer.p1_xy();
    auto& p2_xy     = quad_layer.p2_xy();

    auto& heights   = quad_layer.heights();
    auto& bases     = quad_layer.bases();

    int  n_bases   = quad_layer.n_bases();

    for ( int i = 0; i < n_bases; ++i )
    {
      // Search radius for vertices in the vicinity of the 
      // projected coordinates
      const double r = TQ_QUAD_LAYER_RANGE * heights[i];

      // Create first triangle (b1,b2,p1)
      Edge* base = bases[i];

      if (!base->in_container())
        continue;

      Triangle* t1 = add_quad_layer_triangle(base, p1_xy[i], r);

      if ( t1 ) 
        p1[i] = &(t1->v3());
      else
        continue;

      // Create second triangle (p1,b2,p2)
      base = front_.get_edge( *p1[i], *b2[i] );

      if ( !base ) 
        continue;

      Triangle* t2 = add_quad_layer_triangle(base, p2_xy[i], r);

      if ( t2 )
        p2[i] = &(t2->v3());
      else
        continue;

      // Merge both triangles t1 & t2 to a quad
      // --> First remove the interior edge between these triangles
      Edge* e_rem = intr_edges_.get_edge( *b2[i], *p1[i] );

      if ( e_rem ) 
        check_remove_interior_edge( *e_rem );
      else
        continue;

      // Remove old triangular elements
      check_remove_triangle( *t1 );
      check_remove_triangle( *t2 );
      t1 = nullptr;
      t2 = nullptr;

      // Create new quadrilateral element
      Quad& q_new = quads_.push_back( *b1[i], *b2[i], *p2[i], *p1[i] );
      q_new.is_active( true );

    }

  } // Mesh::create_quad_layer_elements() 
  
  /*------------------------------------------------------------------
  | This is a helper function for the generation of a triangle 
  | during the quad layer generation
  ------------------------------------------------------------------*/
  Triangle* add_quad_layer_triangle(Edge* base, 
                                    const Vec2d& xy, double r)
  {
    TriVector new_tris {};

    Triangle* tri = nullptr;

    Vertex& v1 = base->v1();
    Vertex& v2 = base->v2();

    // Look for vertices in the vicinity of the projected coordinate 
    // xy, which could be used to construct a new triangle
    VertexVector vertex_candidates = find_local_vertices(xy, r);

    // Create potential triangles with all found vertices
    check_vertex_candidates(vertex_candidates, *base, new_tris);

    // If potential triangles have been found, use the best one
    tri = choose_best_triangle(new_tris, *base);

    // If no proper triangle has been found, create a new vertex
    // at the projected position xy and create a new triangle 
    // with this new vertex
    if ( !tri )
    {
      Vertex& v_new = verts_.push_back( xy );
      tri = &( tris_.push_back(v1, v2, v_new) );

      // If the created entities are invalid, clean up
      if ( !check_vertex( v_new ) || !check_triangle( *tri ) )
      {
        check_remove_triangle( *tri );
        check_remove_vertex( v_new );
        tri = nullptr;
      }
      // Otherwise, update the advancing front with the new triangle
      else
      {
        update_front( *base, v_new, *tri );
        v_new.is_fixed( true );
      }
    }

    return tri;

  } // Mesh::add_quad_layer_triangle()

  /*------------------------------------------------------------------
  | In some cases, gaps might be formed during the previous quad layer
  | generation steps. In this function, these gaps are closed with 
  | triangular elements.
  |
  |              p1[i]
  |      v      x 
  |     x       :
  |             :
  |  p2[i-1]    :  
  |   x.........x-------------x
  |             | b1[i]        b2[i]
  |             |           
  |             |
  |             |
  |             x
  |               
  ------------------------------------------------------------------*/
  void finish_quad_layer(QuadLayer& quad_layer)
  {
    auto& b1        = quad_layer.b1();

    auto& p1        = quad_layer.p1();
    auto& p2        = quad_layer.p2();

    int  n_bases   = quad_layer.n_bases();

    for ( int i = 1; i < n_bases; ++i )
    {
      if ( !p1[i] || !p2[i-1] || p1[i] == p2[i-1] )
        continue;

      Vertex& a = *p2[i-1];
      Vertex& b = *b1[i];
      Vertex& c = *p1[i];

      const Vec2d l1 = a.xy()-b.xy();
      const Vec2d l2 = c.xy()-b.xy();
      const double alpha = angle(l1,l2);

      // Don't add a new vertex and instead only connect (a,b,c)
      if ( alpha <= TQ_QUAD_LAYER_ANGLE )
      {
        Triangle* t_new = &( tris_.push_back(a, b, c) );

        if ( !check_triangle( *t_new ) )
        {
          check_remove_triangle( *t_new );
        }
        else
        {
          Edge* base = front_.get_edge( b, c );
          update_front( *base, a, *t_new );
        }
      }
      // Create new vertex and then generate two triangles
      else
      {
        const Vec2d v_xy = b.xy() + l1 + l2;

        Vertex& v_new = verts_.push_back( v_xy );

        Triangle* t1_new = &( tris_.push_back(a, b, v_new) );
        Triangle* t2_new = &( tris_.push_back(b, c, v_new) );

        if (  !check_vertex( v_new ) 
           || !check_triangle( *t1_new ) || !check_triangle( *t2_new ) )
        {
          check_remove_triangle( *t1_new );
          check_remove_triangle( *t2_new );
          check_remove_vertex( v_new );
        }
        else
        {
          Edge* base = nullptr;

          base = front_.get_edge( a, b );
          update_front( *base, v_new, *t1_new );

          base = front_.get_edge( b, c );
          update_front( *base, v_new, *t2_new );

          v_new.is_fixed( true );
        }
      }
    }

  } // Mesh::finish_quad_layer()

  /*------------------------------------------------------------------
  | Triangulate front edges with critical angles
  | 
  |                      p2
  |                      x
  |              v       |
  |             o        |
  |               .      | e_next 
  |                 .    |
  |                   .  |
  |      p1             .|
  |     x----------------x
  |           e_cur       c
  | 
  ------------------------------------------------------------------*/
  void prepare_quad_layer_front(Edge*& e_start, Edge*& e_end, 
                                const double height)
  {
    Edge* e_cur  = e_start;

    Edge* e_last = e_end->get_next_edge();
    Edge* e_prev = e_start->get_prev_edge();

    int edge_count = 0;
    int n_edges = static_cast<int>(front_.size());

    do 
    {
      Edge* e_next = e_cur->get_next_edge();

      ASSERT( e_next, "INVALID DATA STRUCTURE" );
      if ( !e_next )
        break;

      Vertex& p1 = e_cur->v1();
      Vertex&  c = e_cur->v2();
      Vertex& p2 = e_next->v2();

      const Vec2d d1 = p1.xy() - c.xy();
      const Vec2d d2 = p2.xy() - c.xy();

      const double alpha = angle(d1,d2);
      const double delta = ( p2.xy() - p1.xy() ).length();

      const double l = 0.5 * ( e_cur->length() + e_next->length() );
      const double h = MIN( l, height );

      if (   !( TQGeom::is_lefton(p1.xy(), p2.xy(), c.xy()) ) 
          && delta <= 4.0 * h  
          && alpha <= TQ_QUAD_LAYER_ANGLE )
      {
        Edge* e_buf = e_next->get_next_edge();

        if ( e_next == e_last )
          e_last = e_buf;

        if ( delta <= h )
        {
          const Vec2d v_xy = 0.5 * (p1.xy() + p2.xy());
          Vertex& v_new = verts_.push_back( v_xy );

          Triangle& t1_new = tris_.push_back(p1, c, v_new);
          Triangle& t2_new = tris_.push_back(c, p2, v_new);

          if (  !check_vertex( v_new ) 
             || !check_triangle( t1_new )
             || !check_triangle( t2_new ) )
          {
            check_remove_triangle( t1_new );
            check_remove_triangle( t2_new );
            check_remove_vertex( v_new );
          }
          else
          {
            update_front( *e_cur, v_new, t1_new );
            update_front( *e_next, v_new, t2_new );

            v_new.is_fixed( true );
          }
        }
        else
        {
          Triangle& t_new = tris_.push_back(p1, c, p2);

          if ( !check_triangle( t_new ) )
            check_remove_triangle( t_new );
          else
            update_front( *e_cur, p2, t_new );
        }

        e_cur = e_buf;
        ++edge_count;
      }
      else
      {
        e_cur = e_next;
        ++edge_count;
      }

    } while(  ( e_cur )
           && ( edge_count < n_edges )
           && ( e_cur != e_last ) );

    e_start = e_prev->get_next_edge();
    e_end   = e_last->get_prev_edge();

  } // Mesh::prepare_quad_layer_front()

  /*------------------------------------------------------------------
  | Get edge 
  ------------------------------------------------------------------*/
  Edge* get_edge(const Vertex& v1, const Vertex& v2, bool dir=false)
  const 
  {
    Edge* found = nullptr;

    found = intr_edges_.get_edge(v1, v2, dir);

    if ( found != nullptr )
      return found;

    found = bdry_edges_.get_edge(v1, v2, dir);

    if ( found != nullptr )
      return found;

    found = front_.get_edge(v1, v2, dir);

    if ( found != nullptr )
      return found;

    return found;

  } // Mesh::get_edge()


  /*------------------------------------------------------------------
  | This function removes an entity and makes sure, that the 
  | removal succeeded
  ------------------------------------------------------------------*/
  inline void check_remove_vertex(Vertex& v)
  {
    bool removed = verts_.remove( v );
    ASSERT( removed, "Failed to remove vertex.");
    (void) removed;
  }

  inline void check_remove_triangle(Triangle& t)
  {
    bool removed = tris_.remove( t );
    ASSERT( removed, "Failed to remove triangle.");
    (void) removed;
  }

  inline void check_remove_quad(Quad& q)
  {
    bool removed = quads_.remove( q );
    ASSERT( removed, "Failed to remove quad.");
    (void) removed;
  }

  inline void check_remove_interior_edge(Edge& e)
  {
    bool removed = intr_edges_.remove( e );
    ASSERT( removed, "Failed to remove interior edge.");
    (void) removed;
  }

  inline void check_remove_boundary_edge(Edge& e)
  {
    bool removed = bdry_edges_.remove( e );
    ASSERT( removed, "Failed to remove interior edge.");
    (void) removed;
  }


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Domain*    domain_ {nullptr};

  Vertices   verts_;
  Triangles  tris_;
  Quads      quads_;
  Front      front_;

  EdgeList   intr_edges_ { TQGeom::Orientation::NONE };
  EdgeList   bdry_edges_ { TQGeom::Orientation::NONE };

  double     mesh_area_ { 0.0 };

}; // Mesh


/*********************************************************************
* Print out the mesh to std::cout
*********************************************************************/
static inline std::ostream& operator<<(std::ostream& os, 
                                       const Mesh& mesh)
{
  os << "VERTICES " << mesh.vertices().size() << "\n";
  for ( const auto& v_ptr : mesh.vertices() )
  {
    os << std::setprecision(5) << std::fixed 
              << v_ptr->xy().x << "," 
              << v_ptr->xy().y << "\n";
  }

  os << "INTERIOREDGES " << mesh.interior_edges().size() << "\n";
  for ( const auto& e_ptr : mesh.interior_edges() )
    os << std::setprecision(0) << std::fixed 
      << std::setw(4) << e_ptr->v1().index() << "," 
      << std::setw(4) << e_ptr->v2().index() << ","
      << std::setw(4) << e_ptr->facet_l()->index() << ","
      << std::setw(4) << e_ptr->facet_r()->index() << "\n";

  os << "BOUNDARYEDGES " << mesh.boundary_edges().size() << "\n";
  for ( const auto& e_ptr : mesh.boundary_edges() )
    os << std::setprecision(0) << std::fixed 
      << std::setw(4) << e_ptr->v1().index() << "," 
      << std::setw(4) << e_ptr->v2().index() << ","
      << std::setw(4) << e_ptr->facet_l()->index() << ","
      << std::setw(4) << e_ptr->marker() << "\n";

  os << "FRONT " << mesh.front().size() << "\n";
  for ( const auto& e_ptr : mesh.front() )
    os << std::setprecision(0) << std::fixed 
      << std::setw(4) << e_ptr->v1().index() << "," 
      << std::setw(4) << e_ptr->v2().index() << ","
      << std::setw(4) << e_ptr->marker() << "\n";

  os << "QUADS " << mesh.quads().size() << "\n";
  for ( const auto& q_ptr : mesh.quads() )
  {
    os << std::setprecision(0) << std::fixed
      << std::setw(4) << q_ptr->v1().index() << ","
      << std::setw(4) << q_ptr->v2().index() << ","
      << std::setw(4) << q_ptr->v3().index() << ","
      << std::setw(4) << q_ptr->v4().index() << "\n";
  }

  os << "TRIANGLES " << mesh.triangles().size() << "\n";
  for ( const auto& t_ptr : mesh.triangles() )
  {
    os << std::setprecision(0) << std::fixed
      << std::setw(4) << t_ptr->v1().index() << ","
      << std::setw(4) << t_ptr->v2().index() << ","
      << std::setw(4) << t_ptr->v3().index() << "\n";
  }

  os << "QUADNEIGHBORS " << mesh.quads().size() << "\n";
  for ( const auto& q_ptr : mesh.quads() )
  {
    os << std::setprecision(0) << std::fixed << std::setw(4) 
      << ( q_ptr->nbr1() ? q_ptr->nbr1()->index() : -1 ) 
      << "," << std::setw(4) 
      << ( q_ptr->nbr2() ? q_ptr->nbr2()->index() : -1 ) 
      << "," << std::setw(4) 
      << ( q_ptr->nbr3() ? q_ptr->nbr3()->index() : -1 ) 
      << "," << std::setw(4) 
      << ( q_ptr->nbr4() ? q_ptr->nbr4()->index() : -1 ) 
      << "\n";
  }

  os << "TRIANGLENEIGHBORS " << mesh.triangles().size() << "\n";
  for ( const auto& t_ptr : mesh.triangles() )
  {
    os << std::setprecision(0) << std::fixed << std::setw(4) 
      << ( t_ptr->nbr1() ? t_ptr->nbr1()->index() : -1 ) 
      << "," << std::setw(4) 
      << ( t_ptr->nbr2() ? t_ptr->nbr2()->index() : -1 ) 
      << "," << std::setw(4) 
      << ( t_ptr->nbr3() ? t_ptr->nbr3()->index() : -1 ) 
      << "\n";
  }


  return os;
} 

} // namespace TQAlgorithm
} // namespace TQMesh
