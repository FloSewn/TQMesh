/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "STLHeaders.h"

#include "Vertex.h"
#include "Edge.h"
#include "Triangle.h"
#include "Front.h"
#include "Domain.h"
#include "Mesh.h"


namespace TQMesh {

using namespace CppUtils;

/*********************************************************************
* 
*********************************************************************/
class FrontUpdate
{
public:

  using VertexVector   = std::vector<Vertex*>;
  using TriVector      = std::vector<Triangle*>;

  /*------------------------------------------------------------------
  | Constructor / Destructor
  ------------------------------------------------------------------*/
  FrontUpdate(Mesh& mesh, const Domain& domain, Front& front)
  : mesh_ { mesh }, domain_ { domain }, front_ { front } {}

  ~FrontUpdate() {}

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  double min_cell_quality() const { return min_cell_quality_; }
  double max_cell_angle() const { return max_cell_angle_; }

  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  void min_cell_quality(double v) { min_cell_quality_ = v; }
  void max_cell_angle(double v) { max_cell_angle_ = v; }

  /*------------------------------------------------------------------
  | Let the front advance  
  ------------------------------------------------------------------*/
  void advance_front(Edge& base, Vertex& v_new, Triangle& t_new)
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

      if ( e1->is_interior() )
        mesh_.add_interior_edge(e1->v1(), e1->v2());
      if ( e2->is_interior() )
        mesh_.add_interior_edge(e2->v1(), e2->v2());

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

      if ( e1->is_interior() )
        mesh_.add_interior_edge(e1->v1(), e1->v2());

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

      if ( e2->is_interior() )
        mesh_.add_interior_edge(e2->v1(), e2->v2());

      front_.remove( *e2 );
      front_.add_edge(base.v1(), v_new);
    }
    // *** Both edges are not connected to vertex ***
    //     -> Create two new edges
    //     -> Vertex now part of the advancing front
    else
    {
      front_.add_edge(base.v1(), v_new);
      front_.add_edge(v_new, base.v2());
    }

    update_front_state(base.v1());
    update_front_state(base.v2());
    update_front_state(v_new);

    // If current base is not at the boundary, add it to the 
    // interior edge list
    if ( base.is_interior() )
      mesh_.add_interior_edge(base.v1(), base.v2());

    // Remove base edge
    front_.remove( base );

    // Mark new triangle as active
    t_new.is_active( true );

    // Add element area to the total mesh area
    mesh_.add_area( t_new.area() );

  } // advance_front() 

  /*------------------------------------------------------------------
  | Update the advancing front 
  ------------------------------------------------------------------*/
  Triangle* update_front(Edge& base_edge,
                         const Vec2d& new_vertex_position,
                         const Vec2d& search_position,
                         double search_range)
  {
    // Create potential triangles with all found vertices
    TriVector new_triangles = 
      create_possible_triangles(base_edge, search_position, search_range);

    if (new_triangles.size() > 0)
    {
      Triangle& t_new = choose_best_triangle(new_triangles, base_edge);
      advance_front(base_edge, t_new.v3(), t_new);
      return &t_new;
    }

    // Check if a potential triangle can be created with the base edge
    // and a newly created vertex
    Vertex&   b1 = base_edge.v1();
    Vertex&   b2 = base_edge.v2();
    Vertex&   v_new = mesh_.add_vertex(new_vertex_position);
    Triangle& t_new = mesh_.add_triangle(b1, b2, v_new);

    // Algorithm fails if new vertex or new triangle is invalid
    if ( remove_from_mesh_if_invalid(v_new, t_new) )
      return nullptr;

    // Update the advancing front with new vertex
    advance_front(base_edge, v_new, t_new);

    return &t_new;

  } // update_front()

  /*------------------------------------------------------------------
  | Update the advancing front 
  ------------------------------------------------------------------*/
  Triangle* update_front_exhaustive(Edge& base_edge, Vertex& v)
  {
    if ( !v.on_front() )
      return nullptr;

    auto o = orientation(base_edge.v1().xy(), 
                         base_edge.v2().xy(), v.xy());

    if ( o != Orientation::CCW )
      return nullptr;
     
    // Create new potential triangle 
    Triangle& t_new 
      = mesh_.add_triangle(base_edge.v1(), base_edge.v2(), v);

    if ( remove_from_mesh_if_invalid(t_new) )
      return nullptr;

    advance_front(base_edge, v, t_new);

    return &t_new;

  } // update_front_exhaustive()

  /*------------------------------------------------------------------
  | Check if mesh entities are invalid. If yes, remove them from
  | the mesh and return true. Else, simply return false.
  ------------------------------------------------------------------*/
  bool remove_from_mesh_if_invalid(Vertex& v)
  {
    if ( vertex_is_valid(v) )
      return false;
    
    mesh_.remove_vertex(v);
    return true;
  } 

  bool remove_from_mesh_if_invalid(Triangle& t)
  {
    if ( triangle_is_valid(t) )
      return false;

    mesh_.remove_triangle(t);
    return true;
  } 

  bool remove_from_mesh_if_invalid(Vertex& v, Triangle& t)
  {
    if ( vertex_is_valid(v) && triangle_is_valid(t) )
      return false;

    mesh_.remove_triangle(t);
    mesh_.remove_vertex(v);
    return true;
  } 

  bool remove_from_mesh_if_invalid(Vertex& v, Triangle& t1, Triangle& t2)
  {
    if ( vertex_is_valid(v) && triangle_is_valid(t1) && triangle_is_valid(t2) )
      return false;

    mesh_.remove_triangle(t1);
    mesh_.remove_triangle(t2);
    mesh_.remove_vertex(v);
    return true;
  }


private:

  /*------------------------------------------------------------------
  | For a given search location and a respective search range,
  | all vertices that are located in this vicinity are checked,
  | if they might be possible candidates for the generation of new
  | triangles
  ------------------------------------------------------------------*/
  TriVector create_possible_triangles(Edge& base_edge,
                                      const Vec2d& search_position,
                                      double search_range)
  {
    Vertices& vertices = mesh_.vertices();

    // Create potential triangles with all vertices in vicinity of 
    // given search position and search range
    TriVector new_triangles {};

    for ( Vertex* v : vertices.get_items(search_position, search_range) )
    {
      // Skip vertices that are not located on the advancing front
      if ( !v->on_front() )
        continue;

      // Skip vertices that are colinear to the current base edge
      if ( orientation( base_edge.v1().xy(), base_edge.v2().xy(), v->xy() )
          == Orientation::CL )
        continue;

      // Create new potential triangle 
      Triangle& t_new = mesh_.add_triangle( base_edge.v1(), base_edge.v2(), *v );

      // Check if new potential triangle is valid
      if ( !remove_from_mesh_if_invalid(t_new) )
        new_triangles.push_back( &t_new );
    }

    return std::move(new_triangles);

  } // create_possible_triangles()

  /*------------------------------------------------------------------
  | We sort a given vector of <new_triangles> in descending order
  | according to the triangle quality.
  | Finally, the advancing front is updated with the triangle of best
  | quality and all other triangles are removed.
  ------------------------------------------------------------------*/
  Triangle& choose_best_triangle(TriVector& new_triangles,
                                 Edge&      base)
  {
    ASSERT( new_triangles.size() > 0,
        "FrontUpdate::choose_best_triangle(): Invalid input vector");
    DEBUG_LOG("VALID TRIANGLES IN NEIGHBORHOOD: " 
       << new_triangles.size()
    );

    std::sort( new_triangles.begin(), new_triangles.end(),
    [this] ( Triangle* t1, Triangle* t2 )
    {
      const double h1 = domain_.size_function( t1->xy() );
      const double h2 = domain_.size_function( t2->xy() );
      const double q1 = t1->quality(h1);
      const double q2 = t2->quality(h2);

      return ( q1 > q2 );
    });

    for (std::size_t i = 1; i < new_triangles.size(); i++)
      mesh_.triangles().remove( *new_triangles[i] );

    return *new_triangles[0];

  } // choose_best_triangle()

  /*------------------------------------------------------------------
  | Check if a triangle is valid. If yes, return true - 
  | else return false.
  ------------------------------------------------------------------*/
  bool triangle_is_valid(const Triangle& tri)
  {
    Vertices&   vertices = mesh_.vertices();
    Triangles& triangles = mesh_.triangles();
    Quads&         quads = mesh_.quads();

    const double rho   = domain_.size_function( tri.xy() );
    const double range = MAX( 2.0 * rho, tri.max_edge_length() );

    DEBUG_LOG("CHECK NEW TRIANGLE: " << tri);

    if ( !tri.is_valid() )
      return false;

    if ( tri.intersects_front( front_, range ) )
    { DEBUG_LOG("  > FRONT INTERSECTION"); return false; }

    if ( tri.intersects_domain( domain_ ) )
    { DEBUG_LOG("  > DOMAIN INTERSECTION"); return false; }

    if ( tri.intersects_vertex( vertices, range ) )
    { DEBUG_LOG("  > VERTEX INTERSECTION"); return false; }

    if ( tri.intersects_triangle( triangles, range ) )
    { DEBUG_LOG("  > TRIANGLE INTERSECTION"); return false; }

    if ( tri.intersects_quad( quads, range ) )
    { DEBUG_LOG("  > QUAD INTERSECTION"); return false; }

    if ( tri.quality(rho) < min_cell_quality_ )
    { DEBUG_LOG("  > BAD TRIANGLE QUALITY"); return false; }

    if ( tri.max_angle() > max_cell_angle_ )
    { DEBUG_LOG("  > BAD MAXIMUM ANGLE"); return false; }

    DEBUG_LOG("  > VALID");
    return true;

  } // Mesh::triangle_is_valid()

  /*------------------------------------------------------------------
  | Check if a vertex is valid. If yes, return true - 
  | else return false.
  ------------------------------------------------------------------*/
  bool vertex_is_valid(const Vertex& v)
  {
    Triangles& triangles = mesh_.triangles();
    Quads&         quads = mesh_.quads();

    const double rho   = domain_.size_function( v.xy() );
    const double range = 2.0 * rho;

    DEBUG_LOG("CHECK NEW VERTEX: " << v);

    if ( !domain_.is_inside( v ) )
    { DEBUG_LOG("  > OUTSIDE DOMAIN"); return false; }

    if ( v.intersects_facet(triangles, range) )
    { DEBUG_LOG("  > TRIANGLE INTERSECTION"); return false; }

    if ( v.intersects_facet(quads, range) )
    { DEBUG_LOG("  > QUAD INTERSECTION"); return false; }

    DEBUG_LOG("  > VALID");
    return true;

  } // Mesh::vertex_is_valid()


  /*------------------------------------------------------------------
  | Update the "on_front" state of a given vertex
  ------------------------------------------------------------------*/
  void update_front_state(Vertex& v)
  {
    bool on_front = false;

    for ( const auto& e_ptr : v.edges() )
      if ( &e_ptr->edgelist() == &front_ )
      {
        on_front = true;
        break;
      }

    if ( on_front )
      v.add_property( VertexProperty::on_front );
    else
      v.remove_property( VertexProperty::on_front );

    return;
  }
  
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Mesh&           mesh_;
  const Domain&   domain_;
  Front&          front_;

  double          min_cell_quality_ = 0.0;
  double          max_cell_angle_   = M_PI;

}; // FrontUpdate


} // namespace TQMesh
