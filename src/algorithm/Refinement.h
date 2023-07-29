/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <vector>

#include "VecND.h"

#include "Vertex.h"
#include "Edge.h"
#include "Triangle.h"
#include "Quad.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* 
*********************************************************************/
class Refinement
{
public:

  using EdgeVector = std::vector<Edge*>;
  using TriVector  = std::vector<Triangle*>;
  using QuadVector = std::vector<Quad*>;

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
  template <typename Mesh>
  static inline bool refine_to_quads(Mesh& mesh)
  {
    mesh.clear_waste();

    // Check if mesh has any twin edges - only meshes without 
    // twin edges can be refined in order to maintain conformity
    // between neighboring meshes
    for ( auto& e_ptr : mesh.boundary_edges() )
      if ( e_ptr->twin_edge() )
        return false;

    // Gather all coarse edges, quads and tris
    EdgeVector coarse_intr_edges {};
    EdgeVector coarse_bdry_edges {};
    TriVector  coarse_tris {};
    QuadVector coarse_quads {};

    for ( auto& e : mesh.interior_edges() )
      coarse_intr_edges.push_back( e.get() );

    for ( auto& e : mesh.boundary_edges() )
      coarse_bdry_edges.push_back( e.get() );

    for ( auto& q_ptr : mesh.quads() )
      coarse_quads.push_back( q_ptr.get() );

    for ( auto& t_ptr : mesh.triangles() )
      coarse_tris.push_back( t_ptr.get() );

    // Refine interior edges
    for ( auto e : coarse_intr_edges )
    {
      Vertex& v = mesh.add_vertex( e->xy() );
      mesh.interior_edges().add_edge( e->v1(), v );
      mesh.interior_edges().add_edge( v, e->v2() );
      e->sub_vertex( &v );

      if ( e->v1().is_fixed() && e->v2().is_fixed() )
        v.is_fixed( true );
    }

    // Refine boundary edges
    for ( auto e : coarse_bdry_edges )
    {
      Vertex& v = mesh.add_vertex( e->xy() );
      mesh.boundary_edges().add_edge( e->v1(), v, e->marker() );
      mesh.boundary_edges().add_edge( v, e->v2(), e->marker() );
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

      Edge* e13 = mesh.get_edge( v1, v3 );
      Edge* e35 = mesh.get_edge( v3, v5 );
      Edge* e51 = mesh.get_edge( v5, v1 );

      Vertex* v2 = e13->sub_vertex();
      Vertex* v4 = e35->sub_vertex();
      Vertex* v6 = e51->sub_vertex();

      // Create new vertex at center of quad
      Vertex& v7 = mesh.add_vertex( t->xy() );

      // Create new sub-quads 
      Quad& q1 = mesh.add_quad( v1, *v2, v7, *v6 );
      Quad& q2 = mesh.add_quad( v3, *v4, v7, *v2 );
      Quad& q3 = mesh.add_quad( v5, *v6, v7, *v4 );

      // New quads get assigned to colors of old element
      q1.color( t->color() );
      q2.color( t->color() );
      q3.color( t->color() );

      // Create new interior edges
      mesh.interior_edges().add_edge( *v2, v7 );
      mesh.interior_edges().add_edge( *v4, v7 );
      mesh.interior_edges().add_edge( *v6, v7 );
    }

    // Refine all quad elements
    for ( auto q : coarse_quads )
    {
      Vertex& v1 = q->v1();
      Vertex& v3 = q->v2();
      Vertex& v5 = q->v3();
      Vertex& v7 = q->v4();

      Edge* e13 = mesh.get_edge( v1, v3 );
      Edge* e35 = mesh.get_edge( v3, v5 );
      Edge* e57 = mesh.get_edge( v5, v7 );
      Edge* e71 = mesh.get_edge( v7, v1 );

      Vertex* v2 = e13->sub_vertex();
      Vertex* v4 = e35->sub_vertex();
      Vertex* v6 = e57->sub_vertex();
      Vertex* v8 = e71->sub_vertex();

      // Create new vertex at center of quad
      Vertex& v9 = mesh.add_vertex( q->xy() );

      // Create new sub-quads 
      Quad& q1 = mesh.add_quad( v1, *v2, v9, *v8 );
      Quad& q2 = mesh.add_quad( v3, *v4, v9, *v2 );
      Quad& q3 = mesh.add_quad( v5, *v6, v9, *v4 );
      Quad& q4 = mesh.add_quad( v7, *v8, v9, *v6 );

      // New quads get assigned to colors of old element
      q1.color( q->color() );
      q2.color( q->color() );
      q3.color( q->color() );
      q4.color( q->color() );

      // Create new interior edges
      mesh.interior_edges().add_edge( *v2, v9 );
      mesh.interior_edges().add_edge( *v4, v9 );
      mesh.interior_edges().add_edge( *v6, v9 );
      mesh.interior_edges().add_edge( *v8, v9 );
    }

    // Remove old entitires
    for ( auto e : coarse_intr_edges )
      mesh.remove_interior_edge( *e );
      
    for ( auto e : coarse_bdry_edges )
      mesh.remove_boundary_edge( *e );

    for ( auto t : coarse_tris )
      mesh.remove_triangle( *t );

    for ( auto q : coarse_quads )
      mesh.remove_quad( *q );

    return true;
    
  } // refine_to_quads()


private:
  /*------------------------------------------------------------------
  | We hide the constructor, since this class acts only as container
  | for static inline functions
  ------------------------------------------------------------------*/
  Refinement() = default;
  ~Refinement() {};
};


} // namespace TQAlgorithm
} // namespace TQMesh
