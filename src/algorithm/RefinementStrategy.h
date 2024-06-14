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

#include "Vertex.h"
#include "Edge.h"
#include "Triangle.h"
#include "Quad.h"
#include "Mesh.h"
#include "Domain.h"

namespace TQMesh {

using namespace CppUtils;


/*********************************************************************
* 
*********************************************************************/
class RefinementStrategy
{
public:
  /*------------------------------------------------------------------
  |
  ------------------------------------------------------------------*/
  RefinementStrategy(Mesh& mesh, const Domain& domain)
  : mesh_   { &mesh }
  , domain_ { &domain }
  {}

  virtual ~RefinementStrategy() {}

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  Mesh& mesh() { return *mesh_; }

  virtual bool refine() = 0;

protected:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Mesh*         mesh_;
  const Domain* domain_;

}; // RefinementStrategy


/*********************************************************************
* Every triangle is refined into three quads and every quads is 
* refined into four quads. This results in an all-quad mesh.
*
*  qi... sub-quads 
*  vi... vertices
*        
*                         
*                      (v5)                        (v7)         (v6)      (v5) 
*                       o                             o<--------o---------o    
*                     /   \                           |         |         ^
*                   / [q3]  \                         |  [q4]   |   [q3]  |
*                 /           \                       |         |         |
*         (v6)  /      (v7)     \  (v4)               |         |(v9)     |  
*             o---------o---------o               (v8)o---------o---------o(v4)
*           /           |           \                 |         |         |
*         /             |             \               |         |         |
*       /     [q1]      |      [q2]     \             |  [q1]   |   [q2]  |
*     /                 |                 \           v         |         |
*   o-------------------o-------------------o         o---------o-------->o
*  (v1)                (v2)                 (v3)   (v1)        (v2)       (v3)
*
*********************************************************************/
class QuadRefinement : public RefinementStrategy
{
public:

  using EdgeVector = std::vector<Edge*>;
  using TriVector  = std::vector<Triangle*>;
  using QuadVector = std::vector<Quad*>;

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  QuadRefinement(Mesh& mesh, const Domain& domain)
  : RefinementStrategy(mesh, domain) 
  {}

  ~QuadRefinement() {}

  /*------------------------------------------------------------------
  | The actual mesh refinement
  ------------------------------------------------------------------*/
  bool refine() override
  {
    // Empty waste 
    mesh_->clear_waste();
    coarse_intr_edges_.clear();
    coarse_bdry_edges_.clear();
    coarse_quads_.clear();     
    coarse_tris_.clear();      

    // Check if mesh has any twin edges - only meshes without 
    // twin edges can be refined in order to maintain conformity
    // between neighboring meshes
    for ( auto& e_ptr : mesh_->boundary_edges() )
      if ( e_ptr->twin_edge() )
        return false;

    gather_entities_to_refine();

    refine_interior_edges();

    refine_boundary_edges();

    refine_triangles();

    refine_quads();

    remove_old_entities();

    return true;
    
  } // refine_to_quads()

private:

  /*------------------------------------------------------------------
  | Gather all coarse edges, quads and tris
  ------------------------------------------------------------------*/
  void gather_entities_to_refine()
  {
    for ( auto& e : mesh_->interior_edges() )
      coarse_intr_edges_.push_back( e.get() );

    for ( auto& e : mesh_->boundary_edges() )
      coarse_bdry_edges_.push_back( e.get() );

    for ( auto& q_ptr : mesh_->quads() )
      coarse_quads_.push_back( q_ptr.get() );

    for ( auto& t_ptr : mesh_->triangles() )
      coarse_tris_.push_back( t_ptr.get() );
  }

  /*------------------------------------------------------------------
  | Refine interior edges
  ------------------------------------------------------------------*/
  void refine_interior_edges()
  {
    for ( auto e : coarse_intr_edges_ )
    {
      Vertex& v = mesh_->add_vertex( e->xy() );
      mesh_->interior_edges().add_edge( e->v1(), v );
      mesh_->interior_edges().add_edge( v, e->v2() );
      e->sub_vertex( &v );

      v.add_property( e->v1().properties() );
      v.add_property( e->v2().properties() );

      // Interior vertices can not be boundary vertices
      if ( v.has_property( VertexProperty::on_boundary ) )
        v.remove_property( VertexProperty::on_boundary );
    }
  }

  /*------------------------------------------------------------------
  | Refine boundary edges
  ------------------------------------------------------------------*/
  void refine_boundary_edges()
  {
    for ( auto e : coarse_bdry_edges_ )
    {
      Vertex& v = mesh_->add_vertex( e->xy() );
      mesh_->boundary_edges().add_edge( e->v1(), v, e->marker() );
      mesh_->boundary_edges().add_edge( v, e->v2(), e->marker() );
      e->sub_vertex( &v );

      v.add_property( e->v1().properties() );
      v.add_property( e->v2().properties() );

      ASSERT( v.has_property( VertexProperty::on_boundary ),
        "RefinementStrategy::refine_boundary_edges(): Missing "
        "vertex property \"on_boundary\".");
    }
  }

  /*------------------------------------------------------------------
  | Refine triangles
  ------------------------------------------------------------------*/
  void refine_triangles()
  {
    for ( auto t : coarse_tris_ )
    {
      Vertex& v1 = t->v1();
      Vertex& v3 = t->v2();
      Vertex& v5 = t->v3();

      Edge* e13 = mesh_->get_edge( v1, v3 );
      Edge* e35 = mesh_->get_edge( v3, v5 );
      Edge* e51 = mesh_->get_edge( v5, v1 );

      Vertex* v2 = e13->sub_vertex();
      Vertex* v4 = e35->sub_vertex();
      Vertex* v6 = e51->sub_vertex();

      // Create new vertex at center of quad
      Vertex& v7 = mesh_->add_vertex( t->xy() );

      v7.add_property( v1.properties() );
      v7.add_property( v3.properties() );
      v7.add_property( v5.properties() );

      // Interior vertices can not be boundary vertices
      if ( v7.has_property( VertexProperty::on_boundary ) )
        v7.remove_property( VertexProperty::on_boundary );

      // Create new sub-quads 
      Quad& q1 = mesh_->add_quad( v1, *v2, v7, *v6 );
      Quad& q2 = mesh_->add_quad( v3, *v4, v7, *v2 );
      Quad& q3 = mesh_->add_quad( v5, *v6, v7, *v4 );

      // New quads get assigned to colors of old element
      q1.color( t->color() );
      q2.color( t->color() );
      q3.color( t->color() );

      // Create new interior edges
      mesh_->interior_edges().add_edge( *v2, v7 );
      mesh_->interior_edges().add_edge( *v4, v7 );
      mesh_->interior_edges().add_edge( *v6, v7 );
    }
  }

  /*------------------------------------------------------------------
  | Refine quads
  ------------------------------------------------------------------*/
  void refine_quads()
  {
    for ( auto q : coarse_quads_ )
    {
      Vertex& v1 = q->v1();
      Vertex& v3 = q->v2();
      Vertex& v5 = q->v3();
      Vertex& v7 = q->v4();

      Edge* e13 = mesh_->get_edge( v1, v3 );
      Edge* e35 = mesh_->get_edge( v3, v5 );
      Edge* e57 = mesh_->get_edge( v5, v7 );
      Edge* e71 = mesh_->get_edge( v7, v1 );

      Vertex* v2 = e13->sub_vertex();
      Vertex* v4 = e35->sub_vertex();
      Vertex* v6 = e57->sub_vertex();
      Vertex* v8 = e71->sub_vertex();

      // Create new vertex at center of quad
      Vertex& v9 = mesh_->add_vertex( q->xy() );

      v9.add_property( v1.properties() );
      v9.add_property( v3.properties() );
      v9.add_property( v5.properties() );
      v9.add_property( v7.properties() );

      // Interior vertices can not be boundary vertices
      if ( v9.has_property( VertexProperty::on_boundary ) )
        v9.remove_property( VertexProperty::on_boundary );

      // Create new sub-quads 
      Quad& q1 = mesh_->add_quad( v1, *v2, v9, *v8 );
      Quad& q2 = mesh_->add_quad( v3, *v4, v9, *v2 );
      Quad& q3 = mesh_->add_quad( v5, *v6, v9, *v4 );
      Quad& q4 = mesh_->add_quad( v7, *v8, v9, *v6 );

      // New quads get assigned to colors of old element
      q1.color( q->color() );
      q2.color( q->color() );
      q3.color( q->color() );
      q4.color( q->color() );

      // Create new interior edges
      mesh_->interior_edges().add_edge( *v2, v9 );
      mesh_->interior_edges().add_edge( *v4, v9 );
      mesh_->interior_edges().add_edge( *v6, v9 );
      mesh_->interior_edges().add_edge( *v8, v9 );
    }
  }

  /*------------------------------------------------------------------
  | Remove old entities 
  ------------------------------------------------------------------*/
  void remove_old_entities()
  {
    for ( auto e : coarse_intr_edges_ )
      mesh_->remove_interior_edge( *e );
      
    for ( auto e : coarse_bdry_edges_ )
      mesh_->remove_boundary_edge( *e );

    for ( auto t : coarse_tris_ )
      mesh_->remove_triangle( *t );

    for ( auto q : coarse_quads_ )
      mesh_->remove_quad( *q );
  }

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  EdgeVector coarse_intr_edges_ {};
  EdgeVector coarse_bdry_edges_ {};
  QuadVector coarse_quads_      {};
  TriVector  coarse_tris_       {};

};


} // namespace TQMesh
