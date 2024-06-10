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
class ModificationStrategy
{
public:
  /*------------------------------------------------------------------
  |
  ------------------------------------------------------------------*/
  ModificationStrategy(Mesh& mesh, const Domain& domain)
  : mesh_   { &mesh }
  , domain_ { &domain }
  {}

  virtual ~ModificationStrategy() {}

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  Mesh& mesh() { return *mesh_; }

  virtual bool modify() = 0;

protected:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Mesh*         mesh_;
  const Domain* domain_;

}; // ModificationStrategy


/*********************************************************************
* Turn triangular elements to quads 
*********************************************************************/
class Tri2QuadStrategy : public ModificationStrategy
{
public:

  using EdgeList   = std::list<Edge*>;

  /*------------------------------------------------------------------
  |
  ------------------------------------------------------------------*/
  Tri2QuadStrategy(Mesh& mesh, const Domain& domain)
  : ModificationStrategy(mesh, domain) {}

  ~Tri2QuadStrategy() {}

  /*------------------------------------------------------------------
  |
  ------------------------------------------------------------------*/
  bool modify() override
  {
    MeshCleanup::setup_facet_connectivity(*mesh_);

    tri_edges_.clear();

    collect_triangle_edges();

    merge_triangles_to_quads();

    mesh_->clear_waste();

    // Fix bad elements
    MeshCleanup::clear_double_quad_edges(*mesh_, true);
    MeshCleanup::clear_double_triangle_edges(*mesh_, false);

    return true;

  } // Tri2QuadStrategy::modify()

private:

  /*------------------------------------------------------------------
  | This function collects all internal edges that are adjacent to
  | two triangles and sorts them with increasing minimum triangle
  | edge length.
  ------------------------------------------------------------------*/
  void collect_triangle_edges()
  {
    // Pick all internal edges, which are adjacent to two triangles
    for ( const auto& e_ptr : mesh_->interior_edges() )
    {
      Facet* f_l = e_ptr->facet_l();
      Facet* f_r = e_ptr->facet_r();

      if (  ( NullFacet::is_not_null(f_l) && f_l->n_vertices() == 3 )
         && ( NullFacet::is_not_null(f_r) && f_r->n_vertices() == 3 ) )
        tri_edges_.push_back( e_ptr.get() );
    }

    // Sort edge list with increasing minimum edge length of their
    // adjacent triangles
    tri_edges_.sort(
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

  } // Tri2QuadStrategy::collect_triangle_edges()

  /*------------------------------------------------------------------
  |                                  v2               q_r
  | Loop over all sorted edges and     *-------------*
  | merge their adjacent triangles     | \           |
  | to quadrilaterals. It may be,      |   \   f_r   |
  | that the triangle of an            |     \       |
  | upcoming edge has already been     |       \     |
  | merged in this process. Thus       |  f_l    \   |
  | we have to make sure, that the     |           \ |
  | triangles to merge still exist.    *-------------*
  |                                   q_l             v1
  ------------------------------------------------------------------*/
  void merge_triangles_to_quads()
  {
    for ( auto& e : tri_edges_ )
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
      mesh_->remove_interior_edge( *e );

      // Create new quadrilateral element
      Quad& q_new = mesh_->add_quad( q_l, v1, q_r, v2, f_l->color() );
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

        Edge* e_share = mesh_->interior_edges().get_edge(q1,q2);

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
      mesh_->remove_triangle( *(static_cast<Triangle*>(f_l)) );
      mesh_->remove_triangle( *(static_cast<Triangle*>(f_r)) );
    }

  } // Tri2QuadStrategy::merge_triangles_to_quads() 

  /*------------------------------------------------------------------
  |
  ------------------------------------------------------------------*/
  EdgeList tri_edges_ {};

}; // Tri2QuadStrategy


} // namespace TQMesh
