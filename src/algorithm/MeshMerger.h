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
* This class handles the merger between two meshes,
* where the vertices, elements and edges of a donor mesh are 
* copied to a receiver mesh - if both meshes share common interface
* edges.
*********************************************************************/
class MeshMerger
{
public:
  using VertexVector   = std::vector<Vertex*>;
  using EdgeVector     = std::vector<Edge*>;
  using VertexPair     = std::vector<std::pair<Vertex*, Vertex*>>;

  /*------------------------------------------------------------------
  | Constructor / Destructor
  ------------------------------------------------------------------*/
  MeshMerger(Mesh& receiver, Mesh& donor)
  : receiver_     { &receiver } 
  , donor_        { &donor }
  , new_vertices_ { donor.vertices().size(), nullptr }
  {}

  virtual ~MeshMerger() {}

  /*------------------------------------------------------------------
  | Copy donor entities to receiver (if they share boundary edges)
  ------------------------------------------------------------------*/
  bool merge()
  {
    if ( merged_ )
      return false;

    MeshCleanup::setup_facet_connectivity(*receiver_);
    MeshCleanup::setup_facet_connectivity(*donor_);

    // Assign indices to donor mesh
    MeshCleanup::assign_mesh_indices(*receiver_);
    MeshCleanup::assign_mesh_indices(*donor_);

    // Search for vertices that are located on the inteface
    // boundary between donor and receiver
    // -> These vertices will not be copied
    if ( !collect_interface_boundary_vertices() )
      return false;

    // Copy all relevant donor entities 
    copy_donor_vertices();
    copy_donor_quads();
    copy_donor_triangles();
    copy_donor_interior_edges();
    copy_donor_boundary_edges();
    
    // Find boundary inteface edges that should be turned into 
    // interior edges
    convert_boundary_interface_edges();

    // Finally, we need to make sure that all boundary vertices 
    // are still marked as "on_boundary" (this information might have
    // been lost while collecting interface boundary vertices )
    for ( const auto& e_ptr : receiver_->boundary_edges() )
    {
      e_ptr->v1().add_property(VertexProperty::on_boundary);
      e_ptr->v2().add_property(VertexProperty::on_boundary);
    }


    // Forbid to run again
    merged_ = true;

    return true;

  } // merge()

private:

  /*------------------------------------------------------------------
  | Loop over all the receiver's boundary edges and loacte all twin 
  | edges that are connected to the donor mesh. 
  | Return false if both meshes share no boundary edges.
  ------------------------------------------------------------------*/
  bool collect_interface_boundary_vertices()
  {
    std::size_t n_twins = 0;

    for ( const auto& e_ptr : receiver_->boundary_edges() )
    {
      Edge* e_twin = e_ptr->twin_edge();

      if ( !e_twin )
        continue;

      Facet* f_twin = e_twin->facet_l();
      ASSERT( f_twin, "MeshBuilder::merge_meshes(): "
        "Invalid boundary edge connectivity." );

      if ( f_twin->mesh() != donor_ )
        continue;

      auto v1_index = e_twin->v1().index();
      auto v2_index = e_twin->v2().index();

      // Keep in mind, that both edges point in different directions
      new_vertices_[v2_index] = &(e_ptr->v1());
      new_vertices_[v1_index] = &(e_ptr->v2());

      // We remove the "on_boundary" property, but we want to preserve the 
      // location of the interface vertices and thus set them as "fixed"
      new_vertices_[v1_index]->remove_property( VertexProperty::on_boundary );
      new_vertices_[v2_index]->remove_property( VertexProperty::on_boundary );
      new_vertices_[v1_index]->add_property( VertexProperty::is_fixed );
      new_vertices_[v2_index]->add_property( VertexProperty::is_fixed );

      ++n_twins;
    }

    return (n_twins > 0);

  } // collect_interface_boundary_vertices()

  /*------------------------------------------------------------------
  | Copy all donor vertices to the receiver, which have not been
  | flagged as interface boundary vertices
  ------------------------------------------------------------------*/
  void copy_donor_vertices()
  {
    for ( const auto& v_ptr : donor_->vertices() )
    {
      auto v_index = v_ptr->index();

      if ( !new_vertices_[v_index] )
      {
        Vertex& v_new = receiver_->add_vertex( v_ptr->xy() );
        v_new.add_property( v_ptr->properties() );

        new_vertices_[v_index] = &v_new;
      }
    }

  } // copy_donor_vertices()

  /*------------------------------------------------------------------
  | Copy all the donor's quad elements
  ------------------------------------------------------------------*/
  void copy_donor_quads()
  {
    for ( const auto& q_ptr : donor_->quads() )
    {
      auto v1_index = q_ptr->v1().index();
      auto v2_index = q_ptr->v2().index();
      auto v3_index = q_ptr->v3().index();
      auto v4_index = q_ptr->v4().index();

      Vertex* v1 = new_vertices_[v1_index];
      Vertex* v2 = new_vertices_[v2_index];
      Vertex* v3 = new_vertices_[v3_index];
      Vertex* v4 = new_vertices_[v4_index];

      ASSERT(v1, "MeshMerger::copy_donor_quads(): Invalid data structure.");
      ASSERT(v2, "MeshMerger::copy_donor_quads(): Invalid data structure.");
      ASSERT(v3, "MeshMerger::copy_donor_quads(): Invalid data structure.");
      ASSERT(v4, "MeshMerger::copy_donor_quads(): Invalid data structure.");

      receiver_->add_quad( *v1, *v2, *v3, *v4, q_ptr->color() );
    }
  } // copy_donor_quads()

  /*------------------------------------------------------------------
  | Copy all the donor's triangle elements
  ------------------------------------------------------------------*/
  void copy_donor_triangles()
  {
    for ( const auto& t_ptr : donor_->triangles() )
    {
      auto v1_index = t_ptr->v1().index();
      auto v2_index = t_ptr->v2().index();
      auto v3_index = t_ptr->v3().index();

      Vertex* v1 = new_vertices_[v1_index];
      Vertex* v2 = new_vertices_[v2_index];
      Vertex* v3 = new_vertices_[v3_index];

      ASSERT(v1, "MeshMerger::copy_donor_triangles(): Invalid data structure.");
      ASSERT(v2, "MeshMerger::copy_donor_triangles(): Invalid data structure.");
      ASSERT(v3, "MeshMerger::copy_donor_triangles(): Invalid data structure.");

      receiver_->add_triangle( *v1, *v2, *v3, t_ptr->color() );
    }
  } // copy_donor_triangles()

  /*------------------------------------------------------------------
  | Copy all the donor's interior edges
  ------------------------------------------------------------------*/
  void copy_donor_interior_edges()
  {
    for ( const auto& e_ptr : donor_->interior_edges() )
    {
      auto v1_index = e_ptr->v1().index();
      auto v2_index = e_ptr->v2().index();

      Vertex* v1 = new_vertices_[v1_index];
      Vertex* v2 = new_vertices_[v2_index];

      ASSERT(v1, "MeshMerger::copy_donor_interior_edges(): Invalid data structure.");
      ASSERT(v2, "MeshMerger::copy_donor_interior_edges(): Invalid data structure.");

      receiver_->interior_edges().add_edge( *v1, *v2 );
    }
  } // copy_donor_interior_edges()

  /*------------------------------------------------------------------
  | Copy all the donor's boundary edges
  ------------------------------------------------------------------*/
  void copy_donor_boundary_edges()
  {
    EdgeList& boundary_edges = receiver_->boundary_edges();

    for ( const auto& e_ptr : donor_->boundary_edges() )
    {
      Edge* e_twin = e_ptr->twin_edge();

      // Skip if edge is an interface to this mesh
      if ( e_twin && 
           e_twin->facet_l() && 
           e_twin->facet_l()->mesh() == receiver_ )
        continue;

      auto v1_index = e_ptr->v1().index();
      auto v2_index = e_ptr->v2().index();

      Vertex* v1 = new_vertices_[v1_index];
      Vertex* v2 = new_vertices_[v2_index];

      ASSERT(v1, "MeshBuilder::merge_meshes(): Invalid data structure.");
      ASSERT(v2, "MeshBuilder::merge_meshes(): Invalid data structure.");

      Edge& e_new = boundary_edges.add_edge(*v1, *v2, e_ptr->marker());

      // Add interface to possible other receiver meshes
      if ( e_twin )
      {
        e_twin->twin_edge( &e_new );
        e_new.twin_edge( e_twin );
      }
    }
  } // copy_donor_boundary_edges()

  /*------------------------------------------------------------------
  | Turn boundary interface edges into interior edges
  ------------------------------------------------------------------*/
  void convert_boundary_interface_edges()
  {
    EdgeVector    bdry_edges_to_remove;
    VertexPair    new_intr_edges;

    for ( const auto& e_ptr : receiver_->boundary_edges() )
    {
      Edge* e_twin = e_ptr->twin_edge();

      if ( !e_twin )
        continue;

      Facet* f_twin = e_twin->facet_l();
      ASSERT( f_twin, "MeshMerger::convert_boundary_interface_edges(): "
          "Invalid boundary edge connectivity." );

      if ( f_twin->mesh() != donor_ )
        continue;

      auto v1_index = e_twin->v1().index();
      auto v2_index = e_twin->v2().index();

      Vertex* v1 = new_vertices_[v1_index];
      Vertex* v2 = new_vertices_[v2_index];

      ASSERT(v1, "MeshMerger::convert_boundary_interface_edges(): Invalid data structure.");
      ASSERT(v2, "MeshMerger::convert_boundary_interface_edges(): Invalid data structure.");

      new_intr_edges.push_back( {v2, v1} );
      bdry_edges_to_remove.push_back( e_ptr.get() );
    }

    // Remove old edges
    for ( Edge* e : bdry_edges_to_remove )
      receiver_->remove_boundary_edge( *e );

    // Add new interior edges
    for ( auto v : new_intr_edges )
      receiver_->interior_edges().add_edge( *v.first, *v.second );

    // Clean up
    receiver_->clear_waste();


  } // convert_boundary_interface_edges()

  /*------------------------------------------------------------------
  | Attribute
  ------------------------------------------------------------------*/
  Mesh*         donor_;
  Mesh*         receiver_;
  VertexVector  new_vertices_;
  bool          merged_ { false };

}; // MeshMerger

} // namespace TQMesh
