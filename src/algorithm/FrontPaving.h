/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <vector>

#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Mesh.h"
#include "FrontAlgorithm.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* 
*********************************************************************/
class FrontPaving : public FrontAlgorithm
{
public:
  using VertexVector   = std::vector<Vertex*>;
  using TriVector      = std::vector<Triangle*>;

  /*------------------------------------------------------------------
  | Constructor / Destructor 
  ------------------------------------------------------------------*/
  FrontPaving(Mesh& mesh, const Domain& domain)
  : FrontAlgorithm(mesh, domain) {}

  ~FrontPaving() {}

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  double wide_search_factor() const { return wide_search_factor_; }
  double min_cell_quality() const { return validator_.min_cell_quality(); }
  double max_cell_angle() const { return validator_.max_cell_angle(); }

  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  void wide_search_factor(double v) { wide_search_factor_ = v; }
  void min_cell_quality(double v) { validator_.min_cell_quality(v); }
  void max_cell_angle(double v) { validator_.max_cell_angle(v); }

  /*------------------------------------------------------------------
  | Pave a given initialized mesh structure with quads
  ------------------------------------------------------------------*/
  bool generate_elements(int n_elements=0) override
  {
    ASSERT(mesh_.n_boundary_edges() != 0,
      "FrontTriangulation::generate_elements(): "
      "Unable to triangulate mesh that has not been prepared yet.");

    // Prepare the mesh  
    Cleanup::setup_facet_connectivity(mesh_);

    // Initialize the advancing front and its base edge
    Edge* base_edge = init_advancing_front();

    // Remove invalid mesh edges that are no longer needed
    remove_invalid_mesh_edges();

    // Perform the actual mesh generation
    bool success = advancing_front_loop(base_edge, n_elements);

    // Finish mesh structure for output
    add_remaining_front_edges_to_mesh();

    // Remove remaining edges from the front
    front_.clear_edges();

    return success;
  }

private:

  /*------------------------------------------------------------------
  | Let the advancing front create a new triangle 
  ------------------------------------------------------------------*/
  bool advance_front_quad(Edge& base_edge, bool wide_search=false)
  {
    Vertex& b1 = base_edge.v1();
    Vertex& b2 = base_edge.v2();

    double h_b1 = domain_.size_function( b1.xy() );
    double h_b2 = domain_.size_function( b2.xy() );

    Vec2d p1 = b1.xy() + h_b1 * base_edge.normal();
    Vec2d p2 = b2.xy() + h_b1 * base_edge.normal();

    double r1 = domain_.size_function( p1 );
    double r2 = domain_.size_function( p2 );

    const double theta_1 = MAX(-0.25,MIN(0.25,1.0-(h_b1-r1)/h_b1));
    const double theta_2 = MAX(-0.25,MIN(0.25,1.0-(h_b2-r2)/h_b2));

    //p1 += theta_1 * base_edge.length() * base_edge.tangent();
    //p2 -= theta_2 * base_edge.length() * base_edge.tangent();

    if (wide_search)
    {
      r1 *= wide_search_factor_;
      r2 *= wide_search_factor_;
    }

    // ****** Create first triangle *******
    TriVector new_tris_p1 {};

    // Find all vertices in vicinity of p1 
    VertexVector vertex_candidates_p1 = find_vertex_candidates(p1, r1);

    // Create potential triangles with all found vertices
    check_vertex_candidates(vertex_candidates_p1, base_edge, new_tris_p1);

    // If potential triangles have been found, choose the best one
    // --> this function also updates the advancing front
    Triangle* t1 = choose_best_triangle(new_tris_p1, base_edge);

    // No triangle has been found yet -> create new one
    if ( !t1 )
    {
      Vertex&   v_new = mesh_.add_vertex( p1 );
      Triangle& t_new = mesh_.add_triangle( b1, b2, v_new );

      // Algorithm fails if new vertex or new triangle is invalid
      if ( validator_.remove_from_mesh_if_invalid(v_new, t_new) )
        return false;

      // Update the advancing front with new vertex
      update_front( base_edge, v_new, t_new );

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
    VertexVector vertex_candidates_p2 = find_vertex_candidates(p2, r2);

    // Create potential triangles with all found vertices
    check_vertex_candidates(vertex_candidates_p2, *diag, new_tris_p2);

    // If potential triangles have been found, choose the best one
    // --> this function also updates the advancing front
    Triangle* t2 = choose_best_triangle(new_tris_p2, *diag);

    // No triangle has been found yet -> create new one
    if ( !t2 )
    {
      Vertex&   v_new = mesh_.add_vertex( p2 );
      Triangle& t_new = mesh_.add_triangle(d1, d2, v_new);

      // Algorithm fails if new vertex or new triangle is invalid
      if ( validator_.remove_from_mesh_if_invalid(v_new, t_new) )
        return true;

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
    Edge* e_rem = mesh_.get_interior_edge(q2, q4);
    mesh_.remove_interior_edge( *e_rem );

    // Remove old triangular elements
    mesh_.remove_triangle( *t1 );
    mesh_.remove_triangle( *t2 );

    // Create new quadrilateral element
    Quad& q_new = mesh_.add_quad( q1, q2, q3, q4 );
    q_new.is_active( true );

    return true;


  } // FrontPaving::advance_front_quad() */

  /*------------------------------------------------------------------
  | The actual main loop for the advancing front mesh generation
  ------------------------------------------------------------------*/
  bool advancing_front_loop(Edge* base_edge, int n_elements)
  {
    unsigned int iteration   = 0;
    bool         wide_search = false;
    int          n_generated = 0;

    while ( true )
    {
      // Advance current base edge 
      if ( advance_front_quad(*base_edge, wide_search) )
      {
        ++n_generated;

        // Sort front edges after a wide search
        if ( wide_search )
          front_.sort_edges( false );

        // Reset iteration counter and wide search
        iteration = 0;
        wide_search = false;

        // Go to the next base edge
        base_edge = front_.set_base_first();
        mesh_.clear_waste();
      }
      // If it failed, go to the next base edge
      else
      {
        base_edge = front_.set_base_next();
        ++iteration;
      }

      // All front edges failed to create new elements
      // --> Activate wide search for neighboring vertices
      //     and re-run the algorithm
      if ( iteration == front_.size() && !wide_search )
      {
        wide_search = true;
        iteration = 0;
      }

      update_progress_bar();

      // No more edges in the advancing front
      // --> Meshing algorithm succeeded
      if ( front_.size() == 0 )
        return true;
      
      // Maximum number of elements has been generated
      if (n_elements > 0 && n_generated == n_elements)
        return true;

      // All front edges faild to create new elements, even
      // when using the wide search 
      // --> Meshing algorithm failed
      if ( iteration == front_.size() && wide_search )
        return false;
    }

  } // FrontTriangulation::advancing_front_loop()



  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  double wide_search_factor_ = 10.0;

}; // FrontPaving

} // namespace TQAlgorithm
} // namespace TQMesh
