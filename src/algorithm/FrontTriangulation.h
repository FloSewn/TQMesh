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
#include "Facet.h"
#include "Boundary.h"
#include "Front.h"
#include "Domain.h"
#include "Mesh.h"
#include "FrontInitializer.h"
#include "FrontAlgorithm.h"
#include "MeshValidator.h"
#include "MeshingFunctions.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* 
*********************************************************************/
class FrontTriangulation : public FrontAlgorithm
{
public:
  using VertexVector   = std::vector<Vertex*>;
  using TriVector      = std::vector<Triangle*>;

  /*------------------------------------------------------------------
  | Constructor / Destructor 
  ------------------------------------------------------------------*/
  FrontTriangulation(Mesh& mesh, const Domain& domain)
  : FrontAlgorithm(mesh, domain) {}

  ~FrontTriangulation() {}

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  double mesh_range_factor() const { return mesh_range_factor_; }
  double base_height_factor() const { return base_height_factor_; }
  double wide_search_factor() const { return wide_search_factor_; }
  double min_cell_quality() const { return validator_.min_cell_quality(); }
  double max_cell_angle() const { return validator_.max_cell_angle(); }
  double base_vertex_factor() const { return base_vertex_factor_; }

  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  void mesh_range_factor(double v) { mesh_range_factor_ = v; }
  void base_height_factor(double v) { base_height_factor_ = v; }
  void wide_search_factor(double v) { wide_search_factor_ = v; }
  void min_cell_quality(double v) { validator_.min_cell_quality(v); }
  void max_cell_angle(double v) { validator_.max_cell_angle(v); }
  void base_vertex_factor(double v) { base_vertex_factor_ = v; }

  /*------------------------------------------------------------------
  | Triangulate a given initialized mesh structure
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
    finish_mesh_for_output();

    // Remove remaining edges from the front
    front_.clear_edges();

    return success;
  }

private:

  /*------------------------------------------------------------------
  | Let the advancing front create a new triangle 
  ------------------------------------------------------------------*/
  bool advance_front_triangle(Edge& base_edge, bool wide_search=false)
  {
    // Obtain the position of a potential new vertex and the radius 
    // to search around it for potential existing vertices
    const double height = domain_.size_function( base_edge.xy() );
    const Vec2d  v_xy   = base_edge.xy() 
                        + base_height_factor_ * height * base_edge.normal();

    double range = mesh_range_factor_ * domain_.size_function( v_xy );
    if (wide_search)
      range *= wide_search_factor_;

    // Find all vertices in the vicinity of the current base edge
    VertexVector vertex_candidates = find_vertex_candidates(v_xy, range);

    // Create potential triangles with all found vertices
    TriVector new_triangles {};
    MeshingFunctions::check_vertex_candidates(vertex_candidates, 
                                              base_edge, 
                                              new_triangles,
                                              mesh_, validator_);

    // If potential triangles have been found, choose the best one
    if ( choose_best_triangle(new_triangles, base_edge) )
      return true;

    // Check if a potential triangle can be created with the base edge
    // and a newly created vertex
    Vertex& b1 = base_edge.v1();
    Vertex& b2 = base_edge.v2();

    Vertex& v_new 
      = MeshingFunctions::create_base_vertex(base_edge, domain_, mesh_,
                                             base_vertex_factor_ );
    Triangle& t_new = mesh_.add_triangle(b1, b2, v_new);
    
    // Algorithm fails if new vertex or new triangle is invalid
    if ( validator_.remove_from_mesh_if_invalid(v_new, t_new) )
      return false;
    
    // Update the advancing front with new vertex
    update_front(base_edge, v_new, t_new);

    return true;

  } // Mesh::advance_front_triangle() */

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
      if ( advance_front_triangle(*base_edge, wide_search) )
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
  double mesh_range_factor_  = 1.0;
  double base_height_factor_ = sqrt(3.0) / 4.0; 
  double wide_search_factor_ = 10.0;
  double base_vertex_factor_ = 2.00;

}; // FrontTriangulation

} // namespace TQAlgorithm
} // namespace TQMesh
