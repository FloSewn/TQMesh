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
#include "MeshCleanup.h"
#include "MeshingStrategy.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* 
*********************************************************************/
class TriangulationStrategy : public MeshingStrategy
{
public:
  using VertexVector   = std::vector<Vertex*>;
  using TriVector      = std::vector<Triangle*>;

  /*------------------------------------------------------------------
  | Constructor / Destructor 
  ------------------------------------------------------------------*/
  TriangulationStrategy(Mesh& mesh, const Domain& domain)
  : MeshingStrategy(mesh, domain) {}

  ~TriangulationStrategy() {}

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  size_t n_elements() const { return n_elements_; }
  double mesh_range_factor() const { return mesh_range_factor_; }
  double wide_search_factor() const { return wide_search_factor_; }
  double min_cell_quality() const { return front_update_.min_cell_quality(); }
  double max_cell_angle() const { return front_update_.max_cell_angle(); }
  double base_vertex_factor() const { return base_vertex_factor_; }

  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  TriangulationStrategy& n_elements(size_t n) 
  { n_elements_ = n; return *this; }
  TriangulationStrategy& mesh_range_factor(double v) 
  { mesh_range_factor_ = v; return *this; }
  TriangulationStrategy& wide_search_factor(double v) 
  { wide_search_factor_ = v; return *this; }
  TriangulationStrategy& min_cell_quality(double v) 
  { front_update_.min_cell_quality(v); return *this; }
  TriangulationStrategy& max_cell_angle(double v) 
  { front_update_.max_cell_angle(v); return *this; }
  TriangulationStrategy& base_vertex_factor(double v) 
  { base_vertex_factor_ = v; return *this; }

  /*------------------------------------------------------------------
  | Triangulate a given initialized mesh structure
  ------------------------------------------------------------------*/
  bool generate_elements() override
  {
    if (mesh_.n_boundary_edges() < 1)
      return false;

    // Reset counter for generated elements
    n_generated_ = 0;

    // Prepare the mesh  
    MeshCleanup::setup_facet_connectivity(mesh_);

    // Initialize the advancing front and its base edge
    Edge* base_edge = init_advancing_front();

    // Remove invalid mesh edges that are no longer needed
    remove_invalid_mesh_edges();

    // Perform the actual mesh generation
    bool success = advancing_front_loop(base_edge, n_elements_);

    // In case of a failed meshing attempt, use the exhaustive 
    // search approach to fill gaps
    if ( !success )
    {
      int n_remaining = MAX(0, static_cast<int>(n_elements_-n_generated_));
      success = exhaustive_search_loop(base_edge, n_remaining);
    }

    // Finish mesh structure for output
    add_remaining_front_edges_to_mesh();

    // Remove remaining edges from the front
    front_.clear_edges();

    // Improve mesh quality
    if ( success ) 
    {
      MeshCleanup::clear_double_quad_edges(mesh_);
      MeshCleanup::clear_double_triangle_edges(mesh_);
      MeshCleanup::merge_degenerate_triangles(mesh_);
    }

    return success;

  } // TriangulationStrategy::generate_elements()


  /*------------------------------------------------------------------
  | Triangulate a given initialized mesh structure using the 
  | exhaustive search approach 
  ------------------------------------------------------------------*/
  bool generate_elements_exhaustive() 
  {
    if (mesh_.n_boundary_edges() < 1)
      return false;

    // Reset counter for generated elements
    n_generated_ = 0;

    // Prepare the mesh  
    MeshCleanup::setup_facet_connectivity(mesh_);

    // Initialize the advancing front and its base edge
    Edge* base_edge = init_advancing_front();

    // Remove invalid mesh edges that are no longer needed
    remove_invalid_mesh_edges();

    // Perform the actual mesh generation
    bool success = exhaustive_search_loop(base_edge, n_elements_);

    // Finish mesh structure for output
    add_remaining_front_edges_to_mesh();

    // Remove remaining edges from the front
    front_.clear_edges();

    return success;

  } // TriangulationStrategy::generate_elements_exhaustive()

private:

  /*------------------------------------------------------------------
  | Let the advancing front create a new triangle 
  ------------------------------------------------------------------*/
  bool advance_front_triangle(Edge& base_edge, bool wide_search=false)
  {
    // Factor h for height of equlateral triangle
    // h := sqrt(3) / 2 
    constexpr double h  = 0.8660254038; 

    const double l1  = base_edge.length() * h * base_vertex_factor_;
    const double l2  = domain_.size_function( base_edge.xy() );
    const double len = MIN(l1, l2);

    // Coordinate of new vertex 
    const Vec2d v_xy = base_edge.xy() + base_edge.normal() * len;
    double range = mesh_range_factor_ * len;

    if (wide_search)
      range *= wide_search_factor_;

    return front_update_.update_front(base_edge, v_xy, v_xy, range);

  } // TriangulationStrategy::advance_front_triangle() */

  /*------------------------------------------------------------------
  | The actual main loop for the advancing front mesh generation
  ------------------------------------------------------------------*/
  bool advancing_front_loop(Edge* base_edge, int n_elements)
  {
    unsigned int iteration   = 0;
    bool         wide_search = false;

    while ( true )
    {
      // Advance current base edge 
      if ( advance_front_triangle(*base_edge, wide_search) )
      {
        ++n_generated_;

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
      if (n_elements > 0 && n_generated_ == n_elements)
        return true;

      // All front edges faild to create new elements, even
      // when using the wide search 
      // --> Meshing algorithm failed
      if ( iteration == front_.size() && wide_search )
        return false;
    }

  } // TriangulationStrategy::advancing_front_loop()


  /*------------------------------------------------------------------
  | Triangulate the advancing front using an exhaustive search 
  | approach
  ------------------------------------------------------------------*/
  bool exhaustive_search_loop(Edge* base_edge, int n_elements)
  {
    unsigned int iteration = 0;
    front_.sort_edges();
    n_generated_ = 0;

    DEBUG_LOG("USE EXHAUSTIVE SEARCH FOR " << n_elements << " ELEMENTS");

    while ( true )
    {
      if ( exhaustive_search_triangle(*base_edge) )
      {
        ++n_generated_;
        iteration = 0;
        base_edge = front_.set_base_first();
        mesh_.clear_waste();
      }
      else
      {
        base_edge = front_.set_base_next();
        ++iteration;
      }

      update_progress_bar();

      if ( base_edge && iteration == front_.size() )
        return false;

      if ( front_.size() == 0 )
        return true;

      if ( n_elements > 0 && n_generated_ == n_elements )
        return true;
    }

    return true;

  } // TriangulationStrategy::exhaustive_search_loop()


  /*------------------------------------------------------------------
  | Search through all advancing front edges for a possible triangle
  ------------------------------------------------------------------*/
  bool exhaustive_search_triangle(Edge& base_edge)
  {
    std::size_t i = 0;

    Edge* next_edge = base_edge.get_next_edge();

    while ( next_edge && next_edge != &base_edge && i < front_.size() )
    {
      Vertex& v = next_edge->v2();

      if ( v != base_edge.v1() && v != base_edge.v2() )
        if ( front_update_.update_front_exhaustive(base_edge, v) )
          return true;

      ++i;
      next_edge = next_edge->get_next_edge();
    }

    return false;

  } // TriangulationStrategy::exhaustive_search_triangle()

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  size_t n_elements_         = 0;
  double mesh_range_factor_  = 1.0;
  double base_vertex_factor_ = 1.5;
  double wide_search_factor_ = 10.0;
  int    n_generated_        = 0;

}; // TriangulationStrategy

} // namespace TQAlgorithm
} // namespace TQMesh
