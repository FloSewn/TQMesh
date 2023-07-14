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
#include "QuadLayerVertices.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* 
*********************************************************************/
class FrontQuadLayering : public FrontAlgorithm
{
public:
  using VertexVector   = std::vector<Vertex*>;
  using TriVector      = std::vector<Triangle*>;

  /*------------------------------------------------------------------
  | Constructor / Destructor 
  ------------------------------------------------------------------*/
  FrontQuadLayering(Mesh& mesh, const Domain& domain)
  : FrontAlgorithm(mesh, domain) {}

  ~FrontQuadLayering() {}

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  size_t n_layers() const { return n_layers_; }
  double first_height() const { return first_height_; }
  double growth_rate() const { return growth_rate_; }
  const Vec2d& starting_position() const { return xy_start_; }
  const Vec2d& ending_position() const { return xy_end_; }

  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  void n_layers(size_t n) { n_layers_ = n; }
  void first_height(double h) { first_height_ = h; }
  void growth_rate(double r) { growth_rate_ = r; }
  void starting_position(const Vec2d& v) { xy_start_ = v; }
  void starting_position(double x, double y) { xy_start_ = {x,y}; }
  void ending_position(const Vec2d& v) { xy_end_ = v; }
  void ending_position(double x, double y) { xy_end_ = {x,y}; }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  bool generate_elements() override
  {
    if (mesh_.n_boundary_edges() < 1)
      return false;

    // Prepare the mesh  
    Cleanup::setup_facet_connectivity(mesh_);
    
    // Initialize the advancing front and its base edge
    init_advancing_front(false);

    // Remove invalid mesh edges that are no longer needed
    remove_invalid_mesh_edges();

    // Perform the actual mesh generation
    double height = first_height_;
    bool success = true;
    for ( size_t i_layer = 0; i_layer < n_layers_; ++i_layer)
    {
      success = generate_quad_layer(height);

      if (!success) break;

      height *= growth_rate_;
    }

    // Finish mesh structure for output
    finish_mesh_for_output();

    // Remove remaining edges from the front
    front_.clear_edges();

    return success;

  } // generate_elements()

private:

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  bool generate_quad_layer(double height)
  {
    // Determine the edges in the advancing front, that correspond 
    // to the current start and ending positions for the quad layer 
    if ( !find_start_and_ending_edges(xy_start_, xy_end_) )
      return false;

    // Create the quad layer structure, which keeps track of the target
    // vertex coordinates, that are projected from the base vertex 
    // coordinates
    QuadLayerVertices quad_layer_verts { *e_start_, *e_end_, 
                                   closed_layer_, height };
    quad_layer_verts.smooth_heights( domain_ );
    quad_layer_verts.setup_vertex_projection( mesh_, front_ );

    // For each base edge in the quad layer, try to create a quad
    // element with its given projected coordinates
    create_quad_layer_elements( quad_layer_verts );

    // Triangulate the quad layer based edges, where the generation
    // of quads did not succeed
    finish_quad_layer( quad_layer_verts );

    // Remove deleted entities
    mesh_.clear_waste();

    // Set new start and ending vertex coordinates
    return find_next_layer_endings( quad_layer_verts );

  } // generate_quad_layer()

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  bool find_next_layer_endings(QuadLayerVertices& quad_layer_verts)
  {
    // Search for starting vertex
    Vertex* v1 = nullptr;
    int n = quad_layer_verts.n_base_edges();

    for (int i = 0; i < n; ++i)
    {
      v1 = quad_layer_verts.v1_proj()[i];

      if ( v1 && v1->on_front() )
        break;
    }

    // Search for ending vertex
    Vertex* v2 = v1;

    if ( !closed_layer_ )
      for (int i = n-1; i >= 0; --i)
      {
        v2 = quad_layer_verts.v2_proj()[i];

        if ( v2 && v2->on_front() )
          break;
      }

    // Check if vertices were found
    if ( !v1 || !v2 )
      return false;

    xy_start_ = v1->xy();
    xy_end_   = v2->xy();

    return true;

  } // FrontQuadLayering::find_new_layer_endings()

  /*------------------------------------------------------------------
  |   
  |            v1_proj       v2_proj
  |              x-------------x-------------
  |              | \           | \          |
  |              |   \         |   \        |
  |              |     \       |     \      |
  |              |       \     |       \    |
  |              |         \   |         \  |
  |              | base_edge \ |           \|
  |     ---------x-------------x------------x-------
  |            v1_base       v2_base
  |   
  ------------------------------------------------------------------*/
  void create_quad_layer_elements(QuadLayerVertices& quad_layer_verts)
  {
    auto& v1_base    = quad_layer_verts.v1_base();
    auto& v2_base    = quad_layer_verts.v2_base();

    auto& v1_proj    = quad_layer_verts.v1_proj();
    auto& v2_proj    = quad_layer_verts.v2_proj();

    auto& v1_proj_xy = quad_layer_verts.v1_proj_xy();
    auto& v2_proj_xy = quad_layer_verts.v2_proj_xy();

    auto& heights    = quad_layer_verts.heights();
    auto& base_edges = quad_layer_verts.base_edges();

    int n_base_edges = quad_layer_verts.n_base_edges();

    for ( int i = 0; i < n_base_edges; ++i )
    {
      DEBUG_LOG("QUAD LAYER BASE " << i);

      // Search radius for vertices in the vicinity of the 
      // projected coordinates
      const double r = quad_layer_range_ * heights[i];



      // Create first triangle (v1_base,v2_base,v1_proj)
      Edge& base = *base_edges[i];

      // Check if edge has already been put to waste 
      if (!base.in_container())
        continue;

      Triangle* t1 
        = front_update_.update_front(base, v1_proj_xy[i], 
                                     v1_proj_xy[i], r);

      if ( t1 == nullptr ) 
        continue;

      Vertex& v_proj_p1 = t1->v3();



      // Create second triangle (v1_proj,v2_base,v2_proj)
      Edge* diagonal = front_.get_edge( v_proj_p1, *v2_base[i] );

      if ( !diagonal ) 
        continue;

      Triangle* t2 
        = front_update_.update_front(*diagonal, v2_proj_xy[i], 
                                     v2_proj_xy[i], r);

      if ( t2 == nullptr )
        continue;

      Vertex& v_proj_p2 = t2->v3();



      // Merge both triangles t1 & t2 to a quad
      // --> First remove the interior edge between these triangles
      Edge* e_rem = mesh_.interior_edges().get_edge( *v2_base[i], v_proj_p1 );

      if ( e_rem == nullptr ) 
        continue;

      mesh_.remove_interior_edge( *e_rem );

      // Remove old triangular elements
      mesh_.remove_triangle( *t1 );
      mesh_.remove_triangle( *t2 );

      // Create new quadrilateral element
      Quad& q_new = mesh_.add_quad( *v1_base[i], *v2_base[i], v_proj_p2, v_proj_p1 );
      q_new.is_active( true );



      // Store pointers to projection vertices
      v1_proj[i] = &v_proj_p1;
      v2_proj[i] = &v_proj_p2;
    }

  } // FrontQuadLayering::create_quad_layer_elements() 

  /*------------------------------------------------------------------
  | In this step we treat specific quad layer bases, where the former
  | algorithm failed.
  ------------------------------------------------------------------*/
  void finish_quad_layer(QuadLayerVertices& quad_layer_verts)
  {
    auto& v1_base        = quad_layer_verts.v1_base();
    auto& v1_proj        = quad_layer_verts.v1_proj();
    auto& v2_proj        = quad_layer_verts.v2_proj();

    int  n_base_edges = quad_layer_verts.n_base_edges();

    for ( int i = 1; i < n_base_edges; ++i )
    {
      if ( !v1_proj[i] || !v2_proj[i-1] || v1_proj[i] == v2_proj[i-1] )
        continue;
      close_quad_layer_gap(*v2_proj[i-1], *v1_base[i], *v1_proj[i]);
    }

  } // FrontQuadLayering::finish_quad_layer()

  /*------------------------------------------------------------------
  | In some cases, gaps might be formed during the initial quad layer
  | generation step. 
  | Here, these gaps are closed with triangular elements.
  |
  |                            v_c     
  |              v_new       o 
  |                 o        :
  |                          :
  |                          :
  |           v_a            : v_b
  |            o.............o-------------------o
  |                          |      base[i]                         
  |                          |
  |                base[i-1] |
  |                          |
  |                          o
  |
  |              -> v_a = v2_proj[i-1]
  |              -> v_b = v1_base[i]
  |              -> v_c = v1_proj[i]
  ------------------------------------------------------------------*/
  void close_quad_layer_gap(Vertex& v_a, Vertex& v_b, Vertex& v_c)
  {
    const Vec2d l1 = v_a.xy()-v_b.xy();
    const Vec2d l2 = v_c.xy()-v_b.xy();
    const double alpha = angle(l1,l2);

    // Don't add v_a new vertex and instead only connect (v_a,v_b,v_c)
    if ( alpha <= quad_layer_angle_ )
    {
      Triangle* t_new = &( mesh_.add_triangle(v_a, v_b, v_c) );

      if ( front_update_.remove_from_mesh_if_invalid(*t_new) )
        return;
       
      Edge* base = front_.get_edge( v_b, v_c );
      front_update_.advance_front( *base, v_a, *t_new );

      return;
    }

    // Create new vertex and then generate two triangles
    const Vec2d v_xy = v_b.xy() + l1 + l2;

    Vertex& v_new = mesh_.add_vertex( v_xy );

    Triangle& t1_new = mesh_.add_triangle(v_a, v_b, v_new);
    Triangle& t2_new = mesh_.add_triangle(v_b, v_c, v_new);

    if ( front_update_.remove_from_mesh_if_invalid(v_new, t1_new, t2_new) )
      return;

    Edge* base = nullptr;

    base = front_.get_edge( v_a, v_b );
    ASSERT( base, "FrontQuadLayering::close_quad_layer_gabs(): "
      "Front data structure seems to be corrupted.");
    front_update_.advance_front( *base, v_new, t1_new );

    base = front_.get_edge( v_b, v_c );
    ASSERT( base, "FrontQuadLayering::close_quad_layer_gabs(): "
      "Front data structure seems to be corrupted.");
    front_update_.advance_front( *base, v_new, t2_new );

    v_new.is_fixed( true );

    return;

  } // close_quad_layer_gap()
  


  /*------------------------------------------------------------------
  | Find the starting end ending edge in the advancing front that 
  | correspond to the currently set start and ending coordinets 
  ------------------------------------------------------------------*/
  bool find_start_and_ending_edges(const Vec2d& xy_start,
                                   const Vec2d& xy_end)
  {
    // Find closest vertices in current front structure to  
    // current start and ending vertex coordinates
    Vertex& v_start = front_.get_closest_vertex( xy_start );
    Vertex& v_end   = front_.get_closest_vertex( xy_end );

    // Get advancing front edges adjacent to input vertices
    Edge* e_start = front_.get_edge(v_start, 1); 
    Edge* e_end   = front_.get_edge(v_end, 2); 

    ASSERT((e_start && e_end), 
      "FrontQuadLayering::find_start_and_ending_edges(): "
      "Vertex-edge-connectivity seems to be corrupted.");

    if ( !front_.is_traversable(*e_start, *e_end) )
      return false;

    bool is_closed = (v_start == v_end);

    // For closed quad layers, try not to start at sharp angle edges
    if ( is_closed )
    {
      const Vec2d& v1 = e_end->v1().xy();
      const Vec2d& v2 = e_end->v2().xy();
      const Vec2d& v3 = e_start->v2().xy();

      const double ang = angle(v1-v2, v3-v2);

      Edge *e_next = e_start->get_next_edge();

      if ( e_next && ang <= quad_layer_angle_ )
      {
        e_end = e_start;
        e_start = e_next; 
      }
    }

    // Finally, set the found edges 
    e_start_      = e_start;
    e_end_        = e_end;
    closed_layer_ = is_closed;

    return true;

  } // find_start_and_ending_edges()

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  size_t n_layers_     {};
  double first_height_ {};
  double growth_rate_  {};

  // Meshing constants
  double quad_layer_angle_ = 1.57079633; // = 1/2 pi
  double quad_layer_range_ = 0.75;

  // Attributes that change during the layer generation
  Vec2d  xy_start_     {};
  Vec2d  xy_end_       {};
  Edge*  e_start_      {nullptr};
  Edge*  e_end_        {nullptr};
  bool   closed_layer_ {false};



}; // FrontQuadLayering

} // namespace TQAlgorithm
} // namespace TQMesh
