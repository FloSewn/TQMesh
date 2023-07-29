/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <vector>

#include "utils.h"
#include "Vertex.h"
#include "Edge.h"
#include "Triangle.h"
#include "Quad.h"
#include "Domain.h"
#include "Boundary.h"
#include "Mesh.h"
#include "FrontAlgorithm.h"
#include "Cleanup.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* A structure that contains the data for a quad layer generation
*
*                   v1_proj_[0]       v1_proj_[1]       v1_proj_[2]
*  v1_proj_[0]      v1_proj_[1]       v1_proj_[2]
*    ^----------------^-----------------^-----------------^-----...
*    |                |                 |                 |
*    |                |                 |                 |
*    |                |                 |                 |
*    | base_edges_[0] |  base_edges_[1] |  base_edges_[2] |
*    o----------------o-----------------o-----------------o-----...
*  v1_base_[0]      v1_base_[1]       v1_base_[2]            
*                   v2_base_[0]       v2_base_[1]       v2_base_[2]
*
*********************************************************************/
class QuadLayerVertices
{
public:

  using DoubleVector = std::vector<double>;
  using Vec2dVector  = std::vector<Vec2d>;
  using VertexVector = std::vector<Vertex*>;
  using EdgeVector   = std::vector<Edge*>;

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  QuadLayerVertices(Edge& e_start, Edge& e_end, 
                    bool is_closed, double height)
  : e_start_   { &e_start }
  , e_end_     { &e_end   }
  , is_closed_ { is_closed }
  , height_    { height }
  {
    // Initialize edges
    Edge* e_cur = e_start_;

    do 
    {
      ASSERT( e_cur, "QuadLayerVertices::QuadLayerVertices(): "
        "Advancing front data structure seems to be corrupted.");

      add_quadlayer_edge( *e_cur );
      e_cur = e_cur->get_next_edge();

    } while ( e_cur != e_end_ );

    // Add also the ending edge
    add_quadlayer_edge( *e_end_ );

  } // QuadLayerVertices()

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  Edge* e_start() const { return e_start_; }
  Edge* e_end() const { return e_end_; }
  bool is_closed() const { return is_closed_; }
  double height() const { return height_; }
  int n_base_edges() const 
  { 
    // Cast number of edges to int for use of modulo function  
    return static_cast<int>( base_edges_.size() ); 
  }

  const EdgeVector& base_edges() const { return base_edges_; }
  EdgeVector& base_edges() { return base_edges_; }

  const VertexVector& v1_base() const { return v1_base_; }
  const VertexVector& v2_base() const { return v2_base_; }
  const VertexVector& v1_proj() const { return v1_proj_; }
  const VertexVector& v2_proj() const { return v2_proj_; }

  VertexVector& v1_base() { return v1_base_; }
  VertexVector& v2_base() { return v2_base_; }
  VertexVector& v1_proj() { return v1_proj_; }
  VertexVector& v2_proj() { return v2_proj_; }

  const Vec2dVector& v1_proj_xy() const { return v1_proj_xy_; }
  Vec2dVector& v1_proj_xy() { return v1_proj_xy_; }

  const Vec2dVector& v2_proj_xy() const { return v2_proj_xy_; }
  Vec2dVector& v2_proj_xy() { return v2_proj_xy_; }

  const DoubleVector& heights() const { return heights_; }
  DoubleVector& heights() { return heights_; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/

  /*------------------------------------------------------------------
  | Smooth the heights of all quad projections in the layer according
  | to the local size function
  ------------------------------------------------------------------*/
  void smooth_heights(const Domain& domain)
  {
    for ( size_t i = 1; i < heights_.size()-1; ++i )
    {
      const double h1 = heights_[i-1];
      const double h2 = heights_[i];
      const double h3 = heights_[i+1];

      const Vec2d& c = base_edges_[i]->xy();
      const double rho = domain.size_function( c );

      heights_[i] = MIN(rho, (h1+h2+h3) / 3.0);
    }

  } // QuadLayerVertices::smooth_heights()


  /*------------------------------------------------------------------
  | In this function we set up the correct coordinates of all the
  | projected base vertices. In the case of a QuadLayerVertices that 
  | is not closed, adjacent edges are eventually refined. 
  ------------------------------------------------------------------*/
  void setup_vertex_projection(Mesh& mesh, Front& front)
  {
    // Update coordinates of projected base vertices
    for ( size_t i = 1; i < base_edges_.size(); ++i )
      adjust_projected_vertex_coordinates(i-1, i);

    if ( is_closed_ )
    {
      adjust_projected_vertex_coordinates( n_base_edges()-1, 0);
    }
    else
    {
      place_start_vertex(mesh.vertices(), front, mesh.boundary_edges());
      place_end_vertex(mesh.vertices(), front, mesh.boundary_edges());
    }

  } // QuadLayerVertices::setup_vertex_projection()



private:

  /*------------------------------------------------------------------
  | Try to adjust the initial projected base vertex coordinates,
  | such that adjacent base edges have the same coordinates of their
  | projected vertices.
  | However, in some cases this would lead to bad elements. In these
  | cases, the respective base is marked as a "wedge", where we will
  | introduce triangles later on, in order to obtain a better quality
  |
  |                   q                    r
  |                    o------------------o
  |                   /   base_edges_[j]
  |                  /
  |                 /
  |                / base_edges_[i] 
  |               / 
  |              / 
  |             o p 
  |
  ------------------------------------------------------------------*/
  void adjust_projected_vertex_coordinates(int i, int j)
  {
    const Vec2d& p = v1_base_[i]->xy();
    const Vec2d& q = v1_base_[j]->xy(); 
    const Vec2d& r = v2_base_[j]->xy();

    const double alpha = angle( p-q, r-q );

    // If both projected vertices are too far apart, we must 
    // create a wedge in between. In this case we use the default 
    // projection coordinates v1_proj_xy and v2_proj_xy that were already 
    // calculated in the constructor
    if ( is_left(p, r, q) && alpha <= quad_layer_angle_ )
      return;

    // Otherwise, the projected vertex coordinate will be placed
    // in between the originally projected coordinates 
    const Vec2d& n1 = base_edges_[i]->normal();
    const double l1 = heights_[i];

    const Vec2d& n2 = base_edges_[j]->normal();
    const double l2 = heights_[j];

    const Vec2d  normal  = 0.5 * (n1 + n2);
    const double l       = 0.5 * (l1 + l2);
    const Vec2d  nn      = normal / normal.norm();

    Vec2d xy_proj = q + nn * l / sin(0.5*alpha);

    v1_proj_xy_[j] = xy_proj;
    v2_proj_xy_[i] = xy_proj;

    return;

  } // QuadLayerVertices::adjust_projected_vertex_coordinates()


  /*------------------------------------------------------------------
  | To prevent degenerate triangles, we check if the projected 
  | vertex of the first quad layer element is located too close 
  | to the previous advancing front edge. We check for the
  | following three scenarios:
  | 1) The distance between v_prev and v1_proj is less than the 
  |    local quad layer height h. In this case we will simply use 
  |    v_prev as v1_proj
  | 2) v_prev and v1_proj are located very close to each other and 
  |    the previous edge length is greater than the quad layer height.
  |    In this case, the previous edge must be split at the location
  |    of v1_proj projected onto the previous edge. This leads to the 
  |    generation of a new vertex v_new, that will be used as v1_proj 
  | 3) v_prev and v1_proj are located very close to each other and 
  |    and the previous edge length is smaller than the quad layer 
  |    height. In this case, v1_proj is set to v_prev
  | 
  |
  |   Scenario 1              Scenario 2            Scenario 3
  |   ----------              ----------            ----------
  |            h
  |      <------------>                               v1_proj                     
  | v_prev    v1_proj    v_prev                          o           
  |     o      o           o                             |            
  |     |     /             \                    v_prev  |           
  |     |    /               \   v1_proj            o    |        
  |     |   /                 x   o                  \   x       
  |     |  /             v_new \  |                   \  |   
  |     | /                     \ |                    \ |  
  |     |/                       \|                     \| 
  |     o------------o            o-----------o          o-----------o
  |  v_start                  v_start                v_start           
  |
  ------------------------------------------------------------------*/
  void place_start_vertex(Vertices& mesh_vertices, 
                          Front&    front, 
                          EdgeList& mesh_bdry_edges)
  {
    // The edge adjacent to the starting edge
    Edge* e_prev = e_start_->get_prev_edge();

    // The current starting base vertex
    Vertex& v_start = *v1_base_[0];

    ASSERT( (e_prev) && (e_prev->v2() == v_start),
    "During the generation of a QuadLayerVertices, an invalid data\n"
    "structure of the advancing front has been provided.\n"
    "It seems that the given advancing front is not a connected\n"
    "list of edges.\n");

    // The adjacent previous vertex 
    Vertex& v_prev  = e_prev->v1();

    // If the previous vertex is located right to the starting edge,
    // we use the default projection coordinate
    if ( !is_left(v1_base_[0]->xy(), v2_base_[0]->xy(), v_prev.xy()) )
      return;

    const double h = heights_[0];
    const double dist_sqr = (v_prev.xy() - v1_proj_xy_[0]).norm_sqr();

    // Scenario 1
    // ----------
    if ( dist_sqr < h * h )
    {
      // Shift v_prev along the edge e_prev, such that it has the same
      // height as v1_proj_xy, then set v1_proj_xy to this coordinate
      Vec2d delta_1 = v_prev.xy() - v1_base_[0]->xy(); 
      Vec2d delta_2 = v1_proj_xy_[0] - v1_base_[0]->xy(); 
      double s = dot(delta_1, delta_2) / delta_1.norm_sqr();
      Vec2d v_prev_new = v1_base_[0]->xy() + delta_1 * s;
      Cleanup::set_vertex_coordinates(v_prev, v_prev_new);
      v1_proj_xy_[0] = v_prev_new;
      return;
    }

    // Scenario 2
    // ----------
    if ( h < e_prev->length() )
    {
      // Get a pointer to the previous edge in the mesh data structure,
      // since it will be removed later on
      Edge* e_bdry_rem 
        = mesh_bdry_edges.get_edge(e_prev->v1(), e_prev->v2());

      // Split the adjacent front edge in two smaller edges
      const Vec2d d1 = v_prev.xy() - v_start.xy();
      const Vec2d d2 = v1_proj_xy_[0] - v_start.xy();
      const double alpha = angle( d1, d2 );
      const double ang_fac = cos(alpha); 
      double sf = (h * ang_fac) / e_prev->length();

      auto new_edges 
        = front.split_edge(*e_prev, mesh_vertices, sf, false);

      // Remove the splitted edge in the mesh data structure
      if ( e_bdry_rem )
      {
        Edge* tmp = e_bdry_rem->get_next_edge();
        ASSERT( tmp, "QuadLayerVertices::place_start_vertex(): "
          "Front data structure seems to be corrupted.");
        
        mesh_bdry_edges.remove( *e_bdry_rem );
        e_bdry_rem = tmp;

        Edge* e1 = new_edges.first;
        Edge* e2 = new_edges.second;
        mesh_bdry_edges.insert_edge(e_bdry_rem->pos(),
                                    e1->v1(), e1->v2(), e1->marker());
        mesh_bdry_edges.insert_edge(e_bdry_rem->pos(),
                                    e2->v1(), e2->v2(), e2->marker());
      }

      // Add new vertex to quad layer structure
      Vertex& v_new = new_edges.first->v2();
      v1_proj_xy_[0] = v_new.xy();

      return;
    }
      
      
    // Scenario 3
    // ----------
    v1_proj_xy_[0] = v_prev.xy();

    return;

  } // QuadLayerVertices::place_start_vertex()


  /*------------------------------------------------------------------
  | 
  |   Scenario 1           Scenario 2            Scenario 3
  |   ----------           ----------            ----------
  |       h
  |  <---------->           
  |    v2_proj  v_next                 v_next       v2_proj
  |       o      o                       o             o
  |        \     |                      /              |
  |         \    |             v2_proj /               |   v_next
  |          \   |                o   x                x   o 
  |           \  |                |  / v_new           |  /
  |            \ |                | /                  | /
  |             \|                |/                   |/
  |  o-----------o       o--------o           o--------o
  |            v_end            v_end                 v_end
  |
  ------------------------------------------------------------------*/
  void place_end_vertex(Vertices& verts, 
                        Front&    front, 
                        EdgeList& bdry_edges)
  {
    // The edge adjacent to the starting edge
    Edge* e_next = e_end_->get_next_edge();

    // The current ending base vertex
    Vertex& v_end = *v2_base_.back();

    ASSERT( (e_next) && (e_next->v1() == v_end),
    "During the generation of a QuadLayerVertices, an invalid data\n"
    "structure of the advancing front has been provided.\n"
    "It seems that the given advancing front is not a connected\n"
    "list of edges.\n");

    // The adjacent previous vertex 
    Vertex& v_next  = e_next->v2();

    // If the next vertex is located right to the starting edge,
    // we use the default projection coordinate
    const Vec2d& xy_end_1 = v1_base_.back()->xy();
    const Vec2d& xy_end_2 = v2_base_.back()->xy();
    //
    if ( !is_left(xy_end_1, xy_end_2, v_next.xy()) )
      return;

    // Check if the segment between v_start and its projected 
    // coordinate intersects with the previous edge 
    // If yes, merge them
    const double h = heights_.back();
    const double dist_sqr = (v_next.xy() - v2_proj_xy_.back()).norm_sqr();

    // Scenario 1
    // ----------
    if ( dist_sqr < h * h )
    {
      // Shift v_next along the edge e_prev, such that it has the same
      // height as v2_proj_xy, then set v2_proj_xy to this coordinate
      Vec2d delta_1 = v_next.xy() - v_end.xy(); 
      Vec2d delta_2 = v2_proj_xy_.back() - v_end.xy(); 
      double s = dot(delta_1, delta_2) / delta_1.norm_sqr();
      Vec2d v_next_new = v_end.xy() + delta_1 * s;
      Cleanup::set_vertex_coordinates(v_next, v_next_new);
      v2_proj_xy_.back() = v_next_new;
      return;
    }

    // Scenario 2
    // ----------
    if ( h < e_next->length() )
    {
      // Get a pointer to the next edge in the mesh data structure,
      // since it will be removed later on
      Edge* e_bdry_rem = bdry_edges.get_edge(e_next->v1(), e_next->v2());

      // Split the adjacent front edge in two smaller edges
      const Vec2d d1 = v_next.xy() - v_end.xy();
      const Vec2d d2 = v2_proj_xy_.back() - v_end.xy();
      const double alpha = angle( d1, d2 );
      const double ang_fac = cos(alpha); 
      double sf = 1.0 - (h * ang_fac) / e_next->length();

      auto new_edges = front.split_edge(*e_next, verts, sf, false);

      // Remove the splitted edge in the mesh data structure
      if ( e_bdry_rem ) 
      {
        Edge* tmp = e_bdry_rem->get_next_edge();
        ASSERT( tmp, "QuadLayerVertices::place_end_vertex(): "
          "Front data structure seems to be corrupted.");
        bdry_edges.remove( *e_bdry_rem );
        e_bdry_rem = tmp;

        Edge* e1 = new_edges.first;
        Edge* e2 = new_edges.second;
        bdry_edges.insert_edge( e_bdry_rem->pos(),
                                e1->v1(), e1->v2(), e1->marker() );
        bdry_edges.insert_edge( e_bdry_rem->pos(),
                                e2->v1(), e2->v2(), e2->marker() );
      }

      // Add new vertex to quad layer structure
      Vertex& v_new = new_edges.first->v2();
      v2_proj_xy_.back() = v_new.xy();

      return;
    }
      
      
    // Scenario 3
    // ----------
    v2_proj_xy_.back() = v_next.xy();

    return;

  } // QuadLayerVertices::place_end_vertex()

  /*------------------------------------------------------------------
  | Add a new edge to the quad layer structure
  ------------------------------------------------------------------*/
  void add_quadlayer_edge(Edge& e_cur)
  {
    // Pointer to base edge
    base_edges_.push_back( &e_cur );

    // Pointers to vertices of current base edge 
    v1_base_.push_back( &e_cur.v1() );
    v2_base_.push_back( &e_cur.v2() );

    // Adjust height, in order to get good aspect ratios
    double h = MIN( height_, e_cur.length() );
    heights_.push_back( h );

    // Coordinates of initial projected base vertices 
    v1_proj_xy_.push_back( e_cur.v1().xy() + e_cur.normal() * h );
    v2_proj_xy_.push_back( e_cur.v2().xy() + e_cur.normal() * h );

    // Since the actual projected vertices are not generated yet,
    // we introduced them in terms of nullpointers
    v1_proj_.push_back( nullptr );
    v2_proj_.push_back( nullptr );

  } // QuadLayerVertices::add_quadlayer_edge()

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Edge*        e_start_   {nullptr};
  Edge*        e_end_     {nullptr};
  bool         is_closed_ { false };
  double       height_    {  0.0  };

  EdgeVector   base_edges_ {};
  VertexVector v1_base_    {};
  VertexVector v2_base_    {};

  VertexVector v1_proj_    {};
  VertexVector v2_proj_    {};

  Vec2dVector  v1_proj_xy_ {};
  Vec2dVector  v2_proj_xy_ {};

  DoubleVector heights_    {};


  double quad_layer_angle_ = 1.57079633; // = 1/2 pi

}; // QuadLayerVertices




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
    add_remaining_front_edges_to_mesh();

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
    create_quad_layer( quad_layer_verts );

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

  } // FrontQuadLayering::find_next_layer_endings()

  /*------------------------------------------------------------------
  | This function creates a new layer of quad elements   
  ------------------------------------------------------------------*/
  void create_quad_layer(QuadLayerVertices& quad_layer_verts)
  {
    auto& v1_base    = quad_layer_verts.v1_base();
    auto& v2_base    = quad_layer_verts.v2_base();

    auto& v1_proj    = quad_layer_verts.v1_proj();
    auto& v2_proj    = quad_layer_verts.v2_proj();

    auto& v1_proj_xy = quad_layer_verts.v1_proj_xy();
    auto& v2_proj_xy = quad_layer_verts.v2_proj_xy();

    auto& base_edges = quad_layer_verts.base_edges();

    int n_base_edges = quad_layer_verts.n_base_edges();

    for ( int i = 0; i < n_base_edges; ++i )
    {
      DEBUG_LOG("QUAD LAYER BASE " << i);

      // Check if edge has already been put to waste 
      if ( !base_edges[i]->in_container() )
        continue;

      // Create new element
      auto v_proj 
        = create_quad_layer_element(*base_edges[i],
                                    *v1_base[i], v1_proj_xy[i],
                                    *v2_base[i], v2_proj_xy[i]);

      // Store pointer to newly generated vertices
      v1_proj[i] = v_proj.first;
      v2_proj[i] = v_proj.second;
    }

  } // FrontQuadLayering::create_quad_layer() 


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
  std::pair<Vertex*, Vertex*>
  create_quad_layer_element(Edge& base,
                            Vertex& v1_base, const Vec2d& v1_proj_xy, 
                            Vertex& v2_base, const Vec2d& v2_proj_xy)
  {
    std::pair<Vertex*, Vertex*> v_proj { nullptr, nullptr };


    // Create first triangle (v1_base, v2_base, v1_proj)
    double r1 = quad_layer_range_ * (v1_base.xy() - v1_proj_xy).norm();
    Triangle* t1 = front_update_.update_front(base, v1_proj_xy, 
                                              v1_proj_xy, r1);
    if ( !t1 )
      return v_proj;

    Vertex& v_proj_p1 = t1->v3();
    v_proj_p1.is_fixed(true);
    v_proj.first = &v_proj_p1;


    // Use the newly generated triangle edge as base for 
    // the second triangle. However, this edge may not be part of the 
    // advancing front anymore, if the mesh is already degenerated
    Edge* diagonal = front_.get_edge( v_proj_p1, v2_base );
    if ( !diagonal ) 
      return v_proj;


    // Create second triangle (v1_proj, v2_base, v2_proj)
    double r2 = quad_layer_range_ * (v2_base.xy() - v2_proj_xy).norm();
    Triangle* t2 
      = front_update_.update_front(*diagonal, v2_proj_xy, 
                                   v2_proj_xy, r2);
    if ( !t2 )
      return v_proj;

    Vertex& v_proj_p2 = t2->v3();
    v_proj_p2.is_fixed(true);
    v_proj.second = &v_proj_p2;


    // Merge both triangles t1 & t2 to a quad
    // --> First remove the interior edge between these triangles
    Edge* e_rem = mesh_.interior_edges().get_edge( v2_base, v_proj_p1 );

    if ( !e_rem ) 
      return v_proj;

    mesh_.remove_interior_edge( *e_rem );

    // Remove old triangular elements
    mesh_.remove_triangle( *t1 );
    mesh_.remove_triangle( *t2 );

    // Create new quadrilateral element
    Quad& q_new = mesh_.add_quad(v1_base, v2_base, v_proj_p2, v_proj_p1);
    q_new.is_active( true );

    return v_proj;

  } // FrontQuadLayering::create_quad_layer_element()


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
  double quad_layer_range_ = 0.9;

  // Attributes that change during the layer generation
  Vec2d  xy_start_     {};
  Vec2d  xy_end_       {};
  Edge*  e_start_      {nullptr};
  Edge*  e_end_        {nullptr};
  bool   closed_layer_ {false};



}; // FrontQuadLayering

} // namespace TQAlgorithm
} // namespace TQMesh
