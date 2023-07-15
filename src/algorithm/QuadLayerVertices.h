/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <algorithm>

#include "VecND.h"

#include "utils.h"
#include "Vertex.h"
#include "Triangle.h"
#include "Quad.h"
#include "Front.h"
#include "Boundary.h"
#include "Domain.h"
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


} // namespace TQAlgorithm
} // namespace TQMesh
