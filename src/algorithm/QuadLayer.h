/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <algorithm>

#include "Vec2.h"
#include "utils.h"

#include "Vertex.h"
#include "Triangle.h"
#include "Quad.h"
#include "Front.h"
#include "Boundary.h"
#include "Domain.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace TQUtils;

/*********************************************************************
* A structure that contains the data for a quad layer generation
*********************************************************************/
class QuadLayer
{

public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  QuadLayer(Edge* e_start, Edge* e_end, bool is_closed, double height)
  : e_start_   { e_start }
  , e_end_     { e_end   }
  , is_closed_ { is_closed }
  , height_    { height }
  {
    // Every base edge gets an associated QuadProjection, which 
    // is basically a container for the generation of a new 
    // quad element
    Edge* e_cur = e_start_;

    do 
    {
      init_edge( e_cur );
      e_cur = e_cur->get_next_edge();

      ASSERT( e_cur, "INVALID DATA STRUCTURE FOR QUAD LAYER");

    } while ( e_cur != e_end_ );

    // Add also the ending edge
    init_edge( e_end );

    // Cast number of edges to int for use of modulo function  
    n_bases_ = static_cast<int>( bases_.size() );


  } // QuadLayer()



  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  Edge* e_start() const { return e_start_; }
  Edge* e_end() const { return e_end_; }
  bool is_closed() const { return is_closed_; }
  double height() const { return height_; }
  int n_bases() const { return n_bases_; }

  const std::vector<Edge*>& bases() const { return bases_; }
  std::vector<Edge*>& bases() { return bases_; }

  const std::vector<Vertex*>& b1() const { return b1_; }
  const std::vector<Vertex*>& b2() const { return b2_; }
  const std::vector<Vertex*>& p1() const { return p1_; }
  const std::vector<Vertex*>& p2() const { return p2_; }

  std::vector<Vertex*>& b1() { return b1_; }
  std::vector<Vertex*>& b2() { return b2_; }
  std::vector<Vertex*>& p1() { return p1_; }
  std::vector<Vertex*>& p2() { return p2_; }

  const std::vector<Vec2d>& p1_xy() const { return p1_xy_; }
  std::vector<Vec2d>& p1_xy() { return p1_xy_; }

  const std::vector<Vec2d>& p2_xy() const { return p2_xy_; }
  std::vector<Vec2d>& p2_xy() { return p2_xy_; }

  const std::vector<double>& heights() const { return heights_; }
  std::vector<double>& heights() { return heights_; }

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

      const Vec2d& c = bases_[i]->xy();
      const double rho = domain.size_function( c );

      heights_[i] = MIN(rho, (h1+h2+h3) / 3.0);
    }

  } // QuadLayer::smooth_heights()


  /*------------------------------------------------------------------
  | In this function we set up the correct coordinates of all the
  | projected base vertices. In the case of a QuadLayer that is not
  | closed, adjacent edges are eventually refined. 
  ------------------------------------------------------------------*/
  void setup_vertex_projection(Vertices& verts, 
                               Front&    front, 
                               EdgeList& bdry_edges)
  {
    // Update coordinates of projected base vertices
    for ( size_t i = 1; i < bases_.size(); ++i )
      project_base_vertices(i-1, i);

    if ( is_closed_ )
    {
      project_base_vertices(n_bases_-1, 0);
    }
    else
    {
      place_start_vertex(verts, front, bdry_edges);
      place_end_vertex(verts, front, bdry_edges);
    }

  } // QuadLayer::setup_vertex_projection()



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
  |                    *------------------*
  |                   /      bases_[j]
  |                  /
  |                 /
  |                / bases_[i] 
  |               / 
  |              / 
  |             * p 
  |
  ------------------------------------------------------------------*/
  void project_base_vertices(int i, int j)
  {
    const Vec2d& p = b1_[i]->xy();
    const Vec2d& q = b1_[j]->xy(); 
    const Vec2d& r = b2_[j]->xy();

    const double alpha = angle( p-q, r-q );

    // If both projected vertices are too far apart, we must 
    // create a wedge in between. In this case we use the default 
    // projection coordinates p1_xy and p2_xy that were already 
    // calculated in the constructor
    if ( TQGeom::is_left(p, r, q) && alpha <= TQ_QUAD_LAYER_ANGLE )
      return;

    // Otherwise, the projected vertex coordinate will be placed
    // in between the originally projected coordinates 
    const Vec2d& n1 = bases_[i]->normal();
    const double l1 = heights_[i];

    const Vec2d& n2 = bases_[j]->normal();
    const double l2 = heights_[j];

    const Vec2d  norm  = 0.5 * (n1 + n2);
    const double l     = 0.5 * (l1 + l2);
    const Vec2d  nn    = norm / norm.length();

    Vec2d xy_proj = q + nn * l / sin(0.5*alpha);

    p1_xy_[j] = xy_proj;
    p2_xy_[i] = xy_proj;

    return;

  } // QuadLayer::project_base_vertices()


  /*------------------------------------------------------------------
  |
  ------------------------------------------------------------------*/
  void place_start_vertex(Vertices& verts, 
                          Front&    front, 
                          EdgeList& bdry_edges)
  {
    // The edge adjacent to the starting edge
    Edge* e_prv = e_start_->get_prev_edge();

    // The current starting base vertex
    Vertex& v_start = *b1_[0];

    ASSERT( (e_prv) && (e_prv->v2() == v_start),
    "During the generation of a QuadLayer, an invalid data\n"
    "structure of the advancing front has been provided.\n"
    "It seems that the given advancing front is not a connected\n"
    "list of edges.\n");

    // The adjacent previous vertex 
    Vertex& v_prev  = e_prv->v1();

    // If the previous vertex is located right to the starting edge,
    // we use the default projection coordinate
    if ( !TQGeom::is_left(b1_[0]->xy(), b2_[0]->xy(), v_prev.xy()) )
      return;

    // Check if the segment between v_start and its projected 
    // coordinate intersects with the previous edge 
    // If yes, merge them
    const double h     = heights_[0];
    const double d_fac = (v_prev.xy() - p1_xy_[0]).length() / h;

    if ( d_fac < 1.0 )
    {
      p1_[0] = &v_prev;
      return;
    }

    // If the projected vertex coordinate is located within the range 
    // of the previous edge, it must be projected onto it
    //
    //                            v_prev
    //    x                     x 
    //     \                   /
    //      \                /
    //       \     p1      /
    //        \   o      /   
    //         o  |    o
    //          \ |  /      
    //           \|/
    //            x---------------------x
    //            x                 
    //              v_start
    //
    //
    if ( h < e_prv->length() )
    {
      const Vec2d d1 = v_prev.xy() - v_start.xy();
      const Vec2d d2 = p1_xy_[0] - v_start.xy();
      const double alpha = angle( d1, d2 );
      const double ang_fac = cos(alpha); // / sin(0.5*alpha);


      // Remove the old edge from the boundary edge list (if contained)
      Edge* e_bdry_rem = bdry_edges.get_edge(e_prv->v1(), e_prv->v2());

      if ( e_bdry_rem ) 
      {
        Edge* tmp = e_bdry_rem->get_next_edge();
        bdry_edges.remove( *e_bdry_rem );
        e_bdry_rem = tmp;
        ASSERT( e_bdry_rem, "INVALID DATA STRUCTURE" );
      }
    
      // Split the adjacent edge in two small edges
      double sf = (h * ang_fac) / e_prv->length();

      auto new_edges = front.split_edge(*e_prv, verts, sf, false);

      // Add potential new boundary edges
      if ( e_bdry_rem )
      {
        Edge* e1 = new_edges.first;
        Edge* e2 = new_edges.second;
        bdry_edges.insert_edge( e_bdry_rem->pos(),
                                e1->v1(), e1->v2(), e1->marker() );
        bdry_edges.insert_edge( e_bdry_rem->pos(),
                                e2->v1(), e2->v2(), e2->marker() );
      }

      // Add new vertex to quad layer structure
      Vertex& v_new = new_edges.first->v2();

      p1_[0]    = &v_new;
      p1_xy_[0] = v_new.xy();

    }
    // Projected vertex is located outside of adjacent edge range
    // --> Set p to be a
    //
    //            o p1          
    //            |             
    //            |           
    //       x    |       v_prev
    //        \   x     x   
    //         \  |    /
    //          \ |  /
    //           \|/
    //            x---------------------x
    //            v_start           
    else
    {
      p1_[0]    = &v_prev;
      p1_xy_[0] = v_prev.xy();
    }

    return;

  } // QuadLayer::place_start_vertex()


  /*------------------------------------------------------------------
  |
  ------------------------------------------------------------------*/
  void place_end_vertex(Vertices& verts, 
                        Front&    front, 
                        EdgeList& bdry_edges)
  {
    // The edge adjacent to the starting edge
    Edge* e_nxt = e_end_->get_next_edge();

    // The current ending base vertex
    Vertex& v_end = *b2_.back();

    ASSERT( (e_nxt) && (e_nxt->v2() == v_end),
    "During the generation of a QuadLayer, an invalid data\n"
    "structure of the advancing front has been provided.\n"
    "It seems that the given advancing front is not a connected\n"
    "list of edges.\n");

    // The adjacent previous vertex 
    Vertex& v_next  = e_nxt->v2();

    // If the next vertex is located right to the starting edge,
    // we use the default projection coordinate
    const Vec2d& xy_end_1 = b1_.back()->xy();
    const Vec2d& xy_end_2 = b2_.back()->xy();
    //
    if ( !TQGeom::is_left(xy_end_1, xy_end_2, v_next.xy()) )
      return;

    // Check if the segment between v_start and its projected 
    // coordinate intersects with the previous edge 
    // If yes, merge them
    const double h     = heights_.back();
    const double d_fac = (v_next.xy() - p2_xy_.back()).length() / h;

    if ( d_fac < 1.0 )
    {
      p2_.back() = &v_next;
      return;
    }

    // If the projected vertex coordinate is located within the range 
    // of the previous edge, it must be projected onto it
    //
    //            v_next
    //           x                     x 
    //            \                   /
    //             \                /
    //              \     p2      /
    //               \   o      /   
    //                o  |    o
    //                 \ |  /    
    //                  \|/
    //   x---------------x
    //         base        v_end       
    //
    //
    if ( h < e_nxt->length() )
    {
      const Vec2d d1 = v_next.xy() - v_end.xy();
      const Vec2d d2 = p2_xy_.back() - v_end.xy();
      const double alpha = angle( d1, d2 );
      const double ang_fac = cos(alpha); // / sin(0.5*alpha);

      // Remove the old edge from the boundary edge list (if contained)
      Edge* e_bdry_rem = bdry_edges.get_edge(e_nxt->v1(), e_nxt->v2());

      if ( e_bdry_rem ) 
      {
        Edge* tmp = e_bdry_rem->get_next_edge();
        bdry_edges.remove( *e_bdry_rem );
        e_bdry_rem = tmp;
        ASSERT( e_bdry_rem, "INVALID DATA STRUCTURE" );
      }
    
      // Split the adjacent edge in two small edges
      double sf = 1.0 - (h * ang_fac) / e_nxt->length();

      auto new_edges = front.split_edge(*e_nxt, verts, sf, false);

      // Add potential new boundary edges
      if ( e_bdry_rem )
      {
        Edge* e1 = new_edges.first;
        Edge* e2 = new_edges.second;
        bdry_edges.insert_edge( e_bdry_rem->pos(),
                                e1->v1(), e1->v2(), e1->marker() );
        bdry_edges.insert_edge( e_bdry_rem->pos(),
                                e2->v1(), e2->v2(), e2->marker() );
      }

      // Add new vertex to quad layer structure
      Vertex& v_new = new_edges.first->v2();

      p2_.back() = &v_new;
      p2_xy_.back() = v_new.xy();

    }
    // Projected vertex is located outside of adjacent edge range
    // --> Set p to be a
    //
    //             p2
    //            o               
    //            |             
    //            |           
    //       x    |       v_next  
    //        \   x     x   
    //         \  |    /
    //          \ |  /
    //           \|/
    //  x---------x
    //             v_end
    //
    else
    {
      p2_.back()    = &v_next;
      p2_xy_.back() = v_next.xy();
    }

    return;

  } // QuadLayer::place_end_vertex()



  /*------------------------------------------------------------------
  | Add a new edge to the quad layer structure
  ------------------------------------------------------------------*/
  void init_edge(Edge* e_cur)
  {
    // Pointer to base edge
    bases_.push_back( e_cur );

    // Pointers to base vertices
    Vertex* v1 = &e_cur->v1();
    Vertex* v2 = &e_cur->v2();

    b1_.push_back( v1 );
    b2_.push_back( v2 );

    // Adjust height, in order to get good aspect ratios
    double h = MIN( height_, e_cur->length() );
    heights_.push_back( h );

    // Coordinates of initial projected base vertices 
    // --> Vertices are not generated yet
    p1_.push_back( nullptr );
    p2_.push_back( nullptr );
    p1_xy_.push_back( v1->xy() + e_cur->normal() * h );
    p2_xy_.push_back( v2->xy() + e_cur->normal() * h );


  } // QuadLayer::init_edge()



  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Edge*                     e_start_   {nullptr};
  Edge*                     e_end_     {nullptr};
  bool                      is_closed_ { false };
  double                    height_    {  0.0  };
  int                       n_bases_   { 0 };

  std::vector<Edge*>        bases_    {};
  std::vector<Vertex*>      b1_       {};
  std::vector<Vertex*>      b2_       {};

  std::vector<Vertex*>      p1_       {};
  std::vector<Vertex*>      p2_       {};

  std::vector<Vec2d>        p1_xy_    {};
  std::vector<Vec2d>        p2_xy_    {};

  std::vector<double>       heights_  {};

}; // QuadLayer


} // namespace TQAlgorithm
} // namespace TQMesh
