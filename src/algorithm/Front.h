/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <list>
#include <stdexcept>
#include <memory>
#include <vector>
#include <algorithm>
#include <numeric>

#include "VecND.h"
#include "Geometry.h"

#include "utils.h"
#include "EdgeList.h"
#include "Vertex.h"
#include "Boundary.h"
#include "Domain.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* The advancing front - defined by a list of edges
* > Must be defined counter-clockwise
*********************************************************************/
class Front : public EdgeList
{
  using BoolVector   = std::vector<bool>;
  using IntVector    = std::vector<int>;
  using VertexVector = std::vector<Vertex*>;
  using EdgeVector   = std::vector<Edge*>;
  using BdryEdgeConn = std::vector<std::vector<EdgeVector>>;
  using NbrMeshConn  = std::vector<std::vector<EdgeVector>>;
  using FrontData    = std::pair<EdgeVector,BoolVector>;

public:

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Front() : EdgeList( Orientation::NONE ) {}
  ~Front() {}

  /*------------------------------------------------------------------
  |
  ------------------------------------------------------------------*/
  template <typename FrontInitializer>
  void init_front_edges(const Domain&            domain, 
                        const FrontInitializer&  front_data,
                        Vertices&                mesh_vertices)
  {
    BoolVector edges_to_refine {};

    for ( size_t i_bdry = 0; i_bdry < domain.size(); ++i_bdry )
    {
      const EdgeVector& front_edges  = front_data.edges()[i_bdry];
      const BoolVector& is_twin_edge = front_data.is_twin_edge()[i_bdry];
      const IntVector&  markers      = front_data.markers()[i_bdry];

      VertexVector new_vertices {};

      for ( size_t i = 0; i < front_edges.size(); ++i )
      {
        Edge* e = front_edges[i];

        Vertex& v1 = ( !is_twin_edge[i] ) ? e->v1() : e->v2();

        Vertex& v_new = mesh_vertices.push_back( v1.xy(), 
                                                 v1.sizing(), 
                                                 v1.range() );
        v_new.on_front( true );
        v_new.on_boundary( true );
        v_new.is_fixed( true );

        new_vertices.push_back( &v_new );
      }

      int n_verts = static_cast<int>( new_vertices.size() );

      for ( size_t i_edge = 0; i_edge < front_edges.size(); ++i_edge )
      {
        int i1 = static_cast<int>( i_edge );
        int i2 = MOD( i1+1, n_verts );

        Vertex* v1 = new_vertices[i1];
        Vertex* v2 = new_vertices[i2];

        Edge& e_new = this->add_edge( *v1, *v2, markers[i_edge] );

        // Connect the boundary edge of the neighbor mesh
        // and the new front edge
        if ( is_twin_edge[i_edge] )
        {
          Edge* twin_edge = front_edges[i_edge];
          e_new.twin_edge( twin_edge );
          twin_edge->twin_edge( &e_new );
          edges_to_refine.push_back( false );
        }
        else
        {
          edges_to_refine.push_back( true );
        }
      }
    }

    // Refine the front edges, but do not refine sub-edges!
    this->refine(domain, edges_to_refine, mesh_vertices);
    
  } // Front::init_front_edges()

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  Edge& base() { return *base_; }
  const Edge& base() const { return *base_; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  void base(Edge& b) { base_ = &b; }

  /*------------------------------------------------------------------
  | Refine the advancing front 
  | Returns the number of new edges
  ------------------------------------------------------------------*/
  int refine(const Domain& domain, 
             const BoolVector& edges_to_refine,
             Vertices& vertices)
  {
    int n_before = edges_.size();

    std::vector<Edge*> marked {};

    size_t i_edge_to_refine = 0;

    // Loop over all edge segments
    for ( const auto& e_ptr : edges_ )
    {
      // Ignore edges that should not be refined
      if ( !edges_to_refine[i_edge_to_refine++] )
        continue;

      Edge& cur_edge = *e_ptr.get();

      const double rho_1 = domain.size_function( cur_edge.v1().xy() );
      const double rho_2 = domain.size_function( cur_edge.v2().xy() );

      // Define local edge direction from vertex v_b to
      // vertex v_a, such that rho_a < rho_b
      bool dir = (rho_1 < rho_2);

      // Create coordinates of new vertices along the 
      // current edge segment
      std::vector<Vec2d> xy_new 
        = create_sub_vertex_coords(cur_edge, dir, rho_1, rho_2, domain);

      // No new vertices have been added 
      // --> continue with next boundary segment
      // Otherwise mark this edge to be removed later on
      if ( xy_new.size() < 3 )
        continue;
      else
        marked.push_back( &cur_edge );

      // Create new vertices and edges
      create_sub_edges( cur_edge, xy_new, vertices );
    }

    // Remove old boundary segments
    for ( auto cur_edge : marked )
      edges_.remove( *cur_edge );

    // Re-compute area of domain
    compute_area();

    return ( edges_.size() - n_before );

  } // Front::refine()

  /*------------------------------------------------------------------
  | Let base point to first edge in container
  ------------------------------------------------------------------*/
  void set_base_first()
  { 
    if (edges_.size() < 1) return;
    base_ = edges_.begin()->get();  
  } 

  /*------------------------------------------------------------------
  | Let base point to next edge in container
  ------------------------------------------------------------------*/
  void set_base_next()
  {
    if (edges_.size() < 1) return;
    auto iter = base_->pos();
    std::advance( iter, 1 );
    if (iter == edges_.end())
      iter = edges_.begin();
    base_ = iter->get();
  }

  /*------------------------------------------------------------------
  | Sort all edges by length in ascending order
  | and lets the base segment point to the first edge
  ------------------------------------------------------------------*/
  void sort_edges(bool ascending = true)
  {
    // Sort by edge lengths in ascending order
    if ( ascending )
      edges_.sort(
      []( std::unique_ptr<Edge>& a, std::unique_ptr<Edge>& b )
      {
        return a->length() < b->length();
      });
    // Sort by edge lengths in descending order
    else
      edges_.sort(
      []( std::unique_ptr<Edge>& a, std::unique_ptr<Edge>& b )
      {
        return a->length() > b->length();
      });

    // Reset base segment
    set_base_first();

  } // Front::sort_edges()

private:

  /*------------------------------------------------------------------
  | Mark new vertices and edges to be part of the
  | advancing front
  ------------------------------------------------------------------*/
  virtual 
  void mark_objects(Vertex& v1, Vertex& v2, Edge& e) override
  {
    v1.on_front( true );
    v2.on_front( true );
  }

  /*------------------------------------------------------------------
  | Create a vector of new vertex coordinates along a front edge. 
  | The vertices are distributed according to the underlying domain 
  | size function. The first coordinates in the returning vector 
  | corresponds to the starting vertex xy_a and the last coordinates 
  | to the ending vertex xy_b. The start and ending of the edge 
  | segment are set through the flag <dir>, such that
  | rho(x_a) < rho(x_b)
  ------------------------------------------------------------------*/
  std::vector<Vec2d> create_sub_vertex_coords(const Edge& e, 
                                              bool   dir,
                                              double rho_1, 
                                              double rho_2,
                                              const Domain& domain)
  {
    // Define local edge direction from vertex v_b to
    // vertex v_a, such that rho_a < rho_b
    const Vertex& v_a = dir ? e.v1() : e.v2();
    const Vertex& v_b = dir ? e.v2() : e.v1();

    // Edge tangent unit vector
    const Vec2d tang = dir ? e.tangent() : -e.tangent(); 

    // Allocate vectors for storage of new vertex coords
    std::vector<Vec2d> xy_new { v_a.xy() };
    double s_last = 0.0;

    // Compute point on abscissa where no new points 
    // will be generated
    const double rho_b = dir ? rho_2 : rho_1;
    const double s_end = 1.0 - 0.5 * rho_b / e.length();

    // Compute new vertex positions
    Vec2d xy { v_a.xy() };

    while ( true )
    {
      // Predictor
      const double rho = domain.size_function( xy );
      const Vec2d xy_p = xy + rho * tang;

      // Corrector
      const double rho_p = domain.size_function( xy_p );
      const Vec2d  dxy_c = 0.5 * ( rho + rho_p ) * tang;
      const Vec2d  xy_c  = xy + dxy_c;

      const double l = ( xy_c - v_a.xy() ).norm();
      const double s = l / e.length();

      xy_new.push_back( xy_c );
      s_last = s;
      xy = xy_c;

      // Stop if segment length is reached
      if ( s > s_end ) break;
    }

    // Set last vertex to x_b
    xy_new[xy_new.size()-1] = v_b.xy();

    // Compute cropped distance
    const Vec2d d_cr = (1.0 - s_last) * e.length() * tang;

    // Compute size function lengths for every vertex 
    // Neglect start and ending vertices
    std::vector<double> rho_i { 0.0 };
    for ( int i = 1; i < xy_new.size()-1; i++)
      rho_i.push_back( domain.size_function( xy_new[i] ) );
    rho_i.push_back( 0.0 );

    // Compute total length from size function
    const double rho_tot = std::accumulate(
        rho_i.begin(), rho_i.end(), 0.0 );
    // Divide lengths by total length
    std::transform( rho_i.begin(), rho_i.end(), rho_i.begin(),
        [rho_tot](double &v){ return v/rho_tot; });

    // Distribute cropped distance among new vertices
    for ( int i = 1; i < xy_new.size()-1; i++)
      xy_new[i] += rho_i[i] * d_cr;

    // Check that all nodes are ordered ascendingly
#ifndef NDEBUG
    double s_prev = 0.0;
    for ( int i = 1; i < xy_new.size(); i++ )
    {
      const double s = ( xy_new[i] - xy_new[0] ).norm();
      ASSERT( s > s_prev, "ADVANCING FRONT REFINEMENT FAILED." );
      s_prev = s;
    }
#endif
    
    // Adjust order of new vertices 
    if ( !dir ) 
      std::reverse( xy_new.begin(), xy_new.end() );

    return xy_new;

  } // create_sub_vertex_coords()


  /*------------------------------------------------------------------
  | Divide a given edge segment into several sub-segements,
  | providing an array of vertices that are aliged on the 
  | given segment. This function also creates the 
  | recpective vertices for the new edge segments.
  ------------------------------------------------------------------*/
  void create_sub_edges(Edge& e, 
                        std::vector<Vec2d>& xy_new, 
                        Vertices& vertices)
  {
    Vertex* v_cur = &( e.v1() );

    for ( int i = 1; i < xy_new.size()-1; i++ )
    {
      Vertex& v_n = vertices.insert( e.v2().pos(), xy_new[i], 1.0 );

      // We fix the position of all new vertices on the front,
      // such that they won't be shifted during grid smoothing 
      v_n.is_fixed(true);

      Edge& e_new = this->insert_edge(e.pos(), *v_cur, v_n, e.marker());
      e_new.v1().on_boundary( true );
      e_new.v2().on_boundary( true );

      v_cur = &v_n;
    }

    Edge& e_new = this->insert_edge( e.pos(), *v_cur, e.v2(), e.marker() );
    e_new.v1().on_boundary( true );
    e_new.v2().on_boundary( true );

  } // create_sub_edges()

  /*------------------------------------------------------------------
  | Attributes 
  ------------------------------------------------------------------*/
  Edge*   base_ = nullptr;

}; // Front

} // namespace TQAlgorithm
} // namespace TQMesh
