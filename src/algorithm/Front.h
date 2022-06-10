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

#include "Vec2.h"
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
  using EdgeVector = std::vector<Edge*>;

public:

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Front() : EdgeList( Orientation::NONE )
  {}

  /*------------------------------------------------------------------
  | Init the front edges from a domain and refine them
  | Since the advancing front is part if the mesh, this function
  | requires the mesh's vertices to setup the edges
  ------------------------------------------------------------------*/
  void init_front_edges(const Domain& domain, Vertices &mesh_vertices)
  {
    // Check for correct number of vertices of domain and mesh
    ASSERT( domain.vertices().size() == mesh_vertices.size(), 
        "Domain and mesh vertex numbers differ." );

    // Initialize indices of original domain vertices
    size_t v_index = 0; 
    for ( const auto& v_ptr : domain.vertices() )
      v_ptr->index( v_index++ );

    // Initialize array with pointers to the mesh vertices
    std::vector<Vertex*> v_ptrs {};

    for ( const auto& v_ptr : mesh_vertices )
      v_ptrs.push_back( v_ptr.get() );


    // Create list edge-vertex-connectivity given in the domain
    std::vector<std::vector<std::pair<size_t,size_t>>> connectivity {};

    size_t i_bdry = 0;
    for ( const auto& boundary : domain )
    {
      connectivity.push_back( {} );

      for ( const auto& e : boundary->edges() )
      {
        connectivity[i_bdry].push_back( 
            {e->v1().index(), e->v2().index()} 
        );
      }

      // Counter for boundaries
      ++i_bdry;
    }


    // Initialize the advancing front edges
    i_bdry = 0;

    for ( const auto& boundary : domain )
    {
      size_t i_edge = 0;

      for ( const auto& e : boundary->edges() )
      {
        // Get mesh vertices of current edge
        size_t i1 = connectivity[i_bdry][i_edge].first;
        size_t i2 = connectivity[i_bdry][i_edge].second;

        Vertex* v1 = v_ptrs[i1];
        Vertex* v2 = v_ptrs[i2];

        // All advancing front edges are assigned to a specified 
        // edge marker
        this->add_edge( *v1, *v2, e->marker() );

        // We fix the position of all initial boundary vertices,
        // such that they won't be shifted during the grid 
        // smoothing process
        v1->is_fixed(true);
        v2->is_fixed(true);

        // Keep track of the original domain boundary
        v1->on_boundary( true );
        v2->on_boundary( true );

        // Counter for boundary edges
        ++i_edge;
      }

      // Counter for boundaries
      ++i_bdry;
    }


    // Check for correct orientation
    ASSERT( check_orientation(), "Invalid edge list orientation." );

    // Refine the front based on the underlying size function
    this->refine(domain, mesh_vertices);


    /*************************************
     * DEPRECATED
    // Collect all front edges of every domain boundary edge
    for ( const auto& boundary : domain )
    {
      for ( const auto& e : boundary->edges() )
      {
        // Add new edge to the boundary edge's data structure
        BdryEdgeData* bdry_data = e->bdry_data();

        ASSERT( bdry_data != nullptr, 
            "Error: The BdryEdgeData structure of a domain-edge "
            "seems not to be properly initialized.");

        // Find all front edges that are located in the vicinity 
        // of the current domain boundary edge
        const Vec2d& c = e->xy();
        double       r = TQMeshEdgeSearchFactor * e->length();

        const Vec2d& v = e->v1().xy();
        const Vec2d& w = e->v2().xy();

        std::vector<Edge*> found = edges_.get_items(c, r);
        
        // Get all edges, that are actually located on the segment
        // of the current domain boundary edge and add them to
        // its boundary data structure
        for ( Edge* e_found : found )
        {
          const Vec2d& p = e_found->v1().xy();
          const Vec2d& q = e_found->v2().xy();

          if ( in_on_segment(v,w,p) && in_on_segment(v,w,q) )
            bdry_data->sub_edges.push_back( e_found );
        }

        // Sort all edges in the boundary data structure in 
        // ascending order to the starting vertex
        std::sort(bdry_data->sub_edges.begin(), bdry_data->sub_edges.end(), 
        [v](Edge* e1, Edge* e2)
        {
          double delta_1 = (e1->xy() - v).length_squared();
          double delta_2 = (e2->xy() - v).length_squared();
          return (delta_1 < delta_2);
        });
      }
    }
    *******************************************/

  } // init_front_edges()

  /*------------------------------------------------------------------
  | Init the front edges from a given domain, as well as a 
  | corresponding connectivity of edges, that are part of another mesh
  | Since the advancing front is part if the mesh, this function
  | requires the mesh's vertices to setup the edges
  ------------------------------------------------------------------*/
  void init_front_from_edges(const Domain& domain, 
                             const std::vector<std::vector<EdgeVector>>& edge_conn,
                             Vertices &mesh_vertices)
  {
    // Stage 1 
    // ----------
    // Loop over all boundaries and their boundary edges of the domain.
    // |
    // |-If a domain boundary edge contains sub-edges from another 
    // | mesh: loop over all these sub-edges
    // | |
    // | |-Create a new vertex for every sub-edge and store a pointer
    // |   Fort this, take that vertex of the sub-edge, which is closer
    // |   to the starting-vertex of the current domain boundary edge
    // |   (because sub-edge might be aligned in the wrong direction)
    // |
    // |-Else: Use starting vertex of current domain boundary edge
    // |       and initialize a new vertex from it. 
    // |       Store a pointer to it.


    std::vector<Vertex*> new_vertices {};

    size_t i_bdry = 0;

    for ( const auto& boundary : domain )
    {
      size_t i_edge = 0;

      for ( const auto& e : boundary->edges() )
      {
        // This is the starting vertex of the current boundary edge
        Vertex& v1_bdry = e->v1();

        // Check if there are given sub-edges 
        if ( edge_conn[i_bdry][i_edge].size() > 0 )
        {
          for ( Edge* sub_edge : edge_conn[i_bdry][i_edge] )
          {
            const Vec2d& v1_xy = sub_edge->v1().xy();
            const Vec2d& v2_xy = sub_edge->v2().xy();

            // Choose vertex that is closer to starting vertex
            // of current boundary edge
            double delta_1 = ( v1_xy - v1_bdry.xy() ).length_squared();
            double delta_2 = ( v2_xy - v1_bdry.xy() ).length_squared();

            Vertex& v1 = (delta_1 < delta_2) 
                       ? sub_edge->v1()
                       : sub_edge->v2();

            Vertex& v_new = mesh_vertices.push_back( v1.xy(), 
                                                     v1.sizing(), 
                                                     v1.range() );
            v_new.on_front( true );
            v_new.on_boundary( true );
            v_new.is_fixed( true );

            new_vertices.push_back( &v_new );
          }
        }
        // Otherwise, no sub-edges exist for this boundary edge,
        // so it must be initialized from its starting and ending 
        // vertex
        else
        {
          Vertex& v_new = mesh_vertices.push_back( v1_bdry.xy(), 
                                                   v1_bdry.sizing(), 
                                                   v1_bdry.range() );
          v_new.on_front( true );
          v_new.on_boundary( true );
          v_new.is_fixed( true );

          new_vertices.push_back( &v_new );
        }
        ++i_edge;
      }
      ++i_bdry;
    }

    // Stage 2
    // ----------
    // Mark indices of all new vertices


    // Stage 3
    // ----------
    // Create a vertex-edge connectivity (similar to init_front_edges)
    // --> Here we must use the same loop structure as in stage 1,
    //     such that the order of the vertex-pointer-vector and the
    //     vertices in the boundary structure is the same
    //
    // Create list edge-vertex-connectivity given in the domain
    std::vector<std::vector<std::vector<std::pair<size_t,size_t>>>> 
    connectivity {};

    // Stage 4
    // ----------
    // Create edges

    // Stage 5
    // ----------
    // Refine -> But do not refine sub-edges!!!

  } // Front::init_front_from_edges()

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
  int refine(const Domain& domain, Vertices& vertices)
  {
    int n_before = edges_.size();

    std::vector<Edge*> marked {};

    // Loop over all edge segments
    for ( const auto& e_ptr : edges_ )
    {
      Edge* e = e_ptr.get();

      const double rho_1 = domain.size_function( e->v1().xy() );
      const double rho_2 = domain.size_function( e->v2().xy() );

      // Define local edge direction from vertex v_b to
      // vertex v_a, such that rho_a < rho_b
      bool dir = (rho_1 < rho_2);

      // Create coordinates of new vertices along the 
      // current edge segment
      std::vector<Vec2d> xy_new 
        = create_sub_vertex_coords(*e, dir, rho_1, rho_2, domain);

      // No new vertices have been added 
      // --> continue with next boundary segment
      // Otherwise mark this edge to be removed later on
      if ( xy_new.size() < 3 )
        continue;
      else
        marked.push_back( e );

      // Create new vertices and edges
      create_sub_edges( *e, xy_new, vertices );
    }

    // Remove old boundary segments
    for ( auto e : marked )
      edges_.remove( *e );

    // Re-compute area of domain
    compute_area();

    return ( edges_.size() - n_before );

  }

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

      const double l = ( xy_c - v_a.xy() ).length();
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
      const double s = ( xy_new[i] - xy_new[0] ).length();
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
      Vertex& v_n = vertices.insert( e.v2().pos(), 
                                      xy_new[i], 1.0 );

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
