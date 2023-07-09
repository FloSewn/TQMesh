/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "Vertex.h"
#include "Edge.h"
#include "Facet.h"
#include "Boundary.h"
#include "Front.h"
#include "Domain.h"
#include "Mesh.h"
#include "FrontInitializer.h"
#include "MeshValidator.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;



/*********************************************************************
* 
*********************************************************************/
class FrontAlgorithm 
{
public:
  using VertexVector   = std::vector<Vertex*>;
  using TriVector      = std::vector<Triangle*>;

  /*------------------------------------------------------------------
  | Constructor / Destructor 
  ------------------------------------------------------------------*/
  FrontAlgorithm(Mesh& mesh, const Domain& domain)
  : mesh_ { mesh }
  , domain_ { domain }
  , validator_ {mesh, domain, front_} 
  {}

  virtual ~FrontAlgorithm() {}

  /*------------------------------------------------------------------
  | Triangulate a given initialized mesh structure
  ------------------------------------------------------------------*/
  virtual bool generate_elements(int n_elements=0) = 0;

protected:

  /*------------------------------------------------------------------
  | This function is part of the advancing front loop.
  | For a given location <xy> and a respective search range <dist>,
  | all vertices that are located in this vicinity are put into a 
  | vector and the sorted ascendingly towards <xy>.
  ------------------------------------------------------------------*/
  VertexVector find_vertex_candidates(const Vec2d& xy, double dist) 
  {
    Vertices& vertices = mesh_.vertices();

    // Get vertices in vicinity of xy  
    VertexVector vertex_candidates = vertices.get_items(xy, dist);

    // Sort vertices in ascending order towards xy
    std::sort( vertex_candidates.begin(), vertex_candidates.end(), 
    [xy] ( const Vertex* a, const Vertex* b )
    {
      return ( (a->xy()-xy).norm_sqr() 
             < (b->xy()-xy).norm_sqr() );
    });

    return std::move( vertex_candidates ); 

  } // find_vertex_candidates() 

  /*------------------------------------------------------------------
  | This function is part of the advancing front loop.
  | We sort a given vector of <new_triangles> in descending order
  | according to the triangle quality.
  | Finally, the advancing front is updated with the triangle of best
  | quality and all other triangles are removed.
  ------------------------------------------------------------------*/
  Triangle* choose_best_triangle(TriVector& new_triangles,
                                 Edge&      base)
  {
    Triangles& triangles = mesh_.triangles();

    if ( new_triangles.size() < 1 )
      return nullptr;

    DEBUG_LOG("VALID TRIANGLES IN NEIGHBORHOOD: " 
       << new_triangles.size()
    );

    std::sort( new_triangles.begin(), new_triangles.end(),
    [this] ( Triangle* t1, Triangle* t2 )
    {
      const double h1 = domain_.size_function( t1->xy() );
      const double h2 = domain_.size_function( t2->xy() );
      const double q1 = t1->quality(h1);
      const double q2 = t2->quality(h2);

      return ( q1 > q2 );
    });

    Triangle* new_tri = new_triangles[0];
    Vertex&   v_adj   = new_tri->v3();

    update_front( base, v_adj, *new_tri );

    for (int i = 1; i < new_triangles.size(); i++)
      triangles.remove( *new_triangles[i] );

    return new_tri;

  } // choose_best_triangle()

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  void update_front(Edge& base, Vertex& v_new, Triangle& t_new)
  {
    // Get advancing front edges adjacent to vertex
    // -> First two vertices of new triangle tri are always
    //    the base edge vertices
    Edge* e1 = front_.get_edge(v_new, base.v1());
    Edge* e2 = front_.get_edge(v_new, base.v2());

    // *** Both edges are connected to vertex ***
    //     -> No new edge must be created
    //     -> Base vertex v1 no longer on the advancing front
    //     -> Base vertex v2 no longer on the advancing front
    //     -> vertex no longer on the advancing front
    if ( e1 && e2 )
    {
      ASSERT( v_new.on_front(), "Grid structure corrupted." );

      base.v1().on_front( false );
      base.v2().on_front( false );
      v_new.on_front( false );

      if ( e1->is_interior() )
        mesh_.add_interior_edge(e1->v1(), e1->v2());
      if ( e2->is_interior() )
        mesh_.add_interior_edge(e2->v1(), e2->v2());

      front_.remove( *e1 );
      front_.remove( *e2 );
    }
    // *** First edge is connected to vertex ***
    //     -> New edge between second base vertex and vertex
    //     -> Base vertex v1 no longer on the advancing front
    //     -> vertex is part of the advancing front
    else if ( e1 && !e2 )
    {
      ASSERT( v_new.on_front(), "Grid structure corrupted." );

      base.v1().on_front( false );

      if ( e1->is_interior() )
        mesh_.add_interior_edge(e1->v1(), e1->v2());

      front_.remove( *e1 );
      front_.add_edge(v_new, base.v2());
    }
    // *** Second edge is connected to vertex ***
    //     -> New edge between first base vertex and vertex
    //     -> Base vertex v2 no longer on the advancing front
    //     -> vertex is part of the advancing front
    else if ( !e1 && e2 )
    {
      ASSERT( v_new.on_front(), "Grid structure corrupted." );

      base.v2().on_front( false );

      if ( e2->is_interior() )
        mesh_.add_interior_edge(e2->v1(), e2->v2());

      front_.remove( *e2 );
      front_.add_edge(base.v1(), v_new);
    }
    // *** Both edges are not connected to vertex ***
    //     -> Create two new edges
    //     -> Vertex now part of the advancing front
    else
    {
      v_new.on_front( true );
      front_.add_edge(base.v1(), v_new);
      front_.add_edge(v_new, base.v2());
    }

    // If current base is not at the boundary, add it to the 
    // interior edge list
    if ( base.is_interior() )
      mesh_.add_interior_edge(base.v1(), base.v2());

    // Remove base edge
    front_.remove( base );

    // Mark new triangle as active
    t_new.is_active( true );

    // Add element area to the total mesh area
    mesh_.add_area( t_new.area() );

  } // update_front() 

  /*------------------------------------------------------------------
  | This function is part of the advancing front loop.
  | We loop over a given set of <vertex_candidates> and check
  | if we can create a possible triangle with the current base edge
  | (<b1>,<b2>) and the given vertices.
  | If it is possible, new triangles are created and pushed back to 
  | the vector <new_triangles>.
  ------------------------------------------------------------------*/
  void check_vertex_candidates(const VertexVector& vertex_candidates,
                               Edge& base_edge, TriVector& new_triangles)
  {
    for ( Vertex* v : vertex_candidates )
    {
      // Skip vertices that are not located on the advancing front
      if ( !v->on_front() )
        continue;

      // Skip vertices that are colinear to the current base edge
      if ( orientation( base_edge.v1().xy(), base_edge.v2().xy(), v->xy() )
          == Orientation::CL )
        continue;

      // Create new potential triangle 
      Triangle& t_new = mesh_.add_triangle( base_edge.v1(), base_edge.v2(), *v );

      // Check if new potential triangle is valid
      if ( !validator_.remove_from_mesh_if_invalid(t_new) )
        new_triangles.push_back( &t_new );
    }

  } // check_vertex_candidates()

  /*------------------------------------------------------------------
  | Initialize the advancing front structure
  ------------------------------------------------------------------*/
  Edge* init_advancing_front()
  {
    front_.init_front(mesh_);
    Edge* base = front_.set_base_first();
    ASSERT( base, "FrontTriangulation::generate_elements(): "
      "Invalid advancing front structure.");
    front_.sort_edges( false );

    return base;
  }

  /*------------------------------------------------------------------
  | Remove invalid mesh edges (which stem from previous
  | meshing approaches)
  ------------------------------------------------------------------*/
  void remove_invalid_mesh_edges()
  {
    auto invalid_mesh_edges = mesh_.get_invalid_edges();
    for ( auto& e_ptr : invalid_mesh_edges )
      mesh_.remove_interior_edge( *e_ptr );
  }

  /*------------------------------------------------------------------
  | Prepare the mesh for the output - this consists mainly of the 
  | addition of remaining advancing front edges, such that the 
  | mesh process can be continued in a subsequent step
  ------------------------------------------------------------------*/
  void finish_mesh_for_output()
  {
    // Setup mesh connectivity
    Cleanup::setup_vertex_connectivity(mesh_);
    Cleanup::setup_facet_connectivity(mesh_);

    // Add remaining front edges to the mesh 
    for ( auto& e_ptr : front_.edges() )
    {
      const Vertex& v1 = e_ptr->v1();
      const Vertex& v2 = e_ptr->v2();

      if ( mesh_.interior_edges().get_edge(v1, v2) )
        continue;

      if ( mesh_.boundary_edges().get_edge(v1, v2) )
        continue;

      mesh_.add_interior_edge(e_ptr->v1(), e_ptr->v2());
    }
  }

  /*------------------------------------------------------------------
  | Update the progress bar
  ------------------------------------------------------------------*/
  void update_progress_bar()
  {
    double state = std::ceil(100.0 * mesh_.area() / domain_.area());
    progress_bar_.update( static_cast<int>(state) );
    progress_bar_.show( LOG_PROPERTIES.get_ostream(INFO) );
  }

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Mesh&         mesh_;
  const Domain& domain_;
  MeshValidator validator_;

  Front         front_ {};
  ProgressBar   progress_bar_ {};

}; // FrontAlgorithm

} // namespace TQAlgorithm
} // namespace TQMesh
