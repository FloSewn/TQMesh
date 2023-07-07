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
//#include "MeshingAlgorithm.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* 
*********************************************************************/
class FrontTriangulation //: public MeshingAlgorithm
{
public:

  using VertexVector   = std::vector<Vertex*>;
  using TriVector      = std::vector<Triangle*>;

  /*------------------------------------------------------------------
  | Constructor / Destructor 
  ------------------------------------------------------------------*/
  FrontTriangulation(Mesh& mesh, const Domain& domain)
  : mesh_ { mesh }
  , domain_ { domain }
  { 
    front_.init_front(mesh);
  }

  ~FrontTriangulation() {}

  /*------------------------------------------------------------------
  | Triangulate a given initialized mesh structure
  ------------------------------------------------------------------*/
  bool generate_elements(Mesh& mesh, const Domain& domain, int n_elements=0) 
  {
    if ( mesh.n_boundary_edges() == 0 )
    {
      LOG(ERROR) << 
      "Unable to triangulate mesh " << mesh.id() << ". " <<
      "The mesh has not been prepared yet.";
      return false;
    }

    // Initialize base edge
    front_.set_base_first();
    Edge* base = &( front_.base() );

    // Invalid base definition
    if ( !base ) 
    {
      LOG(ERROR) << 
      "Unable to triangulate mesh " << mesh.id() << ". " <<
      "The mesh's advancing front structure seems to corrupted.";
      return false;
    }

    LOG(INFO) << "Start triangulation of mesh " << mesh.id() << ".";

    // Sort the advancing front edges
    front_.sort_edges( false );

    // Start advancing front loop
    ProgressBar progress_bar {};
    unsigned int iter = 0;
    bool wide_search = false;
    int sort_iter = 0;
    int n_generated = 0;

    while ( true )
    {
      // Try to advance the current base edge
      bool success = advance_front_triangle(*base, wide_search);

      // If it worked, reset iteration counter and wide search
      // and go to the next base edge
      if ( success )
      {
        ++n_generated;

        // Sort front edges after a wide search
        if ( wide_search )
          front_.sort_edges( false );

        iter = 0;
        wide_search = false;

        // Sort the advancing front edges
        if ( sort_iter == sort_triangulation_front_ )
        {
          front_.sort_edges( false );
          sort_iter = 0;
        }
        else
        {
          sort_iter = MOD( sort_iter+1, INT_MAX-1 );
        }

        front_.set_base_first();
        base = &( front_.base() );

        mesh.clear_waste();
      }
      // If it failed, go to the next base edge
      else
      {
        front_.set_base_next();
        base = &( front_.base() );
        ++iter;
      }

      // All front edges failed to create new elements
      // --> Activate wide search for neighboring vertices
      //     and re-run the algorithm
      if ( iter == front_.size() && !wide_search )
      {
        wide_search = true;
        iter = 0;
      }

      // Update progress bar
      double state = std::ceil(100.0 * mesh.area() / domain.area());
      progress_bar.update( static_cast<int>(state) );
      progress_bar.show( LOG_PROPERTIES.get_ostream(INFO) );

      // No more edges in the advancing front
      // --> Meshing algorithm succeeded
      if ( front_.size() == 0 )
      {
        LOG(INFO) << "[";
        LOG(INFO) << "The meshing process succeeded.";
        break;
      }
      
      // Maximum number of elements has been generated
      if (n_elements > 0 && n_generated == n_elements)
      {
        LOG(INFO) << "[";
        LOG(INFO) << "The meshing process was aborted.";
        break;
      }

      // All front edges faild to create new elements, even
      // when using the wide search 
      // --> Meshing algorithm failed
      if ( iter == front_.size() && wide_search )
      {
        LOG(INFO) << "[";
        LOG(ERROR) << "The meshing process failed.";
        return false;
      }
    }

    Cleanup::assign_size_function_to_vertices(mesh, domain);
    Cleanup::assign_mesh_indices(mesh);

    return true;
  }

private:

  /*------------------------------------------------------------------
  | Check if a triangle is valid. If yes, return true - 
  | else return false.
  ------------------------------------------------------------------*/
  bool triangle_is_valid(const Triangle& tri)
  {
    Vertices&   vertices = mesh_.vertices();
    Triangles& triangles = mesh_.triangles();
    Quads&         quads = mesh_.quads();

    const double rho   = domain_.size_function( tri.xy() );
    const double range = 2.0 * rho;

    DEBUG_LOG("CHECK NEW TRIANGLE: " << tri);

    if ( !tri.is_valid() )
      return false;

    if ( tri.intersects_front( front_, range ) )
    { DEBUG_LOG("  > FRONT INTERSECTION"); return false; }

    if ( tri.intersects_vertex( vertices, range ) )
    { DEBUG_LOG("  > VERTEX INTERSECTION"); return false; }

    if ( tri.intersects_triangle( triangles, range ) )
    { DEBUG_LOG("  > TRIANGLE INTERSECTION"); return false; }

    if ( tri.intersects_quad( quads, range ) )
    { DEBUG_LOG("  > QUAD INTERSECTION"); return false; }

    if ( tri.quality(rho) < min_cell_quality_ )
    { DEBUG_LOG("  > BAD TRIANGLE QUALITY"); return false; }

    if ( tri.max_angle() > max_cell_angle_ )
    { DEBUG_LOG("  > BAD MAXIMUM ANGLE"); return false; }

    DEBUG_LOG("  > VALID");
    return true;

  } // Mesh::triangle_is_valid()

  /*------------------------------------------------------------------
  | Check if a vertex is valid. If yes, return true - 
  | else return false.
  ------------------------------------------------------------------*/
  bool vertex_is_valid(const Vertex& v)
  {
    Triangles& triangles = mesh_.triangles();
    Quads&         quads = mesh_.quads();

    const double rho   = domain_.size_function( v.xy() );
    const double range = 2.0 * rho;

    DEBUG_LOG("CHECK NEW VERTEX: " << v);

    if ( !domain_.is_inside( v ) )
    { DEBUG_LOG("  > OUTSIDE DOMAIN"); return false; }

    if ( v.intersects_facet(triangles, range) )
    { DEBUG_LOG("  > TRIANGLE INTERSECTION"); return false; }

    if ( v.intersects_facet(quads, range) )
    { DEBUG_LOG("  > QUAD INTERSECTION"); return false; }

    DEBUG_LOG("  > VALID");
    return true;

  } // Mesh::vertex_is_valid()

  /*------------------------------------------------------------------
  | Check if mesh entities are invalid. If yes, remove them and 
  | return true. Else, simply return false.
  ------------------------------------------------------------------*/
  bool remove_if_invalid(Vertex& v)
  {
    if ( vertex_is_valid(v) )
      return false;
    mesh_.remove_vertex(v);
    return true;
  } 

  bool remove_if_invalid(Triangle& t)
  {
    if ( triangle_is_valid(t) )
      return false;
    mesh_.remove_triangle(t);
    return true;
  } 

  bool remove_if_invalid(Vertex& v, Triangle& t)
  {
    if ( !vertex_is_valid(v) || !triangle_is_valid(t) )
    {
      mesh_.remove_triangle(t);
      mesh_.remove_vertex(v);
      return true;
    }
    return false;
  } 

  bool remove_if_invalid(Vertex& v, Triangle& t1, Triangle& t2)
  {
    if (  !vertex_is_valid(v) 
       || !triangle_is_valid(t1) || !triangle_is_valid(t2) )
    {
      mesh_.remove_triangle(t1);
      mesh_.remove_triangle(t2);
      mesh_.remove_vertex(v);
      return true;
    }
    return false;
  }

  /*------------------------------------------------------------------
  | This function is part of the advancing front loop.
  | For a given location <xy> and a respective search range <dist>,
  | all vertices that are located in this vicinity are put into a 
  | vector and the sorted ascendingly towards <xy>.
  | If the flag <wide_search> is activated, the search range is 
  | enlarged by a preset factor, in order to finde more vertex 
  | candidates.
  ------------------------------------------------------------------*/
  VertexVector find_local_vertices(const Vec2d& xy, double dist, 
                                   bool wide_search=false )
  {
    Vertices& vertices = mesh_.vertices();

    if (wide_search)
      dist *= wide_search_factor_;

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

  } // find_local_vertices() 

  /*------------------------------------------------------------------
  | This function is part of the advancing front loop.
  | We loop over a given set of <vertex_candidates> and check
  | if we can create a possible triangle with the current base edge
  | (<b1>,<b2>) and the given vertices.
  | If it is possible, new triangles are created and pushed back to 
  | the vector <new_triangles>.
  ------------------------------------------------------------------*/
  void check_vertex_candidates(const VertexVector& vertex_candidates,
                               Edge& base, TriVector& new_triangles)
  {
    for ( Vertex* v : vertex_candidates )
    {
      // Skip vertices that are not located on the advancing front
      if ( !v->on_front() )
        continue;

      // Skip vertices that are colinear to the current base edge
      if ( orientation( base.v1().xy(), base.v2().xy(), v->xy() )
          == Orientation::CL )
        continue;

      // Create new potential triangle 
      Triangle& t_new = mesh_.add_triangle( base.v1(), base.v2(), *v );

      // Check if new potential triangle is valid
      if ( !remove_if_invalid(t_new) )
        new_triangles.push_back( &t_new );
    }
  } // check_vertex_candidates()

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  Vertex& get_base_vertex(const Edge& base)
  {
    // Half of the factor h for height of equlateral triangle
    // h := sqrt(3) / 2  -   h_fac := h / 2
    constexpr double h_fac = 0.4330127019; 
    const double v_fac = h_fac * base_vertex_factor_;

    // Obtain size function value at the centroid of an equlateral
    // triangle, created from the current base edge
    Vec2d c = base.xy() + base.normal() * base.length() * v_fac;
    const double rho = domain_.size_function(c);

    // Coordinate of new vertex 
    Vec2d xy = base.xy() + base.normal() * rho;

    return mesh_.add_vertex( xy );

  } // get_base_vertex() 


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
  | Let the advancing front create a new triangle 
  ------------------------------------------------------------------*/
  bool advance_front_triangle(Edge& base, bool wide_search=false)
  {
    // Constants
    const double f1 = base_height_factor_;
    const double f2 = mesh_range_factor_;

    // Obtain the position of a potential new vertex and the radius 
    // to search around it for potential existing vertices
    const double height = domain_.size_function( base.xy() );
    const Vec2d  v_xy   = base.xy() + f1 * height * base.normal();
    const double r      = domain_.size_function( v_xy );

    // Find all vertices in the vicinity of the current base edge
    VertexVector vertex_candidates 
      = find_local_vertices(v_xy, f2 * r, wide_search);

    // Create potential triangles with all found vertices
    TriVector new_triangles {};
    check_vertex_candidates(vertex_candidates, base, new_triangles);

    // If potential triangles have been found, choose the best one
    if ( choose_best_triangle(new_triangles, base) )
      return true;

    // Check if a potential triangle can be created with the base edge
    // and a newly created vertex
    Vertex& b1 = base.v1();
    Vertex& b2 = base.v2();

    Vertex&   v_new = get_base_vertex( base );
    Triangle& t_new = mesh_.add_triangle(b1, b2, v_new);
    
    // Algorithm fails if new vertex or new triangle is invalid
    if ( remove_if_invalid(v_new, t_new) )
      return false;
    
    // Update the advancing front with new vertex
    update_front(base, v_new, t_new);

    return true;

  } // Mesh::advance_front_triangle() */


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Mesh&         mesh_;
  const Domain& domain_;
  Front         front_ {};

  int    sort_triangulation_front_ = -1;
  double mesh_range_factor_        = 1.0;
  double base_height_factor_       = 0.43; // ~ sqrt(3) / 4
  double wide_search_factor_       = 10.0;
  double min_cell_quality_         = 0.0;
  double max_cell_angle_           = M_PI;
  double base_vertex_factor_       = 2.00;

}; // FrontTriangulation
 

} // namespace TQAlgorithm
} // namespace TQMesh
