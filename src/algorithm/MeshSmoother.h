/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "VecND.h"
#include "ProgressBar.h"

#include "Vertex.h"
#include "Mesh.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* This class is used to implement smoothing algorithms for a 
* generated mesh
*
* References:
* -----------
*   Torsion spring approach:
*   Zhou and Shimada, An angle-based approach to two-dimensional 
*   mesh smoothing, 
*   IMR 2000 (2000), 373-384
*
*********************************************************************/
class MeshSmoother
{

  using VertexConnectivity = std::vector<std::pair<Vertex*, std::vector<Vertex*>>>;

public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  MeshSmoother() = default;

  /*------------------------------------------------------------------
  | The interface to run smoothing algorithms on a given mesh
  ------------------------------------------------------------------*/
  void smooth(const Domain& domain, Mesh& mesh, int iterations,
              double a_fac=0.5, double eps=0.75, double decay=1.0)
  {
    LOG(INFO) << "Start with smoothing of mesh " << mesh.id();

    ProgressBar progress_bar {};

    unfix_vertices( mesh );

    VertexConnectivity v_conn = init_vertex_connectivity( mesh );

    // Apply smoothing
    for (int i = 0; i < iterations; ++i)
    {
      double progress = std::ceil(100.0 * (i+1.0) / iterations);
      progress_bar.update( static_cast<int>(progress) );
      progress_bar.show( LOG_PROPERTIES.get_ostream(INFO) );

      smooth_torsion(domain, mesh, v_conn, 2, a_fac, eps, decay);
      smooth_laplace(domain, mesh, v_conn, 4, eps, decay);

      eps *= decay;
    }

    LOG(INFO) << "";
    LOG(INFO) << "Smoothing is completed.";

  } // MeshSmoother::smooth()


private:
  /*------------------------------------------------------------------
  | This function checks the element qualities and in case of 
  | bad elements, it sets fixed vertices to be unfixed (e.g. vertices
  | in the quad layer), such that they can be shifted for a better
  | mesh quality
  ------------------------------------------------------------------*/
  void unfix_vertices(Mesh& mesh) 
  {
    Triangles& tris = mesh.triangles();
    Quads& quads = mesh.quads();

    for ( const auto& t_ptr : tris )
    {
      if ( t_ptr->max_angle() > crit_tri_max_angle_ )
      {
        t_ptr->v1().is_fixed( false );
        t_ptr->v2().is_fixed( false );
        t_ptr->v3().is_fixed( false );
      }
    }

    for ( const auto& q_ptr : quads )
    {
      if ( q_ptr->min_angle() < crit_quad_min_angle_ )
      {
        q_ptr->v1().is_fixed( false );
        q_ptr->v2().is_fixed( false );
        q_ptr->v3().is_fixed( false );
        q_ptr->v4().is_fixed( false );
      }
    }

  } // MeshSmoother::unfix_vertices()

  /*------------------------------------------------------------------
  | Apply the grid smoothing with a torsion spring approach 
  ------------------------------------------------------------------*/
  void smooth_torsion(const Domain& domain, Mesh& mesh, 
                      const VertexConnectivity& v_conn, int iterations, 
                      double a_fac=0.5, double eps=0.75, 
                      double decay=1.0)
  {
    Triangles& triangles  = mesh.triangles();
    Quads&     quads      = mesh.quads();
    EdgeList&  intr_edges = mesh.interior_edges();

    for (int iter = 0; iter < iterations; ++iter)
    {
      std::vector<Vec2d> xy_new ( v_conn.size() );

      for (std::size_t i_v = 0; i_v < v_conn.size(); ++i_v) 
      {
        Vertex* v  = v_conn[i_v].first;
        auto& nbrs = v_conn[i_v].second;

        // Fixed vertices keep their location
        if ( v->is_fixed() || v->on_boundary() )
        {
          xy_new[i_v] = v->xy();
          continue;
        }

        // Compute new vertex location
        int n_nbrs = static_cast<int>( nbrs.size() );

        Vec2d xy_m { 0.0, 0.0 };

        for ( int j = 0; j < n_nbrs; ++j ) 
        {
          int i = MOD(j-1, n_nbrs);
          int k = MOD(j+1, n_nbrs);

          const Vec2d vi = v->xy() - nbrs[i]->xy();
          const Vec2d vj = v->xy() - nbrs[j]->xy();
          const Vec2d vk = v->xy() - nbrs[k]->xy();

          const double a1 = angle( vj, vk );
          const double a2 = angle( vj, vi );

          const double b = a_fac * (a2 - a1);

          const double sin_b = sin( b );
          const double cos_b = cos( b );

          const double xj = nbrs[j]->xy().x;
          const double yj = nbrs[j]->xy().y;

          xy_m.x += xj + cos_b * vj.x - sin_b * vj.y;
          xy_m.y += yj + sin_b * vj.x + cos_b * vj.y;
        }

        xy_m /= static_cast<double>( n_nbrs );

        // Compute relaxed new vertex location
        const Vec2d xy_old = v->xy();
        const Vec2d d_xy   = xy_m - xy_old;
        const Vec2d xy_n   = xy_old + eps * d_xy;
        
        // Check validity of new coordinate
        if ( check_coordinate( xy_n, *v, domain, mesh ) )
          xy_new[i_v] = xy_n;
        else
          xy_new[i_v] = xy_old;
      }

      // Apply new vertex positions
      for ( std::size_t i = 0; i < xy_new.size(); ++i )
      {
        Vertex* v = v_conn[i].first;
        v->adjust_xy( xy_new[i] );
      }

      // Update mesh entities
      for ( auto& tri : triangles )
        tri->update_metrics();

      for ( auto& quad : quads )
        quad->update_metrics();

      for ( auto& edge : intr_edges )
        edge->update_metrics();

      eps *= -decay;

    }

  } // MeshSmoother::smooth_torsion()

  /*------------------------------------------------------------------
  | Apply the grid smoothing with a laplace approach 
  ------------------------------------------------------------------*/
  void smooth_laplace(const Domain& domain, Mesh& mesh, 
                      const VertexConnectivity& v_conn, int iterations,
                      double eps=0.75, double decay=1.0)
  {
    Triangles& triangles  = mesh.triangles();
    Quads&     quads      = mesh.quads();
    EdgeList&  intr_edges = mesh.interior_edges();

    for (int iter = 0; iter < iterations; ++iter)
    {
      std::vector<Vec2d> xy_new ( v_conn.size() );

      for (std::size_t i_v = 0; i_v < v_conn.size(); ++i_v) 
      {
        Vertex* v  = v_conn[i_v].first;
        auto& nbrs = v_conn[i_v].second;

        // Fixed vertices keep their location
        if ( v->is_fixed() || v->on_boundary() )
        {
          xy_new[i_v] = v->xy();
          continue;
        }

        // Compute new vertex location
        int n_nbrs = static_cast<int>( nbrs.size() );

        Vec2d xy_m { 0.0, 0.0 };

        for ( int j = 0; j < n_nbrs; ++j )
          xy_m += nbrs[j]->xy();

        xy_m /= static_cast<double>( n_nbrs );

        // Compute relaxed new vertex location
        const Vec2d xy_old = v->xy();
        const Vec2d d_xy   = xy_m - xy_old;
        const Vec2d xy_n   = xy_old + eps * d_xy;
        
        // Check validity of new coordinate
        if ( check_coordinate( xy_n, *v, domain, mesh ) )
          xy_new[i_v] = xy_n;
        else
          xy_new[i_v] = xy_old;
      }

      // Apply new vertex positions
      for ( std::size_t i = 0; i < xy_new.size(); ++i )
      {
        Vertex* v = v_conn[i].first;
        v->adjust_xy( xy_new[i] );
      }

      // Update mesh entities
      for ( auto& tri : triangles )
        tri->update_metrics();

      for ( auto& quad : quads )
        quad->update_metrics();

      for ( auto& edge : intr_edges )
        edge->update_metrics();


      eps *= -decay;

    }

  } // MeshSmoother::smooth_laplace()

  /*------------------------------------------------------------------
  | Sets up the vertex->vertex connectivity that is needed for the 
  | smoothing 
  ------------------------------------------------------------------*/
  VertexConnectivity init_vertex_connectivity(Mesh& mesh)
  {
    VertexConnectivity v_conn {};

    Vertices& vertices = mesh.vertices();

    // For each vertex, gather all vertices that are connected
    // to it via its adjacent facets
    for ( auto& v_ptr : vertices )
    {
      v_conn.push_back( { v_ptr.get() ,{} } );

      auto& v_nbrs = v_conn.back();
      auto& nbrs   = v_nbrs.second;

      // Gather neighbors
      for ( auto f : v_ptr->facets() )
      {
        for ( size_t i = 0; i < f->n_vertices(); ++i )
        {
          Vertex& v_cur = f->vertex(i);

          if (v_cur == *v_ptr)
            continue;

          nbrs.push_back( &v_cur );
        }
      }

      // Sort neighbors by angle
      const Vec2d xy = v_ptr->xy();

      std::sort( nbrs.begin(), nbrs.end(),
      [xy] ( Vertex* v1, Vertex* v2 )
      {
        const Vec2d dxy1 = v1->xy() - xy;
        const Vec2d dxy2 = v2->xy() - xy;
        const double a1 = std::atan2(dxy1.y, dxy1.x);
        const double a2 = std::atan2(dxy2.y, dxy2.x);

        return ( a1 < a2 );
      });
    }

    return std::move( v_conn );

  } // Smoothin::init_vertex_connectivity()

  /*------------------------------------------------------------------
  | Apply the grid smoothing with a laplace approach 
  ------------------------------------------------------------------*/
  bool check_coordinate(const Vec2d& xy_n, const Vertex& v,
                        const Domain& domain, const Mesh& mesh)
  {
    const Triangles& triangles  = mesh.triangles();
    const Quads&     quads      = mesh.quads();

    const double rho = domain.size_function( xy_n );
    const double range = 2.0 * rho;

    for ( const auto& tri : triangles.get_items(xy_n, range) )
    {
      const Vertex& v1 = tri->v1();
      const Vertex& v2 = tri->v2();
      const Vertex& v3 = tri->v3();

      if (  v == v1 || v == v2 || v == v3 )
        continue;

      if ( in_on_triangle(xy_n, v1.xy(), v2.xy(), v3.xy()) )
        return false;
    }

    for ( const auto& quad : quads.get_items(xy_n, range) )
    {
      const Vertex& v1 = quad->v1();
      const Vertex& v2 = quad->v2();
      const Vertex& v3 = quad->v3();
      const Vertex& v4 = quad->v3();

      if (  v == v1 || v == v2 || v == v3 || v == v4 )
        continue;

      if ( in_on_quad(xy_n, v1.xy(), v2.xy(), v3.xy(), v4.xy()) )
        return false;
    }

    return ( domain.is_inside( xy_n ) );

  } // MeshSmoother::check_coordinate_validity()

  /*------------------------------------------------------------------
  | Attributes 
  ------------------------------------------------------------------*/
  double crit_tri_min_angle_   {  10.0 * M_PI / 180.0 };
  double crit_tri_max_angle_   { 150.0 * M_PI / 180.0 };

  double crit_quad_min_angle_  {  10.0 * M_PI / 180.0 };
  double crit_quad_max_angle_  { 150.0 * M_PI / 180.0 };


}; // MeshSmoother

} // namespace TQAlgorithm
} // namespace TQMesh
