/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "Vec2.h"
#include "utils.h"
#include "ProgressBar.h"

#include "Vertex.h"
#include "Triangle.h"
#include "Quad.h"
#include "Mesh.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace TQUtils;

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
class Smoother
{

  using VertexConnectivity = std::vector<std::pair<Vertex*, std::vector<Vertex*>>>;

public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Smoother() = default;

  /*------------------------------------------------------------------
  | The interface to run smoothing algorithms on a given mesh
  ------------------------------------------------------------------*/
  void smooth(const Domain& domain, Mesh& mesh, int iterations,
              double a_fac=0.5, double eps=0.75, double decay=1.0)
  {
    MSG("START MESH SMOOTHING");

    ProgressBar progress_bar {};

    unfix_vertices( mesh );

    VertexConnectivity v_conn = init_vertex_connectivity( mesh );

    // Apply smoothing
    for (int i = 0; i < iterations; ++i)
    {
      double progress = std::ceil(100.0 * (i+1.0) / iterations);
      progress_bar.update( static_cast<int>(progress) );
      progress_bar.show( std::clog );

      smooth_torsion(domain, mesh, v_conn, 2, a_fac, eps, decay);
      smooth_laplace(domain, mesh, v_conn, 4, eps, decay);

      eps *= decay;
    }

    MSG("");
    MSG("DONE!");

  } // Smoother::smooth()


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

  } // Smoother::unfix_vertices()

  /*------------------------------------------------------------------
  | Apply the grid smoothing with a torsion spring approach 
  ------------------------------------------------------------------*/
  void smooth_torsion(const Domain& domain, Mesh& mesh, 
                      const VertexConnectivity& v_conn, int iterations, 
                      double a_fac=0.5, double eps=0.75, 
                      double decay=1.0)
  {
    Triangles& triangles = mesh.triangles();
    Quads& quads = mesh.quads();

    for (int iter = 0; iter < iterations; ++iter)
    {
      for ( auto v_nbrs : v_conn )
      {
        Vertex* v  = v_nbrs.first;
        auto& nbrs = v_nbrs.second;

        if ( v->is_fixed() || v->on_boundary() )
          continue;

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

        const Vec2d xy_old = v->xy();
        const Vec2d d_xy   = xy_m - xy_old;
        const Vec2d xy_n   = xy_old + eps * d_xy;

        v->xy( xy_n );

        const double rho = domain.size_function( xy_n );
        const double range = 2.0 * rho;

        if (  v->intersects_facet(triangles, range)
           || v->intersects_facet(quads, range) 
           || !domain.is_inside( *v ) )
        {
          v->xy( xy_old );
        }

      }

      eps *= -decay;

    }

  } // Smoother::smooth_torsion()

  /*------------------------------------------------------------------
  | Apply the grid smoothing with a laplace approach 
  ------------------------------------------------------------------*/
  void smooth_laplace(const Domain& domain, Mesh& mesh, 
                      const VertexConnectivity& v_conn, int iterations,
                      double eps=0.75, double decay=1.0)
  {
    Triangles& triangles = mesh.triangles();
    Quads& quads = mesh.quads();

    for (int iter = 0; iter < iterations; ++iter)
    {
      for ( auto v_nbrs : v_conn )
      {
        Vertex* v  = v_nbrs.first;
        auto& nbrs = v_nbrs.second;

        if ( v->is_fixed() || v->on_boundary() )
          continue;

        int n_nbrs = static_cast<int>( nbrs.size() );

        Vec2d xy_m { 0.0, 0.0 };

        for ( int j = 0; j < n_nbrs; ++j )
          xy_m += nbrs[j]->xy();

        xy_m /= static_cast<double>( n_nbrs );

        const Vec2d xy_old = v->xy();

        const Vec2d d_xy   = xy_m - xy_old;
        const Vec2d xy_n   = xy_old + eps * d_xy;

        v->xy( xy_n );

        const double rho = domain.size_function( xy_n );
        const double range = 2.0 * rho;

        if (  v->intersects_facet(triangles, range)
           || v->intersects_facet(quads, range) 
           || !domain.is_inside( *v ) )
          v->xy( xy_old );
      }

      eps *= -decay;

    }

  } // Smoother::smooth_laplace()

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
  | Attributes 
  ------------------------------------------------------------------*/
  double crit_tri_min_angle_   {  10.0 * M_PI / 180.0 };
  double crit_tri_max_angle_   { 150.0 * M_PI / 180.0 };

  double crit_quad_min_angle_  {  10.0 * M_PI / 180.0 };
  double crit_quad_max_angle_  { 150.0 * M_PI / 180.0 };


}; // Smoother

} // namespace TQAlgorithm
} // namespace TQMesh
