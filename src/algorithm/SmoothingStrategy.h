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

class MixedSmoothingStrategy;

/*********************************************************************
* This class is used to implement smoothing algorithms for a 
* generated mesh
*
*********************************************************************/
class SmoothingStrategy
{
  using VertexConnectivity = std::vector<std::pair<Vertex*, 
                                         std::vector<Vertex*>>>;

public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  SmoothingStrategy(Mesh& mesh, const Domain& domain)
  : mesh_   { &mesh }
  , domain_ { &domain }
  {}
  
  virtual ~SmoothingStrategy() {}

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  Mesh& mesh() { return *mesh_; }

  /*------------------------------------------------------------------
  | The interface to run smoothing algorithms on a given mesh
  ------------------------------------------------------------------*/
  virtual bool smooth(int iterations) = 0;

protected:
  /*------------------------------------------------------------------
  | This function checks the element qualities and in case of 
  | bad elements, it sets fixed vertices to be unfixed (e.g. vertices
  | in the quad layer), such that they can be shifted for a better
  | mesh quality
  ------------------------------------------------------------------*/
  void unfix_vertices() 
  {
    Triangles& tris  = mesh_->triangles();
    Quads&     quads = mesh_->quads();

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

  } // SmoothingStrategy::unfix_vertices()

  /*------------------------------------------------------------------
  | Sets up the vertex->vertex connectivity that is needed for the 
  | smoothing 
  ------------------------------------------------------------------*/
  void init_vertex_connectivity()
  {
    // Clear the vertex connectivity container
    v_conn_.clear();

    Vertices& vertices = mesh_->vertices();

    // For each vertex, gather all vertices that are connected
    // to it via its adjacent facets
    for ( auto& v_ptr : vertices )
    {
      v_conn_.push_back( { v_ptr.get() ,{} } );

      auto& v_nbrs = v_conn_.back();
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

  } // Smoothin::init_vertex_connectivity()

  /*------------------------------------------------------------------
  | Apply the grid smoothing with a laplace approach 
  ------------------------------------------------------------------*/
  bool check_coordinate(const Vec2d& xy_n, const Vertex& v)
  {
    const Triangles& triangles  = mesh_->triangles();
    const Quads&     quads      = mesh_->quads();

    const double rho = domain_->size_function( xy_n );
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

    return ( domain_->is_inside( xy_n ) );

  } // SmoothingStrategy::check_coordinate_validity()

  /*------------------------------------------------------------------
  | Attributes 
  ------------------------------------------------------------------*/
  Mesh*              mesh_;
  const Domain*      domain_;
  VertexConnectivity v_conn_ {};

  double crit_tri_min_angle_   {  10.0 * M_PI / 180.0 };
  double crit_tri_max_angle_   { 150.0 * M_PI / 180.0 };

  double crit_quad_min_angle_  {  10.0 * M_PI / 180.0 };
  double crit_quad_max_angle_  { 150.0 * M_PI / 180.0 };


}; // SmoothingStrategy


/*********************************************************************
* This class implements a simple Laplacian smoothing  
*********************************************************************/
class LaplaceSmoothingStrategy : public SmoothingStrategy
{
  friend class MixedSmoothingStrategy;

public:

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  LaplaceSmoothingStrategy(Mesh& mesh, const Domain& domain) 
  : SmoothingStrategy(mesh, domain) {}

  ~LaplaceSmoothingStrategy() {}

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  double epsilon() const { return eps_; }
  double decay() const { return decay_; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  LaplaceSmoothingStrategy& epsilon(double e) { eps_ = e; return *this;} 
  LaplaceSmoothingStrategy& decay(double d) { decay_ = d; return *this;} 

  /*------------------------------------------------------------------
  | Apply mesh smoothing
  ------------------------------------------------------------------*/
  bool smooth(int iterations) override
  {
    unfix_vertices();

    init_vertex_connectivity();

    smooth_laplace(iterations);
    
    return true;

  } // LaplaceSmoothingStrategy::smooth()


private:
  /*------------------------------------------------------------------
  | Apply the grid smoothing with a laplace approach 
  ------------------------------------------------------------------*/
  void smooth_laplace(int iterations)
  {
    Triangles& triangles  = mesh_->triangles();
    Quads&     quads      = mesh_->quads();
    EdgeList&  intr_edges = mesh_->interior_edges();

    for (int iter = 0; iter < iterations; ++iter)
    {
      std::vector<Vec2d> xy_new ( v_conn_.size() );

      for (std::size_t i_v = 0; i_v < v_conn_.size(); ++i_v) 
      {
        Vertex* v  = v_conn_[i_v].first;
        auto& nbrs = v_conn_[i_v].second;

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
        const Vec2d xy_n   = xy_old + eps_ * d_xy;
        
        // Check validity of new coordinate
        if ( check_coordinate( xy_n, *v ) )
          xy_new[i_v] = xy_n;
        else
          xy_new[i_v] = xy_old;
      }

      // Apply new vertex positions
      for ( std::size_t i = 0; i < xy_new.size(); ++i )
      {
        Vertex* v = v_conn_[i].first;
        v->adjust_xy( xy_new[i] );
      }

      // Update mesh entities
      for ( auto& tri : triangles )
        tri->update_metrics();

      for ( auto& quad : quads )
        quad->update_metrics();

      for ( auto& edge : intr_edges )
        edge->update_metrics();

      eps_ *= -decay_;
    }

  } // SmoothingStrategy::smooth_laplace()


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  double eps_   = 0.75;
  double decay_ = 1.00;

}; // LaplaceSmoothingStrategy


/*********************************************************************
* This class implement the torsion spring smoothing algorithm
*
* References:
* -----------
*   Torsion spring approach:
*   Zhou and Shimada, An angle-based approach to two-dimensional 
*   mesh smoothing, 
*   IMR 2000 (2000), 373-384
*
*********************************************************************/
class TorsionSmoothingStrategy : public SmoothingStrategy
{
  friend class MixedSmoothingStrategy;

public:

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  TorsionSmoothingStrategy(Mesh& mesh, const Domain& domain) 
  : SmoothingStrategy(mesh, domain) {}

  ~TorsionSmoothingStrategy() {}

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  double epsilon() const { return eps_; }
  double decay() const { return decay_; }
  double angle_factor() const { return angle_factor_; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  TorsionSmoothingStrategy& epsilon(double e) { eps_ = e; return *this;} 
  TorsionSmoothingStrategy& decay(double d) { decay_ = d; return *this;} 
  TorsionSmoothingStrategy& angle_factor(double a) 
  { angle_factor_ = a; return *this; } 

  /*------------------------------------------------------------------
  | Apply mesh smoothing
  ------------------------------------------------------------------*/
  bool smooth(int iterations) override
  {
    unfix_vertices();

    init_vertex_connectivity();

    smooth_torsion(iterations);
    
    return true;

  } // TorsionSmoothingStrategy::smooth()

private:
  /*------------------------------------------------------------------
  | Apply the grid smoothing with a torsion spring approach 
  ------------------------------------------------------------------*/
  void smooth_torsion(int iterations)
  {
    Triangles& triangles  = mesh_->triangles();
    Quads&     quads      = mesh_->quads();
    EdgeList&  intr_edges = mesh_->interior_edges();

    for (int iter = 0; iter < iterations; ++iter)
    {
      std::vector<Vec2d> xy_new ( v_conn_.size() );

      for (std::size_t i_v = 0; i_v < v_conn_.size(); ++i_v) 
      {
        Vertex* v  = v_conn_[i_v].first;
        auto& nbrs = v_conn_[i_v].second;

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

          const double b = angle_factor_ * (a2 - a1);

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
        const Vec2d xy_n   = xy_old + eps_ * d_xy;
        
        // Check validity of new coordinate
        if ( check_coordinate( xy_n, *v ) )
          xy_new[i_v] = xy_n;
        else
          xy_new[i_v] = xy_old;
      }

      // Apply new vertex positions
      for ( std::size_t i = 0; i < xy_new.size(); ++i )
      {
        Vertex* v = v_conn_[i].first;
        v->adjust_xy( xy_new[i] );
      }

      // Update mesh entities
      for ( auto& tri : triangles )
        tri->update_metrics();

      for ( auto& quad : quads )
        quad->update_metrics();

      for ( auto& edge : intr_edges )
        edge->update_metrics();

      eps_ *= -decay_;
    }
  } // SmoothingStrategy::smooth_torsion() */

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  double eps_          = 0.75;
  double decay_        = 1.00;
  double angle_factor_ = 0.5;

}; // TorsionSmoothingStrategy


/*********************************************************************
* This class implements a mixed smoothing approach that utilized
* both laplacian and torsion smoothing 
*********************************************************************/
class MixedSmoothingStrategy : public SmoothingStrategy
{
public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  MixedSmoothingStrategy(Mesh& mesh, const Domain& domain) 
  : SmoothingStrategy(mesh, domain) {}

  ~MixedSmoothingStrategy() {}

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  double epsilon() const { return eps_; }
  double decay() const { return decay_; }
  double angle_factor() const { return angle_factor_; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  MixedSmoothingStrategy& epsilon(double e) { eps_ = e; return *this;} 
  MixedSmoothingStrategy& decay(double d) { decay_ = d; return *this;} 
  MixedSmoothingStrategy& angle_factor(double a) 
  { angle_factor_ = a; return *this; } 

  /*------------------------------------------------------------------
  | Apply mesh smoothing
  ------------------------------------------------------------------*/
  bool smooth(int iterations) override
  {
    LaplaceSmoothingStrategy laplace { *mesh_, *domain_ };
    laplace.epsilon( eps_ );
    laplace.decay( decay_ );
    laplace.unfix_vertices();
    laplace.init_vertex_connectivity();

    TorsionSmoothingStrategy torsion { *mesh_, *domain_ };
    torsion.epsilon( eps_ );
    torsion.decay( decay_ );
    torsion.angle_factor( angle_factor_ );
    torsion.unfix_vertices();
    torsion.init_vertex_connectivity();

    for (int i = 0; i < iterations; ++i)
    {
      torsion.smooth_torsion(2);
      laplace.smooth_laplace(2);

      eps_ *= decay_;

      laplace.epsilon( eps_ );
      torsion.epsilon( eps_ );
    }
    
    return true;

  } // MixedSmoothingStrategy::smooth()

private:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  double eps_          = 0.75;
  double decay_        = 1.00;
  double angle_factor_ = 0.5;

}; // MixedSmoothingStrategy

} // namespace TQAlgorithm
} // namespace TQMesh
