/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "VecND.h"

#include "Vertex.h"
#include "Mesh.h"
#include "SmoothingStrategy.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* This class implements a simple Laplacian smoothing  
*********************************************************************/
class LaplaceSmoothingStrategy : public SmoothingStrategy
{

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
  void epsilon(double e) { eps_ = e; } 
  void decay(double d) { decay_ = d; } 

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


private:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  double eps_   = 0.75;
  double decay_ = 1.00;

}; // LaplaceSmoothingStrategy

} // namespace TQAlgorithm
} // namespace TQMesh
