/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "TQMesh.h"

namespace TQMesh {

using namespace CppUtils;

class MixedSmoothingStrategy;

/*********************************************************************
* This class is used to implement smoothing algorithms for a 
* generated mesh
*
*********************************************************************/
class SmoothingStrategy
{
public:
  using VConn          = std::pair<Vertex*,std::vector<Vertex*>>;
  using Connectivities = std::vector<VConn>;
  using Vec2dVector    = std::vector<Vec2d>;

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
  | This is the general loop for smoothing strategies
  ------------------------------------------------------------------*/
  virtual Vec2d compute_displacement(const VConn& v_conn) const = 0;

  /*------------------------------------------------------------------
  | For each vertex in a quad layer, we collect its direction to
  | the nearest point at a boundary edge. This information is used
  | later on, in order to preserve the prescribed quad layer height 
  | of these vertices during the smoothing process
  ------------------------------------------------------------------*/
  void collect_dispalcement_directions() 
  {
    Vertices& vertices = mesh_->vertices();
    EdgeList& boundary_edges = mesh_->boundary_edges();

    // Clear vector that contains the direction to the nearest
    // boundary edge for all quad layer vertices
    bdry_direction_.clear();
    bdry_direction_.assign(vertices.size(), {0.0, 0.0} );

    if ( !quad_layer_smoothing_ )
      return;

    for ( auto& v_ptr : vertices )
    {
      // Skip all vertices except those located in quad layers
      if ( !v_ptr->has_property(VertexProperty::in_quad_layer) )
        continue;

      // Obtain nearest boundary edge
      const Vec2d& v_xy = v_ptr->xy();
      Edge* e_ptr = boundary_edges.get_nearest_edge( v_xy );

      if ( !e_ptr )
        continue;

      // Locate closest point to current vertex on the given edge
      const Vec2d& e1_xy = e_ptr->v1().xy(); 
      const Vec2d& e2_xy = e_ptr->v2().xy(); 

      const Vec2d delta_e = e2_xy - e1_xy;
      const Vec2d delta_v = v_xy - e1_xy;

      const double dot_p = dot(delta_v, delta_e) / delta_e.norm_sqr();
      const double t = std::clamp(dot_p, 0.0, 1.0);
      const Vec2d p_xy = e1_xy + t * delta_e;
      const Vec2d delta_p = p_xy - v_xy;

      bdry_direction_.push_back( delta_p / delta_p.norm() );
    }

  } // collect_dispalcement_directions()

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
  | Check if a given coordinate is valid for a vertex
  ------------------------------------------------------------------*/
  bool check_new_vertex_coordinate(const Vec2d& xy_n, const Vertex& v) const
  {
    if ( !domain_->is_inside( xy_n ) )
        return false;

    return true;

  } // Smoothing::check_new_vertex_coordinate()

  /*------------------------------------------------------------------
  | Check if a new coordinate for a given vertex is valid
  ------------------------------------------------------------------*/
  bool check_coordinate(const Vec2d& xy_n, const Vertex& v) const
  {
    const Triangles& triangles  = mesh_->triangles();
    const Quads&     quads      = mesh_->quads();

    const double range = 2.0 * (xy_n - v.xy()).norm();

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

  } // SmoothingStrategy::check_coordinate()

  /*------------------------------------------------------------------
  | Check if a new coordinate for a given vertex is valid
  ------------------------------------------------------------------*/
  bool new_vertex_position_is_valid(const Vertex& v, double range) const
  {
    if ( !domain_->is_inside( v ) )
      return false;

    // Check facet orientation
    for ( const auto& f_ptr : v.facets() )
    {
      if ( f_ptr->n_vertices() == 3 )
      {
        const Triangle& tri = *(static_cast<const Triangle*>(f_ptr));
        const Vec2d& v1 = tri.v1().xy();
        const Vec2d& v2 = tri.v2().xy();
        const Vec2d& v3 = tri.v3().xy();
        if ( orientation(v1, v2, v3) != Orientation::CCW )
          return false;
      }

      if ( f_ptr->n_vertices() == 4 )
      {
        const Quad& quad = *(static_cast<const Quad*>(f_ptr));
        const Vec2d& v1 = quad.v1().xy();
        const Vec2d& v2 = quad.v2().xy();
        const Vec2d& v3 = quad.v3().xy();
        const Vec2d& v4 = quad.v4().xy();
        if ( orientation(v1, v2, v3) != Orientation::CCW )
          return false;
        if ( orientation(v2, v3, v4) != Orientation::CCW )
          return false;
      }
    }

    // Check edge crossings between local facet edges
    double search_radius = 0.0;
    for ( const auto& f_ptr : v.facets() )
      search_radius = MAX(search_radius, f_ptr->max_edge_length());
    
    const EdgeList& interior_edges = mesh_->interior_edges();
    const EdgeList& boundary_edges = mesh_->boundary_edges();

    auto i_edges = interior_edges.get_edges(v.xy(), search_radius);
    auto b_edges = boundary_edges.get_edges(v.xy(), search_radius);

    for (std::size_t i = 0; i < i_edges.size(); ++i)
    {
      const Vec2d& v1 = i_edges[i]->v1().xy();
      const Vec2d& v2 = i_edges[i]->v2().xy();

      for (std::size_t j = 0; j < i_edges.size(); ++j)
      {
        if (i == j) 
          continue;

        const Vec2d& w1 = i_edges[j]->v1().xy();
        const Vec2d& w2 = i_edges[j]->v2().xy();

        if ( line_line_crossing( v1, v2, w1, w2 ) )
          return false;
      }

      for (std::size_t j = 0; j < b_edges.size(); ++j)
      {
        const Vec2d& w1 = b_edges[j]->v1().xy();
        const Vec2d& w2 = b_edges[j]->v2().xy();

        if ( line_line_crossing( v1, v2, w1, w2 ) )
          return false;
      }
    }

    return true;

    //auto Edges = mesh_->get_interior_edges(

    /*

    for ( const auto& f_ptr : v.facets() )
    {
      // Check correct facet orientation!!

      if ( f_ptr->n_vertices() == 3 )
      {
        const Triangle& tri = *(static_cast<const Triangle*>(f_ptr));
        if ( tri.intersects_domain(*domain_) ) 
          return false;
      }

      if ( f_ptr->n_vertices() == 4 )
      {
        const Quad& quad = *(static_cast<const Quad*>(f_ptr));
        if ( quad.intersects_domain(*domain_) ) 
          return false;
      }

    }

    if ( v.intersects_facet(mesh_->triangles(), 2.0*range) )
      return false;

    if ( v.intersects_facet(mesh_->quads(), 2.0*range) )
      return false;

    if ( v.intersects_mesh_edges(*mesh_, 2.0*range, edge_factor * rho) )
      return false;
    */

    return true;

  } // SmoothingStrategy::new_vertex_position_is_valid()

  /*------------------------------------------------------------------
  | This is the general loop for smoothing strategies
  ------------------------------------------------------------------*/
  void smoothing_iteration() const
  {
    for (std::size_t i_v = 0; i_v < v_conn_.size(); ++i_v) 
    {
      Vertex& v  = *(v_conn_[i_v].first);

      // Fixed vertices keep their location
      if ( v.has_property( VertexProperty::is_fixed )    || 
           v.has_property( VertexProperty::on_boundary )  )
        continue;

      // If quad layer smoothing is not enabled, 
      // fix quad layer vertices
      if ( !quad_layer_smoothing_ && 
           v.has_property( VertexProperty::in_quad_layer ) ) 
        continue;

      // Compute vertex displacement (handeled differently by 
      // each smoothing strategy)
      Vec2d delta = compute_displacement( v_conn_[i_v] );

      // For quad layer vertices, we substract the projection of the
      // displacement onto the directional vector towards the nearest
      // boundary, in order to preserve quad layer heights
      if ( v.has_property( VertexProperty::in_quad_layer ) )
      {
        const Vec2d& bdry_dir = bdry_direction_[i_v];
        delta -= bdry_dir * dot(delta, bdry_dir);
      }

      // Check the validity of the new coordinate and possibly change 
      // it back, if it violates the mesh structure
      // --> xy_old must be a copy of v.xy(), in order to store its
      //     state!
      const Vec2d xy_old = v.xy();
      const Vec2d xy_n   = xy_old + eps_ * delta;

      if ( check_new_vertex_coordinate( xy_n, v ) )
      {
        MeshCleanup::set_vertex_coordinates(v, xy_n);

        if ( !new_vertex_position_is_valid(v, (xy_old-xy_n).norm()) )
          MeshCleanup::set_vertex_coordinates(v, xy_old);
      }
    }

  } // SmoothingStrategy::smoothing_iteration()

  /*------------------------------------------------------------------
  | This is the general loop for smoothing strategies
  ------------------------------------------------------------------*/
  void smoothing_loop(int iterations) 
  {
    for (int iter = 0; iter < iterations; ++iter)
    {
      smoothing_iteration();
      eps_ *= -decay_;
    }

  } // SmoothingStrategy::smoothing_loop() 

  /*------------------------------------------------------------------
  | Attributes 
  ------------------------------------------------------------------*/
  Mesh*              mesh_;
  const Domain*      domain_;
  Connectivities     v_conn_ {};
  Vec2dVector        bdry_direction_ {};

  double             eps_                  = 0.75;
  double             decay_                = 1.00;
  bool               quad_layer_smoothing_ = false;
  double             angle_factor_         = 0.5;

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
  bool quad_layer_smoothing() const { return quad_layer_smoothing_; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  LaplaceSmoothingStrategy& epsilon(double e) { eps_ = e; return *this;} 
  LaplaceSmoothingStrategy& decay(double d) { decay_ = d; return *this;} 
  LaplaceSmoothingStrategy& quad_layer_smoothing(bool b) 
  { quad_layer_smoothing_ = b; return *this;} 

  /*------------------------------------------------------------------
  | Apply mesh smoothing
  ------------------------------------------------------------------*/
  bool smooth(int iterations) override
  {
    MeshCleanup::setup_facet_connectivity(*mesh_);

    init_vertex_connectivity();

    collect_dispalcement_directions();

    smoothing_loop(iterations);
    
    return true;

  } // LaplaceSmoothingStrategy::smooth()


protected:
  /*------------------------------------------------------------------
  | This is the general loop for smoothing strategies
  ------------------------------------------------------------------*/
  Vec2d compute_displacement(const VConn& v_conn) const override
  {
    Vertex* v  = v_conn.first;
    auto& nbrs = v_conn.second;
    int n_nbrs = static_cast<int>( nbrs.size() );

    Vec2d xy_m { 0.0, 0.0 };

    for ( int j = 0; j < n_nbrs; ++j )
      xy_m += nbrs[j]->xy();

    xy_m /= static_cast<double>( n_nbrs );

    return xy_m - v->xy();

  } // LaplaceSmoothingStrategy::compute_displacement()

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
  bool quad_layer_smoothing() const { return quad_layer_smoothing_; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  TorsionSmoothingStrategy& epsilon(double e) { eps_ = e; return *this;} 
  TorsionSmoothingStrategy& decay(double d) { decay_ = d; return *this;} 
  TorsionSmoothingStrategy& angle_factor(double a) 
  { angle_factor_ = a; return *this; } 
  TorsionSmoothingStrategy& quad_layer_smoothing(bool b) 
  { quad_layer_smoothing_ = b; return *this;} 

  /*------------------------------------------------------------------
  | Apply mesh smoothing
  ------------------------------------------------------------------*/
  bool smooth(int iterations) override
  {
    MeshCleanup::setup_facet_connectivity(*mesh_);

    init_vertex_connectivity();

    collect_dispalcement_directions();

    smoothing_loop(iterations);
    
    return true;

  } // TorsionSmoothingStrategy::smooth()

protected:
  /*------------------------------------------------------------------
  | This is the general loop for smoothing strategies
  ------------------------------------------------------------------*/
  Vec2d compute_displacement(const VConn& v_conn) const override
  {
    Vertex* v  = v_conn.first;
    auto& nbrs = v_conn.second;
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

    return xy_m - v->xy();

  } // TorsionSmoothingStrategy::compute_displacement()

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
  bool quad_layer_smoothing() const { return quad_layer_smoothing_; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  MixedSmoothingStrategy& epsilon(double e) { eps_ = e; return *this;} 
  MixedSmoothingStrategy& decay(double d) { decay_ = d; return *this;} 
  MixedSmoothingStrategy& angle_factor(double a) 
  { angle_factor_ = a; return *this; } 
  MixedSmoothingStrategy& quad_layer_smoothing(bool b) 
  { quad_layer_smoothing_ = b; return *this;} 

  /*------------------------------------------------------------------
  | Apply mesh smoothing
  ------------------------------------------------------------------*/
  bool smooth(int iterations) override
  {
    MeshCleanup::setup_facet_connectivity(*mesh_);

    LaplaceSmoothingStrategy laplace { *mesh_, *domain_ };
    laplace.epsilon( eps_ );
    laplace.decay( decay_ );
    laplace.quad_layer_smoothing( quad_layer_smoothing_ );
    laplace.init_vertex_connectivity();
    laplace.collect_dispalcement_directions();

    TorsionSmoothingStrategy torsion { *mesh_, *domain_ };
    torsion.epsilon( eps_ );
    torsion.decay( decay_ );
    torsion.angle_factor( angle_factor_ );
    torsion.quad_layer_smoothing( quad_layer_smoothing_ );
    torsion.init_vertex_connectivity();
    torsion.collect_dispalcement_directions();

    for (int i = 0; i < iterations; ++i)
    {
      torsion.smoothing_iteration();
      laplace.smoothing_iteration();

      eps_ *= decay_;

      laplace.epsilon( eps_ );
      torsion.epsilon( eps_ );
    }
    
    return true;

  } // MixedSmoothingStrategy::smooth()

protected:
  /*------------------------------------------------------------------
  | Dummy function
  ------------------------------------------------------------------*/
  Vec2d compute_displacement(const VConn& v_conn) const override
  { return {0.0, 0.0}; }

}; // MixedSmoothingStrategy

} // namespace TQMesh
