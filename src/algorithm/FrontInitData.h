/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <TQMeshConfig.h>
#include "STLHeaders.h"
#include "CppUtils.h"

#include "EdgeList.h"
#include "Vertex.h"
#include "Domain.h"
#include "Mesh.h"

namespace TQMesh {

using namespace CppUtils;

/*********************************************************************
* This class is used to initialize a new advancing front structure
* from mesh boundary edges of a given domain.
* If a bunch of additional meshes is provided, which might be 
* adjacent to the given domain, FrontInitData takes their boundary
* edges into account, in order to preserve mesh adjacency.
*********************************************************************/
class FrontInitData__OLD
{
public:
  using IntVector    = std::vector<int>;
  using BoolVector   = std::vector<bool>;
  using EdgeVector   = std::vector<Edge*>;
  using MeshVector   = std::vector<Mesh*>;

  /*------------------------------------------------------------------
  | Constructor / Destructor
  ------------------------------------------------------------------*/
  FrontInitData__OLD(const Domain& domain, const MeshVector& meshes) 
  { 
    collect_front_edges(domain, meshes); 
  }

  FrontInitData__OLD(const Domain& domain) 
  { 
    MeshVector dummy {};  
    collect_front_edges(domain, dummy); 
  }

  ~FrontInitData__OLD() {};

  /*------------------------------------------------------------------
  | Move operator
  ------------------------------------------------------------------*/
  FrontInitData__OLD(FrontInitData__OLD&& f)
  : edges_ { std::move(f.edges_) }
  , colors_ { std::move(f.colors_) }
  , is_twin_edge_ { std::move(f.is_twin_edge_) }
  {}
  
  FrontInitData__OLD& operator=(FrontInitData__OLD&& f)
  {
    edges_ = std::move(f.edges_);
    colors_ = std::move(f.colors_);
    is_twin_edge_ = std::move(f.is_twin_edge_);
    return *this;
  }

  /*------------------------------------------------------------------
  | Getter
  ------------------------------------------------------------------*/
  std::vector<EdgeVector>& edges() { return edges_; }
  const std::vector<EdgeVector>& edges() const { return edges_; }

  std::vector<BoolVector>& is_twin_edge() { return is_twin_edge_; }
  const std::vector<BoolVector>& is_twin_edge() const { return is_twin_edge_; }

  std::vector<IntVector>& colors() { return colors_; }
  const std::vector<IntVector>& colors() const { return colors_; }

  std::size_t size() const { return edges_.size(); }

private:
  /*------------------------------------------------------------------
  | Traverses all boundaries of a given domain and collects the 
  | respective boundary edges, as well as their color and orientation.
  | Additionally, the boundary edges of given meshes are also 
  | considered, if they overlap with the domain boundary. This 
  | is required to construct conformal grids.
  | Edges which result from adjacent meshes are oriented in the 
  | opposite direction. This can be used to identify these edges
  | from "normal" boundary edges (e.g. for the advancing front 
  | refinement).
  ------------------------------------------------------------------*/
  void collect_front_edges(const Domain& domain,
                           const MeshVector& meshes)
  {
    // Collect advancing front edges that stem from exterior & interior
    // domain boundary edges.
    // If an exterior domain edge is adjacent to a neighbouring mesh,
    // take that neighbor mesh edges to preserve mesh adjacency 
    for ( const auto& boundary : domain )
    {
      EdgeVector edges {};
      IntVector  colors {};
      BoolVector is_twin_edge {};

      for ( const auto& e : boundary->edges() )
      {
        EdgeVector nbr_edges = get_neighbor_mesh_edges(*e, meshes);

        // Add edges from adjacent mesh
        if ( nbr_edges.size() > 0 )
        {
          for ( Edge* nbr_e : nbr_edges )
          {
            edges.push_back( nbr_e );
            is_twin_edge.push_back( true );
            colors.push_back( nbr_e->color() );
          }
        }
        // Add edges from domain boundary
        else
        {
          edges.push_back( e.get() ) ;
          is_twin_edge.push_back( false );
          colors.push_back( e->color() );
        }
      }

      edges_.push_back( edges );
      is_twin_edge_.push_back( is_twin_edge );
      colors_.push_back( colors );
    }

    // Collect advancing front edges that stem from fixed interior
    // edges of the given domain
    // For every fixed edge, we a also create an additional ghost
    // edge that is oriented in the opposite direction
    //
    // EdgeVector edges {};
    // IntVector  colors {};
    // BoolVector is_twin_edge {};

    // for ( const auto& e : domain.fixed_edges() )
    // {
    //   edges.push_back( e.get() );
    //   is_twin_edge.push_back( false );
    //   colors.push_back( e->color() );
    // }

    // edges_.push_back( edges );
    // is_twin_edge_.push_back( is_twin_edge );
    // colors_.push_back( colors );

    // **** ----------------------------------------------------------
    // **** Feature implementation: Fixed interior edges 
    // **** ----------------------------------------------------------
    // - Here we could add additional fixed interor edges that have been
    //   defined by the user
    //
    // - These edges must be added twice - in forward and backward 
    //   direction
    //
    // - The problem with the former is, that this will cause these two 
    //   interior edges to be added in the mesh structure 
    //   -> FrontUpdate.advance_front() - FrontUpdate.h, line 126
    //   although it should only be one single edge
    //
    // - Possible approaches to overcome this:
    //   > Track these double-edges and include only one of them.
    //     However, this will cause major changes in 
    //     FrontUpdate.advance_front()
    //   > How to do this? Add further if-statement in 
    //     FrontUpdate.h, line 126, where it is checked if the base 
    //     is interior
    //   > We could add another function, that checks if the base is 
    //     doubly defined 
    //   > This would result in the definition of a new edge color,
    //     Additionally to INTERIOR_EDGE_COLOR, we could use 
    //     INTERIOR_TWIN_EDGE, or something
    //   > Thus, when adding a fixed interior edge here,
    //     we should already give it a special color that can later
    //     on be used to indicate it (maybe the first one is a default
    //     interior front edge and the second one is a special 
    //     INTERIOR_TWIN_EDGE
    //
    //  - It is important, that we also use some kind of closed edgelist 
    //    for these interior fixed edges, so that the upcoming methods 
    //    of Front for adding vertices and edges can be used
    //  - Maybe the user passes a list of interior edges and these 
    //    will be converted into some kind of an interior boundary 
    //  - However, this would require some more thoughts about the
    //    user interface for interior fixed edges, in order to define
    //    multiple line segements
    //    > Instead of an "add_fixed_edge()" function, we should provide
    //      something like "add_fixed_edge_semgents()", where multiple
    //      fixed edges are passed at once by a list of vertices
    //    > Similar to boundaries, the Domain then requires an additional
    //      list of "interior_edge_segments_"
    //
    //     

  } // FrontInitData__OLD::collect_front_edges()

  /*------------------------------------------------------------------
  | For a given domain boundary edge <e>, 
  | search all boundary edges of a 
  | neighboring mesh, that are contained within <e>.
  | This function is used upon the initialization of the advancing
  | front.
  ------------------------------------------------------------------*/
  EdgeVector get_neighbor_mesh_edges(const Edge& e, 
                                     const MeshVector& meshes) 
  const
  {
    EdgeVector   nbr_edges {};

    const Vec2d& c = e.xy();
    double       r = edge_overlap_range_ * e.length();

    const Vec2d& v = e.v1().xy();
    const Vec2d& w = e.v2().xy();

    for ( Mesh* m : meshes )
    {
      ASSERT( m, "Neighbor mesh data structure is invalid." );

      // Get edges in the vicinity of the given edge
      EdgeVector edges_in_vicinity = m->get_bdry_edges(c, r);

      // Get all edges, that are actually located on the segment
      // of the current domain boundary edge 
      for ( Edge* e_vicinity : edges_in_vicinity )
      {
        const Vec2d& p = e_vicinity->v1().xy();
        const Vec2d& q = e_vicinity->v2().xy();

        if ( in_on_segment(v,w,p) && in_on_segment(v,w,q) )
          nbr_edges.push_back( e_vicinity );
      }

      // In case edges have been found, sort them in 
      // ascending direction to the stating vertex
      std::sort(nbr_edges.begin(), nbr_edges.end(), 
      [v](Edge* e1, Edge* e2)
      {
        double delta_1 = (e1->xy() - v).norm_sqr();
        double delta_2 = (e2->xy() - v).norm_sqr();
        return (delta_1 < delta_2);
      });
    }

    return std::move(nbr_edges);

  } // FrontInitData__OLD::get_neighbor_mesh_edges()


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  std::vector<EdgeVector> edges_        {};
  std::vector<IntVector>  colors_       {};
  std::vector<BoolVector> is_twin_edge_ {};

  double edge_overlap_range_ { 1.5 };

}; // FrontInitData__OLD



/*********************************************************************
*
*********************************************************************/
class FrontInitData
{
public:
  using IntVector     = std::vector<int>;
  using BoolVector    = std::vector<bool>;
  using IDPairVector  = std::vector<std::pair<std::size_t,std::size_t>>;
  using VertexVector  = std::vector<Vertex*>;
  using EdgeVector    = std::vector<Edge*>;
  using MeshVector    = std::vector<Mesh*>;

  /*------------------------------------------------------------------
  | Constructor / Destructor
  ------------------------------------------------------------------*/
  FrontInitData(const Domain& domain, const MeshVector& meshes) 
  { 
    collect_front_edges(domain, meshes); 
  }

  FrontInitData(const Domain& domain) 
  { 
    MeshVector dummy {};  
    collect_front_edges(domain, dummy); 
  }

  ~FrontInitData() {};

  /*------------------------------------------------------------------
  | Move operator
  ------------------------------------------------------------------*/
  FrontInitData(FrontInitData&& f)
  : edges_ { std::move(f.edges_) }
  , colors_ { std::move(f.colors_) }
  , is_twin_edge_ { std::move(f.is_twin_edge_) }
  {}
  
  FrontInitData& operator=(FrontInitData&& f)
  {
    edges_ = std::move(f.edges_);
    colors_ = std::move(f.colors_);
    is_twin_edge_ = std::move(f.is_twin_edge_);
    return *this;
  }

  /*------------------------------------------------------------------
  | Getter
  ------------------------------------------------------------------*/
  EdgeVector& edges() { return edges_; }
  const EdgeVector& edges() const { return edges_; }

  VertexVector& vertices() { return vertices_; }
  const VertexVector& vertices() const { return vertices_; }

  IDPairVector& vertex_ids() { return vertex_ids_; }
  const IDPairVector& vertex_ids() const { return vertex_ids_; }

  BoolVector& is_twin_edge() { return is_twin_edge_; }
  const BoolVector& is_twin_edge() const { return is_twin_edge_; }

  IntVector& colors() { return colors_; }
  const IntVector& colors() const { return colors_; }

  std::size_t size() const { return edges_.size(); }

private:
  /*------------------------------------------------------------------
  | For a given domain boundary edge <e>, 
  | search all boundary edges of a 
  | neighboring mesh, that are contained within <e>.
  | This function is used upon the initialization of the advancing
  | front.
  ------------------------------------------------------------------*/
  EdgeVector get_neighbor_mesh_edges(const Edge& e, 
                                     const MeshVector& meshes) 
  const
  {
    EdgeVector   nbr_edges {};

    const Vec2d& c = e.xy();
    double       r = edge_overlap_range_ * e.length();

    const Vec2d& v = e.v1().xy();
    const Vec2d& w = e.v2().xy();

    for ( Mesh* m : meshes )
    {
      ASSERT( m, "Neighbor mesh data structure is invalid." );

      // Get edges in the vicinity of the given edge
      EdgeVector edges_in_vicinity = m->get_bdry_edges(c, r);

      // Get all edges, that are actually located on the segment
      // of the current domain boundary edge 
      for ( Edge* e_vicinity : edges_in_vicinity )
      {
        const Vec2d& p = e_vicinity->v1().xy();
        const Vec2d& q = e_vicinity->v2().xy();

        if ( in_on_segment(v,w,p) && in_on_segment(v,w,q) )
          nbr_edges.push_back( e_vicinity );
      }

      // In case edges have been found, sort them in 
      // ascending direction to the stating vertex
      std::sort(nbr_edges.begin(), nbr_edges.end(), 
      [v](Edge* e1, Edge* e2)
      {
        double delta_1 = (e1->xy() - v).norm_sqr();
        double delta_2 = (e2->xy() - v).norm_sqr();
        return (delta_1 < delta_2);
      });
    }

    return std::move(nbr_edges);

  } // FrontInitData::get_neighbor_mesh_edges()


  /*------------------------------------------------------------------
  | Traverses all boundaries of a given domain and collects the 
  | respective boundary edges, as well as their color and orientation.
  | Additionally, the boundary edges of given meshes are also 
  | considered, if they overlap with the domain boundary. This 
  | is required to construct conformal grids.
  | Edges which result from adjacent meshes are oriented in the 
  | opposite direction. This can be used to identify these edges
  | from "normal" boundary edges (e.g. for the advancing front 
  | refinement).
  ------------------------------------------------------------------*/
  void collect_front_edges(const Domain& domain,
                           const MeshVector& meshes)
  {
    // Collect advancing front edges that stem from exterior & interior
    // domain boundary edges.
    // If an exterior domain edge is adjacent to a neighbouring mesh,
    // take that neighbor mesh edges to preserve mesh adjacency 
    for ( const auto& boundary : domain )
    {
      for ( const auto& e : boundary->edges() )
      {
        EdgeVector nbr_edges = get_neighbor_mesh_edges(*e, meshes);

        // Add edges from adjacent mesh
        if ( nbr_edges.size() > 0 )
          for ( Edge* nbr_e : nbr_edges )
            add_edge_data( *nbr_e, true );
        // Add edges from domain boundary
        else
          add_edge_data( *e, false );
      }
    }

    // Collect advancing front edges that stem from fixed interior
    // edges of the given domain
    // For every fixed edge, we a also create an additional ghost
    // edge that is oriented in the opposite direction
    for ( const auto& e : domain.fixed_edges() )
      add_edge_data( *e, false );


  } // FrontInitData::collect_front_edges()

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  void add_edge_data(const Edge& e, bool is_twin)
  {
    const std::size_t i1 = add_edge_vertex( e.v1() );
    const std::size_t i2 = add_edge_vertex( e.v2() );

    vertex_ids_.push_back( {i1, i2} );
    edges_.push_back( const_cast<Edge*>( &e ) );
    is_twin_edge_.push_back( is_twin );
    colors_.push_back( e.color() );
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  std::size_t add_edge_vertex(const Vertex& v)
  {
    std::size_t vertex_location = 0;

    for ( Vertex* v_ptr : vertices_ )
    {
      if ( &v == v_ptr )
        return vertex_location;

      ++vertex_location;
    }

    vertices_.push_back( const_cast<Vertex*>( &v ) );

    return vertex_location;
  }

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  VertexVector  vertices_     {};
  IDPairVector  vertex_ids_   {};
  EdgeVector    edges_        {};
  IntVector     colors_       {};
  BoolVector    is_twin_edge_ {};

  double edge_overlap_range_ { 1.5 };

}; // FrontInitData

} // namespace TQMesh
