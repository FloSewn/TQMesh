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

// https://stackoverflow.com/questions/68321552/how-to-define-less-than-operator-or-stdless-struct-for-eigenvector3f
// https://en.cppreference.com/w/cpp/algorithm/lexicographical_compare

struct LexicographicLess 
{
  template<class T>
  bool operator()(T const& lhs, T const& rhs) const 
  {
    return std::lexicographical_compare(lhs.begin(), lhs.end(), 
                                        rhs.begin(), rhs.end());
  }
};

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
  using Vec2dMap      = std::map<Vec2d,std::size_t,LexicographicLess>;
  using VecVector     = std::vector<Vec2d>;

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
    for ( const auto& e : domain.fixed_edges() )
      add_edge_data( *e, false );

  } // FrontInitData::collect_front_edges()

  /*------------------------------------------------------------------
  |
  ------------------------------------------------------------------*/
  void collect_front_vertices(const Domain& domain,
                              const MeshVector& meshes)
  {


  } // FrontInitData::collect_front_vertices()

  

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
  | 
  ------------------------------------------------------------------*/
  void add_edge_data(const Edge& e, bool is_twin)
  {
    const std::size_t i1 = add_edge_vertex( e.v1() );
    const std::size_t i2 = add_edge_vertex( e.v2() );

    if ( is_twin )
      vertex_ids_.push_back( {i2, i1} );
    else
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
    // If vertex is already included, return its index
    if ( vertex_map_.count(v.xy()) )
      return vertex_map_.at( v.xy() );

    // Add new vertex and store its index to the vertex-map
    const auto i = vertices_.size();

    vertex_map_[v.xy()] = i;

    vertices_.push_back( const_cast<Vertex*>( &v ) );

    return i;
  }

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  VertexVector  vertices_     {};
  IDPairVector  vertex_ids_   {};
  Vec2dMap      vertex_map_   {};

  EdgeVector    edges_        {};
  IntVector     colors_       {};
  BoolVector    is_twin_edge_ {};
  

  double edge_overlap_range_ { 1.5 };

}; // FrontInitData

} // namespace TQMesh
