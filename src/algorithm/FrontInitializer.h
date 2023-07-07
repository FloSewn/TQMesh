/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <vector>

#include "VecND.h"
#include "Geometry.h"

#include "utils.h"
#include "Vertex.h"
#include "Edge.h"
#include "EdgeList.h"
#include "Domain.h"
#include "Boundary.h"
#include "Mesh.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* The advancing front - defined by a list of edges
* > Must be defined counter-clockwise
*********************************************************************/
class FrontInitializer
{
public:
  using IntVector    = std::vector<int>;
  using BoolVector   = std::vector<bool>;
  using EdgeVector   = std::vector<Edge*>;
  using MeshVector   = std::vector<Mesh*>;

  /*------------------------------------------------------------------
  | Constructor / Destructor
  ------------------------------------------------------------------*/
  FrontInitializer(const Domain& domain, const MeshVector& meshes) 
  { collect_front_edges(domain, meshes); }

  FrontInitializer(const Domain& domain) 
  { MeshVector dummy {};  collect_front_edges(domain, dummy); }

  ~FrontInitializer() {};

  /*------------------------------------------------------------------
  | Move operator
  ------------------------------------------------------------------*/
  FrontInitializer(FrontInitializer&& f)
  : edges_ { std::move(f.edges_) }
  , markers_ { std::move(f.markers_) }
  , is_twin_edge_ { std::move(f.is_twin_edge_) }
  {}
  
  FrontInitializer& operator=(FrontInitializer&& f)
  {
    edges_ = std::move(f.edges_);
    markers_ = std::move(f.markers_);
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

  std::vector<IntVector>& markers() { return markers_; }
  const std::vector<IntVector>& markers() const { return markers_; }

private:
  /*------------------------------------------------------------------
  | Traverses all boundaries of a given domain and collects the 
  | respective boundary edges, as well as their marker and orientation.
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
    for ( const auto& boundary : domain )
    {
      EdgeVector edges {};
      IntVector  markers {};
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
            markers.push_back( nbr_e->marker() );
          }
        }
        // Add edges from domain boundary
        else
        {
          edges.push_back( e.get() ) ;
          is_twin_edge.push_back( false );
          markers.push_back( e->marker() );
        }
      }

      edges_.push_back( edges );
      is_twin_edge_.push_back( is_twin_edge );
      markers_.push_back( markers );
    }

  } // FrontInitializer::collect_front_edges()

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
    double       r = CONSTANTS.edge_search_factor() * e.length();

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

      // If edges have been located, terminate the search
      // --> First come, first serve
      //if ( nbr_edges.size() > 0 )
      //  break;
    }

    return std::move(nbr_edges);

  } // FrontInitializer::get_neighbor_mesh_edges()


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  std::vector<EdgeVector> edges_        {};
  std::vector<IntVector>  markers_      {};
  std::vector<BoolVector> is_twin_edge_ {};

}; // FrontInitializer


} // namespace TQAlgorithm
} // namespace TQMesh
