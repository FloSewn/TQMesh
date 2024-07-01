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
* from an already existing mesh 
*********************************************************************/
class FrontInitData
{
public:
  using IntVector    = std::vector<int>;
  using BoolVector   = std::vector<bool>;
  using EdgeVector   = std::vector<Edge*>;
  using MeshVector   = std::vector<Mesh*>;

  /*------------------------------------------------------------------
  | Constructor / Destructor
  ------------------------------------------------------------------*/
  FrontInitData(const Domain& domain, const MeshVector& meshes) 
  { collect_front_edges(domain, meshes); }

  FrontInitData(const Domain& domain) 
  { MeshVector dummy {};  collect_front_edges(domain, dummy); }

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

  } // FrontInitData::collect_front_edges()

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

      // If edges have been located, terminate the search
      // --> First come, first serve
      //if ( nbr_edges.size() > 0 )
      //  break;
    }

    return std::move(nbr_edges);

  } // FrontInitData::get_neighbor_mesh_edges()


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  std::vector<EdgeVector> edges_        {};
  std::vector<IntVector>  colors_       {};
  std::vector<BoolVector> is_twin_edge_ {};

  double edge_overlap_range_ { 1.5 };

}; // FrontInitData


} // namespace TQMesh
