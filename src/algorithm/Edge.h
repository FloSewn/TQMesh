/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <iomanip>   
#include <algorithm>   

#include "utils.h"
#include "Geometry.h"
#include "VecND.h"
#include "Container.h"

#include "Vertex.h"
#include "Facet.h"
#include "NullFacet.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* Forward declarations
*********************************************************************/
class EdgeList;


/*********************************************************************
* A simple edge class that is stored in the Container
* The edge is defined by two vertices v1 and v2. It is directed from
* v1 to v2. The edge normal is aligned to the left of the edge 
* direction.
*
*            n
*            ^  
*            |
*     v1 ---------> v2
*
*
*********************************************************************/
class Edge : public ContainerEntry<Edge>
{
public:
  /*------------------------------------------------------------------
  | Constructor 
  ------------------------------------------------------------------*/
  Edge(Vertex& v1, Vertex& v2, EdgeList& edgelist, int m) 
  : ContainerEntry<Edge>( 0.5 * ( v1.xy()+v2.xy() ) )
  , v1_       {&v1}
  , v2_       {&v2}
  , edgelist_ {&edgelist} 
  , marker_   {m}
  {
    ASSERT((v1_ && v2_),
        "Failed to create edge structure due to given nullptr." );

    update_metrics(false);

    v1_->add_edge( *this );
    v2_->add_edge( *this );
  }

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  int marker() const { return marker_; }

  EdgeList& edgelist() { return *edgelist_; }
  const EdgeList& edgelist() const { return *edgelist_; }

  const Vertex& v1() const { return *v1_; };
  const Vertex& v2() const { return *v2_; };
  Vertex& v1() { return *v1_; };
  Vertex& v2() { return *v2_; };

  double length() const { return length_; }
  const Vec2d& normal() const { return norm_;}
  const Vec2d& tangent() const { return tang_;}

  const Facet* facet_l() const { return face_l_; }
  const Facet* facet_r() const { return face_r_; }
  Facet* facet_l() { return face_l_; }
  Facet* facet_r() { return face_r_; }

  Vertex* sub_vertex() const { return sub_vertex_; }
  Edge* twin_edge() const { return twin_edge_; }

  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  void facet_l(Facet* f) { face_l_ = f; }
  void facet_r(Facet* f) { face_r_ = f; }

  void sub_vertex(Vertex* v) { sub_vertex_ = v; }
  void twin_edge(Edge* e) { twin_edge_ = e; }

  /*------------------------------------------------------------------
  | Function returns, if edges is located on a boundary
  | or if it is in the interior of the domain
  ------------------------------------------------------------------*/
  bool on_boundary() const 
  { return ( marker_ != CONSTANTS.interior_edge_marker() ); }
  bool is_interior() const
  { return !on_boundary(); }

  /*------------------------------------------------------------------
  | Get the next edge, that is connected to the ending vertex of  
  | this edge. The next edge must also have the same marker
  | as this edge.
  | The first vertex of the next edge must be the last vertex of the 
  | current edge.
  | This function can be used to traverse a connected list of edges
  ------------------------------------------------------------------*/
  Edge* get_next_edge() const
  {
    for ( auto e : v2_->edges() )
    {
      if ( e == this )
        continue;

      if ( &(e->edgelist()) != edgelist_ )
        continue;

      if ( &(e->v1()) != v2_ )
        continue; 

      return e;
    }

    return nullptr;

  } // get_next_edge()

  /*------------------------------------------------------------------
  | Get the previous edge, that is connected to the starting vertex of  
  | this edge. The previous edge must also have the same marker
  | as this edge.
  | The last vertex of the found edge must be the first vertex of the
  | current edge.
  | This function can be used to traverse a connected list of edges
  ------------------------------------------------------------------*/
  Edge* get_prev_edge() const
  {
    for ( auto e : v1_->edges() )
    {
      if ( e == this )
        continue;

      if ( &(e->edgelist()) != edgelist_ )
        continue;

      if ( &(e->v2()) != v1_ )
        continue;

      return e;
    }

    return nullptr;

  } // get_prev_edge()

  /*------------------------------------------------------------------
  | Update the edge if its vertices changed
  ------------------------------------------------------------------*/
  void update_metrics(bool update_centroid=true) 
  {

    if ( update_centroid )
    {
      Vec2d xy_new = 0.5 * ( v1_->xy() + v2_->xy() );
      bool success = container_->update( *this, xy_new );
      ASSERT( success, "Edge::update_metrics(): "
          "Failed to update edge centroid.");
      (void) success;
    }

    const Vec2d d_xy = v2_->xy() - v1_->xy();

    length_ = d_xy.norm();
    tang_   = d_xy / length_;

    norm_.x = -tang_.y;
    norm_.y =  tang_.x;
  }

  /*------------------------------------------------------------------
  | Destructor function for container garbage collection
  ------------------------------------------------------------------*/
  void container_destructor() override
  { 
    if (v1_) v1_->remove_edge( *this );
    if (v2_) v2_->remove_edge( *this );

    v1_ = nullptr;
    v2_ = nullptr;
  }


private:
  /*------------------------------------------------------------------
  | Edge attributes 
  ------------------------------------------------------------------*/
  Vertex*             v1_       { nullptr };
  Vertex*             v2_       { nullptr };
  EdgeList*           edgelist_ { nullptr };
  int                 marker_   { -1 };

  // Edge properties
  double              length_      { 0.0 };
  Vec2d               tang_        { 0.0, 0.0 };
  Vec2d               norm_        { 0.0, 0.0 };

  // Pointer to adjacent facets
  Facet*              face_l_ { &NullFacet::get_instance() };
  Facet*              face_r_ { &NullFacet::get_instance() };

  // Sub vertex for quad refinement of the mesh
  Vertex*             sub_vertex_ {nullptr};

  // Twin edge of a neighbor mesh
  Edge*               twin_edge_ {nullptr};

}; // Edge 

/*********************************************************************
* Define general edge container declaration
*********************************************************************/
using Edges = Container<Edge>;

/*********************************************************************
* Edge ostream operator overload
*********************************************************************/
static std::ostream& operator<<(std::ostream& os, 
                                const Edge& e)
{ return os << e.v1() << " --> "  << e.v2(); }

/***********************************************************
* Edge equality operator 
***********************************************************/
static bool operator==(const Edge& e1, const Edge& e2)
{ return &e1 == &e2; }
static bool operator!=(const Edge& e1, const Edge& e2)
{ return !(e1 == e2); }

} // namespace TQAlgorithm
} // namespace TQMesh
