/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <list>
#include <iomanip>   
#include <algorithm>   

#include "utils.h"
#include "geometry.h"
#include "Vec2.h"
#include "Container.h"

#include "Vertex.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace TQUtils;

/*********************************************************************
* Forward declarations
*********************************************************************/
class Facet;
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
class Edge
{
public:

  friend Container<Edge>;
  using List = typename Container<Edge>::List;
  using ContainerIterator = typename Container<Edge>::List::iterator;
  using EdgeArray = std::array<Edge*,2>;

  /*------------------------------------------------------------------
  | Constructor 
  ------------------------------------------------------------------*/
  Edge(Vertex& v1, Vertex& v2, EdgeList& edgelist, int m) 
  : v1_       {&v1}
  , v2_       {&v2}
  , edgelist_ {&edgelist} 
  , marker_   {m}
  {
    ASSERT((v1_ && v2_),
        "Failed to create edge structure due to given nullptr." );
      
    const Vec2d d_xy = v2_->xy() - v1_->xy();

    xy_     = 0.5 * ( v1_->xy() + v2_->xy() );
    length_ = d_xy.length();
    tang_   = d_xy / length_;

    norm_.x = -tang_.y;
    norm_.y =  tang_.x;

    v1_->add_edge( *this );
    v2_->add_edge( *this );
  }

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  const Vec2d& xy() const { return xy_; }
  int marker() const { return marker_; }

  EdgeList& edgelist() { return *edgelist_; }
  const EdgeList& edgelist() const { return *edgelist_; }

  const ContainerIterator& pos() const { return pos_; }

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

  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  void facet_l(Facet* f) { face_l_ = f; }
  void facet_r(Facet* f) { face_r_ = f; }

  void sub_vertex(Vertex* v) { sub_vertex_ = v; }

  /*------------------------------------------------------------------
  | Function returns, if edges is located on a boundary
  | or if it is in the interior of the domain
  ------------------------------------------------------------------*/
  bool on_boundary() const 
  { return ( marker_ != TQ_INTR_EDGE_MARKER ); }
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
  | Mandatory container functions 
  ------------------------------------------------------------------*/
  bool in_container() const { return in_container_; }

private:
  /*------------------------------------------------------------------
  | Mandatory container functions 
  ------------------------------------------------------------------*/
  void container_destructor() 
  { 
    if (v1_) v1_->remove_edge( *this );
    if (v2_) v2_->remove_edge( *this );

    v1_ = nullptr;
    v2_ = nullptr;
  }

  /*------------------------------------------------------------------
  | Edge attributes 
  ------------------------------------------------------------------*/
  Vertex*             v1_       { nullptr };
  Vertex*             v2_       { nullptr };
  EdgeList*           edgelist_ { nullptr };
  int                 marker_   { -1 };

  Vec2d               xy_          { 0.0, 0.0 };
  double              length_      { 0.0 };
  Vec2d               tang_        { 0.0, 0.0 };
  Vec2d               norm_        { 0.0, 0.0 };

  Facet*              face_l_ {nullptr};
  Facet*              face_r_ {nullptr};

  Vertex*             sub_vertex_    {nullptr};

  // Mandatory container attributes
  ContainerIterator   pos_;
  bool                in_container_;

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
