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

#include "utils.h"
#include "Vertex.h"
#include "Facet.h"

namespace TQMesh {

using namespace CppUtils;

/*********************************************************************
* Forward declarations
*********************************************************************/
class EdgeList;

/*********************************************************************
* Different properties for mesh edges 
* Edges can be classified as: 
* - Being located on the domain boundary ("on_boundary")
*   If they are not located on the boundary, they are classified as 
*   internal edges (which is also the default condition)
* - Being "ghost"-edges ("is_ghost") - these are special internal 
*   edges, that are required for the placement of interior edges 
*   by the user
*********************************************************************/
enum class EdgeProperty : uint8_t {
  no_property    = 0b00000000,
  on_boundary    = 0b00000001,
  is_fixed       = 0b00000010,
  is_ghost       = 0b00000100,
};

// Overloaded bitwise NOT operator
static inline EdgeProperty operator~(EdgeProperty prop) 
{ return static_cast<EdgeProperty>(~static_cast<uint8_t>(prop)); }

// Overloaded compound bitwise OR operator for adding properties
static inline EdgeProperty&
operator|=(EdgeProperty& lhs, EdgeProperty rhs) 
{ 
  lhs = static_cast<EdgeProperty>(
    static_cast<uint8_t>(lhs) | static_cast<uint8_t>(rhs)
  );
  return lhs;
}

// Overloaded compound bitwise AND operator for removing properties
static inline EdgeProperty& 
operator&=(EdgeProperty& lhs, EdgeProperty rhs) 
{
  lhs = static_cast<EdgeProperty>(
    static_cast<uint8_t>(lhs) & static_cast<uint8_t>(rhs)
  );
  return lhs;
}

// Overloaded bitwise AND operator for checking edge properties
static inline bool 
operator&(EdgeProperty lhs, EdgeProperty rhs) 
{ return (static_cast<uint8_t>(lhs) & static_cast<uint8_t>(rhs)) != 0; }

// Overloaded bitwise OR operator for checking edge properties
static inline bool 
operator|(EdgeProperty lhs, EdgeProperty rhs) 
{ return (static_cast<uint8_t>(lhs) | static_cast<uint8_t>(rhs)) != 0; }


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
  Edge(Vertex& v1, Vertex& v2, EdgeList& edgelist, int color) 
  : ContainerEntry<Edge>( 0.5 * ( v1.xy()+v2.xy() ) )
  , v1_       {&v1}
  , v2_       {&v2}
  , edgelist_ {&edgelist} 
  , color_    {color}
  {
    ASSERT((v1_ && v2_),
        "Failed to create edge structure due to given nullptr." );

    update_metrics(false);

    ASSERT( length_ > 0.0, 
        "Edge: Invalid edge definition - edge length is zero.");

    v1_->add_edge( *this );
    v2_->add_edge( *this );
  }

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  int color() const { return color_; }

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
  | Handle edge properties 
  ------------------------------------------------------------------*/
  const EdgeProperty& properties() const { return properties_; }
  void add_property(EdgeProperty p) { properties_ |= p; }
  void set_property(EdgeProperty p) { properties_ = p; }
  void remove_property(EdgeProperty p) { properties_ &= ~p; }
  bool has_property(EdgeProperty p) const { return (properties_ & p); }
  bool has_no_property() const 
  { return (properties_ == EdgeProperty::no_property); }

  bool on_boundary() const { return has_property(EdgeProperty::on_boundary); }
  bool is_interior() const { return !on_boundary(); }
  bool is_ghost() const { return has_property(EdgeProperty::is_ghost); }
  bool is_fixed() const { return has_property(EdgeProperty::is_fixed); }

  /*------------------------------------------------------------------
  | Get the next edge, that is connected to the ending vertex of  
  | this edge. The next edge must also have the same color
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
  | this edge. The previous edge must also have the same color
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
  int                 color_    { -1 };

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

  EdgeProperty        properties_ { EdgeProperty::no_property};

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

} // namespace TQMesh
