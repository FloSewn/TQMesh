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

#include "Geometry.h"
#include "Vec2.h"
#include "Container.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* Forward declaration
*********************************************************************/
class Edge;
class Facet;

/*********************************************************************
* This class describes a Vertex in a 2 dimensional domain
*********************************************************************/
class Vertex
{
public:

  friend Container<Vertex>;
  using ContainerIterator = Container<Vertex>::List::iterator;
  using EdgeList = std::list<Edge*>;
  using FacetList = std::list<Facet*>;
  using VertexVector = std::vector<Vertex*>;

  /*------------------------------------------------------------------
  | Constructor 
  ------------------------------------------------------------------*/
  Vertex(double x, double y, double s=1.0, double r=1.0) 
  : xy_ {x, y}, sizing_ {s}, range_ {r}
  {}
  Vertex(const Vec2d& c, double s=1.0, double r=1.0) 
  : xy_ {c}, sizing_ {s}, range_ {r}
  {}

  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  void index (unsigned int i) { index_ = i; }
  void on_front( bool f ) { on_front_ = f; }
  void on_boundary( bool b ) { on_bdry_ = b; }
  void xy(const Vec2d& c) { xy_ = c; }
  void is_fixed( bool f ) { is_fixed_ = f; }

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  const Vec2d xy() const { return xy_; }
  double sizing() const { return sizing_; }
  double range() const { return range_; }
  const ContainerIterator& pos() const { return pos_; }
  unsigned int index() const { return index_; }
  bool on_front() const { return on_front_; }
  bool on_boundary() const { return on_bdry_; }
  bool is_fixed() const { return is_fixed_; }

  EdgeList& edges() { return edges_;}
  const EdgeList& edges() const { return edges_;}

  const Edge& edges(size_t i) const 
  { 
    ASSERT(!( i < 0 || i >= edges_.size() ),
            "Invalid access to vertex edge list." );
    auto iter = edges_.begin();
    std::advance( iter, i );
    ASSERT( *iter, "Vertex edge is nullptr.");
    return *(*iter);
  }

  const FacetList& facets() const { return facets_;}
  const Facet& facets(size_t i) const
  {
    ASSERT(!( i < 0 || i >= facets_.size() ),
            "Invalid access to vertex triangle list." );
    auto iter = facets_.begin();
    std::advance( iter, i );
    ASSERT( *iter, "Vertex triagle is nullptr.");
    return *(*iter);
  }

  const VertexVector& vertices() const { return verts_; }
  VertexVector& vertices() { return verts_; }
  const Vertex* adjacent_vertex(size_t i) const { return verts_[i]; }
  Vertex* adjacent_vertex(size_t i) { return verts_[i]; }

  /*------------------------------------------------------------------
  | Add / remove adjacent simplices 
  ------------------------------------------------------------------*/
  void add_edge(Edge& e) { edges_.push_back(&e); }
  void remove_edge(Edge& e) { edges_.remove(&e); }

  void add_facet(Facet& t) { facets_.push_back(&t); }
  void remove_facet(Facet& t) { facets_.remove(&t); }

  /*------------------------------------------------------------------
  | Functions for adjacency checks
  ------------------------------------------------------------------*/
  bool is_adjacent(const Edge& q) 
  {
    for ( auto e : edges_ )
      if ( &q == e ) return true;
    return false;
  }

  bool is_adjacent(const Facet& q) 
  {
    for ( auto f : facets_ )
      if ( &q == f ) return true;
    return false;
  }

  /*------------------------------------------------------------------
  | Check if vertex intersects with
  ------------------------------------------------------------------*/
  template <typename T>
  bool intersects_facet(const Container<T>& facets,
                        const double range) const
  {
    for ( const auto& f : facets.get_items(xy_, range) )
      if ( f->intersects_vertex( *this ) )
        return true;

    return false;

  } // intersects_facet()

  /*------------------------------------------------------------------
  | Check if the facets that are adjacent to the vertex do interesect
  | with each other 
  ------------------------------------------------------------------*/
  template <typename T>
  bool adjacent_facet_intersection(const std::list<T*>& facets) const
  {
    for ( const auto& f : facets )
    {
      for ( const auto& n : facets )
      {
        if ( f == n ) 
          continue;

        for (size_t i = 0; i < f->n_vertices(); ++i)
          if ( n->intersects_vertex( f->vertex(i) ) )
            return true;
      }
    }
    return false;
  }

  /*------------------------------------------------------------------
  | Mandatory container functions 
  ------------------------------------------------------------------*/
  bool in_container() const { return in_container_; }

private:
  /*------------------------------------------------------------------
  | Mandatory container functions 
  ------------------------------------------------------------------*/
  void container_destructor() {}

  /*------------------------------------------------------------------
  | Vertex attributes 
  ------------------------------------------------------------------*/
  Vec2d               xy_;
  EdgeList            edges_    { };
  FacetList           facets_   { };
  VertexVector        verts_    { };

  double              sizing_   { 1.0 };
  double              range_    { 1.0 };
  unsigned int        index_    {  0  };
  bool                on_front_ { false };
  bool                on_bdry_  { false };
  bool                is_fixed_ { false };

  // Mandatory container attributes
  ContainerIterator   pos_;
  bool                in_container_;

}; // Vertex 

/*********************************************************************
* Define general vertex container declaration
*********************************************************************/
using Vertices = Container<Vertex>;

/*********************************************************************
* Vertex ostream operator overload
*********************************************************************/
static std::ostream& operator<<(std::ostream& os, 
                                const Vertex& v)
{ return os << v.xy(); }

/***********************************************************
* Vertex equality operator 
***********************************************************/
static bool operator==(const Vertex& v1, const Vertex& v2)
{ return &v1 == &v2; }
static bool operator!=(const Vertex& v1, const Vertex& v2)
{ return !(v1 == v2); }

} // namespace TQAlgorithm
} // namespace TQMesh
