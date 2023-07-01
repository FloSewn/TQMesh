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
#include "VecND.h"
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
class Vertex : public ContainerEntry<Vertex>
{
public:
  using EdgeList = std::list<Edge*>;
  using FacetList = std::list<Facet*>;
  using VertexVector = std::vector<Vertex*>;

  /*------------------------------------------------------------------
  | Constructor 
  ------------------------------------------------------------------*/
  Vertex(double x, double y, double s=1.0, double r=1.0) 
  : ContainerEntry(x,y), sizing_ {s}, range_ {r}
  {}
  Vertex(const Vec2d& c, double s=1.0, double r=1.0) 
  : ContainerEntry(c), sizing_ {s}, range_ {r}
  {}

  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  void index (unsigned int i) { index_ = i; }
  void on_front( bool f ) { on_front_ = f; }
  void on_boundary( bool b ) { on_bdry_ = b; }
  void is_fixed( bool f ) { is_fixed_ = f; }

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  double       sizing() const { return sizing_; }
  double       range() const { return range_; }
  unsigned int index() const { return index_; }
  bool         on_front() const { return on_front_; }
  bool         on_boundary() const { return on_bdry_; }
  bool         is_fixed() const { return is_fixed_; }

  /*------------------------------------------------------------------
  | Access edges that are adjacent to this vertex 
  ------------------------------------------------------------------*/
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

  /*------------------------------------------------------------------
  | Access facets that are adjacent to this vertex 
  ------------------------------------------------------------------*/
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

  /*------------------------------------------------------------------
  | Access vertices that are adjacent to this vertex 
  ------------------------------------------------------------------*/
  const VertexVector& vertices() const { return verts_; }
  VertexVector& vertices() { return verts_; }
  const Vertex* adjacent_vertex(size_t i) const { return verts_[i]; }
  Vertex* adjacent_vertex(size_t i) { return verts_[i]; }

  /*------------------------------------------------------------------
  | Add / remove adjacent simplices 
  ------------------------------------------------------------------*/
  void add_vertex(Vertex& v) { verts_.push_back(&v); }
  void remove_vertex(Vertex& v) 
  { verts_.erase(std::remove(verts_.begin(), verts_.end(), &v), verts_.end()); }

  void add_edge(Edge& e) { edges_.push_back(&e); }
  void remove_edge(Edge& e) { edges_.remove(&e); }

  void add_facet(Facet& t) { facets_.push_back(&t); }
  void remove_facet(Facet& t) { facets_.remove(&t); }

  /*------------------------------------------------------------------
  | Functions for adjacency checks
  ------------------------------------------------------------------*/
  bool is_adjacent(const Vertex& q) 
  {
    for ( auto v : verts_ )
      if ( &q == v ) return true;
    return false;
  }

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
    for ( const auto& f : facets.get_items(this->xy_, range) )
      if ( f->intersects_vertex( *this ) )
        return true;

    return false;

  } // intersects_facet()

private:
  /*------------------------------------------------------------------
  | Vertex attributes 
  ------------------------------------------------------------------*/
  EdgeList            edges_    { };
  FacetList           facets_   { };
  VertexVector        verts_    { };

  double              sizing_   { 1.0 };
  double              range_    { 1.0 };
  unsigned int        index_    {  0  };
  bool                on_front_ { false };
  bool                on_bdry_  { false };
  bool                is_fixed_ { false };

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

/*********************************************************************
* Vertex equality operator 
*********************************************************************/
static bool operator==(const Vertex& v1, const Vertex& v2)
{ return &v1 == &v2; }
static bool operator!=(const Vertex& v1, const Vertex& v2)
{ return !(v1 == v2); }

} // namespace TQAlgorithm
} // namespace TQMesh
