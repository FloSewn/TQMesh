/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <stdexcept>

#include "Container.h"
#include "utils.h"
#include "Vec2.h"
#include "geometry.h"

#include "Edge.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace TQUtils;

/*********************************************************************
* This class represents a directed list of edges
*
*
*********************************************************************/
class EdgeList
{
public:

  using iterator       = Container<Edge>::iterator;
  using const_iterator = Container<Edge>::const_iterator;

  iterator begin() { return edges_.begin(); }
  iterator end() { return edges_.end(); }

  const_iterator begin() const { return edges_.begin(); }
  const_iterator end() const { return edges_.end(); }

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  EdgeList( TQGeom::Orientation orient ) : orient_ { orient }
  {
    ASSERT( ( orient != TQGeom::Orientation::CL ),
        "Invalid edge list orientation.");
  }

  /*------------------------------------------------------------------
  | Copy Constructor
  ------------------------------------------------------------------*/
  EdgeList(const EdgeList& el) : orient_ { el.orient_ }
  {
    for ( auto& e : el.edges_ )
      add_edge( e->v1(), e->v2(), e->marker() );
  }

  /*------------------------------------------------------------------
  | Move Constructor
  ------------------------------------------------------------------*/
  EdgeList(EdgeList&& el) 
  : orient_ { el.orient_ }
  , edges_ { std::move( el.edges_ ) }
  , area_ { el.area_ }
  {}

  /*------------------------------------------------------------------
  | This function takes care of the garbage collection 
  ------------------------------------------------------------------*/
  void clear_waste()
  {
    edges_.clear_waste();
  }

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  size_t size() const { return edges_.size(); }
  double area() const { return area_; }
  bool is_ccw() const { return (orient_ == TQGeom::Orientation::CCW); }
  TQGeom::Orientation orient() const { return orient_; }

  const auto& edges() const { return edges_; }
  auto& edges() { return edges_; }

  /*------------------------------------------------------------------
  | Access operator
  ------------------------------------------------------------------*/
  const Edge& operator[](int i) const { return edges_[i]; }
  Edge& operator[](int i) { return edges_[i]; }

  /*------------------------------------------------------------------
  | Insert a new edge to the list at a specified position
  ------------------------------------------------------------------*/
  Edge& insert_edge(auto pos, Vertex& v1, Vertex& v2, 
                    int marker=TQ_INTR_EDGE_MARKER)
  {
    Edge& e = edges_.insert(pos, v1, v2, *this, marker);
    if ( orient_ != TQGeom::Orientation::NONE )
      compute_area();
    
    // Mark the added objects somehow
    mark_objects(v1, v2, e);

    return e;

  } // EdgeList::insert_edge() 

  /*------------------------------------------------------------------
  | Push a new edge back to the edge list
  | By default, all edges get colored by the interior edge marker.
  |
  | > EdgeLists with clockwise / counter-clockwise orientation must 
  |   be closed and must be made up of a single edge path (each 
  |   vertex is connected to a maximum of two edges from such a list) 
  | 
  | > EdgeList with NONE orientation can be open, vertices can be 
  |   connected to more than two edges of this list type
  ------------------------------------------------------------------*/
  Edge& add_edge(Vertex& v1, Vertex& v2, 
                 int marker=TQ_INTR_EDGE_MARKER)
  { 
    if ( orient_ == TQGeom::Orientation::NONE || edges_.size() < 1 )
      return insert_edge( edges_.end(), v1, v2, marker ); 

    // Check that the new edge is connected to the last edge
    const Vertex& v_last = edges_.back().v2();

    if ( v1 != v_last )
      throw std::runtime_error(
          "Failed to add an Edge to EdgeList. The new edge to add "
          "is not connected to the last Edge in the EdgeList.");

    return insert_edge( edges_.end(), v1, v2, marker ); 
  } 

  /*------------------------------------------------------------------
  | Remove an edge from the edge list
  ------------------------------------------------------------------*/
  bool remove(Edge& edge) { return edges_.remove( edge ); }

  /*------------------------------------------------------------------
  | Clear all edges and eventually the associated vertices
  ------------------------------------------------------------------*/
  void clear_edges()
  {
    size_t n_edges = edges_.size();

    for ( size_t i = 0; i < n_edges; ++i )
      remove( edges_[0] );

    edges_.clear_waste();

  } // EdgeList::clear_edges();

  /*------------------------------------------------------------------
  | Split a given edge. The value <s> must be in the interval (0,1)
  | and controls where the edge will be split
  | * s -> 0 ---> split near vertex v1
  | * s -> 1 ---> split near vertex v2
  | The function returns a pair of pointers to the generated edges.
  | The old edge is removed from the edgelist.
  ------------------------------------------------------------------*/
  std::pair<Edge*, Edge*>
  split_edge(Edge& edge, Vertices& vertices, const double s,
             bool fix_new_vertex=false) 
  {
    if ( &edge.edgelist() != this )
      return std::move( std::pair<Edge*,Edge*>(nullptr,nullptr) );

    if ( s >= 1.0 || s <= 0.0 )
      return std::move( std::pair<Edge*,Edge*>(nullptr,nullptr) );

    const double q = (1.0-s);

    Vertex& v1 = edge.v1();
    Vertex& v2 = edge.v2();

    const Vec2d xy_new  = s * v1.xy()     + q * v2.xy();
    const double sizing = s * v1.sizing() + q * v2.sizing();
    const double range  = s * v1.range()  + q * v2.range();

    // Place new vertex between v1 and v2
    Vertex& v_new = vertices.insert(v2.pos(), xy_new, sizing, range);

    if ( fix_new_vertex )
      v_new.is_fixed( true );

    int marker = edge.marker();
    bool on_boundary = edge.on_boundary();

    Edge& e1_new = this->insert_edge(edge.pos(), v1, v_new, marker);
    Edge& e2_new = this->insert_edge(edge.pos(), v_new, v2, marker);

    v1.on_boundary( on_boundary );
    v_new.on_boundary( on_boundary );
    v2.on_boundary( on_boundary );

    // Remove old edge
    this->remove(edge);


    return std::move( std::pair<Edge*,Edge*>( &e1_new, &e2_new) );

  } // EdgeList::split_edge()

  /*------------------------------------------------------------------
  | Check if a simplex is inside the area that is surrounded by the 
  | edges. If the object is located on the edge semgents, it is 
  | treated as lying inside.
  | Source: http://alienryderflex.com/polygon/
  ------------------------------------------------------------------*/
  template <typename T>
  bool is_inside(const T& s) const 
  {
    if (edges_.size() < 3)
      return false;

    const Vec2d obj = s.xy();
    
    int count = 0;

    for ( const auto& e_ptr : edges_ )
    {
      const Vec2d v1 = e_ptr->v1().xy();
      const Vec2d v2 = e_ptr->v2().xy();

      // Object equals edge vertex
      if ( obj == v1  ||  obj == v2 )
        return true;

      // Object on edge
      if ( EQ( obj.y, v2.y ) && EQ( obj.y, v1.y ) )
        if ( TQGeom::in_on_segment(v1, v2, obj) )
          return true;

      // Crossing lines
      if (  ( obj.y > v2.y  &&  obj.y <= v1.y )
         || ( obj.y > v1.y  &&  obj.y <= v2.y ) )
      {
        double dx = (v1.x - v2.x);
        double m  = (obj.y - v2.y) / (v1.y - v2.y);
        double x  = v2.x + m * dx;
        if ( x < obj.x ) ++count;
      }
    }

    return ( (count&1) == 1 ); // := (count%2 == 1)

  } // EdgeList::is_inside()

  /*------------------------------------------------------------------
  | Check if the edgelist is traversable for given start and ending
  | vertices
  ------------------------------------------------------------------*/
  bool is_traversable(Edge& e_start, Edge& e_end)
  {
    Edge* e_cur      = &e_start;
    int   edge_count = 0;

    int   n_edges    = static_cast<int>(edges_.size());

    do
    {
      e_cur = e_cur->get_next_edge();
      ++edge_count;

    } while (  ( e_cur ) 
            && ( edge_count < n_edges ) 
            && ( e_cur != &e_end ) );

    if ( e_cur != &e_end ) 
      return false;

    return true;

  } // EdgeList::is_traversable()

  /*------------------------------------------------------------------
  | Search for an edge that is part of this edgelist and which 
  | connects the given vertices v1 and v2
  | If <dir> is set to true, the edge direction v1->v2 is also 
  | considered in the search
  ------------------------------------------------------------------*/
  Edge* 
  get_edge(const Vertex& v1, const Vertex& v2, bool dir=false) const
  {
    // Consider edge direction
    if (dir)
    {
      for ( const auto& e : v1.edges() )
      {
        if ( &e->edgelist() != this )
          continue;

        if ( (e->v1() == v1 && e->v2() == v2) )
          return e;
      }
    }
    // Do not consider edge direction
    else
    {
      for ( const auto& e : v1.edges() )
      {
        if ( &e->edgelist() != this )
          continue;

        if ( (e->v1() == v1 && e->v2() == v2) ||
             (e->v2() == v1 && e->v1() == v2)  )
          return e;
      }
    }

    return nullptr;

  } // get_edge()

  /*------------------------------------------------------------------
  | Search for an edge that is part of this edgelist and which 
  | shares the given vertex v 
  | * If <pos> is set to 1, v must be the first vertex 
  |   of the target edge
  | * If <pos> is set to 2, v must be the second vertex 
  |   of the target edge
  | * If <pos> is set to 0, v can be first or second vertex 
  |   of the target edge
  ------------------------------------------------------------------*/
  Edge* 
  get_edge(const Vertex& v, int pos) const
  {
    Edge* found = nullptr;

    // Return the first edge with the prescribed marker, irregarding
    // of the vertex position in the found edge
    if ( pos == 0 )
    {
      for ( auto e : v.edges() )
        if (  &e->edgelist() == this )
        { found = e; break; }
    }
    // Return the first edge with the prescribed marker, where the 
    // vertex must be defined as the first vertex v1 of the edge
    else if ( pos == 1 )
    {
      for ( auto e : v.edges() )
        if (  &e->edgelist() == this && e->v1() == v )
        { found = e; break; }
    }
    // Return the first edge with the prescribed marker, where the 
    // vertex must be defined as the second vertex v2 of the edge
    else if ( pos == 2 )
    {
      for ( auto e : v.edges() )
        if (  &e->edgelist() == this && e->v2() == v )
        { found = e; break; }
    }
    else
    {
      MSG("Invalid argument for function Vertex::get_edge_from_vertex()");
    }

    return found;

  } // get_edge()


protected:
  /*------------------------------------------------------------------
  | Compute the area enclosed by all edges. 
  | Splits the edge list into triangles and sums up
  | the triangle areas.
  ------------------------------------------------------------------*/
  void compute_area()
  {
    if ( edges_.size() < 3 ) return;

    // Take first node of first edge as reference point
    auto e_ptr = edges_.begin();
    const Vec2d& ref = e_ptr->get()->v1().xy();

    // Define triangles through upcoming edges and reference 
    // node and sum up all triangle areas 
    double area {0.0};

    for ( ++e_ptr; e_ptr != edges_.end(); ++e_ptr)
    {
      const Vec2d& d1 = e_ptr->get()->v1().xy() - ref;
      const Vec2d& d2 = e_ptr->get()->v2().xy() - ref;

      area += cross(d1, d2);
    }

    // CW defined edge lists will have negative areas
    area_ = 0.5 * area;

  } // EdgeList::compute_area()

  /*------------------------------------------------------------------
  | Check the orientation of the edgelist for correctness
  ------------------------------------------------------------------*/
  bool check_orientation()
  {
    if ( orient_ == TQGeom::Orientation::NONE )
      return true;

    if (edges_.size() < 2)
      return true;

    if (edges_.size() < 3)
    {
      auto e1_ptr = edges_.begin();
      auto e2_ptr = std::next(e1_ptr);

      const Edge* e1 = e1_ptr->get();
      const Edge* e2 = e2_ptr->get();

      // Choose both nodes from e1
      const auto p = e1->v1().xy();
      const auto q = e1->v2().xy();

      const auto r = (  e2->v1() != e1->v1() 
                     && e2->v1() != e1->v2() ) 
                     ? e2->v1().xy() 
                     : e2->v2().xy();

      auto ori = TQGeom::orientation( p, q, r );

      // Allow colinear edges
      return ( ori == orient_ || ori == TQGeom::Orientation::CL );
    }
    else
    {
      const bool ccw = is_ccw();
      return (   ( ccw && area_ >= 0.0) 
              || (!ccw && area_ <= 0.0)  );
    }

  } // EdgeList::check_orientation()

  /*------------------------------------------------------------------
  | Used for derived classes to adjust new objects 
  ------------------------------------------------------------------*/
  virtual void mark_objects(Vertex& v1, Vertex& v2, Edge& e) 
  {}

  /*------------------------------------------------------------------
  | EdgeList attributes
  ------------------------------------------------------------------*/
  TQGeom::Orientation orient_;
  Container<Edge>     edges_ {};
  double              area_  {0.0};


}; // EdgeList


/*********************************************************************
* Print out edge list to std::cout
*********************************************************************/
static inline std::ostream& operator<<(std::ostream& os, 
                                       const EdgeList& el)
{
  for ( const auto& e : el.edges() )
    os << *e << "\n";

  return os;
} 

} // namespace TQAlgorithm
} // namespace TQMesh
