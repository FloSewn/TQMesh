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

using namespace CppUtils;

/*********************************************************************
* Forward declaration
*********************************************************************/
class Edge;
class Facet;

/*********************************************************************
* Different properties for mesh vertices 
*********************************************************************/
enum class VertexProperty : uint8_t {
  no_property   = 0b00000000,
  on_front      = 0b00000001,
  in_quad_layer = 0b00000010,
  is_fixed      = 0b00000100,
  on_boundary   = 0b00001000,
};

// Overloaded bitwise NOT operator
static inline VertexProperty operator~(VertexProperty prop) 
{ return static_cast<VertexProperty>(~static_cast<uint8_t>(prop)); }

// Overloaded compound bitwise OR operator for adding properties
static inline VertexProperty&
operator|=(VertexProperty& lhs, VertexProperty rhs) 
{ 
  lhs = static_cast<VertexProperty>(
    static_cast<uint8_t>(lhs) | static_cast<uint8_t>(rhs)
  );
  return lhs;
}

// Overloaded compound bitwise AND operator for removing properties
static inline VertexProperty& 
operator&=(VertexProperty& lhs, VertexProperty rhs) 
{
  lhs = static_cast<VertexProperty>(
    static_cast<uint8_t>(lhs) & static_cast<uint8_t>(rhs)
  );
  return lhs;
}

// Overloaded bitwise AND operator for checking vertex properties
static inline bool 
operator&(VertexProperty lhs, VertexProperty rhs) 
{ return (static_cast<uint8_t>(lhs) & static_cast<uint8_t>(rhs)) != 0; }



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
  Vertex(double x, double y, double s=0.0, double r=0.0) 
  : ContainerEntry(x,y), mesh_size_ {s}, size_range_ {r}
  {}

  Vertex(const Vec2d& c, double s=0.0, double r=0.0) 
  : ContainerEntry(c), mesh_size_ {s}, size_range_ {r}
  {}

  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  void mesh_size(double s) { mesh_size_ = s; }
  void size_range(double r) { size_range_ = r; }
  void index (unsigned int i) { index_ = i; }

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  double       mesh_size() const { return mesh_size_; }
  double       size_range() const { return size_range_; }
  unsigned int index() const { return index_; }


  /*------------------------------------------------------------------
  | Handle vertex properties 
  ------------------------------------------------------------------*/
  const VertexProperty& properties() const { return properties_; }
  void add_property(VertexProperty p) { properties_ |= p; }
  void set_property(VertexProperty p) { properties_ = p; }
  void remove_property(VertexProperty p) { properties_ &= ~p; }
  bool has_property(VertexProperty p) const { return (properties_ & p); }
  bool has_no_property() const 
  { return (properties_ == VertexProperty::no_property); }

  bool on_front() const { return has_property(VertexProperty::on_front); }
  bool on_boundary() const { return has_property(VertexProperty::on_boundary); }
  bool is_fixed() const { return has_property(VertexProperty::is_fixed); }
  bool in_quad_layer() const { return has_property(VertexProperty::in_quad_layer); }

  /*------------------------------------------------------------------
  | Change vertex coordinate 
  | ATTENTION: This function does not update the edges / facets that
  | are attached to this vertex!
  | -> Use MeshCleanup::set_vertex_coordinates() instead!!!
  ------------------------------------------------------------------*/
  void adjust_xy(const Vec2d& xy)
  {
    bool success = container_->update( *this, xy );
    ASSERT( success, "Vertex::adjust_xy(): "
        "Failed to update vertex coordinate.");
    (void) success;
  }

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

  void add_edge(Edge& e) { edges_.push_back(&e); }
  void remove_edge(Edge& e) { edges_.remove(&e); }

  void add_facet(Facet& t) { facets_.push_back(&t); }
  void remove_facet(Facet& t) { facets_.remove(&t); }

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


  /*------------------------------------------------------------------
  | Check if the vertex is located within a minimum distance to 
  | a mesh's edges
  ------------------------------------------------------------------*/
  template <typename Mesh>
  bool intersects_mesh_edges(const Mesh& mesh, 
                             const double search_range,
                             const double limit_range) const
  {
    const double limit_sqr = limit_range * limit_range;

    for ( const auto& e_ptr : mesh.get_intr_edges(this->xy_, search_range) )
    {
      const Vec2d& xy1 = e_ptr->v1().xy();
      const Vec2d& xy2 = e_ptr->v2().xy();
      const double dist_sqr = distance_point_edge_sqr(this->xy_, xy1, xy2);

      if (dist_sqr <= limit_sqr)
        return true;
    }

    for ( const auto& e_ptr : mesh.get_bdry_edges(this->xy_, search_range) )
    {
      const Vec2d& xy1 = e_ptr->v1().xy();
      const Vec2d& xy2 = e_ptr->v2().xy();
      const double dist_sqr = distance_point_edge_sqr(this->xy_, xy1, xy2);

      if (dist_sqr <= limit_sqr)
        return true;
    }

    return false;

  } // intersects_mesh_edges()

private:
  /*------------------------------------------------------------------
  | Vertex attributes 
  ------------------------------------------------------------------*/
  EdgeList            edges_    { };
  FacetList           facets_   { };

  double              mesh_size_  { 0.0 };
  double              size_range_ { 0.0 };
  unsigned int        index_      { 0   };
  VertexProperty      properties_ { VertexProperty::no_property};

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

} // namespace TQMesh
