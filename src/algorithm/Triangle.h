/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <iostream>
#include <list>
#include <array>

#include "VecND.h"
#include "Geometry.h"

#include "utils.h"
#include "Vertex.h"
#include "Facet.h"

#include "Domain.h"
#include "FacetGeometry.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

class Mesh;

/*********************************************************************
* A simple triangle - Must be defined CCW
* =======================================
*
*  Arrangement of tri neighbors:
*  -----------------------------
*
*                       v2                  v1... vertices
*                       o                   ei... edges
*                     /   \                 fi... neighboring 
*                   /       \                     facets
*            f1   /           \   f0
*               / e1         e0 \
*             /                   \
*           /                       \
*         /            e2             \
*       o-------------------------------o
*     v0                                 v1
*                      f2
*
*
*  Splitting into sub-quads:
*  -------------------------
*
*                       v2                  
*                       o                   
*                     /   \
*                   /  q2   \
*            e1   /           \   e0
*               o-------o-------o
*             /         |         \
*           /    q0     |    q1     \
*         /             |             \
*       o---------------o---------------o
*     v0                                 v1
*                      e2
*
*
*********************************************************************/
class Triangle : public Facet, public ContainerEntry<Triangle>
{
public:

  using FacetArray = std::array<Facet*,3>;
  using VertexArray = std::array<Vertex*,3>;

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Triangle(Vertex& v1, Vertex& v2, Vertex& v3)
  : ContainerEntry<Triangle> { TriangleGeometry::calc_centroid(v1, v2, v3) }
  , vertices_ {&v1, &v2, &v3}
  {
    std::fill(facets_.begin(), facets_.end(), &NullFacet::get_instance()); 

    update_metrics(false);

    vertices_[0]->add_facet( *this );
    vertices_[1]->add_facet( *this );
    vertices_[2]->add_facet( *this );

    ASSERT( (edge_lengths_[0] > 0.0), "Invalid triangle: Vertices collapse.");
    ASSERT( (edge_lengths_[1] > 0.0), "Invalid triangle: Vertices collapse.");
    ASSERT( (edge_lengths_[2] > 0.0), "Invalid triangle: Vertices collapse.");
  }
    
  /*------------------------------------------------------------------
  | Container destructor function
  ------------------------------------------------------------------*/
  void container_destructor() override
  {
    if (vertices_[0]) vertices_[0]->remove_facet( *this );
    if (vertices_[1]) vertices_[1]->remove_facet( *this );
    if (vertices_[2]) vertices_[2]->remove_facet( *this );

    vertices_[0] = nullptr;
    vertices_[1] = nullptr;
    vertices_[2] = nullptr;
  }

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  const Vertex& vertex(size_t i) const { return *vertices_[i]; }
  Vertex&       vertex(size_t i)       { return *vertices_[i]; }
  const Vertex& v1() const { return *vertices_[0]; }
  Vertex&       v1()       { return *vertices_[0]; }
  const Vertex& v2() const { return *vertices_[1]; }
  Vertex&       v2()       { return *vertices_[1]; }
  const Vertex& v3() const { return *vertices_[2]; }
  Vertex&       v3()       { return *vertices_[2]; }

  const Facet*  neighbor(size_t i) const { return facets_[i]; }
  Facet*        neighbor(size_t i)       { return facets_[i]; }
  const Facet*  nbr1() const { return facets_[0]; }
  Facet*        nbr1()       { return facets_[0]; }
  const Facet*  nbr2() const { return facets_[1]; }
  Facet*        nbr2()       { return facets_[1]; }
  const Facet*  nbr3() const { return facets_[2]; }
  Facet*        nbr3()       { return facets_[2]; }

  size_t        n_vertices() const override { return 3; }
  const  Vec2d& xy() const override { return ContainerEntry<Triangle>::xy_; }
  const  Vec2d& circumcenter() const { return circumcenter_; }
  Mesh*         mesh() const { return mesh_; }
  int           color() const override { return color_; }
  int           index() const override { return index_; }
  bool          is_active() const { return active_; }
  //bool          marker() const { return marker_; }
  double        area() const { return area_; }
  double        circumradius() const { return circumradius_; }
  double        min_angle() const { return min_angle_; }
  double        max_angle() const { return max_angle_; }
  double        edgelength(unsigned int i) const { return edge_lengths_[i]; }
  double        angle(unsigned int i) const { return angles_[i]; }
  double        min_edge_length() const { return min_edge_length_; }
  double        max_edge_length() const { return max_edge_length_; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  void neighbor(size_t i, Facet* f) override { facets_[i] = f; }
  void nbr1(Facet* f) { facets_[0] = f; }
  void nbr2(Facet* f) { facets_[1] = f; }
  void nbr3(Facet* f) { facets_[2] = f; }

  void mesh(Mesh* m) { mesh_ = m; }
  void color(int i) override { color_ = i; }
  void index(int i) override { index_ = i; }
  void is_active(bool a) { active_ = a; }
  //void marker(bool c){ marker_ = c; }

  /*------------------------------------------------------------------
  | Returns the index of a triangle vertex for a given input vertex
  | Returns -1 if no vertex is found
  ------------------------------------------------------------------*/
  int get_vertex_index(const Vertex& v) const override
  {
    if ( &v == vertices_[0] ) return 0;
    if ( &v == vertices_[1] ) return 1;
    if ( &v == vertices_[2] ) return 2;
    return -1;
  } 

  /*------------------------------------------------------------------
  | Returns the index of a triangle edge for two given input vertices
  | Returns -1 if no edge is found
  |
  |                   vertices_[2]
  |                  x
  |                 / \
  |                /   \
  |               /     \
  |          e1  /       \ e0
  |             /         \
  |            /           \
  |           /             \
  |          /               \
  |         x-----------------x
  |   vertices_[0]   e2    vertices_[1]
  |
  ------------------------------------------------------------------*/
  int get_edge_index(const Vertex& v1, const Vertex& v2) const override
  {
    if ( (&v1==vertices_[0] && &v2==vertices_[1]) || 
         (&v1==vertices_[1] && &v2==vertices_[0]) )
      return 2;

    if ( (&v1==vertices_[1] && &v2==vertices_[2]) || 
         (&v1==vertices_[2] && &v2==vertices_[1]) )
      return 0;

    if ( (&v1==vertices_[2] && &v2==vertices_[0]) || 
         (&v1==vertices_[0] && &v2==vertices_[2]) )
      return 1;

    return -1;
  } 

  /*------------------------------------------------------------------
  | Returns true if the triangle intersects with a vertex 
  | --> Vertex is located within the triangle or on its
  |     edges
  ------------------------------------------------------------------*/
  bool intersects_vertex(const Vertex& v) const override
  { return TriangleGeometry::check_intersection(*this, v); }

  /*------------------------------------------------------------------
  | Returns true, if any triangle vertex is not within a given domain.
  | Triangle vertices are allowed to be located on domain edges.
  | Domain vertices are allowed to be placed on triangle edges.
  ------------------------------------------------------------------*/
  bool intersects_domain(const Domain& domain) const 
  { return TriangleGeometry::check_intersection(*this, domain); }

  /*------------------------------------------------------------------
  | Returns true if the triangle intersects with a triangle 
  | in a given Container of Triangles 
  | The factor range scales the vicinity range from which 
  | triangles to pick from
  ------------------------------------------------------------------*/
  template <typename T>
  bool intersects_triangle(const Container<T>& tris,
                           const double range) const
  { return TriangleGeometry::check_intersection(*this, tris, range); }

  /*------------------------------------------------------------------
  | Returns true if the triangle intersects with a quad 
  | in a given Container of Quadrilaterals 
  | The factor range scales the vicinity range from which 
  | triangles to pick from
  ------------------------------------------------------------------*/
  template <typename T>
  bool intersects_quad(const Container<T>& quads,
                       const double range) const
  { return TriangleGeometry::check_intersection(*this, quads, range); }

  /*------------------------------------------------------------------
  | Returns true if a triangle edge intersects with an edge of 
  | the advancing front.  
  | The factor range scales the vicinity range from which 
  | front edges to pick from
  ------------------------------------------------------------------*/
  template <typename Front>
  bool intersects_front(const Front& front, const double range) const 
  { return TriangleGeometry::check_intersection(*this, front, range); }

  /*------------------------------------------------------------------
  | Returns true if the triangle encloses an advancing front vertex.
  | The factor range scales the vicinity range from which 
  | front edges to pick from
  ------------------------------------------------------------------*/
  bool intersects_vertex(const Vertices& verts,
                         const double range) const 
  { return TriangleGeometry::check_intersection(*this, verts, range); }

  /*------------------------------------------------------------------
  | Compute the triangle quality based on the local mesh scale h
  ------------------------------------------------------------------*/
  double quality(const double h) const 
  { return TriangleGeometry::calc_quality(edge_lengths_, shape_factor_, h); }

  /*------------------------------------------------------------------
  | Returns true if the triangle is valid
  ------------------------------------------------------------------*/
  bool is_valid() const
  { return TriangleGeometry::check_validity(area_, edge_lengths_); }

  /*------------------------------------------------------------------
  | Update the triangle if its vertices changed
  ------------------------------------------------------------------*/
  void update_metrics(bool update_centroid=true) override
  {
    const Vertex& v1 = *vertices_[0];
    const Vertex& v2 = *vertices_[1];
    const Vertex& v3 = *vertices_[2];

    if ( update_centroid )
    {
      Vec2d xy_new = TriangleGeometry::calc_centroid(v1, v2, v3);
      bool success = container_->update( *this, xy_new );
      ASSERT( success, "Triangle::update_metrics(): "
          "Failed to update triangle centroid.");
      (void) success;
    }

    area_            = TriangleGeometry::calc_area( v1, v2, v3 );
    circumcenter_    = TriangleGeometry::calc_circumcenter( v1, v2, v3 );
    circumradius_    = TriangleGeometry::calc_circumradius( v1, circumcenter_ );
    edge_lengths_[0] = TriangleGeometry::calc_edge_length( v1, v2 );
    edge_lengths_[1] = TriangleGeometry::calc_edge_length( v2, v3 );
    edge_lengths_[2] = TriangleGeometry::calc_edge_length( v3, v1 );
    min_edge_length_ = edge_lengths_.min(); 
    max_edge_length_ = edge_lengths_.max(); 
    angles_          = TriangleGeometry::calc_angles( v1, v2, v3, edge_lengths_ );
    min_angle_       = angles_.min(); 
    max_angle_       = angles_.max(); 
    shape_factor_    = TriangleGeometry::calc_shape_factor(edge_lengths_, area_);

  }

private:
  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  VertexArray          vertices_;
  FacetArray           facets_;

  int                  color_           {DEFAULT_ELEMENT_COLOR};
  int                  index_           {-1};
  bool                 active_          {false};

  Mesh*                mesh_            {nullptr};

  Vec2d                circumcenter_    {0.0,0.0};
  double               area_            {0.0};
  double               circumradius_    {0.0};
  double               min_angle_       {0.0};
  double               max_angle_       {0.0};
  double               shape_factor_    {0.0};
  double               quality_         {0.0};
  double               min_edge_length_ {0.0};
  double               max_edge_length_ {0.0};
 
  Vec3d                edge_lengths_    {0.0};
  Vec3d                angles_          {0.0};

}; // Triangle

/*********************************************************************
* Define general triangle container declaration
*********************************************************************/
using Triangles = Container<Triangle>;

/*********************************************************************
* Triangle ostream overload 
*********************************************************************/
static std::ostream& operator<<(std::ostream& os, const Triangle& t)
{ return os << t.v1() << " -> " << t.v2() << " -> " << t.v3(); }

/*********************************************************************
* Edge equality operator 
*********************************************************************/
static bool operator==(const Triangle& t1, const Triangle& t2)
{ return (&t1 == &t2); }
static bool operator!=(const Triangle& t1, const Triangle& t2)
{ return !(t1 == t2); }

} // namespace TQAlgorithm
} // namespace TQMesh
