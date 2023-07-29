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
* A simple quadrilateral - Must be defined CCW
* ============================================
*
*  Arrangement of quad neighbors:
*  ------------------------------
*
*                      f1
*           v3                     v2       vi... vertices
*             o<------------------o         ei... edges
*             |        e1         ^         fi... neighboring 
*             |                   |               facets
*             |                   |
*             |                   |
*         f2  | e2             e0 |  f0
*             |                   |
*             |                   |
*             |                   |
*             v         e3        |
*             o------------------>o
*           v0                     v1
*                       f3
*
*
*  Splitting into sub-quads:
*  -------------------------
*
*                      e1
*           v3                     v2       qi... sub-quads 
*             o<--------o---------o
*             |         |         ^
*             |         |         |
*             |   q3    |   q2    |
*             |         |         |
*         e2  o---------o---------o  e0
*             |         |         |
*             |   q0    |   q1    |
*             |         |         |
*             v         |         |
*             o---------o-------->o
*           v0                     v1
*                       e3
*
*
*********************************************************************/
class Quad : public Facet, public ContainerEntry<Quad>
{
public:

  using FacetArray = std::array<Facet*,4>;
  using VertexArray = std::array<Vertex*,4>;

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Quad(Vertex& v1, Vertex& v2, Vertex& v3, Vertex& v4)
  : ContainerEntry<Quad> { QuadGeometry::calc_centroid(v1, v2, v3, v4) }
  , vertices_ {&v1, &v2, &v3, &v4}
  {
    std::fill(facets_.begin(), facets_.end(), &NullFacet::get_instance()); 

    update_metrics(false);

    vertices_[0]->add_facet( *this );
    vertices_[1]->add_facet( *this );
    vertices_[2]->add_facet( *this );
    vertices_[3]->add_facet( *this );

    ASSERT( (edge_lengths_[0] > 0.0), "Invalid quad: Vertices collapse.");
    ASSERT( (edge_lengths_[1] > 0.0), "Invalid quad: Vertices collapse.");
    ASSERT( (edge_lengths_[2] > 0.0), "Invalid quad: Vertices collapse.");
    ASSERT( (edge_lengths_[3] > 0.0), "Invalid quad: Vertices collapse.");
  }

  /*------------------------------------------------------------------
  | Container destructor function
  ------------------------------------------------------------------*/
  void container_destructor() override
  {
    if (vertices_[0]) vertices_[0]->remove_facet( *this );
    if (vertices_[1]) vertices_[1]->remove_facet( *this );
    if (vertices_[2]) vertices_[2]->remove_facet( *this );
    if (vertices_[3]) vertices_[3]->remove_facet( *this );

    vertices_[0] = nullptr;
    vertices_[1] = nullptr;
    vertices_[2] = nullptr;
    vertices_[3] = nullptr;
  }

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  const Vertex& vertex(size_t i) const { return *vertices_[i]; }
  Vertex&       vertex(size_t i) { return *vertices_[i]; }
  const Vertex& v1() const { return *vertices_[0]; }
  Vertex&       v1() { return *vertices_[0]; }
  const Vertex& v2() const { return *vertices_[1]; }
  Vertex&       v2() { return *vertices_[1]; }
  const Vertex& v3() const { return *vertices_[2]; }
  Vertex&       v3() { return *vertices_[2]; }
  const Vertex& v4() const { return *vertices_[3]; }
  Vertex&       v4() { return *vertices_[3]; }

  const Facet*  neighbor(size_t i) const { return facets_[i]; }
  Facet*        neighbor(size_t i) { return facets_[i]; }
  const Facet*  nbr1() const { return facets_[0]; }
  Facet*        nbr1() { return facets_[0]; }
  const Facet*  nbr2() const { return facets_[1]; }
  Facet*        nbr2() { return facets_[1]; }
  const Facet*  nbr3() const { return facets_[2]; }
  Facet*        nbr3() { return facets_[2]; }
  const Facet*  nbr4() const { return facets_[3]; }
  Facet*        nbr4() { return facets_[3]; }

  size_t        n_vertices() const override { return 4; }
  const Vec2d&  xy() const override { return ContainerEntry<Quad>::xy_; }
  const Vec2d&  circumcenter() const { return circumcenter_; }
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
  void nbr4(Facet* f) { facets_[3] = f; }

  void mesh(Mesh* m) { mesh_ = m; }
  void color(int i) override { color_ = i; }
  void index(int i) override { index_ = i; }
  void is_active(bool a) { active_ = a; }
  //void marker(bool c){ marker_ = c; }

  /*------------------------------------------------------------------
  | Returns the index of a quad vertex for a given input vertex
  | Returns -1 if no vertex is found
  ------------------------------------------------------------------*/
  int get_vertex_index(const Vertex& v) const override
  {
    if ( &v == vertices_[0] ) return 0;
    if ( &v == vertices_[1] ) return 1;
    if ( &v == vertices_[2] ) return 2;
    if ( &v == vertices_[3] ) return 3;
    return -1;
  } 

  /*------------------------------------------------------------------
  | Returns the index of a quad edge for two given input vertices
  | Returns -1 if no edge is found
  |
  |
  | vertices_[3]      e1     vertices_[2]
  |        x--------------------x
  |        |                    |
  |        |                    |
  |        |                    |
  |     e2 |                    | e0
  |        |                    |
  |        |                    |
  |        |                    |
  |        x--------------------x
  |  vertices_[0]    e3      vertices_[1]
  |
  ------------------------------------------------------------------*/
  int get_edge_index(const Vertex& v1, const Vertex& v2) const override
  {
    if ( (&v1==vertices_[0] && &v2==vertices_[1]) || 
         (&v1==vertices_[1] && &v2==vertices_[0]) )
      return 3;

    if ( (&v1==vertices_[1] && &v2==vertices_[2]) || 
         (&v1==vertices_[2] && &v2==vertices_[1]) )
      return 0;

    if ( (&v1==vertices_[2] && &v2==vertices_[3]) || 
         (&v1==vertices_[3] && &v2==vertices_[2]) )
      return 1;

    if ( (&v1==vertices_[3] && &v2==vertices_[0]) || 
         (&v1==vertices_[0] && &v2==vertices_[3]) )
      return 2;

    return -1;
  } 


  /*------------------------------------------------------------------
  | Returns true if the quad intersects with a vertex 
  | --> Vertex is located within the quad or on its
  |     edges
  ------------------------------------------------------------------*/
  bool intersects_vertex(const Vertex& v) const override
  { return QuadGeometry::check_intersection(*this, v); }

  /*------------------------------------------------------------------
  | Returns true, if any quad vertex is not within a given domain.
  | Triangle vertices are allowed to be located on domain edges.
  ------------------------------------------------------------------*/
  bool intersects_domain(const Domain& domain) const
  { return QuadGeometry::check_intersection(*this, domain); }

  /*------------------------------------------------------------------
  | Returns true if the quad intersects with a triangle 
  | in a given Container of Triangles 
  | The factor range scales the vicinity range from which 
  | triangles to pick from
  ------------------------------------------------------------------*/
  template <typename T>
  bool intersects_triangle(const Container<T>& tris,
                           const double range) const
  { return QuadGeometry::check_intersection(*this, tris, range); }

  /*------------------------------------------------------------------
  | Returns true if the quad intersects with a quad 
  | in a given Container of Quadrilaterals 
  | The factor range scales the vicinity range from which 
  | triangles to pick from
  ------------------------------------------------------------------*/
  template <typename Q>
  bool intersects_quad(const Container<Q>& quads,
                       const double range) const
  { return QuadGeometry::check_intersection(*this, quads, range); }

  /*------------------------------------------------------------------
  | Returns true if a quad edge is too close to a vertex in a given 
  | advancing front 
  | The factor range scales the vicinity range from which 
  | vertices to pick from
  | The factor min_dist_sqr defines the minimum squared 
  | distance that an advancing front vertex must be located 
  | from a quad edge
  ------------------------------------------------------------------*/
  template <typename Front>
  bool intersects_front(const Front& front,
                        const double range) const
  { return QuadGeometry::check_intersection(*this, front, range); }

  /*------------------------------------------------------------------
  | Returns true if the quad encloses an advancing front vertex.
  | The factor range scales the vicinity range from which 
  | front edges to pick from
  ------------------------------------------------------------------*/
  bool intersects_vertex(const Vertices& verts,
                         const double range) const 
  { return QuadGeometry::check_intersection(*this, verts, range); }

  /*------------------------------------------------------------------
  | Compute the quad quality based on the local mesh scale h
  ------------------------------------------------------------------*/
  double quality(const double h) const 
  { return QuadGeometry::calc_quality(edge_lengths_, shape_factor_, h); }

  /*------------------------------------------------------------------
  | Returns true if the quad is valid
  ------------------------------------------------------------------*/
  bool is_valid() const
  { return QuadGeometry::check_validity(area_, edge_lengths_); }

  /*------------------------------------------------------------------
  | Update the quad if its vertices changed
  ------------------------------------------------------------------*/
  void update_metrics(bool update_centroid=true) override
  {
    const Vertex& v1 = *vertices_[0];
    const Vertex& v2 = *vertices_[1];
    const Vertex& v3 = *vertices_[2];
    const Vertex& v4 = *vertices_[3];

    if ( update_centroid )
    {
      Vec2d xy_new = QuadGeometry::calc_centroid(v1, v2, v3, v4);
      bool success = container_->update( *this, xy_new );
      ASSERT( success, "Quad::update_metrics(): "
          "Failed to update quad centroid.");
      (void) success;
    }

    area_            = QuadGeometry::calc_area( v1, v2, v3, v4 );
    circumcenter_    = QuadGeometry::calc_circumcenter( v1, v2, v3, v4 );
    circumradius_    = QuadGeometry::calc_circumradius( v1, v2, v3, v4, circumcenter_ );
    edge_lengths_[0] = QuadGeometry::calc_edge_length( v1, v2 );
    edge_lengths_[1] = QuadGeometry::calc_edge_length( v2, v3 );
    edge_lengths_[2] = QuadGeometry::calc_edge_length( v3, v4 );
    edge_lengths_[3] = QuadGeometry::calc_edge_length( v4, v1 );
    min_edge_length_ = edge_lengths_.min(); 
    max_edge_length_ = edge_lengths_.max(); 
    angles_          = QuadGeometry::calc_angles( v1, v2, v3, v4, edge_lengths_ );
    min_angle_       = angles_.min(); 
    max_angle_       = angles_.max(); 
    shape_factor_    = QuadGeometry::calc_shape_factor(edge_lengths_, area_);
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

  Vec4d                edge_lengths_    {0.0};
  Vec4d                angles_          {0.0};

}; // Quad

/*********************************************************************
* Define general quad container declaration
*********************************************************************/
using Quads = Container<Quad>;

/*********************************************************************
* Quad ostream overload 
*********************************************************************/
static std::ostream& operator<<(std::ostream& os, 
                                const Quad& q)
{
  return os << q.v1() << " -> " << q.v2() << " -> "
            << q.v3() << " -> " << q.v4();
}


/*********************************************************************
* Edge equality operator 
*********************************************************************/
static bool operator==(const Quad& q1, const Quad& q2)
{ return (&q1 == &q2); }
static bool operator!=(const Quad& q1, const Quad& q2)
{ return !(q1 == q2); }

} // namespace TQAlgorithm
} // namespace TQMesh
