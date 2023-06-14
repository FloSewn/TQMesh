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

#include "Vec2.h"
#include "Geometry.h"

#include "utils.h"
#include "Vertex.h"
#include "Facet.h"

#include "Front.h"
#include "Domain.h"

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

  using DoubleArray = std::array<double,4>;
  using FacetArray = std::array<Facet*,4>;
  using VertexArray = std::array<Vertex*,4>;

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Quad(Vertex& v1, Vertex& v2, Vertex& v3, Vertex& v4)
  : ContainerEntry<Quad> { calc_centroid(v1, v2, v3, v4) }
  , v_ {&v1, &v2, &v3, &v4}
  {
    area_            = calc_area( v1, v2, v3, v4 );
    circumcenter_    = calc_circumcenter( v1, v2, v3, v4 );
    circumradius_    = calc_circumradius( v1, v2, v3, v4, circumcenter_ );
    edge_lengths_[0] = calc_edge_length( v1, v2 );
    edge_lengths_[1] = calc_edge_length( v2, v3 );
    edge_lengths_[2] = calc_edge_length( v3, v4 );
    edge_lengths_[3] = calc_edge_length( v4, v1 );
    min_edge_length_ = calc_minimum_edge_length( edge_lengths_ );
    max_edge_length_ = calc_maximum_edge_length( edge_lengths_ );
    angles_          = calc_angles( v1, v2, v3, v4, edge_lengths_ );
    min_angle_       = calc_minimum_angle(angles_);
    max_angle_       = calc_maximum_angle(angles_);
    shape_factor_    = calc_shape_factor(edge_lengths_, area_);

    v_[0]->add_facet( *this );
    v_[1]->add_facet( *this );
    v_[2]->add_facet( *this );
    v_[3]->add_facet( *this );

    ASSERT( (edge_lengths_[0] > 0.0), "Invalid quad: Vertices collapse.");
    ASSERT( (edge_lengths_[1] > 0.0), "Invalid quad: Vertices collapse.");
    ASSERT( (edge_lengths_[2] > 0.0), "Invalid quad: Vertices collapse.");
    ASSERT( (edge_lengths_[3] > 0.0), "Invalid quad: Vertices collapse.");
  }

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  const Vertex& vertex(size_t i) const { return *v_[i]; }
  const Vertex& v1() const { return *v_[0]; }
  const Vertex& v2() const { return *v_[1]; }
  const Vertex& v3() const { return *v_[2]; }
  const Vertex& v4() const { return *v_[3]; }
  size_t n_vertices() const { return 4; }

  Vertex& vertex(size_t i) { return *v_[i]; }
  Vertex& v1() { return *v_[0]; }
  Vertex& v2() { return *v_[1]; }
  Vertex& v3() { return *v_[2]; }
  Vertex& v4() { return *v_[3]; }

  const Facet* neighbor(size_t i) const { return f_[i]; }
  const Facet* nbr1() const { return f_[0]; }
  const Facet* nbr2() const { return f_[1]; }
  const Facet* nbr3() const { return f_[2]; }
  const Facet* nbr4() const { return f_[3]; }

  Facet* neighbor(size_t i) { return f_[i]; }
  Facet* nbr1() { return f_[0]; }
  Facet* nbr2() { return f_[1]; }
  Facet* nbr3() { return f_[2]; }
  Facet* nbr4() { return f_[3]; }

  const Vec2d& xy() const override { return ContainerEntry<Quad>::xy_; }
  const Vec2d& circumcenter() const { return circumcenter_; }

  Mesh* mesh() const { return mesh_; }
  int color() const { return color_; }
  int index() const { return index_; }
  bool is_active() const { return active_; }
  bool marker() const { return marker_; }

  double area() const { return area_; }
  double circumradius() const { return circumradius_; }
  double min_angle() const { return min_angle_; }
  double max_angle() const { return max_angle_; }

  double edgelength(unsigned int i) const { return edge_lengths_[i]; }
  double angle(unsigned int i) const { return angles_[i]; }

  double min_edge_length() const { return min_edge_length_; }
  double max_edge_length() const { return max_edge_length_; }

  double quality(const double h) const 
  { 
    const double e1 = edge_lengths_[0];
    const double e2 = edge_lengths_[1];
    const double e3 = edge_lengths_[2];
    const double e4 = edge_lengths_[3];

    const double f1_1 = e1 / h;
    const double f1_2 = e2 / h;
    const double f1_3 = e3 / h;
    const double f1_4 = e4 / h;

    const double f2_1 = h / e1;
    const double f2_2 = h / e2;
    const double f2_3 = h / e3;
    const double f2_4 = h / e4;

    const double q1 = MIN(f1_1, f2_1);
    const double q2 = MIN(f1_2, f2_2);
    const double q3 = MIN(f1_3, f2_3);
    const double q4 = MIN(f1_4, f2_4);
    
    return shape_factor_ * q1 * q2 * q3 * q4; 
  }


  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  void neighbor(size_t i, Facet* f) { f_[i] = f; }
  void nbr1(Facet* f) { f_[0] = f; }
  void nbr2(Facet* f) { f_[1] = f; }
  void nbr3(Facet* f) { f_[2] = f; }
  void nbr4(Facet* f) { f_[3] = f; }

  void mesh(Mesh* m) { mesh_ = m; }
  void color(int i) { color_ = i; }
  void index(int i) { index_ = i; }
  void is_active(bool a) { active_ = a; }
  void marker(bool c){ marker_ = c; }

  /*------------------------------------------------------------------
  | Returns true if the quad is valid
  ------------------------------------------------------------------*/
  bool is_valid() const
  {
    if ( area_ <= 0.0 )
    {
      DEBUG_LOG("  | NON-POSITIVE QUAD AREA " << area_);
      return false;
    }
    return true;
  }

  /*------------------------------------------------------------------
  | Returns the index of a quad vertex for a given input vertex
  | Returns -1 if no vertex is found
  ------------------------------------------------------------------*/
  int get_vertex_index(const Vertex& v) const
  {
    if ( &v == v_[0] )
      return 0;
    if ( &v == v_[1] )
      return 1;
    if ( &v == v_[2] )
      return 2;
    if ( &v == v_[3] )
      return 3;

    return -1;

  } //Quad::get_vertex_index()

  /*------------------------------------------------------------------
  | Returns the index of a quad edge for two given input vertices
  | Returns -1 if no edge is found
  |
  |
  |      v_[3]       e1        v_[2]
  |        x--------------------x
  |        |                    |
  |        |                    |
  |        |                    |
  |     e2 |                    | e0
  |        |                    |
  |        |                    |
  |        |                    |
  |        x--------------------x
  |       v_[0]     e3         v_[1]
  |
  ------------------------------------------------------------------*/
  int get_edge_index(const Vertex& v1, const Vertex& v2) const
  {
    if ( (&v1==v_[0] && &v2==v_[1]) || (&v1==v_[1] && &v2==v_[0]) )
      return 3;

    if ( (&v1==v_[1] && &v2==v_[2]) || (&v1==v_[2] && &v2==v_[1]) )
      return 0;

    if ( (&v1==v_[2] && &v2==v_[3]) || (&v1==v_[3] && &v2==v_[2]) )
      return 1;

    if ( (&v1==v_[3] && &v2==v_[0]) || (&v1==v_[0] && &v2==v_[3]) )
      return 2;

    return -1;

  } // get_edge_index()


  /*------------------------------------------------------------------
  | Returns true if the quad intersects with a vertex 
  | --> Vertex is located within the quad or on its
  |     edges
  ------------------------------------------------------------------*/
  bool intersects_vertex(const Vertex& v) const
  {
    if (  v == *v_[0] || v == *v_[1] || v == *v_[2] || v == *v_[3] )
      return false;

    return in_on_quad(v.xy(), 
        v_[0]->xy(), v_[1]->xy(), v_[2]->xy(), v_[3]->xy());

  } // Quad::intersects_vertex() 

  /*------------------------------------------------------------------
  | Returns true, if any quad vertex is not within a given domain.
  | Triangle vertices are allowed to be located on domain edges.
  ------------------------------------------------------------------*/
  bool intersects_domain(const Domain& domain) const
  {
    return !( domain.is_inside( *v_[0] ) 
           && domain.is_inside( *v_[1] ) 
           && domain.is_inside( *v_[2] ) 
           && domain.is_inside( *v_[3] ) );

  } // Quad::intersects_domain()

  /*------------------------------------------------------------------
  | Returns true if the quad intersects with a triangle 
  | in a given Container of Triangles 
  | The factor range scales the vicinity range from which 
  | triangles to pick from
  ------------------------------------------------------------------*/
  template <typename T>
  bool intersects_triangle(const Container<T>& tris,
                           const double range) const
  {
    for ( const auto& t : tris.get_items(xy_, range) )
    {
      // Ignore inactive elements
      if ( !t->is_active() ) continue;

      const Vec2d& p1 = t->v1().xy();
      const Vec2d& q1 = t->v2().xy();
      const Vec2d& r1 = t->v3().xy();

      const Vec2d& p2 = v_[0]->xy();
      const Vec2d& q2 = v_[1]->xy();
      const Vec2d& r2 = v_[2]->xy();
      const Vec2d& s2 = v_[3]->xy();

      if ( tri_quad_intersection( p1,q1,r1, p2,q2,r2,s2 ) )
        return true;
    }

    return false;

  } // Quad::intersects_tri() 

  /*------------------------------------------------------------------
  | Returns true if the quad intersects with a quad 
  | in a given Container of Quadrilaterals 
  | The factor range scales the vicinity range from which 
  | triangles to pick from
  ------------------------------------------------------------------*/
  template <typename Q>
  bool intersects_quad(const Container<Q>& quads,
                       const double range) const
  {
    for ( const auto& q : quads.get_items(xy_, range) )
    {
      // Ignore inactive elements
      if ( !q->is_active() ) continue;

      const Vec2d& p1 = q->v1().xy();
      const Vec2d& q1 = q->v2().xy();
      const Vec2d& r1 = q->v3().xy();
      const Vec2d& s1 = q->v4().xy();

      const Vec2d& p2 = v_[0]->xy();
      const Vec2d& q2 = v_[1]->xy();
      const Vec2d& r2 = v_[2]->xy();
      const Vec2d& s2 = v_[2]->xy();

      if ( quad_quad_intersection( p2,q2,r2,s2, p1,q1,r1,s1 ) )
        return true;
    }

    return false;

  } // Quad::intersects_quad() 

  /*------------------------------------------------------------------
  | Returns true if a quad edge is too close to a vertex in a given 
  | advancing front 
  | The factor range scales the vicinity range from which 
  | vertices to pick from
  | The factor min_dist_sqr defines the minimum squared 
  | distance that an advancing front vertex must be located 
  | from a quad edge
  ------------------------------------------------------------------*/
  bool intersects_front(const Front& front,
                        const double range,
                        const double min_dist_sqr) const
  {
    for (const auto& e : front.edges().get_items(xy_, range))
    {
      const Vertex& v   = e->v1();

      const Vec2d& v_xy = v.xy();
      const Vec2d& q1   = v_[0]->xy();
      const Vec2d& q2   = v_[1]->xy();
      const Vec2d& q3   = v_[2]->xy();
      const Vec2d& q4   = v_[3]->xy();

      if ( v == *v_[0] || v == *v_[1] || v == *v_[2] || v == *v_[3] )
        continue;

      if (vertex_edge_dist_sqr(v_xy, q1,q2) < min_dist_sqr ||
          vertex_edge_dist_sqr(v_xy, q2,q3) < min_dist_sqr ||
          vertex_edge_dist_sqr(v_xy, q3,q4) < min_dist_sqr ||
          vertex_edge_dist_sqr(v_xy, q4,q1) < min_dist_sqr  )
        return true;
    }

    return false;

  } // Quad::intersects_front()


  /*------------------------------------------------------------------
  | Container destructor function
  ------------------------------------------------------------------*/
  void container_destructor() override
  {
    if (v_[0]) v_[0]->remove_facet( *this );
    if (v_[1]) v_[1]->remove_facet( *this );
    if (v_[2]) v_[2]->remove_facet( *this );
    if (v_[3]) v_[3]->remove_facet( *this );

    v_[0] = nullptr;
    v_[1] = nullptr;
    v_[2] = nullptr;
    v_[3] = nullptr;
  }

  /*------------------------------------------------------------------
  | Compute quad centroid
  ------------------------------------------------------------------*/
  static inline Vec2d
  calc_centroid(const Vertex& v1, const Vertex& v2, 
                const Vertex& v3, const Vertex& v4) 
  { return 0.25 * (v1.xy() + v2.xy() + v3.xy() + v4.xy()); }

  /*------------------------------------------------------------------
  | Compute quad area
  ------------------------------------------------------------------*/
  static inline double 
  calc_area(const Vertex& v1, const Vertex& v2, 
            const Vertex& v3, const Vertex& v4)
  {
    const Vec2d& e1 = v2.xy() - v1.xy();
    const Vec2d& e2 = v3.xy() - v1.xy();
    const Vec2d& e3 = v4.xy() - v1.xy();
    const double a1 = 0.5 * cross(e1, e2); 
    const double a2 = 0.5 * cross(e2, e3); 
    return a1 + a2;
  }

  /*------------------------------------------------------------------
  | Compute quad circumcenter
  ------------------------------------------------------------------*/
  static inline Vec2d
  calc_circumcenter(const Vertex& v1, const Vertex& v2, 
                    const Vertex& v3, const Vertex& v4)
  { return 0.25 * (v1.xy() + v2.xy() + v3.xy() + v4.xy()); }

  /*------------------------------------------------------------------
  | Compute quad circumference radius
  ------------------------------------------------------------------*/
  static inline double
  calc_circumradius(const Vertex& v1, const Vertex& v2, 
                    const Vertex& v3, const Vertex& v4,
                    const Vec2d& circumcenter)
  {
    DoubleArray r = { ( circumcenter - v1.xy() ).length(),
                      ( circumcenter - v2.xy() ).length(),
                      ( circumcenter - v3.xy() ).length(),
                      ( circumcenter - v4.xy() ).length() };

    return *std::max_element( r.begin(), r.end() );
  }

  /*------------------------------------------------------------------
  | Compute quad edge lengths
  ------------------------------------------------------------------*/
  static inline double
  calc_edge_length(const Vertex& v1, const Vertex& v2)
  { return ( v2.xy() - v1.xy() ).length(); }

  static inline double 
  calc_minimum_edge_length(const DoubleArray& edge_lengths)
  { return *std::min_element(edge_lengths.begin(), edge_lengths.end()); }

  static inline double 
  calc_maximum_edge_length(const DoubleArray& edge_lengths)
  { return *std::max_element(edge_lengths.begin(), edge_lengths.end()); }

  /*------------------------------------------------------------------
  | Compute quad angles
  ------------------------------------------------------------------*/
  static inline DoubleArray
  calc_angles(const Vertex& v1, const Vertex& v2, 
              const Vertex& v3, const Vertex& v4,
              const DoubleArray& edge_lengths)
  {
    const Vec2d& p = v1.xy();
    const Vec2d& q = v2.xy();
    const Vec2d& r = v3.xy();
    const Vec2d& s = v4.xy();

    double l1 = edge_lengths[0];
    double l2 = edge_lengths[1];
    double l3 = edge_lengths[2];
    double l4 = edge_lengths[3];

    double a1 = acos( dot( q-p,  s-p ) / (l4*l1) );
    double a2 = acos( dot( p-q,  r-q ) / (l1*l2) );
    double a3 = acos( dot( s-r,  q-r ) / (l2*l3) );
    double a4 = acos( dot( r-s,  p-s ) / (l3*l4) );

    return { a1, a2, a3, a4 };
  }

  static inline double 
  calc_minimum_angle(const DoubleArray& angles)
  { return *std::min_element(angles.begin(), angles.end()); }

  static inline double 
  calc_maximum_angle(const DoubleArray& angles)
  { return *std::max_element(angles.begin(), angles.end()); }

  /*------------------------------------------------------------------
  | Compute quad shape factor
  ------------------------------------------------------------------*/
  static inline double calc_shape_factor(DoubleArray& edge_lengths,
                                         double area)
  {
    // --> WARNING: Use same formulation as for triangles
    //--------------------------------------------------
    // The norm_factor is used, in order to get:
    // Shape factor -> 1 for equilateral triangles
    // Shape factor -> 0 for bad triangles
    const double norm_factor  = 3.4641016151377544;
    const double edge_sum_sqr = edge_lengths[0] * edge_lengths[0]
                              + edge_lengths[1] * edge_lengths[1]
                              + edge_lengths[2] * edge_lengths[2]
                              + edge_lengths[3] * edge_lengths[3];

    if ( edge_sum_sqr > 0.0 )
      return norm_factor * area / edge_sum_sqr;

    return 0.0;
  } 

private:

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  VertexArray          v_               { nullptr };
  FacetArray           f_               { nullptr };


  int                  color_           {CONSTANTS.default_element_color()};
  int                  index_           {-1};
  bool                 active_          {false};
  bool                 marker_          {false};

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

  DoubleArray          edge_lengths_    {0.0};
  DoubleArray          angles_          {0.0};


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
