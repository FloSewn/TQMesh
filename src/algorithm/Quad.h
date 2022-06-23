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
class Quad : public Facet
{
public:

  friend Container<Quad>;
  using ContainerIterator = Container<Quad>::List::iterator;
  using DoubleArray = std::array<double,4>;
  using FacetArray = std::array<Facet*,4>;
  using VertexArray = std::array<Vertex*,4>;

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Quad(Vertex& v1, Vertex& v2, Vertex& v3, Vertex& v4)
  : v_ {&v1, &v2, &v3, &v4}
  {
    calc_centroid();
    calc_area();
    calc_circumcenter();
    calc_edgelengths();
    calc_angles();
    calc_shape_factor();
    calc_quality();

    v_[0]->add_facet( *this );
    v_[1]->add_facet( *this );
    v_[2]->add_facet( *this );
    v_[3]->add_facet( *this );
  }

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  const ContainerIterator& pos() const { return pos_; }

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

  const Vec2d& xy() const { return xy_; }
  const Vec2d& circumcenter() const { return circ_centr_; }

  int mesh_id() const { return mesh_id_; }
  int index() const { return index_; }
  bool is_active() const { return active_; }
  bool marker() const { return marker_; }

  double area() const { return area_; }
  double circumradius() const { return circ_radius_; }
  double min_angle() const { return min_angle_; }
  double max_angle() const { return max_angle_; }

  double edgelength(unsigned int i) const { return edge_len_[i]; }
  double angle(unsigned int i) const { return angles_[i]; }

  double min_edge_length() const { return min_edge_len_; }
  double max_edge_length() const { return max_edge_len_; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  void neighbor(size_t i, Facet* f) { f_[i] = f; }
  void nbr1(Facet* f) { f_[0] = f; }
  void nbr2(Facet* f) { f_[1] = f; }
  void nbr3(Facet* f) { f_[2] = f; }
  void nbr4(Facet* f) { f_[3] = f; }

  void mesh_id(int i) { mesh_id_ = i; }
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
      DBG_MSG("  | NON-POSITIVE QUAD AREA " << area_);
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

  } //Mesh::get_vertex_index()

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
  | Mandatory container functions 
  ------------------------------------------------------------------*/
  bool in_container() const { return in_container_; }


private:

  /*------------------------------------------------------------------
  | Compute quad centroid
  ------------------------------------------------------------------*/
  void calc_centroid()
  { xy_ = 0.25 * (v_[0]->xy()+v_[1]->xy()+v_[2]->xy()+v_[3]->xy());  } 

  /*------------------------------------------------------------------
  | Compute quad area
  ------------------------------------------------------------------*/
  void calc_area()
  {
    const Vec2d& e1 = v_[1]->xy() - v_[0]->xy();
    const Vec2d& e2 = v_[2]->xy() - v_[0]->xy();
    const Vec2d& e3 = v_[3]->xy() - v_[0]->xy();
    const double a1 = 0.5 * cross(e1, e2); 
    const double a2 = 0.5 * cross(e2, e3); 
    area_ = a1 + a2;
  }

  /*------------------------------------------------------------------
  | Compute quad circumcenter
  ------------------------------------------------------------------*/
  void calc_circumcenter()
  {
    circ_centr_  = xy_;

    DoubleArray r = { ( xy_ - v_[0]->xy() ).length(),
                      ( xy_ - v_[1]->xy() ).length(),
                      ( xy_ - v_[2]->xy() ).length(),
                      ( xy_ - v_[3]->xy() ).length() };

    circ_radius_ = *std::max_element( r.begin(), r.end() );
  }

  /*------------------------------------------------------------------
  | Compute quad edge lengths
  ------------------------------------------------------------------*/
  void calc_edgelengths()
  {
    edge_len_[0] = ( v_[1]->xy() - v_[0]->xy() ).length();
    edge_len_[1] = ( v_[2]->xy() - v_[1]->xy() ).length();
    edge_len_[2] = ( v_[3]->xy() - v_[2]->xy() ).length();
    edge_len_[3] = ( v_[0]->xy() - v_[3]->xy() ).length();

    ASSERT( (edge_len_[0] > 0.0), "Invalid quad definition.");
    ASSERT( (edge_len_[1] > 0.0), "Invalid quad definition.");
    ASSERT( (edge_len_[2] > 0.0), "Invalid quad definition.");
    ASSERT( (edge_len_[3] > 0.0), "Invalid quad definition.");

    min_edge_len_ = *std::min_element(edge_len_.begin(),
                                      edge_len_.end()  );
    max_edge_len_ = *std::max_element(edge_len_.begin(),
                                      edge_len_.end()  );
  }

  /*------------------------------------------------------------------
  | Compute quad angles
  ------------------------------------------------------------------*/
  void calc_angles()
  {
    const Vec2d& p = v_[0]->xy();
    const Vec2d& q = v_[1]->xy();
    const Vec2d& r = v_[2]->xy();
    const Vec2d& s = v_[3]->xy();

    double l1 = edge_len_[0];
    double l2 = edge_len_[1];
    double l3 = edge_len_[2];
    double l4 = edge_len_[3];

    angles_[0] = acos( dot( q-p,  s-p ) / (l4*l1) );
    angles_[1] = acos( dot( p-q,  r-q ) / (l1*l2) );
    angles_[2] = acos( dot( s-r,  q-r ) / (l2*l3) );
    angles_[3] = acos( dot( r-s,  p-s ) / (l3*l4) );

    min_angle_ = *std::min_element(angles_.begin(), 
                                   angles_.end()  );

    max_angle_ = *std::max_element(angles_.begin(), 
                                   angles_.end()  );
  }

  /*------------------------------------------------------------------
  | Compute quad shape factor
  ------------------------------------------------------------------*/
  void calc_shape_factor()
  {
  } 

  /*------------------------------------------------------------------
  | Compute quad quality factor
  ------------------------------------------------------------------*/
  void calc_quality()
  { 
  } 

  /*------------------------------------------------------------------
  | Mandatory container functions 
  ------------------------------------------------------------------*/
  void container_destructor() 
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
  | 
  ------------------------------------------------------------------*/
  VertexArray          v_             { nullptr };
  FacetArray           f_             { nullptr };

  Vec2d                xy_            {0.0, 0.0};
  Vec2d                circ_centr_    {0.0,0.0};

  int                  mesh_id_       {TQMeshDefaultMeshId};
  int                  index_         {-1};
  bool                 active_        {false};
  bool                 marker_        {false};

  double               area_          {0.0};
  double               circ_radius_   {0.0};
  double               min_angle_     {0.0};
  double               max_angle_     {0.0};
  double               shape_fac_     {0.0};
  double               quality_       {0.0};
  double               min_edge_len_  {0.0};
  double               max_edge_len_  {0.0};

  DoubleArray          edge_len_      {0.0};
  DoubleArray          angles_        {0.0};

  ContainerIterator    pos_;
  bool                 in_container_;

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
