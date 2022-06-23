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
class Triangle : public Facet
{
public:

  friend Container<Triangle>;
  using ContainerIterator = Container<Triangle>::List::iterator;
  using DoubleArray = std::array<double,3>;
  using FacetArray = std::array<Facet*,3>;
  using VertexArray = std::array<Vertex*,3>;

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Triangle(Vertex& v1, Vertex& v2, Vertex& v3)
  : v_ {&v1, &v2, &v3}
  {
    calc_centroid();
    calc_area();
    calc_circumcenter();
    calc_edgelengths();
    calc_angles();
    calc_shape_factor();
    //calc_quality();

    v_[0]->add_facet( *this );
    v_[1]->add_facet( *this );
    v_[2]->add_facet( *this );
  }

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  const ContainerIterator& pos() const { return pos_; }

  const Vertex& vertex(size_t i) const { return *v_[i]; }
  const Vertex& v1() const { return *v_[0]; }
  const Vertex& v2() const { return *v_[1]; }
  const Vertex& v3() const { return *v_[2]; }
  size_t n_vertices() const { return 3; }

  Vertex& vertex(size_t i) { return *v_[i]; }
  Vertex& v1() { return *v_[0]; }
  Vertex& v2() { return *v_[1]; }
  Vertex& v3() { return *v_[2]; }

  const Facet* neighbor(size_t i) const { return f_[i]; }
  const Facet* nbr1() const { return f_[0]; }
  const Facet* nbr2() const { return f_[1]; }
  const Facet* nbr3() const { return f_[2]; }

  Facet* neighbor(size_t i) { return f_[i]; }
  Facet* nbr1() { return f_[0]; }
  Facet* nbr2() { return f_[1]; }
  Facet* nbr3() { return f_[2]; }

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

  double quality(const double h) const 
  { 
    const double e1 = edge_len_[0];
    const double e2 = edge_len_[1];
    const double e3 = edge_len_[2];

    const double f1_1 = e1 / h;
    const double f1_2 = e2 / h;
    const double f1_3 = e3 / h;

    const double f2_1 = h / e1;
    const double f2_2 = h / e2;
    const double f2_3 = h / e3;

    const double q1 = MIN(f1_1, f2_1);
    const double q2 = MIN(f1_2, f2_2);
    const double q3 = MIN(f1_3, f2_3);
    
    return shape_fac_ * q1 * q2 * q3; 
  }

  double edgelength(unsigned int i) const { return edge_len_[i]; }
  double angle(unsigned int i) const { return angles_[i]; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  void neighbor(size_t i, Facet* f) { f_[i] = f; }
  void nbr1(Facet* f) { f_[0] = f; }
  void nbr2(Facet* f) { f_[1] = f; }
  void nbr3(Facet* f) { f_[2] = f; }

  void mesh_id(int i) { mesh_id_ = i; }
  void index(int i) { index_ = i; }
  void is_active(bool a) { active_ = a; }
  void marker(bool c){ marker_ = c; }

  /*------------------------------------------------------------------
  | Returns true if the triangle is valid
  ------------------------------------------------------------------*/
  bool is_valid() const
  {
    if ( area_ <= 0.0 )
    {
      DBG_MSG("  > NON-POSITIVE TRIANGLE AREA " << area_);
      return false;
    }

    if ( edge_len_[0]<=0.0 || edge_len_[1]<=0.0 || edge_len_[2]<=0.0 )
    {
      DBG_MSG("  > NON-POSITIVE TRIANGLE EDGE LENGTH ");
      return false;
    }

    return true;
  }

  /*------------------------------------------------------------------
  | Returns the index of a triangle vertex for a given input vertex
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

    return -1;

  } //Mesh::get_vertex_index()

  /*------------------------------------------------------------------
  | Returns the index of a triangle edge for two given input vertices
  | Returns -1 if no edge is found
  |
  |                   v_[2]
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
  |        v_[0]     e2       v_[1]
  |
  ------------------------------------------------------------------*/
  int get_edge_index(const Vertex& v1, const Vertex& v2) const
  {
    if ( (&v1==v_[0] && &v2==v_[1]) || (&v1==v_[1] && &v2==v_[0]) )
      return 2;

    if ( (&v1==v_[1] && &v2==v_[2]) || (&v1==v_[2] && &v2==v_[1]) )
      return 0;

    if ( (&v1==v_[2] && &v2==v_[0]) || (&v1==v_[0] && &v2==v_[2]) )
      return 1;

    return -1;

  } // get_edge_index()


  /*------------------------------------------------------------------
  | Returns true if the triangle intersects with a vertex 
  | --> Vertex is located within the triangle or on its
  |     edges
  ------------------------------------------------------------------*/
  bool intersects_vertex(const Vertex& v) const
  {
    if (  v == *v_[0] || v == *v_[1] || v == *v_[2] )
      return false;

    return in_on_triangle(v.xy(), 
        v_[0]->xy(), v_[1]->xy(), v_[2]->xy());

  } // Triangle::intersects_vertex() 

  /*------------------------------------------------------------------
  | Returns true, if any triangle vertex is not within a given domain.
  | Triangle vertices are allowed to be located on domain edges.
  | Domain vertices are allowed to be placed on triangle edges.
  ------------------------------------------------------------------*/
  bool intersects_domain(const Domain& domain) const 
  {
    return !( domain.is_inside( *v_[0] ) 
           && domain.is_inside( *v_[1] ) 
           && domain.is_inside( *v_[2] ) );

  } // Triangle::intersects_domain()

  /*------------------------------------------------------------------
  | Returns true if the triangle intersects with a triangle 
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
      // Ignore same triangles
      if (t == this) continue;

      // Ignore inactive elements
      if ( !t->is_active() ) continue;

      const Vec2d& p1 = t->v1().xy();
      const Vec2d& q1 = t->v2().xy();
      const Vec2d& r1 = t->v3().xy();
      const Vec2d& c1 = t->xy();

      const Vec2d& p2 = v_[0]->xy();
      const Vec2d& q2 = v_[1]->xy();
      const Vec2d& r2 = v_[2]->xy();
      const Vec2d& c2 = this->xy();

      // Check for triangle edge intersection
      if ( tri_tri_intersection( p1,q1,r1, p2,q2,r2 ) )
        return true;

      // Check if one triangle contains the other
      if (  in_triangle( c1, p2,q2,r2 ) 
         || in_triangle( c2, p1,q1,r1 ) )
        return true;
    }

    return false;

  } // Triangle::intersects_tri() 

  /*------------------------------------------------------------------
  | Returns true if the triangle intersects with a quad 
  | in a given Container of Quadrilaterals 
  | The factor range scales the vicinity range from which 
  | triangles to pick from
  ------------------------------------------------------------------*/
  template <typename T>
  bool intersects_quad(const Container<T>& quads,
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

      if ( tri_quad_intersection( p2,q2,r2, p1,q1,r1,s1 ) )
        return true;
    }

    return false;

  } // Triangle::intersects_quad() 

  /*------------------------------------------------------------------
  | Returns true if a triangle edge intersects with an edge of 
  | the advancing front.  
  | The factor range scales the vicinity range from which 
  | front edges to pick from
  ------------------------------------------------------------------*/
  bool intersects_front(const Front& front,
                        const double range) const 
  {
    for (const auto& e : front.edges().get_items(xy_, range))
    {
      const Vec2d& e1 = e->v1().xy();
      const Vec2d& e2 = e->v2().xy();

      const Vec2d& t1   = v_[0]->xy();
      const Vec2d& t2   = v_[1]->xy();
      const Vec2d& t3   = v_[2]->xy();

      if ( line_tri_intersection( e1,e2, t1,t2,t3 ) )
        return true;
    }

    return false;

  } // Triangle::intersects_front()

  /*------------------------------------------------------------------
  | Returns true if the triangle encloses an advancing front vertex.
  | The factor range scales the vicinity range from which 
  | front edges to pick from
  ------------------------------------------------------------------*/
  bool intersects_vertex(const Vertices& verts,
                         const double range) const 
  {
    for (const auto& v : verts.get_items(xy_, range))
    {
      if (  v == v_[0] || v == v_[1] || v == v_[2] )
        continue;

      const Vec2d& t1   = v_[0]->xy();
      const Vec2d& t2   = v_[1]->xy();
      const Vec2d& t3   = v_[2]->xy();

      const Vec2d& xy   = v->xy();

      if ( in_on_triangle( xy, t1,t2,t3 ) )
        return true;
    }

    return false;

  } // intersects_vertex()

  /*------------------------------------------------------------------
  | Mandatory container functions 
  ------------------------------------------------------------------*/
  bool in_container() const { return in_container_; }

private:

  /*------------------------------------------------------------------
  | Compute triangle centroid
  ------------------------------------------------------------------*/
  void calc_centroid()
  { 
    xy_ = ( v_[0]->xy() + v_[1]->xy() + v_[2]->xy() ) / 3.0; 
  } 

  /*------------------------------------------------------------------
  | Compute triangle area
  ------------------------------------------------------------------*/
  void calc_area()
  {
    const Vec2d& e1 = v_[1]->xy() - v_[0]->xy();
    const Vec2d& e2 = v_[2]->xy() - v_[0]->xy();
    area_ = 0.5 * cross(e1, e2); 
  }

  /*------------------------------------------------------------------
  | Compute triangle circumcenter
  ------------------------------------------------------------------*/
  void calc_circumcenter()
  {
    const Vec2d& B = v_[1]->xy() - v_[0]->xy();
    const Vec2d& C = v_[2]->xy() - v_[0]->xy();
    double D = 2.0 * cross(B, C);
    double Ux = ( C.y * (B.x*B.x + B.y*B.y) 
                - B.y * (C.x*C.x + C.y*C.y) ) / D;
    double Uy = ( B.x * (C.x*C.x + C.y*C.y) 
                - C.x * (B.x*B.x + B.y*B.y) ) / D;

    circ_centr_.x = Ux + v_[0]->xy().x;
    circ_centr_.y = Uy + v_[0]->xy().y;

    circ_radius_ = sqrt( Ux * Ux + Uy * Uy );
  }

  /*------------------------------------------------------------------
  | Compute triangle edge lengths
  ------------------------------------------------------------------*/
  void calc_edgelengths()
  {
    edge_len_[0] = ( v_[1]->xy() - v_[0]->xy() ).length();
    edge_len_[1] = ( v_[2]->xy() - v_[1]->xy() ).length();
    edge_len_[2] = ( v_[0]->xy() - v_[2]->xy() ).length();

    min_edge_len_ = *std::min_element(edge_len_.begin(),
                                      edge_len_.end()  );
    max_edge_len_ = *std::max_element(edge_len_.begin(),
                                      edge_len_.end()  );
  }

  /*------------------------------------------------------------------
  | Compute triangle angles
  ------------------------------------------------------------------*/
  void calc_angles()
  {
    const Vec2d& p = v_[0]->xy();
    const Vec2d& q = v_[1]->xy();
    const Vec2d& r = v_[2]->xy();

    double l1 = edge_len_[0];
    double l2 = edge_len_[1];
    double l3 = edge_len_[2];

    double a1 = acos( dot( q-p,  r-p ) / (l1*l3) );
    double a2 = acos( dot( p-q,  r-q ) / (l1*l2) );
    double a3 = acos( dot( p-r,  q-r ) / (l2*l3) );

    angles_[0] = a1;
    angles_[1] = a2;
    angles_[2] = a3;

    min_angle_ = *std::min_element(angles_.begin(), 
                                   angles_.end()  );

    max_angle_ = *std::max_element(angles_.begin(), 
                                   angles_.end()  );
  }

  /*------------------------------------------------------------------
  | Compute triangle shape factor
  ------------------------------------------------------------------*/
  void calc_shape_factor()
  {
    // The norm_factor is used, in order to get:
    // Shape factor -> 1 for equilateral triangles
    // Shape factor -> 0 for bad triangles

    const double norm_factor  = 3.4641016151377544;
    const double edge_sum_sqr = edge_len_[0] * edge_len_[0]
                              + edge_len_[1] * edge_len_[1]
                              + edge_len_[2] * edge_len_[2];

    if ( edge_sum_sqr > 0.0 )
      shape_fac_ = norm_factor * area_ / edge_sum_sqr;
    else
      shape_fac_ = 0.0;

  } // Triangle::calc_shape_factor()

  /*------------------------------------------------------------------
  | Mandatory container functions 
  ------------------------------------------------------------------*/
  void container_destructor() 
  {
    if (v_[0]) v_[0]->remove_facet( *this );
    if (v_[1]) v_[1]->remove_facet( *this );
    if (v_[2]) v_[2]->remove_facet( *this );

    v_[0] = nullptr;
    v_[1] = nullptr;
    v_[2] = nullptr;
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  VertexArray          v_ { nullptr };
  FacetArray           f_ { nullptr };


  Vec2d                xy_           {0.0, 0.0};
  Vec2d                circ_centr_   {0.0,0.0};

  int                  mesh_id_      {TQMeshDefaultMeshId};
  int                  index_        {-1};
  bool                 active_       {false};
  bool                 marker_       {false};

  double               area_         {0.0};
  double               circ_radius_  {0.0};
  double               min_angle_    {0.0};
  double               max_angle_    {0.0};
  double               shape_fac_    {0.0};
  double               quality_      {0.0};
  double               min_edge_len_ {0.0};
  double               max_edge_len_ {0.0};
 
  DoubleArray          edge_len_     {0.0};
  DoubleArray          angles_       {0.0};

  ContainerIterator    pos_;
  bool                 in_container_;

}; // Triangle

/*********************************************************************
* Define general triangle container declaration
*********************************************************************/
using Triangles = Container<Triangle>;

/*********************************************************************
* Triangle ostream overload 
*********************************************************************/
static std::ostream& operator<<(std::ostream& os, 
                                const Triangle& t)
{
  return os << t.v1() << " -> " << t.v2() << " -> "
            << t.v3();
}


/*********************************************************************
* Edge equality operator 
*********************************************************************/
static bool operator==(const Triangle& t1, const Triangle& t2)
{ return (&t1 == &t2); }
static bool operator!=(const Triangle& t1, const Triangle& t2)
{ return !(t1 == t2); }

} // namespace TQAlgorithm
} // namespace TQMesh
