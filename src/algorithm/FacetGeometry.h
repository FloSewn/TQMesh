/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "VecND.h"
#include "Geometry.h"

#include "Vertex.h"

namespace TQMesh {

using namespace CppUtils;

/*********************************************************************
* Utility class for the calculation of triangle-related geometry
*********************************************************************/
class TriangleGeometry 
{
public:

  using VertexArray = std::array<Vertex*,3>;

  /*------------------------------------------------------------------
  | Compute triangle centroid
  ------------------------------------------------------------------*/
  static inline Vec2d
  calc_centroid(const Vertex& v1, const Vertex& v2, const Vertex& v3) 
  { return (v1.xy() + v2.xy() + v3.xy()) / 3.0; }

  /*------------------------------------------------------------------
  | Compute triangle area
  ------------------------------------------------------------------*/
  static inline double 
  calc_area(const Vertex& v1, const Vertex& v2, const Vertex& v3)
  {
    const Vec2d& e1 = v2.xy() - v1.xy();
    const Vec2d& e2 = v3.xy() - v1.xy();
    return 0.5 * cross(e1, e2); 
  }

  /*------------------------------------------------------------------
  | Compute triangle circumcenter
  ------------------------------------------------------------------*/
  static inline Vec2d
  calc_circumcenter(const Vertex& v1, const Vertex& v2, const Vertex& v3)
  {
    const Vec2d& B = v2.xy() - v1.xy();
    const Vec2d& C = v3.xy() - v1.xy();

    double D = 2.0 * cross(B, C);
    double Ux = ( C.y * (B.x*B.x + B.y*B.y) 
                - B.y * (C.x*C.x + C.y*C.y) ) / D;
    double Uy = ( B.x * (C.x*C.x + C.y*C.y) 
                - C.x * (B.x*B.x + B.y*B.y) ) / D;

    return { Ux + v1.xy().x, 
             Uy + v1.xy().y };
  }

  /*------------------------------------------------------------------
  | Compute triangle circumference radius
  ------------------------------------------------------------------*/
  static inline double
  calc_circumradius(const Vertex& v1, const Vec2d& circumcenter)
  {
    double ux = circumcenter.x - v1.xy().x;
    double uy = circumcenter.y - v1.xy().y;
    return sqrt(ux * ux + uy * uy);
  }

  /*------------------------------------------------------------------
  | Compute triangle edge lengths
  ------------------------------------------------------------------*/
  static inline double
  calc_edge_length(const Vertex& v1, const Vertex& v2) 
  { return ( v2.xy() - v1.xy() ).norm(); }

  /*------------------------------------------------------------------
  | Compute triangle angles
  ------------------------------------------------------------------*/
  static inline Vec3d 
  calc_angles(const Vertex& v1, const Vertex& v2, const Vertex& v3,
              const Vec3d& edge_lengths)
  {
    const Vec2d& p = v1.xy();
    const Vec2d& q = v2.xy();
    const Vec2d& r = v3.xy();

    double l1 = edge_lengths[0];
    double l2 = edge_lengths[1];
    double l3 = edge_lengths[2];

    double a1 = acos( dot( q-p,  r-p ) / (l1*l3) );
    double a2 = acos( dot( p-q,  r-q ) / (l1*l2) );
    double a3 = acos( dot( p-r,  q-r ) / (l2*l3) );

    return { a1, a2, a3 };
  }

  /*------------------------------------------------------------------
  | Compute triangle shape factor
  ------------------------------------------------------------------*/
  static inline double calc_shape_factor(Vec3d& edge_lengths,
                                         double area)
  {
    // The norm_factor is used, in order to get:
    // Shape factor -> 1 for equilateral triangles
    // Shape factor -> 0 for bad triangles

    const double norm_factor  = 3.4641016151377544;
    const double edge_sum_sqr = edge_lengths[0] * edge_lengths[0]
                              + edge_lengths[1] * edge_lengths[1]
                              + edge_lengths[2] * edge_lengths[2];

    if ( edge_sum_sqr > 0.0 )
      return norm_factor * area / edge_sum_sqr;

    return 0.0;
  } 

  /*------------------------------------------------------------------
  | Compute triangle quality value
  ------------------------------------------------------------------*/
  static inline double calc_quality(const Vec3d& edge_lengths, 
                                    const double shape_factor,
                                    const double mesh_size)
  {
    const double f1_1 = edge_lengths[0] / mesh_size;
    const double f1_2 = edge_lengths[1] / mesh_size;
    const double f1_3 = edge_lengths[2] / mesh_size;

    const double f2_1 = mesh_size / edge_lengths[0];
    const double f2_2 = mesh_size / edge_lengths[1];
    const double f2_3 = mesh_size / edge_lengths[2];

    const double q1 = MIN(f1_1, f2_1);
    const double q2 = MIN(f1_2, f2_2);
    const double q3 = MIN(f1_3, f2_3);
    
    return shape_factor * q1 * q2 * q3; 
  }

  /*------------------------------------------------------------------
  | Check if a triangle is valid
  ------------------------------------------------------------------*/
  static inline bool check_validity(const double area, 
                                    const Vec3d& edge_lengths)
  {
    if ( area <= 0.0 )
    {
      DEBUG_LOG("  > NON-POSITIVE TRIANGLE AREA " << area);
      return false;
    }

    if ( edge_lengths[0]<=0.0 || edge_lengths[1]<=0.0 || 
         edge_lengths[2]<=0.0  )
    {
      DEBUG_LOG("  > NON-POSITIVE TRIANGLE EDGE LENGTH ");
      return false;
    }

    return true;
  }

  /*------------------------------------------------------------------
  | Check intersection between a triangle <tri> and a vertex <v>
  ------------------------------------------------------------------*/
  template <typename T>
  static inline bool check_intersection(const T& tri,
                                        const Vertex& v)
  {
    if (  v == tri.v1() || v == tri.v2() || v == tri.v3() )
      return false;
    return in_on_triangle(v.xy(), tri.v1().xy(), tri.v2().xy(), tri.v3().xy());
  }

  /*------------------------------------------------------------------
  | Check intersection between a triangle <tri> and a <domain>
  ------------------------------------------------------------------*/
  template <typename T, typename D>
  static inline bool check_intersection(const T& tri,
                                        const D& domain)
  { 
    const Vec2d& xy = tri.xy();
    const double r  = tri.max_edge_length();

    const Vec2d& t1 = tri.v1().xy();
    const Vec2d& t2 = tri.v2().xy();
    const Vec2d& t3 = tri.v3().xy();

    for ( const auto& e_ptr : domain.get_edges(xy, r) )
    {
      const Vec2d& e1 = e_ptr->v1().xy();
      const Vec2d& e2 = e_ptr->v2().xy();

      if ( line_tri_crossing( e1,e2, t1,t2,t3 ) )
        return true;
    }

    return false;
  }

  /*------------------------------------------------------------------
  | Check intersection between a triangle <tri> and all triangles of a 
  | given container <container> that are located within <range> 
  ------------------------------------------------------------------*/
  template <typename T>
  static inline bool check_intersection(const T& tri,
                                        const Container<T>& container,
                                        const double range)
  {
    for ( const auto& t : container.get_items(tri.xy(), range) )
    {
      // Ignore same triangles
      if (t == &tri) continue;

      // Ignore inactive elements
      if ( !t->is_active() ) continue;

      const Vec2d& p1 = t->v1().xy();
      const Vec2d& q1 = t->v2().xy();
      const Vec2d& r1 = t->v3().xy();
      const Vec2d& c1 = t->xy();

      const Vec2d& p2 = tri.v1().xy();
      const Vec2d& q2 = tri.v2().xy();
      const Vec2d& r2 = tri.v3().xy();
      const Vec2d& c2 = tri.xy();

      // Check for triangle edge intersection
      if ( tri_tri_intersection( p1,q1,r1, p2,q2,r2 ) )
        return true;

      // Check if one triangle contains the other
      if ( in_triangle( c1, p2,q2,r2 ) || in_triangle( c2, p1,q1,r1 ) )
        return true;
    }

    return false;
  }

  /*------------------------------------------------------------------
  | Check intersection between a triangle <tri> and all quads of a 
  | given container <container> that are located within <range> 
  ------------------------------------------------------------------*/
  template <typename T, typename Q>
  static inline bool check_intersection(const T& tri,
                                        const Container<Q>& container,
                                        const double range)
  {
    for ( const auto& q : container.get_items(tri.xy(), range) )
    {
      // Ignore inactive elements
      if ( !q->is_active() ) continue;

      const Vec2d& p1 = q->v1().xy();
      const Vec2d& q1 = q->v2().xy();
      const Vec2d& r1 = q->v3().xy();
      const Vec2d& s1 = q->v4().xy();
      const Vec2d& c1 = q->xy();

      const Vec2d& p2 = tri.v1().xy();
      const Vec2d& q2 = tri.v2().xy();
      const Vec2d& r2 = tri.v3().xy();
      const Vec2d& c2 = tri.xy();

      // Check for edge intersection
      if ( tri_quad_intersection( p2,q2,r2, p1,q1,r1,s1 ) )
        return true;

      // Check if one element contains the other
      if ( in_triangle( c1, p2,q2,r2 ) || in_quad( c2, p1,q1,r1,s1 ) )
        return true;
    }

    return false;
  }

  /*------------------------------------------------------------------
  | Check intersection between a triangle <tri> and edges of a 
  | given advancing front structure that are located within <range> 
  ------------------------------------------------------------------*/
  template <typename T, typename F>
  static inline bool check_intersection(const T& tri, const F& front,
                                        const double range)
  {
    const Vec2d& t1 = tri.v1().xy();
    const Vec2d& t2 = tri.v2().xy();
    const Vec2d& t3 = tri.v3().xy();

    for (const auto& e : front.edges().get_items(tri.xy(), range))
    {
      const Vec2d& e1 = e->v1().xy();
      const Vec2d& e2 = e->v2().xy();

      if ( line_tri_intersection( e1,e2, t1,t2,t3 ) )
        return true;
    }

    return false;
  } 

  /*------------------------------------------------------------------
  | Check intersection between a triangle <tri> and all vertices of a 
  | given container <container> that are located within <range> 
  ------------------------------------------------------------------*/
  template <typename T>
  static inline bool check_intersection(const T& tri,
                                        const Vertices& container,
                                        const double range)
  {
    for (const auto& v : container.get_items(tri.xy(), range))
    {
      if (  v == &tri.v1() || v == &tri.v2() || v == &tri.v3() )
        continue;

      const Vec2d& t1 = tri.v1().xy();
      const Vec2d& t2 = tri.v2().xy();
      const Vec2d& t3 = tri.v3().xy();

      const Vec2d& xy = v->xy();

      if ( in_on_triangle( xy, t1,t2,t3 ) )
        return true;
    }

    return false;
  }

private:
  /*------------------------------------------------------------------
  | We hide the constructor, since this class acts only as container
  | for static inline functions
  ------------------------------------------------------------------*/
  TriangleGeometry() = default;
  ~TriangleGeometry() {};

}; // TriangleGeometry


/*********************************************************************
* Utility class for the calculation of quad-related geometry
*********************************************************************/
class QuadGeometry 
{
public:

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
    Vec4d r = { ( circumcenter - v1.xy() ).norm(),
                ( circumcenter - v2.xy() ).norm(),
                ( circumcenter - v3.xy() ).norm(),
                ( circumcenter - v4.xy() ).norm() };

    return *std::max_element( r.begin(), r.end() );
  }

  /*------------------------------------------------------------------
  | Compute quad edge lengths
  ------------------------------------------------------------------*/
  static inline double
  calc_edge_length(const Vertex& v1, const Vertex& v2)
  { return ( v2.xy() - v1.xy() ).norm(); }

  /*------------------------------------------------------------------
  | Compute quad angles
  ------------------------------------------------------------------*/
  static inline Vec4d
  calc_angles(const Vertex& v1, const Vertex& v2, 
              const Vertex& v3, const Vertex& v4,
              const Vec4d& edge_lengths)
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

  /*------------------------------------------------------------------
  | Compute quad shape factor
  ------------------------------------------------------------------*/
  static inline double calc_shape_factor(Vec4d& edge_lengths,
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

  /*------------------------------------------------------------------
  | Compute quad quality value
  ------------------------------------------------------------------*/
  static inline double calc_quality(const Vec4d& edge_lengths,
                                    const double shape_factor,
                                    const double mesh_size)
  {
    const double f1_1 = edge_lengths[0] / mesh_size;
    const double f1_2 = edge_lengths[1] / mesh_size;
    const double f1_3 = edge_lengths[2] / mesh_size;
    const double f1_4 = edge_lengths[3] / mesh_size;

    const double f2_1 = mesh_size / edge_lengths[0];
    const double f2_2 = mesh_size / edge_lengths[1];
    const double f2_3 = mesh_size / edge_lengths[2];
    const double f2_4 = mesh_size / edge_lengths[3];

    const double q1 = MIN(f1_1, f2_1);
    const double q2 = MIN(f1_2, f2_2);
    const double q3 = MIN(f1_3, f2_3);
    const double q4 = MIN(f1_4, f2_4);
    
    return shape_factor * q1 * q2 * q3 * q4; 
  }

  /*------------------------------------------------------------------
  | Check if a quad is valid
  ------------------------------------------------------------------*/
  static inline bool check_validity(const double area, 
                                    const Vec4d& edge_lengths)
  {
    if ( area <= 0.0 )
    {
      DEBUG_LOG("  | NON-POSITIVE QUAD AREA " << area);
      return false;
    }

    if ( edge_lengths[0]<=0.0 || edge_lengths[1]<=0.0 || 
         edge_lengths[2]<=0.0 || edge_lengths[3]<=0.0  )
    {
      DEBUG_LOG("  > NON-POSITIVE QUAD EDGE LENGTH ");
      return false;
    }

    return true;
  }

  /*------------------------------------------------------------------
  | Check intersection between a quad <quad> and a vertex <v>
  ------------------------------------------------------------------*/
  template <typename Q>
  static inline bool check_intersection(const Q& quad,
                                        const Vertex& v)
  {
    if (  v == quad.v1() || v == quad.v2() || 
          v == quad.v3() || v == quad.v4()  )
      return false;

    return in_on_quad(v.xy(), quad.v1().xy(), quad.v2().xy(),
                              quad.v3().xy(), quad.v4().xy());
  }

  /*------------------------------------------------------------------
  | Check intersection between a quad <quad> and a <domain>
  ------------------------------------------------------------------*/
  template <typename Q, typename D>
  static inline bool check_intersection(const Q& quad,
                                        const D& domain)
  { return !( domain.is_inside( quad.v1() ) 
           && domain.is_inside( quad.v2() ) 
           && domain.is_inside( quad.v3() ) 
           && domain.is_inside( quad.v4() ) ); }

  /*------------------------------------------------------------------
  | Check intersection between a quad <quad> and all triangles of a 
  | given container <container> that are located within <range> 
  ------------------------------------------------------------------*/
  template <typename Q, typename T>
  static inline bool check_intersection(const Q& quad,
                                        const Container<T>& container,
                                        const double range)
  {
    for ( const auto& t : container.get_items(quad.xy(), range) )
    {
      // Ignore inactive elements
      if ( !t->is_active() ) continue;

      const Vec2d& p1 = t->v1().xy();
      const Vec2d& q1 = t->v2().xy();
      const Vec2d& r1 = t->v3().xy();
      const Vec2d& c1 = t->xy();

      const Vec2d& p2 = quad.v1().xy();
      const Vec2d& q2 = quad.v2().xy();
      const Vec2d& r2 = quad.v3().xy();
      const Vec2d& s2 = quad.v4().xy();
      const Vec2d& c2 = quad.xy();

      // Check for edge intersection
      if ( tri_quad_intersection( p1,q1,r1, p2,q2,r2,s2 ) )
        return true;

      // Check if one element contains the other
      if ( in_triangle( c2, p1,q1,r1 ) || in_quad( c1, p2,q2,r2,s2 ) )
        return true;
    }

    return false;
  }

  /*------------------------------------------------------------------
  | Check intersection between a <quad> and all quads of a 
  | given container <container> that are located within <range> 
  ------------------------------------------------------------------*/
  template <typename Q>
  static inline bool check_intersection(const Q& quad,
                                        const Container<Q>& container,
                                        const double range)
  {
    for ( const auto& q : container.get_items(quad.xy(), range) )
    {
      // Ignore same quads
      if (q == &quad) continue;

      // Ignore inactive elements
      if ( !q->is_active() ) continue;

      const Vec2d& p1 = q->v1().xy();
      const Vec2d& q1 = q->v2().xy();
      const Vec2d& r1 = q->v3().xy();
      const Vec2d& s1 = q->v4().xy();
      const Vec2d& c1 = q->xy();

      const Vec2d& p2 = quad.v1().xy();
      const Vec2d& q2 = quad.v2().xy();
      const Vec2d& r2 = quad.v3().xy();
      const Vec2d& s2 = quad.v4().xy();
      const Vec2d& c2 = quad.xy();

      // Check for edge intersection
      if ( quad_quad_intersection( p2,q2,r2,s2, p1,q1,r1,s1 ) )
        return true;

      // Check if one quad contains the other
      if ( in_quad( c1, p2,q2,r2,s2 ) || in_quad( c2, p1,q1,r1,s1 ) )
        return true;
    }

    return false;
  }

  /*------------------------------------------------------------------
  | Check intersection between a <quad> and edges of a 
  | given advancing front structure that are located within <range> 
  ------------------------------------------------------------------*/
  template <typename Q, typename F>
  static inline bool check_intersection(const Q& quad, const F& front,
                                        const double range)
  {
    for (const auto& e : front.edges().get_items(quad.xy(), range))
    {
      const Vec2d& e1 = e->v1().xy();
      const Vec2d& e2 = e->v2().xy();

      const Vec2d& q1 = quad.v1().xy();
      const Vec2d& q2 = quad.v2().xy();
      const Vec2d& q3 = quad.v3().xy();
      const Vec2d& q4 = quad.v4().xy();

      if ( line_quad_intersection( e1,e2, q1,q2,q3,q4 ) )
        return true;
    }

    return false;
  }

  /*------------------------------------------------------------------
  | Check intersection between a <quad> and all vertices of a 
  | given container <container> that are located within <range> 
  ------------------------------------------------------------------*/
  template <typename Q>
  static inline bool check_intersection(const Q& quad,
                                        const Vertices& container,
                                        const double range)
  {
    for (const auto& v : container.get_items(quad.xy(), range))
    {
      if (  v == &quad.v1() || v == &quad.v2() || 
            v == &quad.v3() || v == &quad.v4()  )
        continue;

      const Vec2d& q1 = quad.v1().xy();
      const Vec2d& q2 = quad.v2().xy();
      const Vec2d& q3 = quad.v3().xy();
      const Vec2d& q4 = quad.v4().xy();

      const Vec2d& xy = v->xy();

      if ( in_on_quad( xy, q1,q2,q3,q4 ) )
        return true;
    }

    return false;
  }

  /*------------------------------------------------------------------
  | Check if a given quad is located too close to advancing
  | front edges that are located within <range> 
  ------------------------------------------------------------------*/
  template <typename Q, typename F>
  static inline bool is_too_close(const Q& quad, const F& front,
                                  const double range,
                                  const double min_dist_sqr)
  {
    for (const auto& e : front.edges().get_items(quad.xy(), range))
    {
      const Vertex& v1  = e->v1();
      const Vertex& v2  = e->v2();

      if ( v1 == quad.v1() || v1 == quad.v2() || 
           v1 == quad.v3() || v1 == quad.v4() )
        continue;

      if ( v2 == quad.v1() || v2 == quad.v2() || 
           v2 == quad.v3() || v2 == quad.v4() )
        continue;

      const Vec2d& v1_xy = v1.xy();
      const Vec2d& v2_xy = v2.xy();

      const Vec2d& q1 = quad.v1().xy();
      const Vec2d& q2 = quad.v2().xy();
      const Vec2d& q3 = quad.v3().xy();
      const Vec2d& q4 = quad.v4().xy();

      if (distance_point_edge_sqr(v1_xy, q1,q2) < min_dist_sqr ||
          distance_point_edge_sqr(v1_xy, q2,q3) < min_dist_sqr ||
          distance_point_edge_sqr(v1_xy, q3,q4) < min_dist_sqr ||
          distance_point_edge_sqr(v1_xy, q4,q1) < min_dist_sqr  )
        return true;

      if (distance_point_edge_sqr(v2_xy, q1,q2) < min_dist_sqr ||
          distance_point_edge_sqr(v2_xy, q2,q3) < min_dist_sqr ||
          distance_point_edge_sqr(v2_xy, q3,q4) < min_dist_sqr ||
          distance_point_edge_sqr(v2_xy, q4,q1) < min_dist_sqr  )
        return true;
    }

    return false;
  }

private:
  /*------------------------------------------------------------------
  | We hide the constructor, since this class acts only as container
  | for static inline functions
  ------------------------------------------------------------------*/
  QuadGeometry() = default;
  ~QuadGeometry() {};

}; // QuadGeometry

} // namespace TQMesh
