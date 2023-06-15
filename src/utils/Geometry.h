/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "MathUtility.h"
#include "VecND.h"

namespace CppUtils {

/*--------------------------------------------------------------------
| Geometrical orientation
--------------------------------------------------------------------*/
enum class Orientation
{
  CW,     // Clockwise
  CCW,    // Counter-Clockwise
  CL,     // Colinear
  NONE    // No specified orientation
};

/*--------------------------------------------------------------------
| Min() / Max() functions
--------------------------------------------------------------------*/
template <typename T>
static inline Vec2<T> bbox_min(const Vec2<T>& a, const Vec2<T>& b)
{ return Vec2<T> { MIN(a[0],b[0]), MIN(a[1], b[1]) }; }
template <typename T>
static inline Vec2<T> bbox_max(const Vec2<T>& a, const Vec2<T>& b)
{ return Vec2<T> { MAX(a[0],b[0]), MAX(a[1], b[1]) }; }

/*--------------------------------------------------------------------
| Check for orientation of three points (p, q, r)
--------------------------------------------------------------------*/
template <typename T>
static inline Orientation orientation(const Vec2<T>& p,
                                      const Vec2<T>& q,
                                      const Vec2<T>& r)
{
  T area2 = (p.x-r.x) * (q.y-r.y)
          - (q.x-r.x) * (p.y-r.y);
  
  if ( ( area2*area2 ) < CPPUTILS_SMALL )
    return Orientation::CL;

  if ( area2 > 0)
    return Orientation::CCW;
  
  return Orientation::CW;
}

/*--------------------------------------------------------------------
| Check if point r lies to the left of segment (p,q) 
--------------------------------------------------------------------*/
template <typename T>
static inline bool is_left(const Vec2<T>& p,
                           const Vec2<T>& q,
                           const Vec2<T>& r)
{
  if (orientation(p,q,r) == Orientation::CCW)
    return true;
  return false;
}

/*--------------------------------------------------------------------
| Check if point r lies to the left of segment (p,q) o
| or on the segment
--------------------------------------------------------------------*/
template <typename T>
static inline bool is_lefton(const Vec2<T>& p,
                             const Vec2<T>& q,
                             const Vec2<T>& r)
{
  if (orientation(p,q,r) == Orientation::CW)
    return false;

  return true;
}

/*--------------------------------------------------------------------
| Check if point r lies within a segment (p,q) 
--------------------------------------------------------------------*/
template <typename T>
static inline bool in_segment(const Vec2<T>& p,
                              const Vec2<T>& q,
                              const Vec2<T>& r)
{
  if (orientation(p,q,r) != Orientation::CL)
    return false;

  const Vec2<T> d_qp  = q-p;
  const Vec2<T> d_rp  = r-p;
  const T t = dot(d_rp, d_qp); 
  const T l2 = d_qp.norm_sqr();

  if ( t > 0 && t < l2)
    return true;
  
  return false;
}

/*--------------------------------------------------------------------
| Check if point r lies within a segment (p,q) or on
| its endings 
--------------------------------------------------------------------*/
template <typename T>
static inline bool in_on_segment(const Vec2<T>& p,
                                 const Vec2<T>& q,
                                 const Vec2<T>& r)
{
  if (orientation(p,q,r) != Orientation::CL)
    return false;

  const Vec2<T> d_qp  = q-p;
  const Vec2<T> d_rp  = r-p;
  const T t = dot(d_rp, d_qp); 
  const T l2 = d_qp.norm_sqr();

  if ( t >= 0 && t <= l2)
    return true;
  
  return false;
}

/*--------------------------------------------------------------------
| Check if two lines (p1,q1) and (p2,q2) intersect
| 
| * Returns true, if segments intersect at any point but
|   their edges
| * Returns true, if one line contains a part of the other
| * Returns false, if both lines share both end points
| * Returns false in all other cases
--------------------------------------------------------------------*/
template <typename T>
static inline bool line_line_intersection(const Vec2<T>& p1,
                                          const Vec2<T>& q1,
                                          const Vec2<T>& p2,
                                          const Vec2<T>& q2)
{
  Orientation o1 = orientation(p1, q1, p2);
  Orientation o2 = orientation(p1, q1, q2);
  Orientation o3 = orientation(p2, q2, p1);
  Orientation o4 = orientation(p2, q2, q1);

  if (  ( (o1 == Orientation::CCW && o2 == Orientation::CW ) ||
          (o1 == Orientation::CW  && o2 == Orientation::CCW) ) 
     && ( (o3 == Orientation::CCW && o4 == Orientation::CW ) ||
          (o3 == Orientation::CW  && o4 == Orientation::CCW) ) )
  {
    return true;
  }

  // (p1,q1) and p2 are collinear and p2 lies on segment (p1,q1)
  if ( (o1 == Orientation::CL) && in_segment(p1,q1,p2) )
    return true;

  // (p1,q1) and q2 are collinear and q2 lies on segment (p1,q1)
  if ( (o2 == Orientation::CL) && in_segment(p1,q1,q2) )
    return true;

  // (p2,q2) and p1 are collinear and p1 lies on segment (p2,q2)
  if ( (o3 == Orientation::CL) && in_segment(p2,q2,p1) )
    return true;

  // (p2,q2) and q1 are collinear and q1 lies on segment (p2,q2)
  if ( (o4 == Orientation::CL) && in_segment(p2,q2,q1) )
    return true;

  return false;
}

/*--------------------------------------------------------------------
| Check if a line segment (a,b) intersects with a 
| triangle (p,q,r) 
| 
| * Returns true, if segment intersects at any 
|   point but the triangle edges
| * Returns true, if one segment edge contains a part of 
|   the triangle
| * Returns false, if both the edge segment and the triangle
|   share both end points
|   -> proper adjacent triangles (share edge points)
| * Returns false in all other cases
--------------------------------------------------------------------*/
template <typename T>
static inline bool line_tri_intersection(const Vec2<T>& a,
                                         const Vec2<T>& b,
                                         const Vec2<T>& p,
                                         const Vec2<T>& q,
                                         const Vec2<T>& r)
{
  return ( line_line_intersection(a,b, p,q) || 
           line_line_intersection(a,b, q,r) || 
           line_line_intersection(a,b, r,p)  );
}

/*--------------------------------------------------------------------
| Check if a line segment (a,b) intersects with a 
| quad (p,q,r,s) 
| 
| * Returns true, if segment intersects at any 
|   point but the quad edges
| * Returns true, if one segment edge contains a part of 
|   the quad
| * Returns false, if both the edge segment and the quad
|   share both end points
|   -> proper adjacent quads (share edge points)
| * Returns false in all other cases
--------------------------------------------------------------------*/
template <typename T>
static inline bool line_quad_intersection(const Vec2<T>& a,
                                          const Vec2<T>& b,
                                          const Vec2<T>& p,
                                          const Vec2<T>& q,
                                          const Vec2<T>& r,
                                          const Vec2<T>& s)
{
  return ( line_line_intersection(a,b, p,q) || 
           line_line_intersection(a,b, q,r) || 
           line_line_intersection(a,b, r,s) ||
           line_line_intersection(a,b, s,p)  );
}


/*--------------------------------------------------------------------
| Check if two triangles (p1,q1,r1) and (p2,q2,r2) intersect
| 
| * Returns true, if triangle segments intersect at any 
|   point but their edges
| * Returns true, if one triangle edge contains a part of 
|   the other
| * Returns false, if both triangles share both end points
|   -> proper adjacent triangles (share edge points)
| * Returns false in all other cases
--------------------------------------------------------------------*/
template <typename T>
static inline bool tri_tri_intersection(const Vec2<T>& p1,
                                        const Vec2<T>& q1,
                                        const Vec2<T>& r1,
                                        const Vec2<T>& p2,
                                        const Vec2<T>& q2,
                                        const Vec2<T>& r2)
{

  return ( line_tri_intersection(p1,q1, p2,q2,r2) ||
           line_tri_intersection(q1,r1, p2,q2,r2) ||
           line_tri_intersection(r1,p1, p2,q2,r2)  );
}

/*--------------------------------------------------------------------
| Check if two quads (p1,q1,r1,s1) and (p2,q2,r2,s2) do 
| intersect
| 
| * Returns true, if quad segments intersect at any 
|   point but their edges
| * Returns true, if one quad edge contains a part of 
|   the other
| * Returns false, if both quads share both end points
|   -> proper adjacent quads (share edge points)
| * Returns false in all other cases
--------------------------------------------------------------------*/
template <typename T>
static inline bool quad_quad_intersection(const Vec2<T>& p1,
                                          const Vec2<T>& q1,
                                          const Vec2<T>& r1,
                                          const Vec2<T>& s1,
                                          const Vec2<T>& p2,
                                          const Vec2<T>& q2,
                                          const Vec2<T>& r2,
                                          const Vec2<T>& s2)
{
  return ( line_quad_intersection(p1,q1, p2,q2,r2,s2) ||
           line_quad_intersection(q1,r1, p2,q2,r2,s2) ||
           line_quad_intersection(r1,s1, p2,q2,r2,s2) ||
           line_quad_intersection(s1,p1, p2,q2,r2,s2)  );
}

/*--------------------------------------------------------------------
| Check if a triangle (p1,q1,r1) and a quad (p2,q2,r2,s2) do 
| intersect
| 
| * Returns true, if facet segments intersect at any 
|   point but their edges
| * Returns true, if one facet edge contains a part of 
|   the other
| * Returns false, if both facets share both end points
|   -> proper adjacent facets (share edge points)
| * Returns false in all other cases
--------------------------------------------------------------------*/
template <typename T>
static inline bool tri_quad_intersection(const Vec2<T>& p1,
                                         const Vec2<T>& q1,
                                         const Vec2<T>& r1,
                                         const Vec2<T>& p2,
                                         const Vec2<T>& q2,
                                         const Vec2<T>& r2,
                                         const Vec2<T>& s2)
{
  return ( line_quad_intersection(p1,q1, p2,q2,r2,s2) ||
           line_quad_intersection(q1,r1, p2,q2,r2,s2) ||
           line_quad_intersection(r1,p1, p2,q2,r2,s2)  );
}


/*--------------------------------------------------------------------
| Check if two rectangles a & b overlap. The rectangles
| are defined by their lower left and upper right 
| corners.
--------------------------------------------------------------------*/
template <typename T>
inline bool rect_overlap(const Vec2<T>& a_lowleft, 
                         const Vec2<T>& a_upright,
                         const Vec2<T>& b_lowleft,
                         const Vec2<T>& b_upright)
{
  return (  (a_lowleft.x <= b_upright.x)
         && (b_lowleft.x <= a_upright.x)
         && (a_lowleft.y <= b_upright.y)
         && (b_lowleft.y <= a_upright.y) );

} /* rect_overlap() */

/*--------------------------------------------------------------------
| Check if a vertex v is inside or on the edge of a 
| rectangle. The rectangle is defined by its lower 
| left and upper right corners.
--------------------------------------------------------------------*/
template <typename T>
inline bool in_on_rect(const Vec2<T>& v,
                       const Vec2<T>& lowleft,
                       const Vec2<T>& upright)
{
  return (  (v.x >= lowleft.x) && (v.y >= lowleft.y)
         && (v.x <= upright.x) && (v.y <= upright.y) );
}

/*--------------------------------------------------------------------
| Check if a vertex v is inside of a 
| rectangle. The rectangle is defined by its lower 
| left and upper right corners.
--------------------------------------------------------------------*/
template <typename T>
inline bool in_rect(const Vec2<T>& v,
                    const Vec2<T>& lowleft,
                    const Vec2<T>& upright)
{
  return (  (v.x > lowleft.x) && (v.y > lowleft.y)
         && (v.x < upright.x) && (v.y < upright.y) );
}

/*--------------------------------------------------------------------
| Check if a vertex v is inside of a triangle (p,q,r).
--------------------------------------------------------------------*/
template <typename T>
inline bool in_on_triangle(const Vec2<T>& v,
                           const Vec2<T>& p,
                           const Vec2<T>& q,
                           const Vec2<T>& r)
{
  return (  is_lefton(p,q,v) 
         && is_lefton(q,r,v) 
         && is_lefton(r,p,v) );
}

/*--------------------------------------------------------------------
| Check if a vertex v is inside of a triangle (p,q,r).
| Triangle edges are not included. 
--------------------------------------------------------------------*/
template <typename T>
inline bool in_triangle(const Vec2<T>& v,
                        const Vec2<T>& p,
                        const Vec2<T>& q,
                        const Vec2<T>& r)
{
  return (  is_left(p,q,v) 
         && is_left(q,r,v) 
         && is_left(r,p,v) );
}

/*--------------------------------------------------------------------
| Check if a vertex v is inside of a quad (p,q,r,s).
--------------------------------------------------------------------*/
template <typename T>
inline bool in_on_quad(const Vec2<T>& v, 
                       const Vec2<T>& p, const Vec2<T>& q,
                       const Vec2<T>& r, const Vec2<T>& s)
{
  return (  is_lefton(p,q,v) 
         && is_lefton(q,r,v) 
         && is_lefton(r,s,v) 
         && is_lefton(s,p,v) );
}

/*--------------------------------------------------------------------
| Check if a vertex v is inside of a quad (p,q,r,s).
| Quad edges are not included
--------------------------------------------------------------------*/
template <typename T>
inline bool in_quad(const Vec2<T>& v, 
                    const Vec2<T>& p, const Vec2<T>& q,
                    const Vec2<T>& r, const Vec2<T>& s)
{
  return (  is_left(p,q,v) 
         && is_left(q,r,v) 
         && is_left(r,s,v) 
         && is_left(s,p,v) );
}

/*--------------------------------------------------------------------
| This function returns the squared normal distance between a point p
| and a line segment defined by its vertices (e1,e2).
| 
| Reference:
| https://stackoverflow.com/questions/849211/shortest-\
| distance-between-a-point-and-a-line-segment
--------------------------------------------------------------------*/
template <typename CoordType, std::size_t Dim>
inline double distance_point_edge_sqr(const VecND<CoordType,Dim>& p,
                                      const VecND<CoordType,Dim>& e1,
                                      const VecND<CoordType,Dim>& e2)
{
  const VecND<CoordType,Dim> delta_edge = e2 - e1;
  const VecND<CoordType,Dim> delta_p = p - e1;

  const CoordType zero {};
  const CoordType one {static_cast<CoordType>(1)};

  // Project p onto segment e1->e2
  const CoordType dot_p = dot(delta_p, delta_edge) / delta_edge.norm_sqr();

  // Clip projection onto the segment 
  const CoordType t = std::clamp(dot_p, zero, one); 

  // Projected p 
  const VecND<CoordType,Dim> proj_p = e1 + t * delta_edge;

  // Squared distance between projected and actual vertex
  return (p - proj_p).norm_sqr();
}



/*--------------------------------------------------------------------
| Check if two segments (which are defined by their ending vertices 
| (v,w) overlap
--------------------------------------------------------------------*/
template <typename T>
inline bool segment_overlap(const Vec2<T>& v1, const Vec2<T> w1,
                            const Vec2<T>& v2, const Vec2<T> w2)
{
  if ( ( v1 == v2 ) && ( w1 == w2 ) ) 
    return true;

  if ( ( v1 == w2 ) && ( w1 == v2 ) ) 
    return true;

  return false;
}

} // namespace CppUtils
