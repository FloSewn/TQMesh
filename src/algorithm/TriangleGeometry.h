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
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* Utility class for the calculation of triangle-related geometry
*********************************************************************/
class TriangleGeometry 
{
public:

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

  static inline double 
  calc_minimum_edge_length(const Vec3d& edge_lengths)
  { return *std::min_element(edge_lengths.begin(), edge_lengths.end()); }

  static inline double 
  calc_maximum_edge_length(const Vec3d& edge_lengths)
  { return *std::max_element(edge_lengths.begin(), edge_lengths.end()); }

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

  static inline double 
  calc_minimum_angle(const Vec3d& angles)
  { return *std::min_element(angles.begin(), angles.end()); }

  static inline double 
  calc_maximum_angle(const Vec3d& angles)
  { return *std::max_element(angles.begin(), angles.end()); }

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

}; // TriangleGeometry

} // namespace TQAlgorithm
} // namespace TQMesh
