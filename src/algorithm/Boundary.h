/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "Vec2.h"
#include "utils.h"
#include "Geometry.h"

#include "EdgeList.h"
#include "Edge.h"
#include "Vertex.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/********************************************************************* 
* Boundary type orientation
*********************************************************************/
enum class BdryType
{
  EXTERIOR,
  INTERIOR 
}; 


/*********************************************************************
* A mesh boundary - defined by a list of edges
* > Interior boundaries are defined clockwise (CW)
* > Exterior boundaries are defined counter-clockwise (CCW)
*********************************************************************/
class Boundary : public EdgeList
{
public:

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Boundary(BdryType btype)
  : EdgeList( (btype == BdryType::EXTERIOR 
              ? Orientation::CCW : Orientation::CW) ) 
  , btype_ { btype }
  { }

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  bool is_exterior() const 
  { return (btype_ == BdryType::EXTERIOR); }
  bool is_interior() const
  { return (btype_ == BdryType::INTERIOR); }

  /*------------------------------------------------------------------
  | Override insert_edge() method of parent EdgeList, since all 
  | boundary edges must be defined with an appropriate boundary 
  | marker (integers > 0)
  ------------------------------------------------------------------*/
  Edge& insert_edge(const_iterator pos, Vertex& v1, Vertex& v2, 
                    int marker)
  { 
    ASSERT( (marker >= 0), 
            "Boundary markers must be greater than zero");

    Edge& new_edge = EdgeList::insert_edge(pos, v1, v2, marker); 

    // Initialize the boundary data structure of the new edge,
    // which keeps track of all mesh edges that will be created 
    // on it
    //new_edge.init_bdry_data();

    return new_edge;
  }

  /*------------------------------------------------------------------
  | Override add_edge() method of parent EdgeList, since all 
  | boundary edges must be defined with an appropriate boundary 
  | marker (integers > 0)
  ------------------------------------------------------------------*/
  Edge& add_edge(Vertex& v1, Vertex& v2, int marker)
  { return this->insert_edge( edges_.end(), v1, v2, marker ); }

  /*------------------------------------------------------------------
  | Set the boundary to square shape
  ------------------------------------------------------------------*/
  void set_shape_square(Vertices& vertices, int marker, 
                        const Vec2d& xy, double w,
                        double sizing=1.0, double range=1.0)
  {
    set_shape_rectangle(vertices, marker, xy, w, w, sizing, range);

  } // Boundary::set_shape_square()

  /*------------------------------------------------------------------
  | Set the boundary to equilateral triangle
  ------------------------------------------------------------------*/
  void set_shape_triangle(Vertices& vertices, int marker,
                          const Vec2d& xy, double a,
                          double sizing=1.0, double range=1.0)
  {
    set_shape_circle(vertices, marker, xy, a/sqrt(3), 3, sizing, range);

  } // Boundary::set_shape_square()

  /*------------------------------------------------------------------
  | Set the boundary to rectangular shape
  ------------------------------------------------------------------*/
  void set_shape_rectangle(Vertices& vertices, int marker,
                           const Vec2d& xy, double w, double h,
                           double sizing=1.0, double range=1.0)
  {
    const double half_w = 0.5 * w;
    const double half_h = 0.5 * h;

    std::vector<Vec2d> v_shape {
      {xy[0] - half_w, xy[1] - half_h},
      {xy[0] + half_w, xy[1] - half_h},
      {xy[0] + half_w, xy[1] + half_h},
      {xy[0] - half_w, xy[1] + half_h},
    };

    create_boundary_shape(vertices, v_shape, marker, sizing, range);

  } // Boundary::set_shape_rectangle()

  /*------------------------------------------------------------------
  | Set the boundary to circular shape
  ------------------------------------------------------------------*/
  void set_shape_circle(Vertices& vertices, int marker,
                        const Vec2d& xy, double r, size_t n=30,
                        double sizing=1.0, double range=1.0)
  {
    if ( n < 3 ) 
      return;

    std::vector<Vec2d> v_shape;

    double delta_a = 2.0 * M_PI / static_cast<double>(n);

    for ( size_t i = 0; i < n; ++i )
    {
      double ang_i = i * delta_a;
      Vec2d d_xy = { cos(ang_i), sin(ang_i) }; 

      v_shape.push_back( xy + r * d_xy );
    }

    create_boundary_shape(vertices, v_shape, marker, sizing, range);

  } // Boundary::set_shape_Circle()


private:

  /*------------------------------------------------------------------
  | Create the boundary from a given shape
  ------------------------------------------------------------------*/
  void create_boundary_shape(Vertices&           vertices, 
                             std::vector<Vec2d>& v_shape,
                             int                 marker,
                             double              sizing,
                             double              range) 
  {
    // Do nothinng if the boundary contains edges
    if ( edges_.size() > 0 )
      return;

    // Create new vertices
    std::vector<Vertex*> new_verts {};

    int N = static_cast<int>(v_shape.size());

    for (int i = 0; i < N; ++i)
    {
      Vertex& v_new = vertices.push_back( v_shape[i], sizing, range );
      new_verts.push_back( &v_new );
    }

    // CCW orientation for exteriror boundary
    if ( btype_ == BdryType::EXTERIOR )
    {
      for (int i = 0; i < N; ++i)
      {
        int j = MOD(i+1, N);
        add_edge( *new_verts[i], *new_verts[j], marker );
      }
    }
    // CW orientation for interior boundary
    else if ( btype_ == BdryType::INTERIOR )
    {
      for (int i = N-1; i >= 0; --i)
      {
        int j = MOD(i+1, N);
        add_edge( *new_verts[j], *new_verts[i], marker );
      }
    }

  };

  /*------------------------------------------------------------------
  | Boundary attributes
  ------------------------------------------------------------------*/
  BdryType   btype_;

}; // Boundary

} // namespace TQAlgorithm
} // namespace TQMesh
