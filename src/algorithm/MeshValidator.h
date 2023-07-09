/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "Vertex.h"
#include "Edge.h"
#include "Facet.h"
#include "Boundary.h"
#include "Front.h"
#include "Domain.h"
#include "Mesh.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* This class manages the removal of mesh entities
*********************************************************************/
class MeshValidator 
{
public:
  /*------------------------------------------------------------------
  | Constructor / Destructor 
  ------------------------------------------------------------------*/
  MeshValidator(Mesh& mesh, const Domain& domain, Front& front)
  : mesh_ { mesh }, domain_ { domain }, front_ { front } {}

  ~MeshValidator() {}

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  double min_cell_quality() const { return min_cell_quality_; }
  double max_cell_angle() const { return max_cell_angle_; }

  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  void min_cell_quality(double v) { min_cell_quality_ = v; }
  void max_cell_angle(double v) { max_cell_angle_ = v; }

  /*------------------------------------------------------------------
  | Check if mesh entities are invalid. If yes, remove them from
  | the mesh and return true. Else, simply return false.
  ------------------------------------------------------------------*/
  bool remove_from_mesh_if_invalid(Vertex& v)
  {
    if ( vertex_is_valid(v) )
      return false;
    
    mesh_.remove_vertex(v);
    return true;
  } 

  bool remove_from_mesh_if_invalid(Triangle& t)
  {
    if ( triangle_is_valid(t) )
      return false;

    mesh_.remove_triangle(t);
    return true;
  } 

  bool remove_from_mesh_if_invalid(Vertex& v, Triangle& t)
  {
    if ( vertex_is_valid(v) && triangle_is_valid(t) )
      return false;

    mesh_.remove_triangle(t);
    mesh_.remove_vertex(v);
    return true;
  } 

  bool remove_from_mesh_if_invalid(Vertex& v, Triangle& t1, Triangle& t2)
  {
    if ( vertex_is_valid(v) && triangle_is_valid(t1) && triangle_is_valid(t2) )
      return false;

    mesh_.remove_triangle(t1);
    mesh_.remove_triangle(t2);
    mesh_.remove_vertex(v);
    return true;
  }

private:

  /*------------------------------------------------------------------
  | Check if a triangle is valid. If yes, return true - 
  | else return false.
  ------------------------------------------------------------------*/
  bool triangle_is_valid(const Triangle& tri)
  {
    Vertices&   vertices = mesh_.vertices();
    Triangles& triangles = mesh_.triangles();
    Quads&         quads = mesh_.quads();

    const double rho   = domain_.size_function( tri.xy() );
    const double range = 2.0 * rho;

    DEBUG_LOG("CHECK NEW TRIANGLE: " << tri);

    if ( !tri.is_valid() )
      return false;

    if ( tri.intersects_front( front_, range ) )
    { DEBUG_LOG("  > FRONT INTERSECTION"); return false; }

    if ( tri.intersects_vertex( vertices, range ) )
    { DEBUG_LOG("  > VERTEX INTERSECTION"); return false; }

    if ( tri.intersects_triangle( triangles, range ) )
    { DEBUG_LOG("  > TRIANGLE INTERSECTION"); return false; }

    if ( tri.intersects_quad( quads, range ) )
    { DEBUG_LOG("  > QUAD INTERSECTION"); return false; }

    if ( tri.quality(rho) < min_cell_quality_ )
    { DEBUG_LOG("  > BAD TRIANGLE QUALITY"); return false; }

    if ( tri.max_angle() > max_cell_angle_ )
    { DEBUG_LOG("  > BAD MAXIMUM ANGLE"); return false; }

    DEBUG_LOG("  > VALID");
    return true;

  } // Mesh::triangle_is_valid()

  /*------------------------------------------------------------------
  | Check if a vertex is valid. If yes, return true - 
  | else return false.
  ------------------------------------------------------------------*/
  bool vertex_is_valid(const Vertex& v)
  {
    Triangles& triangles = mesh_.triangles();
    Quads&         quads = mesh_.quads();

    const double rho   = domain_.size_function( v.xy() );
    const double range = 2.0 * rho;

    DEBUG_LOG("CHECK NEW VERTEX: " << v);

    if ( !domain_.is_inside( v ) )
    { DEBUG_LOG("  > OUTSIDE DOMAIN"); return false; }

    if ( v.intersects_facet(triangles, range) )
    { DEBUG_LOG("  > TRIANGLE INTERSECTION"); return false; }

    if ( v.intersects_facet(quads, range) )
    { DEBUG_LOG("  > QUAD INTERSECTION"); return false; }

    DEBUG_LOG("  > VALID");
    return true;

  } // Mesh::vertex_is_valid()


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Mesh&           mesh_;
  const Domain&   domain_;
  Front&          front_;

  double          min_cell_quality_   = 0.0;
  double          max_cell_angle_     = M_PI;

}; // MeshValidator

} // namespace TQAlgorithm
} // namespace TQMesh
