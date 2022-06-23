/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "Vec2.h"
#include "Geometry.h"
#include "Container.h"
#include "utils.h"

namespace TQMesh {
namespace TQAlgorithm {

class Vertex;

/*********************************************************************
* This class defines any two dimensional facet
*********************************************************************/
class Facet
{
  public:

  /*------------------------------------------------------------------
  | Constructor 
  ------------------------------------------------------------------*/
  Facet() = default;

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  virtual double area() const = 0;
  virtual double min_angle() const = 0;
  virtual double max_angle() const = 0;

  virtual size_t n_vertices() const { return 0; }
  virtual const Vec2d& xy() const = 0; 
  virtual const Vertex& vertex(size_t i) const = 0;
  virtual Vertex& vertex(size_t i) = 0;
  virtual int index() const = 0;
  virtual bool is_active() const = 0;
  virtual bool marker() const = 0;
  virtual int mesh_id() const = 0;

  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  virtual void mesh_id(int i) = 0;
  virtual void index(int i) = 0;
  virtual void is_active(bool a) = 0;
  virtual void marker(bool c) = 0;

  /*------------------------------------------------------------------
  | Functions to return indices of vertices, edges...
  ------------------------------------------------------------------*/
  virtual int get_vertex_index(const Vertex& v1) const = 0;
  virtual 
  int get_edge_index(const Vertex& v1, const Vertex& v2) const = 0;


  /*------------------------------------------------------------------
  | Set Facet neighbor 
  ------------------------------------------------------------------*/
  virtual void neighbor(size_t i, Facet* f) = 0;

  /*------------------------------------------------------------------
  | Facet intersection
  ------------------------------------------------------------------*/
  virtual bool intersects_vertex(const Vertex& v) const
  { ASSERT(0,"THIS SHOULD NOT BE CALLED"); return false; }


}; // Facet


} // namespace TQAlgorithm
} // namespace TQMesh
