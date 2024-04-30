/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "TQMesh.h"

namespace TQMesh {

class Mesh;

/*********************************************************************
* This class defines any two dimensional facet
*********************************************************************/
class Facet
{
  public:

  /*------------------------------------------------------------------
  | Constructor / Destructor
  ------------------------------------------------------------------*/
  Facet() = default;
  virtual ~Facet() {};

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  virtual const Vertex& vertex(size_t i) const = 0;
  virtual Vertex&       vertex(size_t i) = 0;
  virtual size_t        n_vertices() const { return 0; }
  virtual const Vec2d&  xy() const = 0; 
  virtual Mesh*         mesh() const = 0;
  virtual int           color() const = 0;
  virtual int           index() const = 0;
  // virtual bool          is_active() const = 0;
  // virtual bool          marker() const = 0;
  virtual double        area() const = 0;
  virtual double        min_angle() const = 0;
  virtual double        max_angle() const = 0;
  virtual double        min_edge_length() const = 0;
  virtual double        max_edge_length() const = 0;


  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  virtual void color(int i) = 0;
  virtual void index(int i) = 0;
  // virtual void is_active(bool a) = 0;
  // virtual void marker(bool c) = 0;
  virtual void mesh(Mesh* m) = 0;

  /*------------------------------------------------------------------
  | Functions to return indices of vertices, edges...
  ------------------------------------------------------------------*/
  virtual int get_vertex_index(const Vertex& v1) const = 0;
  virtual int get_edge_index(const Vertex& v1, const Vertex& v2) const = 0;

  /*------------------------------------------------------------------
  | Set Facet neighbor 
  ------------------------------------------------------------------*/
  virtual void neighbor(size_t i, Facet* f) = 0;

  /*------------------------------------------------------------------
  | Facet intersection
  ------------------------------------------------------------------*/
  virtual bool intersects_vertex(const Vertex& v) const
  { ASSERT(0,"THIS SHOULD NOT BE CALLED"); return false; }

  /*------------------------------------------------------------------
  | Update the facet metrics if its vertices changed
  ------------------------------------------------------------------*/
  virtual void update_metrics(bool update_centroid=true) = 0;

}; // Facet


/*********************************************************************
*  
*********************************************************************/
class NullFacet : public Facet
{
public:
  /*------------------------------------------------------------------
  | Singleton: Get access to instance 
  ------------------------------------------------------------------*/
  static NullFacet& get_instance() 
  {
    static NullFacet null_facet; 
    return null_facet;
  }

  /*------------------------------------------------------------------
  | Singleton: Check if a given facet is the NullFacet
  ------------------------------------------------------------------*/
  static inline bool is_null(Facet* f) 
  { return f == &NullFacet::get_instance(); }

  static inline bool is_not_null(Facet* f) 
  { return f != &NullFacet::get_instance(); }

  /*------------------------------------------------------------------
  | Singleton: Delete copy & assignment operator
  ------------------------------------------------------------------*/
  NullFacet(const NullFacet&) = delete;
  NullFacet& operator=(const NullFacet&) = delete; 

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  const Vertex& vertex(size_t i) const override
  { 
    ASSERT( false, "Invalid NullFacet access: vertex().");
    return dummy_; 
  }

  Vertex& vertex(size_t i) override
  { 
    ASSERT( false, "Invalid NullFacet access: vertex().");
    return dummy_; 
  }

  const Vec2d& xy() const override 
  { 
    ASSERT( false, "Invalid NullFacet access: xy()."); 
    return dummy_.xy(); 
  }

  size_t n_vertices() const override { return 0; }
  int    color() const override { return -1; }
  int    index() const override { return -1; }
  double area() const override { return 0; }
  double min_angle() const override { return 0; }
  double max_angle() const override { return 0; }
  double min_edge_length() const override { return 0; }
  double max_edge_length() const override { return 0; }
  Mesh*  mesh() const override { return nullptr; }

  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  void color(int i) override
  { ASSERT( false, "Invalid NullFacet access: color()."); }

  void index(int i) override
  { ASSERT( false, "Invalid NullFacet access: index()."); }

  void mesh(Mesh* m) override
  { ASSERT( false, "Invalid NullFacet access: mesh()."); }


  /*------------------------------------------------------------------
  | Functions to return indices of vertices, edges...
  ------------------------------------------------------------------*/
  int get_vertex_index(const Vertex& v) const override
  {
    ASSERT( false, "Invalid NullFacet access: get_vertex_index().");
    return -1;
  } 

  int get_edge_index(const Vertex& v1, const Vertex& v2) const override
  {
    ASSERT( false, "Invalid NullFacet access: get_edge_index().");
    return -1;
  }

  /*------------------------------------------------------------------
  | Set Facet neighbor 
  ------------------------------------------------------------------*/
  void neighbor(size_t i, Facet* f) override 
  { ASSERT( false, "Invalid NullFacet access: neighbor()."); }

  /*------------------------------------------------------------------
  | Facet intersection
  ------------------------------------------------------------------*/
  bool intersects_vertex(const Vertex& v) const override
  { 
    ASSERT( false, "Invalid NullFacet access: intersects_vertex()."); 
    return false;
  }

  /*------------------------------------------------------------------
  | Facet metric update
  ------------------------------------------------------------------*/
  void update_metrics(bool update_centroid=true) override {}

private: 
  /*------------------------------------------------------------------
  | Constructor / Destructor
  ------------------------------------------------------------------*/
  NullFacet() = default;
  ~NullFacet() {}

  /*------------------------------------------------------------------
  | Dummy attributes
  ------------------------------------------------------------------*/
  Vertex dummy_ {0.0, 0.0, 0.0, 0.0};

}; // NullFacet

} // namespace TQMesh
