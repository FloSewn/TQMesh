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
#include "Container.h"
#include "utils.h"

#include "Vertex.h"
#include "Edge.h"
#include "NullFacet.h"

namespace TQMesh {
namespace TQAlgorithm {

/*********************************************************************
* 
*********************************************************************/
class NullEdge : public Edge
{
public:
  /*------------------------------------------------------------------
  | Singleton: Get access to instance 
  ------------------------------------------------------------------*/
  static NullEdge& get_instance() 
  {
    static NullEdge null_edge; 
    return null_edge;
  }

  /*------------------------------------------------------------------
  | Singleton: Check if a given edge is the NullEdge
  ------------------------------------------------------------------*/
  static inline bool is_null(Edge* e) 
  { return e == &NullEdge::get_instance(); }

  static inline bool is_not_null(Edge* e) 
  { return e != &NullEdge::get_instance(); }

  /*------------------------------------------------------------------
  | Singleton: Delete copy & assignment operator
  ------------------------------------------------------------------*/
  NullEdge(const NullEdge&) = delete;
  NullEdge& operator=(const NullEdge&) = delete; 

  /*------------------------------------------------------------------
  | Override getter functions
  ------------------------------------------------------------------*/
  int marker() const 
  { 
    ASSERT(false, "Invalid NullEdge access: marker()");
    return -1; 
  }

  double length() const 
  { 
    ASSERT(false, "Invalid NullEdge access: length()");
    return 0.0; 
  }

  const Vec2d& normal() const 
  { 
    ASSERT(false, "Invalid NullEdge access: normal()");
    return vec_dummy_;
  }

  const Vec2d& tangent() const 
  { 
    ASSERT(false, "Invalid NullEdge access: tangent()");
    return vec_dummy_;
  }

  const Vertex& v1() const 
  { 
    ASSERT(false, "Invalid NullEdge access: v1()");
    return v1_dummy_; 
  };

  const Vertex& v2() const 
  { 
    ASSERT(false, "Invalid NullEdge access: v2()");
    return v2_dummy_; 
  };

  Vertex& v1() 
  { 
    ASSERT(false, "Invalid NullEdge access: v1()");
    return v1_dummy_; 
  };

  Vertex& v2() 
  { 
    ASSERT(false, "Invalid NullEdge access: v2()");
    return v2_dummy_; 
  };

  EdgeList& edgelist()  
  { 
    ASSERT(false, "Invalid NullEdge access: edgelist()");
    return *el_dummy_; 
  }

  const EdgeList& edgelist() const  
  { 
    ASSERT(false, "Invalid NullEdge access: edgelist()");
    return *el_dummy_; 
  }

  Edge* get_next_edge() const 
  {
    ASSERT(false, "Invalid NullEdge access: get_next_edge()");
    return &NullEdge::get_instance();
  }

  Edge* get_prev_edge() const 
  {
    ASSERT(false, "Invalid NullEdge access: get_prev_edge()");
    return &NullEdge::get_instance();
  }

private: 
  /*------------------------------------------------------------------
  | Constructor / Destructor
  ------------------------------------------------------------------*/
  NullEdge() : Edge(v1_dummy_, v2_dummy_, *el_dummy_, -1) {}

  /*------------------------------------------------------------------
  | Dummy attributes
  ------------------------------------------------------------------*/
  Vertex    v1_dummy_  {0.0, 0.0, 0.0, 0.0};
  Vertex    v2_dummy_  {0.0, 0.0, 0.0, 0.0};
  EdgeList* el_dummy_  { nullptr };
  Vec2d     vec_dummy_ {0.0, 0.0};

}; // NullFacet


} // namespace TQAlgorithm
} // namespace TQMesh
