/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <algorithm>
#include <limits.h>

#include "VecND.h"
#include "VtkIO.h"

#include "utils.h"
#include "Vertex.h"
#include "Triangle.h"
#include "Quad.h"
#include "NullFacet.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* The mesh 
*********************************************************************/
class Mesh
{
public:

  using EdgeVector     = std::vector<Edge*>;
  using VertexVector   = std::vector<Vertex*>;
  using TriVector      = std::vector<Triangle*>;
  using QuadVector     = std::vector<Quad*>;

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Mesh(int       mesh_id=DEFAULT_MESH_ID,
       int       element_color=DEFAULT_ELEMENT_COLOR,
       double    qtree_scale=ContainerQuadTreeScale,
       size_t    qtree_items=ContainerQuadTreeItems, 
       size_t    qtree_depth=ContainerQuadTreeDepth)
  : mesh_id_    { ABS(mesh_id) }
  , elem_color_ { ABS(element_color) }
  , verts_      { qtree_scale, qtree_items, qtree_depth }
  , quads_      { qtree_scale, qtree_items, qtree_depth }
  , tris_       { qtree_scale, qtree_items, qtree_depth }
  { }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  friend std::ostream& operator<<(std::ostream& os, const Mesh& mesh);

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  void add_area(double a) { mesh_area_ += a; }

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  double area()             const { return mesh_area_; }
  int    id()               const { return mesh_id_; }
  int    element_color()    const { return elem_color_; }
  bool   completed()        const { return mesh_completed_; }
  size_t n_vertices()       const { return verts_.size(); }
  size_t n_quads()          const { return quads_.size(); }
  size_t n_triangles()      const { return tris_.size(); }
  size_t n_elements()       const { return quads_.size() + tris_.size(); }
  size_t n_interior_edges() const { return intr_edges_.size(); }
  size_t n_boundary_edges() const { return bdry_edges_.size(); }
  size_t n_edges()          const { return intr_edges_.size() + bdry_edges_.size(); }

  const Vertices& vertices() const { return verts_; }
  Vertices& vertices() { return verts_; }

  const Quads& quads() const { return quads_; }
  Quads& quads() { return quads_; }

  const Triangles& triangles() const { return tris_; }
  Triangles& triangles() { return tris_; }

  const EdgeList& interior_edges() const { return intr_edges_; }
  EdgeList& interior_edges() { return intr_edges_; }

  const EdgeList& boundary_edges() const { return bdry_edges_; }
  EdgeList& boundary_edges() { return bdry_edges_; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  void element_color(int c) { elem_color_ = c; }

  /*------------------------------------------------------------------
  | This function takes care of the garbage collection 
  ------------------------------------------------------------------*/
  void clear_waste()
  {
    verts_.clear_waste();
    quads_.clear_waste();
    tris_.clear_waste();
    intr_edges_.clear_waste();
    bdry_edges_.clear_waste();
  }

  /*------------------------------------------------------------------
  | Check if the mesh is empty
  ------------------------------------------------------------------*/
  bool is_empty() const 
  { return ( n_vertices()==0 && n_elements()==0 && n_edges()==0 ); }

  /*------------------------------------------------------------------
  | Return mesh entities within a given position and radius
  ------------------------------------------------------------------*/
  VertexVector 
  get_vertices(const Vec2d& center, double radius) const
  { return std::move( verts_.get_items(center, radius ) ); }

  QuadVector
  get_quads(const Vec2d& center, double radius) const
  { return std::move( quads_.get_items(center, radius ) ); }

  TriVector
  get_triangles(const Vec2d& center, double radius) const
  { return std::move( tris_.get_items(center, radius ) ); }

  EdgeVector
  get_intr_edges(const Vec2d& center, double radius) const
  { return std::move( intr_edges_.get_edges(center, radius ) ); }

  EdgeVector
  get_bdry_edges(const Vec2d& center, double radius) const
  { return std::move( bdry_edges_.get_edges(center, radius ) ); }

  /*------------------------------------------------------------------
  | Return all valid interior mesh edges.
  | An interior edge is valid, if it is connected to its left and 
  | its right facet.
  | -> Assumes initialized facet connectivity
  ------------------------------------------------------------------*/
  EdgeVector get_valid_interior_edges() const
  {
    EdgeVector valid_edges {};

    for ( const auto& e_ptr : intr_edges_ )
      if ( NullFacet::is_not_null( e_ptr->facet_l() ) && 
           NullFacet::is_not_null( e_ptr->facet_r() )  )
        valid_edges.push_back( e_ptr.get() );

    return std::move(valid_edges);
  }

  /*------------------------------------------------------------------
  | Return all valid boundary mesh edges.
  | A boundary edge is valid, if it is connected to its left facet.
  | -> Assumes initialized facet connectivity
  ------------------------------------------------------------------*/
  EdgeVector get_valid_boundary_edges() const 
  {
    EdgeVector valid_edges {};

    for ( const auto& e_ptr : bdry_edges_ )
      if ( NullFacet::is_not_null( e_ptr->facet_l() ) )
        valid_edges.push_back( e_ptr.get() );

    return std::move(valid_edges);
  }

  /*------------------------------------------------------------------
  | Return boundary edges that act as interfaces to neighboring 
  | meshes (= twin edges)
  | -> Assumes initialized facet connectivity
  ------------------------------------------------------------------*/
  EdgeVector get_interface_edges() const
  {
    EdgeVector interface_edges {};

    for ( const auto& e_ptr : bdry_edges_ )
      if (e_ptr->twin_edge())
        interface_edges.push_back( e_ptr.get() );

    return std::move( interface_edges );
  }

  /*------------------------------------------------------------------
  | Return all mesh edges that are part of the advancing front:
  | -> Either interior edges that are not connected to two facets 
  |    or boundary edges that are not connected to their left facet.
  | -> These edges are used for the advancing front generation 
  |    of the mesh
  | -> This function assumes an initialized facet connectivity
  ------------------------------------------------------------------*/
  EdgeVector get_front_edges() const
  {
    EdgeVector invalid_edges {};

    for ( const auto& e_ptr : intr_edges_ )
      if ( NullFacet::is_null( e_ptr->facet_l() ) || 
           NullFacet::is_null( e_ptr->facet_r() )  )
        invalid_edges.push_back( e_ptr.get() );

    for ( const auto& e_ptr : bdry_edges_ )
      if ( NullFacet::is_null( e_ptr->facet_l() ) )
        invalid_edges.push_back( e_ptr.get() );

    return std::move(invalid_edges);
  }

  /*------------------------------------------------------------------
  | For a given pair of vertices (v1,v2) return a corresponding
  | interior edge from the mesh that connects them
  ------------------------------------------------------------------*/
  Edge* 
  get_interior_edge(const Vertex& v1, const Vertex& v2, bool dir=false) 
  const 
  { return intr_edges_.get_edge(v1, v2, dir); }

  Edge* 
  get_boundary_edge(const Vertex& v1, const Vertex& v2, bool dir=false) 
  const 
  { return bdry_edges_.get_edge(v1, v2, dir); }

  /*------------------------------------------------------------------
  | For a given pair of vertices (v1,v2) return a corresponding
  | edge from the mesh that connects them
  ------------------------------------------------------------------*/
  Edge* get_edge(const Vertex& v1, const Vertex& v2, bool dir=false)
  const 
  {
    Edge* found = get_interior_edge(v1, v2, dir);

    if ( found != nullptr )
      return found;

    found = get_boundary_edge(v1, v2, dir);

    return found;

  } // Mesh::get_edge()


  /*------------------------------------------------------------------
  | Add new mesh vertices  
  ------------------------------------------------------------------*/
  Vertex& add_vertex(const Vec2d& xy)
  {
    Vertex& v_new = verts_.push_back( xy );
    return v_new;
  } 

  /*------------------------------------------------------------------
  | Add new mesh quads  
  ------------------------------------------------------------------*/
  Quad& add_quad(Vertex& v1, Vertex& v2, Vertex& v3, Vertex& v4,
                 int color=-1)
  {
    Quad& q_new = quads_.push_back(v1, v2, v3, v4);

    if ( color < 0 )
      q_new.color( elem_color_ );
    else
      q_new.color( color );

    q_new.mesh( this );
    return q_new;
  } 

  /*------------------------------------------------------------------
  | Add new mesh triangles  
  ------------------------------------------------------------------*/
  Triangle& add_triangle(Vertex& v1, Vertex& v2, Vertex& v3, 
                         int color=-1)
  {
    Triangle& t_new = tris_.push_back(v1, v2, v3);

    if ( color < 0 )
      t_new.color( elem_color_ );
    else 
      t_new.color( color );

    t_new.mesh( this );
    return t_new;
  } 

  /*------------------------------------------------------------------
  | Add new interior mesh edge   
  ------------------------------------------------------------------*/
  Edge& add_interior_edge(Vertex& v1, Vertex& v2)
  { return intr_edges_.add_edge(v1, v2, INTERIOR_EDGE_MARKER); }

  /*------------------------------------------------------------------
  | Add new boundary mesh edge   
  ------------------------------------------------------------------*/
  Edge& add_boundary_edge(Vertex& v1, Vertex& v2, int marker)
  { return bdry_edges_.add_edge(v1, v2, marker); }

  /*------------------------------------------------------------------
  | These functions remove mesh entities and makes sure, that the 
  | removal succeeded
  ------------------------------------------------------------------*/
  void remove_vertex(Vertex& v)
  {
    bool removed = verts_.remove( v );
    ASSERT( removed, "Failed to remove vertex."); (void) removed;
  }

  void remove_quad(Quad& q)
  {
    bool removed = quads_.remove( q );
    ASSERT( removed, "Failed to remove quad."); (void) removed;
  }

  void remove_triangle(Triangle& t)
  {
    bool removed = tris_.remove( t );
    ASSERT( removed, "Failed to remove triangle."); (void) removed;
  }

  void remove_interior_edge(Edge& e)
  {
    if ( &e.edgelist() != &intr_edges_ )
      return;

    bool removed = intr_edges_.remove( e );
    ASSERT( removed, "Failed to remove interior edge."); (void) removed;
  }

  void remove_boundary_edge(Edge& e)
  {
    if ( &e.edgelist() != &bdry_edges_ )
      return;

    bool removed = bdry_edges_.remove( e );
    ASSERT( removed, "Failed to remove interior edge."); (void) removed;
  }


private:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  int        mesh_id_        { 0 };
  int        elem_color_     { DEFAULT_ELEMENT_COLOR };
  bool       mesh_completed_ { false };
  double     mesh_area_      { 0.0 };

  Vertices   verts_;
  Quads      quads_;
  Triangles  tris_;
  EdgeList   intr_edges_ { Orientation::NONE };
  EdgeList   bdry_edges_ { Orientation::NONE };

}; // Mesh


/*********************************************************************
* Print out the mesh to std::cout
*********************************************************************/
inline std::ostream& operator<<(std::ostream& os, const Mesh& mesh)
{
  auto interior_edges  = mesh.get_valid_interior_edges();
  auto boundary_edges  = mesh.get_valid_boundary_edges();
  auto interface_edges = mesh.get_interface_edges();
  auto front_edges     = mesh.get_front_edges();


  os << "MESH " << mesh.id() << "\n";

  // Print out vertex coordinates
  os << "VERTICES " << mesh.vertices().size() << "\n";
  for ( const auto& v_ptr : mesh.vertices() )
    os << std::setprecision(5) << std::fixed 
       << v_ptr->xy().x << "," 
       << v_ptr->xy().y << "\n";

  // Print out all valid interior edges
  os << "INTERIOREDGES " << interior_edges.size() << "\n";
  for ( auto& e_ptr : interior_edges )
    os << std::setprecision(0) << std::fixed 
      << std::setw(4) << e_ptr->v1().index() << "," 
      << std::setw(4) << e_ptr->v2().index() << ","
      << std::setw(4) << e_ptr->facet_l()->index() << ","
      << std::setw(4) << e_ptr->facet_r()->index() << "\n";

  // Print out all valid boundary edges
  os << "BOUNDARYEDGES " << boundary_edges.size() << "\n";
  for ( auto& e_ptr : boundary_edges )
    os << std::setprecision(0) << std::fixed 
      << std::setw(4) << e_ptr->v1().index() << "," 
      << std::setw(4) << e_ptr->v2().index() << ","
      << std::setw(4) << e_ptr->facet_l()->index() << ","
      << std::setw(4) << e_ptr->marker() << "\n";

  // Print out all interface edges to other meshes
  os << "INTERFACEEDGES " << interface_edges.size() << "\n";
  for ( auto& e_ptr : interface_edges )
    os << std::setprecision(0) << std::fixed 
      << std::setw(4) << e_ptr->v1().index() << "," 
      << std::setw(4) << e_ptr->v2().index() << ","
      << std::setw(4) << e_ptr->facet_l()->index() << ","
      << std::setw(4) << e_ptr->twin_edge()->facet_l()->index() << ","
      << std::setw(4) << e_ptr->twin_edge()->facet_l()->color() << "\n";

  // Print out all advancing front edges
  os << "FRONT " << front_edges.size() << "\n";
  for ( auto& e_ptr : front_edges )
    os << std::setprecision(0) << std::fixed 
      << std::setw(4) << e_ptr->v1().index() << "," 
      << std::setw(4) << e_ptr->v2().index() << ","
      << std::setw(4) << -1 << "\n";

  // Print out all quads
  os << "QUADS " << mesh.quads().size() << "\n";
  for ( const auto& q_ptr : mesh.quads() )
    os << std::setprecision(0) << std::fixed
      << std::setw(4) << q_ptr->v1().index() << ","
      << std::setw(4) << q_ptr->v2().index() << ","
      << std::setw(4) << q_ptr->v3().index() << ","
      << std::setw(4) << q_ptr->v4().index() << ","
      << std::setw(4) << q_ptr->color() << "\n";

  // Print out all triangles
  os << "TRIANGLES " << mesh.triangles().size() << "\n";
  for ( const auto& t_ptr : mesh.triangles() )
    os << std::setprecision(0) << std::fixed
      << std::setw(4) << t_ptr->v1().index() << ","
      << std::setw(4) << t_ptr->v2().index() << ","
      << std::setw(4) << t_ptr->v3().index() << ","
      << std::setw(4) << t_ptr->color() << "\n";

  // Print out all quad neighbors
  os << "QUADNEIGHBORS " << mesh.quads().size() << "\n";
  for ( const auto& q_ptr : mesh.quads() )
  {
    auto nbr1_index = q_ptr->nbr1() ? q_ptr->nbr1()->index() : -1;
    auto nbr2_index = q_ptr->nbr2() ? q_ptr->nbr2()->index() : -1;
    auto nbr3_index = q_ptr->nbr3() ? q_ptr->nbr3()->index() : -1;
    auto nbr4_index = q_ptr->nbr4() ? q_ptr->nbr4()->index() : -1;

    os << std::setprecision(0) << std::fixed << std::setw(4) 
      << nbr1_index << "," << std::setw(4) 
      << nbr2_index << "," << std::setw(4) 
      << nbr3_index << "," << std::setw(4) 
      << nbr4_index << "\n";
  }

  // Print out all triangle neighbors
  os << "TRIANGLENEIGHBORS " << mesh.triangles().size() << "\n";
  for ( const auto& t_ptr : mesh.triangles() )
  {
    auto nbr1_index = t_ptr->nbr1() ? t_ptr->nbr1()->index() : -1;
    auto nbr2_index = t_ptr->nbr2() ? t_ptr->nbr2()->index() : -1;
    auto nbr3_index = t_ptr->nbr3() ? t_ptr->nbr3()->index() : -1;

    os << std::setprecision(0) << std::fixed << std::setw(4) 
      << nbr1_index << "," << std::setw(4) 
      << nbr2_index << "," << std::setw(4) 
      << nbr3_index << "\n";
  }

  os << "SIZEFUNCTION " << mesh.vertices().size() << "\n";
  for ( const auto& v_ptr : mesh.vertices() )
  {
    os << std::setprecision(5) << std::fixed 
       << v_ptr->mesh_size() << "\n";
  }

  return os;
} 



} // namespace TQAlgorithm
} // namespace TQMesh
