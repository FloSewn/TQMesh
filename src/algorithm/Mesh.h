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
//#include "VtkIO.h"

#include "utils.h"
#include "Vertex.h"
#include "Triangle.h"
#include "Quad.h"
#include "Facet.h"

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
<<<<<<< HEAD
  , tris_       { qtree_scale, qtree_items, qtree_depth }
  , quads_      { qtree_scale, qtree_items, qtree_depth }
  , front_      { }
  { }
  
  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  friend std::ostream& operator<<(std::ostream& os, const Mesh& mesh);

  /*------------------------------------------------------------------
  | Return mesh entities within a given position and radius
  ------------------------------------------------------------------*/
  VertexVector 
  get_vertices(const Vec2d& center, double radius) const
  { return std::move( verts_.get_items(center, radius ) ); }

  EdgeVector
  get_intr_edges(const Vec2d& center, double radius) const
  { return std::move( intr_edges_.get_edges(center, radius ) ); }

  EdgeVector
  get_bdry_edges(const Vec2d& center, double radius) const
  { return std::move( bdry_edges_.get_edges(center, radius ) ); }

  TriVector
  get_triangles(const Vec2d& center, double radius) const
  { return std::move( tris_.get_items(center, radius ) ); }

  QuadVector
  get_quads(const Vec2d& center, double radius) const
  { return std::move( quads_.get_items(center, radius ) ); }

  /*------------------------------------------------------------------
  | Add new mesh vertices  
  ------------------------------------------------------------------*/
  Vertex& add_vertex(const Vec2d& xy)
  {
    Vertex& v_new = verts_.push_back( xy );
    return v_new;
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
  | Add another mesh to the neighbor-mesh-list. This is only 
  | possible prior to the initialization of the advancing front.
  ------------------------------------------------------------------*/
  bool add_neighbor_mesh(Mesh& mesh)
  { 
    if ( mesh_initialized_ )
    {
      LOG(ERROR) << 
      "The advancing front of mesh " << mesh_id_ << " " << 
      "is already initialized. "
      "Adding other meshes as neighbors is only possible prior to "
      "the advancing front initialization.";
      return false;
    }

    neighbor_meshes_.push_back( &mesh ); 

    return true;
  }

  /*------------------------------------------------------------------
  | Remove another mesh to the neighbor-mesh-list
  ------------------------------------------------------------------*/
  void remove_neighbor_mesh(Mesh& mesh)
  {
    auto it = std::find(neighbor_meshes_.begin(), 
                        neighbor_meshes_.end(), 
                        &mesh);
    if ( it != neighbor_meshes_.end() )
      neighbor_meshes_.erase( it );
  }

  /*------------------------------------------------------------------
  | Initialize the advancing front structure of the mesh
  ------------------------------------------------------------------*/
  bool init_advancing_front()
  {
    // Check if boundaries are defined 
    if ( domain_->size() < 1 )
    {
      LOG(ERROR) << 
      "Can not initialize the mesh's advancing front, "
      "since no domain boundaries are defined.";
      return false;
    }

    // Check if domain boundaries are traversable and 
    // if an exterior boundary exists
    bool extr_bdry_found = false;

    for ( const auto& boundary : *domain_ )
    {
      if ( boundary->is_exterior() )
        extr_bdry_found = true;

      Edge& e_start = boundary->edges()[0];
      if ( !boundary->is_traversable(e_start, e_start) )
      {
        LOG(ERROR) << 
        "Can not initialize advancing front. "
        "Mesh boundaries are not traversable.";
        return false;
      }
    }

    if ( !extr_bdry_found )
    {
      LOG(ERROR) << 
      "Can not initialize advancing front. "
      "No exterior mesh boundary found.";
      return false;
    }

    // Count the number of overlaps with other meshes
    size_t n_overlaps = 0;

    for ( auto neighbor_mesh : neighbor_meshes_ )
    {
      Domain& neighbor_domain = neighbor_mesh->domain();

      n_overlaps += domain_->get_overlaps( neighbor_domain );
    }

    if ( n_overlaps > 0 )
    {
      LOG(INFO) <<
      "Mesh " << mesh_id_ << " overlaps with " << n_overlaps << 
      " edges of other meshes. These edges will be considered "
      "in the meshing process.";
    }

    // Get all edges that will define the advancing front
    FrontInitData front_data = collect_front_edges();

    // Initialize edges in the advancing front
    front_.init_front_edges(*domain_, front_data, verts_);
     
    // Setup the mesh boundary edges from the initial advancing front
    for ( auto& e : front_ )
    {
      Vertex& v1 = e->v1();
      Vertex& v2 = e->v2();
      int marker = e->marker();
      Edge& e_new = bdry_edges_.add_edge( v1, v2, marker );

      // Connect boundary edges of this mesh and its parner mesh
      Edge* e_twin = e->twin_edge();

      if ( e_twin )
      {
        e_twin->twin_edge( &e_new );
        e_new.twin_edge( e_twin );
        e->twin_edge( nullptr );
      }
    }

    // Set the mesh to be initialized
    mesh_initialized_ = true;

    return true;

  } // Mesh::init_advancing_front()

  /*------------------------------------------------------------------
  | This function takes care of the garbage collection 
  ------------------------------------------------------------------*/
  void clear_waste()
  {
    verts_.clear_waste();
    tris_.clear_waste();
    quads_.clear_waste();
    front_.clear_waste();
    intr_edges_.clear_waste();
    bdry_edges_.clear_waste();
  }

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  const Domain& domain() const 
  { ASSERT(domain_,"Invalid mesh domain"); return *domain_; }
  Domain& domain() 
  { ASSERT(domain_,"Invalid mesh domain"); return *domain_; }

  const Front& front() const { return front_; }
  Front& front() { return front_; }

  const Triangles& triangles() const { return tris_; }
  Triangles& triangles() { return tris_; }

  const Quads& quads() const { return quads_; }
  Quads& quads() { return quads_; }

  const Vertices& vertices() const { return verts_; }
  Vertices& vertices() { return verts_; }

  const EdgeList& interior_edges() const { return intr_edges_; }
  EdgeList& interior_edges() { return intr_edges_; }

  const EdgeList& boundary_edges() const { return bdry_edges_; }
  EdgeList& boundary_edges() { return bdry_edges_; }

  double area() const { return mesh_area_; }
  int id() const { return mesh_id_; }
  int element_color() const { return elem_color_; }
  size_t neighbor_meshes() const { return neighbor_meshes_.size(); }
  bool completed() const { return mesh_completed_; }
  bool initialized() const { return mesh_initialized_; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  void element_color(int c) { elem_color_ = c; }

  /*------------------------------------------------------------------
  | Export the mesh to a file
  ------------------------------------------------------------------*/
  void write_to_file(const std::string& path, ExportType type )
  {
    // Prepare all indices of mesh entities
    assign_mesh_indices();

    if ( type == ExportType::cout )
      std::cout << (*this);

    if ( type == ExportType::txt )
      write_to_txt_file( path );

    if ( type == ExportType::vtu )
      write_to_vtu_file( path );

    
  } // Mesh::write_to_file()

  /*------------------------------------------------------------------
  | This function assigns a corresponding global index to each entity
  | of the current mesh
  | --> Quads are stored prior to triangles
  ------------------------------------------------------------------*/
  void assign_mesh_indices()
  {
    unsigned int v_index = 0; 
    unsigned int e_index = 0;

    // Vertices 
    for ( const auto& v_ptr : verts_ )
      v_ptr->index( v_index++ );

    // Quads
    for ( const auto& q_ptr : quads_ )
      q_ptr->index( e_index++ );

    // Triangles 
    for ( const auto& t_ptr : tris_ )
      t_ptr->index( e_index++ );

  } // assign_mesh_indices()

  /*------------------------------------------------------------------
  | Every triangle is refined into three quads and every 
  | quads is refined into four quads.
  | This results in an all-quad mesh.
  |
  |                        
  |          (v7)         (v6)      (v5)      qi... sub-quads 
  |             o<--------o---------o         vi... vertices
  |             |         |         ^
  |             |  [q4]   |   [q3]  |
  |             |         |         |
  |             |         |(v9)     |  
  |         (v8)o---------o---------o(v4) 
  |             |         |         |
  |             |         |         |
  |             |  [q1]   |   [q2]  |
  |             v         |         |
  |             o---------o-------->o
  |          (v1)        (v2)       (v3)
  |                         
  |                         
  |                        (v5)                 
  |                         o                   
  |                       /   \
  |                     / [q3]  \
  |                   /           \
  |           (v6)  /      (v7)     \  (v4)
  |               o---------o---------o
  |             /           |           \
  |           /             |             \
  |         /     [q1]      |      [q2]     \
  |       /                 |                 \
  |     o-------------------o-------------------o
  |    (v1)                (v2)                 (v3)
  |
  ------------------------------------------------------------------*/
  bool refine_to_quads()
  {
    // Check if the mesh was properly generated
    if ( !mesh_completed_ )
    {
      LOG(ERROR) <<
      "Unable to refine mesh " << mesh_id_ << ", " <<
      "since the mesh elements have not yet been generated.";
      return false;
    }

    // Check if this mesh is connected to any other meshes
    // If yes, do not refine
    if ( neighbor_meshes() > 0 )
    {
      LOG(ERROR) << 
      "Mesh " << mesh_id_ << " can not be refined to quad elements, "
      "since it is connected to other meshes. ";
      return false;
    }

    LOG(INFO) << "Start with quad refinement of mesh " << mesh_id_;
    clear_waste();

    // Gather all coarse edges, quads and tris
    EdgeVector coarse_intr_edges {};
    EdgeVector coarse_bdry_edges {};
    TriVector  coarse_tris {};
    QuadVector coarse_quads {};

    for ( auto& e : intr_edges_ )
      coarse_intr_edges.push_back( e.get() );

    for ( auto& e : bdry_edges_ )
      coarse_bdry_edges.push_back( e.get() );

    for ( auto& t_ptr : tris_ )
      coarse_tris.push_back( t_ptr.get() );

    for ( auto& q_ptr : quads_ )
      coarse_quads.push_back( q_ptr.get() );

    // Refine interior edges
    for ( auto e : coarse_intr_edges )
    {
      Vertex& v = add_vertex( e->xy() );
      intr_edges_.add_edge( e->v1(), v );
      intr_edges_.add_edge( v, e->v2() );
      e->sub_vertex( &v );

      if ( e->v1().is_fixed() && e->v2().is_fixed() )
        v.is_fixed( true );
    }

    // Refine boundary edges
    for ( auto e : coarse_bdry_edges )
    {
      Vertex& v = add_vertex( e->xy() );
      bdry_edges_.add_edge( e->v1(), v, e->marker() );
      bdry_edges_.add_edge( v, e->v2(), e->marker() );
      e->sub_vertex( &v );
      v.on_boundary( true );

      if ( e->v1().is_fixed() && e->v2().is_fixed() )
        v.is_fixed( true );
    }

    // Refine all triangle elements
    for ( auto t : coarse_tris )
    {
      Vertex& v1 = t->v1();
      Vertex& v3 = t->v2();
      Vertex& v5 = t->v3();

      Edge* e13 = get_edge( v1, v3 );
      Edge* e35 = get_edge( v3, v5 );
      Edge* e51 = get_edge( v5, v1 );

      Vertex* v2 = e13->sub_vertex();
      Vertex* v4 = e35->sub_vertex();
      Vertex* v6 = e51->sub_vertex();

      // Create new vertex at center of quad
      Vertex& v7 = add_vertex( t->xy() );

      // Create new sub-quads 
      Quad& q1 = add_quad( v1, *v2, v7, *v6 );
      Quad& q2 = add_quad( v3, *v4, v7, *v2 );
      Quad& q3 = add_quad( v5, *v6, v7, *v4 );

      // New quads get assigned to colors of old element
      q1.color( t->color() );
      q2.color( t->color() );
      q3.color( t->color() );

      // Create new interior edges
      intr_edges_.add_edge( *v2, v7 );
      intr_edges_.add_edge( *v4, v7 );
      intr_edges_.add_edge( *v6, v7 );

    }

    // Refine all quad elements
    for ( auto q : coarse_quads )
    {
      Vertex& v1 = q->v1();
      Vertex& v3 = q->v2();
      Vertex& v5 = q->v3();
      Vertex& v7 = q->v4();

      Edge* e13 = get_edge( v1, v3 );
      Edge* e35 = get_edge( v3, v5 );
      Edge* e57 = get_edge( v5, v7 );
      Edge* e71 = get_edge( v7, v1 );

      Vertex* v2 = e13->sub_vertex();
      Vertex* v4 = e35->sub_vertex();
      Vertex* v6 = e57->sub_vertex();
      Vertex* v8 = e71->sub_vertex();

      // Create new vertex at center of quad
      Vertex& v9 = add_vertex( q->xy() );

      // Create new sub-quads 
      Quad& q1 = add_quad( v1, *v2, v9, *v8 );
      Quad& q2 = add_quad( v3, *v4, v9, *v2 );
      Quad& q3 = add_quad( v5, *v6, v9, *v4 );
      Quad& q4 = add_quad( v7, *v8, v9, *v6 );

      // New quads get assigned to colors of old element
      q1.color( q->color() );
      q2.color( q->color() );
      q3.color( q->color() );
      q4.color( q->color() );

      // Create new interior edges
      intr_edges_.add_edge( *v2, v9 );
      intr_edges_.add_edge( *v4, v9 );
      intr_edges_.add_edge( *v6, v9 );
      intr_edges_.add_edge( *v8, v9 );

    }

    // Remove old entitires
    for ( auto e : coarse_intr_edges )
      remove_interior_edge( *e );
      
    for ( auto e : coarse_bdry_edges )
      remove_boundary_edge( *e );

    for ( auto t : coarse_tris )
      remove_triangle( *t );

    for ( auto q : coarse_quads )
      remove_quad( *q );


    // Re-initialize vertex-to-vertex connectivity
    setup_vertex_connectivity();

    // Update connectivity between elements and edges
    setup_facet_connectivity();

    LOG(INFO) << "Mesh refinement is completed.";
      
    return true;
      
  } // Mesh::refine_to_quads()

  /*------------------------------------------------------------------
  | Algorithm to create a single quad layer for a connected list of
  | advancing front edges that start with v_start and end with v_end
  ------------------------------------------------------------------*/
  bool create_quad_layers(Vertex& v_start, Vertex& v_end, 
                          size_t n_layers, double first_height,
                          double growth_ratio)
  {
    if ( !mesh_initialized_ )
    {
      LOG(ERROR) << 
      "Unable to create quad layers for mesh " << mesh_id_ << ". " <<
      "The mesh's advancing front has not yet been initialized.";
      return false;
    }

    if ( mesh_completed_ )
    {
      LOG(ERROR) <<
      "Unable to create quad layers for mesh " << mesh_id_ << ". " <<
      "Quad layers must be created prior to the meshing process. ";
    }

    bool quads_generated;

    double height = first_height;

    Vertex* v1 = &v_start;
    Vertex* v2 = &v_end;

    for ( size_t i = 0; i < n_layers; ++i )
    {
      LOG(INFO) << 
      "Generate quad layer (" << i+1 << 
      "/" << n_layers << ") for mesh " << mesh_id_;
      quads_generated = add_quad_layer( v1, v2, height );

      if ( !quads_generated )
      {
        LOG(ERROR) <<
        "Generation of quad layer (" << i+1 << 
        "/" << n_layers << ") failed. " <<
        "Quad layer generation is aborted.";

        return false;
      }

      height *= growth_ratio;
    }

    setup_facet_connectivity();
    merge_triangles_to_quads();

    return true;

  } // Mesh::create_quad_layers()

  /*------------------------------------------------------------------
  | Create a triangular mesh using the advancing front algorithm
  ------------------------------------------------------------------*/
  bool triangulate()
  {
    if ( !mesh_initialized_ )
    {
      LOG(ERROR) << 
      "Unable to triangulate mesh " << mesh_id_ << ". " <<
      "The mesh's advancing front has not yet been initialized.";
      return false;
    }

    unsigned int iter = 0;
    bool wide_search = false;
    int sort_iter = 0;

    ProgressBar progress_bar {};

    // Initialize base edge
    front_.set_base_first();
    Edge* base = &( front_.base() );

    // Invalid base definition
    if ( !base ) 
    {
      LOG(ERROR) << 
      "Unable to triangulate mesh " << mesh_id_ << ". " <<
      "The mesh's advancing front structure seems to corrupted.";
      return false;
    }

    LOG(INFO) << "Start triangulation of mesh " << mesh_id_ << ".";

    // Sort the advancing front edges
    front_.sort_edges( false );

    // Start advancing front loop
    while ( true )
    {
      // Try to advance the current base edge
      bool success = advance_front_triangle(*base, wide_search);

      // If it worked, reset iteration counter and wide search
      // and go to the next base edge
      if ( success )
      {
        // Sort front edges after a wide search
        if ( wide_search )
          front_.sort_edges( false );

        iter = 0;
        wide_search = false;

        // Sort the advancing front edges
        if ( sort_iter == CONSTANTS.sort_triangulation_front() )
        {
          front_.sort_edges( false );
          sort_iter = 0;
        }
        else
        {
          sort_iter = MOD( sort_iter+1, INT_MAX-1 );
        }

        front_.set_base_first();
        base = &( front_.base() );

        clear_waste();
      }
      // If it failed, go to the next base edge
      else
      {
        front_.set_base_next();
        base = &( front_.base() );
        ++iter;
      }

      // All front edges failed to create new elements
      // --> Activate wide search for neighboring vertices
      //     and re-run the algorithm
      if ( iter == front_.size() && !wide_search )
      {
        wide_search = true;
        iter = 0;
      }

      // Update progress bar
      double state = std::ceil(100.0 * area() / domain_->area());
      progress_bar.update( static_cast<int>(state) );
      progress_bar.show( LOG_PROPERTIES.get_ostream(INFO) );

      // No more edges in the advancing front
      // --> Meshing algorithm succeeded
      if ( front_.size() == 0 )
      {
        LOG(INFO) << "[";
        LOG(INFO) << "The meshing process succeeded.";
        break;
      }

      // All front edges faild to create new elements, even
      // when using the wide search 
      // --> Meshing algorithm failed
      if ( iter == front_.size() && wide_search )
      {
        LOG(INFO) << "[";
        LOG(ERROR) << "The meshing process failed.";
        return false;
      }
    }

    // Initialize vertex-to-vertex connectivity
    setup_vertex_connectivity();

    // Initialize facet-to-facet connectivity
    // as well as connectivity between edges and facets
    setup_facet_connectivity();

    // Set the mesh generation to be completed
    mesh_completed_ = true;

    return true;

  } // Mesh::triangulate()

  /*------------------------------------------------------------------
  | Create a quadrilateral mesh using the advancing front algorithm
  ------------------------------------------------------------------*/
  bool pave()
  {
    if ( !mesh_initialized_ )
    {
      LOG(ERROR) << 
      "Unable to pave mesh " << mesh_id_ << ". " <<
      "The mesh's advancing front has not yet been initialized.";
      return false;
    }

    unsigned int iter = 0;
    bool wide_search = false;

    ProgressBar progress_bar {};

    // Initialize base edge
    front_.set_base_first();
    Edge* base = &( front_.base() );

    // Invalid base definition
    if ( !base ) 
    {
      LOG(ERROR) << 
      "Unable to triangulate mesh " << mesh_id_ << ". " <<
      "The mesh's advancing front structure seems to corrupted.";
      return false;
    }

    LOG(INFO) << "Start paving of mesh " << mesh_id_ << ".";

    // Start advancing front loop
    while ( true )
    {
      // Try to advance the current base edge
      bool success = advance_front_quad(*base, wide_search, -1.0);

      // If it worked, reset iteration counter and wide search
      // and go to the next base edge
      if ( success )
      {
        iter = 0;
        wide_search = false;

        front_.set_base_first();
        base = &( front_.base() );

        clear_waste();
      }
      // If it failed, go to the next base edge
      else
      {
        front_.set_base_next();
        base = &( front_.base() );
        ++iter;
      }

      // All front edges failed to create new elements
      // --> Activate wide search for neighboring vertices
      //     and re-run the algorithm
      if ( iter == front_.size() && !wide_search )
      {
        wide_search = true;
        iter = 0;
      }

      // Update progress bar
      double state = std::ceil(100.0 * area() / domain_->area());
      progress_bar.update( static_cast<int>(state) );
      progress_bar.show( LOG_PROPERTIES.get_ostream(INFO) );

      // No more edges in the advancing front
      // --> Meshing algorithm succeeded
      if ( front_.size() == 0 )
      {
        LOG(INFO) << "";
        LOG(INFO) << "The meshing process succeeded.";
        break;
      }

      // All front edges faild to create new elements, even
      // when using the wide search 
      // --> Meshing algorithm failed
      if ( iter == front_.size() && wide_search )
      {
        LOG(INFO) << "";
        LOG(ERROR) << "The meshing process failed.";
        return false;
      }
    }

    // Initialize vertex-to-vertex connectivity
    setup_vertex_connectivity();

    // Initialize facet-to-facet connectivity
    // as well as connectivity between edges and facets
    setup_facet_connectivity();

    // Set the mesh generation to be completed
    mesh_completed_ = true;

    // Merge remaining triangles to quads
    merge_triangles_to_quads();

    return true;

  } // Mesh::pave()


  /*------------------------------------------------------------------
  |    
  ------------------------------------------------------------------*/
  void merge_degenerate_triangles()
  {
    // Collect all vertices, that are adjacent to exactly three 
    // triangles
    std::list<Vertex*> bad_vertices {};

    for ( const auto& vertex_ptr : verts_ )
    {
      if ( vertex_ptr->facets().size() != 3 )
        continue;
      if ( vertex_ptr->edges().size() != 3 )
        continue;
      if ( vertex_ptr->facets(0).n_vertices() != 3 )
        continue;
      if ( vertex_ptr->facets(1).n_vertices() != 3 )
        continue;
      if ( vertex_ptr->facets(2).n_vertices() != 3 )
        continue;

      auto color = vertex_ptr->facets(0).color();

      if ( vertex_ptr->facets(1).color() != color )
        continue;
      if ( vertex_ptr->facets(2).color() != color )
        continue;

      bad_vertices.push_back( vertex_ptr.get() );
    }


    for ( auto& bad_vertex : bad_vertices )
    {
      const Triangle* t0 = static_cast<const Triangle*>(
          &bad_vertex->facets(0));

      const Triangle* t1 = static_cast<const Triangle*>(
          &bad_vertex->facets(1));

      const Triangle* t2 = static_cast<const Triangle*>(
          &bad_vertex->facets(2));

      auto color = t0->color();

      // Get surrounding vertices
      int i0 = t0->get_vertex_index( *bad_vertex );
      int i1 = t1->get_vertex_index( *bad_vertex );
      int i2 = t2->get_vertex_index( *bad_vertex );

      Vertex& v0 = const_cast<Vertex&>( t0->vertex(MOD(i0+1,3)) );
      Vertex& v1 = const_cast<Vertex&>( t1->vertex(MOD(i1+1,3)) );
      Vertex& v2 = const_cast<Vertex&>( t2->vertex(MOD(i2+1,3)) );

      // Remove interior edges 
      Edge* e0 = intr_edges_.get_edge(*bad_vertex, v0);
      Edge* e1 = intr_edges_.get_edge(*bad_vertex, v1);
      Edge* e2 = intr_edges_.get_edge(*bad_vertex, v2);

      ASSERT( e0, "Interior edge not defined.\n" );
      remove_interior_edge( *e0 );

      ASSERT( e1, "Interior edge not defined.\n" );
      remove_interior_edge( *e1 );

      ASSERT( e2, "Interior edge not defined.\n" );
      remove_interior_edge( *e2 );

      // Remove triangles
      remove_triangle(const_cast<Triangle&>(*t0));
      remove_triangle(const_cast<Triangle&>(*t1));
      remove_triangle(const_cast<Triangle&>(*t2));

      // Remove vertex
      remove_vertex( *bad_vertex );

      // Add new triangle
      auto o = orientation(v0.xy(), v1.xy(), v2.xy());

      if ( o == Orientation::CCW ) 
        add_triangle( v0, v1, v2, color );
      else
        add_triangle( v0, v2, v1, color );
    }

    // Remove deleted entities
    clear_waste();

    // Re-initialize vertex-to-vertex connectivity
    setup_vertex_connectivity();

    // Re-initialize facet-to-facet connectivity
    setup_facet_connectivity();

  } // Mesh::merge_degenerate_triangles() 

  /*------------------------------------------------------------------
  | This function loops over all internal edges and, if possible,
  | merges two adjacent triangles to one quad element.
  |
  | -> This function requires a finalized mesh and that the functions
  |    setup_vertex_connectivity() and setup_facet_connectivity()
  |    have been called before.
  ------------------------------------------------------------------*/
  void merge_triangles_to_quads()
  {
    // Pick all internal edges, which are adjacent to two triangles
    std::list<Edge*> tri_edges {};

    for ( const auto& e_ptr : intr_edges_ )
    {
      Facet* f_l = e_ptr->facet_l();
      Facet* f_r = e_ptr->facet_r();

      if (  (f_l && f_l->n_vertices() == 3) 
         && (f_r && f_r->n_vertices() == 3) )
        tri_edges.push_back( e_ptr.get() );
    }

    // Sort edge list with increasing minimum angles of the 
    // adjacent triangles
    tri_edges.sort(
    []( Edge* a, Edge* b )
    {
      const double a_l = ( a->facet_l() )
                       ? a->facet_l()->min_angle() 
                       : TQ_MAX;
      const double a_r = ( a->facet_r() )
                       ? a->facet_r()->min_angle() 
                       : TQ_MAX;
      const double a_ang = MIN(a_l, a_r);

      const double b_l = ( b->facet_l() )
                       ? b->facet_l()->min_angle() 
                       : TQ_MAX;
      const double b_r = ( b->facet_r() )
                       ? b->facet_r()->min_angle() 
                       : TQ_MAX;
      const double b_ang = MIN(b_l, b_r);

      return a_ang < b_ang;

    });

    // Loop over all chosen edges and merge the adjacent triangles
    // to quadrilaterals. It may be, that the triangle of an 
    // upcoming edge has already been merged in this process. Thus
    // we have to make sure, that the triangles to merge still
    // exist.
    //
    //   v2               q_r
    //     *-------------*
    //     | \           |
    //     |   \   f_r   |
    //     |     \       |
    //     |       \     |
    //     |  f_l    \   |
    //     |           \ |
    //     *-------------*
    //    q_l             v1
    //
    for ( auto& e : tri_edges )
    {
      Facet* f_l = e->facet_l();
      Facet* f_r = e->facet_r();

      if ( f_l->n_vertices() != 3 || f_r->n_vertices() != 3 )
        continue;

      Vertex& v1 = e->v1();
      Vertex& v2 = e->v2();


      int i_l = f_l->get_edge_index(v1, v2);
      int i_r = f_r->get_edge_index(v1, v2);

      Vertex& q_l = f_l->vertex(i_l);
      Vertex& q_r = f_r->vertex(i_r);

      // Remove internal edge
      remove_interior_edge( *e );
      e = nullptr;

      // Create new quadrilateral element
      Quad& q_new = add_quad( q_l, v1, q_r, v2 );
      q_new.is_active( true );

      // Update internal edge connectivity to prevent bad memory 
      // access for upcoming internal edges that are no longer 
      // adjacent to two triangles due to the merge
      for (int i = 0; i < 4; ++i)
      {
        int i1 = i;
        int i2 = MOD(i+1, 4);

        Vertex& q1 = q_new.vertex(i1);
        Vertex& q2 = q_new.vertex(i2);

        Edge* e_share = intr_edges_.get_edge(q1,q2);

        if ( !e_share ) continue;

        Facet* t_l = e_share->facet_l();
        Facet* t_r = e_share->facet_r();

        if ( t_l && (t_l == f_l || t_l == f_r) )
          e_share->facet_l( &q_new );

        if ( t_r && (t_r == f_l || t_r == f_r) )
          e_share->facet_r( &q_new );
      }

      // Remove triangles
      remove_triangle( *(static_cast<Triangle*>(f_l)) );
      remove_triangle( *(static_cast<Triangle*>(f_r)) );

    }

    // Re-initialize facet-to-facet connectivity
    setup_facet_connectivity();

    // Bad elements may have been created up to this point
    // due to the merging of triangles to quads
    // The two upcoming function fix these bad elements
    clean_double_quad_edges();
    clean_double_triangle_edges();
    
    // Remove deleted entities
    clear_waste();

    // Re-initialize vertex-to-vertex connectivity
    setup_vertex_connectivity();

    // Re-initialize facet-to-facet connectivity
    setup_facet_connectivity();

  } // Mesh::merge_triangles_to_quads()

  /*------------------------------------------------------------------
  | Merge this mesh with one of its neighbor meshes. All elements, 
  | edges and vertices of the neighbor mesh will be added to this 
  | mesh. The other mesh will then be removed from this mesh's 
  | neighbor list and instead will be placed in its list of merged
  | meshes.
  ------------------------------------------------------------------*/
  void merge_neighbor_mesh(Mesh& neighbor)
  {
    if ( !mesh_completed_ || !mesh_initialized_ )
    {
      LOG(ERROR) <<
      "Unable to merge mesh " << mesh_id_ << 
      " and mesh " << neighbor.id() <<
      ", since mesh " << mesh_id_ << " is not yet completed.";
    }

    if ( !neighbor.completed() || !neighbor.initialized() )
    {
      LOG(ERROR) <<
      "Unable to merge mesh " << mesh_id_ << 
      " and mesh " << neighbor.id() <<
      ", since mesh " << neighbor.id() << " is not yet completed.";
    }

    LOG(INFO) <<
    "Merging mesh " << mesh_id_ << " with mesh " << neighbor.id();

    // Assign indices to neighbor mesh
    neighbor.assign_mesh_indices();

    // Search for vertices that are located on the inteface
    // boundary between this mesh and the neighbor mesh
    VertexVector new_vertices( neighbor.vertices().size(), nullptr );

    for ( const auto& e_ptr : bdry_edges_ )
    {
      Edge* e_twin = e_ptr->twin_edge();

      if ( !e_twin )
        continue;

      Facet* f_twin = e_twin->facet_l();
      ASSERT( f_twin, "Invalid boundary edge connectivity." );

      if ( f_twin->mesh() != &neighbor )
        continue;

      auto v1_index = e_twin->v1().index();
      auto v2_index = e_twin->v2().index();

      // Keep in mind, that both edges point in different directions
      new_vertices[v2_index] = &(e_ptr->v1());
      new_vertices[v1_index] = &(e_ptr->v2());
    }

    // Add the neighbor's vertices which are not located
    // on the interface boundary
    for ( const auto& v_ptr : neighbor.vertices() )
    {
      auto v_index = v_ptr->index();

      if ( !new_vertices[v_index]  )
      {
        Vertex& v_new = add_vertex( v_ptr->xy() );
        new_vertices[v_index] = &v_new;
      }
    }
    
    // Add the neighbor's triangle elements
    for ( const auto& t_ptr : neighbor.triangles() )
    {
      auto v1_index = t_ptr->v1().index();
      auto v2_index = t_ptr->v2().index();
      auto v3_index = t_ptr->v3().index();

      Vertex* v1 = new_vertices[v1_index];
      Vertex* v2 = new_vertices[v2_index];
      Vertex* v3 = new_vertices[v3_index];

      ASSERT( v1, "Invalid data structure.");
      ASSERT( v2, "Invalid data structure.");
      ASSERT( v3, "Invalid data structure.");

      add_triangle( *v1, *v2, *v3, t_ptr->color() );
    }

    // Add the neighbor's quad elements
    for ( const auto& q_ptr : neighbor.quads() )
    {
      auto v1_index = q_ptr->v1().index();
      auto v2_index = q_ptr->v2().index();
      auto v3_index = q_ptr->v3().index();
      auto v4_index = q_ptr->v4().index();

      Vertex* v1 = new_vertices[v1_index];
      Vertex* v2 = new_vertices[v2_index];
      Vertex* v3 = new_vertices[v3_index];
      Vertex* v4 = new_vertices[v4_index];

      ASSERT( v1, "Invalid data structure.");
      ASSERT( v2, "Invalid data structure.");
      ASSERT( v3, "Invalid data structure.");
      ASSERT( v4, "Invalid data structure.");

      add_quad( *v1, *v2, *v3, *v4, q_ptr->color() );
    }

    // Add the neighbor's interior edges 
    for ( const auto& e_ptr : neighbor.interior_edges() )
    {
      auto v1_index = e_ptr->v1().index();
      auto v2_index = e_ptr->v2().index();

      Vertex* v1 = new_vertices[v1_index];
      Vertex* v2 = new_vertices[v2_index];

      ASSERT( v1, "Invalid data structure.");
      ASSERT( v2, "Invalid data structure.");

      intr_edges_.add_edge( *v1, *v2 );
    }

    // Add the neighbor's boundary edges 
    for ( const auto& e_ptr : neighbor.boundary_edges() )
    {
      Edge* e_twin = e_ptr->twin_edge();

      // Skip if edge is an interface to this mesh
      if ( e_twin && 
           e_twin->facet_l() && 
           e_twin->facet_l()->mesh() == this )
        continue;

      auto v1_index = e_ptr->v1().index();
      auto v2_index = e_ptr->v2().index();

      Vertex* v1 = new_vertices[v1_index];
      Vertex* v2 = new_vertices[v2_index];

      ASSERT( v1, "Invalid data structure.");
      ASSERT( v2, "Invalid data structure.");

      Edge& e_new = bdry_edges_.add_edge( *v1, *v2, e_ptr->marker() );

      // Add interface to possible other neighbor meshes
      if ( e_twin )
      {
        e_twin->twin_edge( &e_new );
        e_new.twin_edge( e_twin );
      }
    }

    // Find boundary inteface edges that should be turned into 
    // interior edges
    std::vector<std::pair<Vertex*, Vertex*>> new_intr_edges {};
    EdgeVector bdry_edges_to_remove {};
    
    for ( const auto& e_ptr : bdry_edges_ )
    {
      Edge* e_twin = e_ptr->twin_edge();

      if ( !e_twin )
        continue;

      Facet* f_twin = e_twin->facet_l();
      ASSERT( f_twin, "Invalid boundary edge connectivity." );

      if ( f_twin->mesh() != &neighbor )
        continue;

      auto v1_index = e_twin->v1().index();
      auto v2_index = e_twin->v2().index();

      Vertex* v1 = new_vertices[v1_index];
      Vertex* v2 = new_vertices[v2_index];

      ASSERT( v1, "Invalid data structure.");
      ASSERT( v2, "Invalid data structure.");

      new_intr_edges.push_back( {v2, v1} );
      bdry_edges_to_remove.push_back( e_ptr.get() );
    }

    // Remove old edges
    for ( Edge* e : bdry_edges_to_remove )
      remove_boundary_edge( *e );

    // Add new interior edges
    for ( auto v : new_intr_edges )
      intr_edges_.add_edge( *v.first, *v.second );

    // Add neighbors of neighbor mesh
    // --> Check if these neighbors are possibly already merged 
    //     with this mesh
    for ( Mesh* nbr : neighbor.neighbor_meshes_ )
    {
      if ( nbr == this )
        continue;

      // Check if neigbhor of neighbor is also adjacent to 
      // this mesh
      auto nbr_found = std::find(neighbor_meshes_.begin(), 
                                 neighbor_meshes_.end(), nbr);

      if( nbr_found != neighbor_meshes_.end() ) 
        continue;

      // Check if neigbhor of neighbor is already merged to 
      // this mesh
      nbr_found = std::find(merged_meshes_.begin(), 
                            merged_meshes_.end(), nbr);

      if( nbr_found != merged_meshes_.end() ) 
        continue;

      // Add neighbor mesh
      neighbor_meshes_.push_back( nbr );
    }

    // Remove neighbor connectivity
    neighbor.remove_neighbor_mesh(*this);
    remove_neighbor_mesh(neighbor);
    merged_meshes_.push_back( &neighbor );
    clear_waste();

    // Re-initialize vertex-to-vertex connectivity
    setup_vertex_connectivity();

    // Update connectivity between elements and edges
    setup_facet_connectivity();

    LOG(INFO) << "Meshes have been merged successfully";

  } // Mesh::merge_neighbor_mesh()

private:
  /*------------------------------------------------------------------
  | Let the advancing front create a new quad element
  ------------------------------------------------------------------*/
  bool advance_front_quad(Edge& base, bool wide_search=false,
                          double height = -1.0)
  {
    // Constants
    const double fq = CONSTANTS.quad_range_factor();

    Vertex& b1 = base.v1();
    Vertex& b2 = base.v2();

    double h_b1 = height;
    double h_b2 = height;

    Vec2d p1 = b1.xy() + h_b1 * base.normal();
    Vec2d p2 = b2.xy() + h_b2 * base.normal();

    if ( height <= 0 )
    {
      h_b1 = domain_->size_function( b1.xy() );
      h_b2 = domain_->size_function( b2.xy() );
      
      p1 = b1.xy() + h_b1 * base.normal();
      p2 = b2.xy() + h_b2 * base.normal();

      const double rho_q1 = domain_->size_function( p1 );
      const double rho_q2 = domain_->size_function( p2 );

      const double theta_1 = MAX(-0.25,MIN(0.25,1.0-(h_b1-rho_q1)/h_b1));
      const double theta_2 = MAX(-0.25,MIN(0.25,1.0-(h_b2-rho_q2)/h_b2));

      p1 += theta_1 * base.length() * base.tangent();
      p2 -= theta_2 * base.length() * base.tangent();
    }

    // Search radii
    double r1 = domain_->size_function( p1 );
    double r2 = domain_->size_function( p2 );

    if ( height > 0 )
    {
      r1 = fq * height;
      r2 = fq * height;
    }

    // ****** Create first triangle *******
    TriVector new_tris_p1 {};

    // Find all vertices in vicinity of p1 
    VertexVector vertex_candidates_p1 
      = find_local_vertices(p1, r1, wide_search );

    // Create potential triangles with all found vertices
    check_vertex_candidates(vertex_candidates_p1, base, new_tris_p1);

    // If potential triangles have been found, choose the best one
    // --> this function also updates the advancing front
    Triangle* t1 = choose_best_triangle(new_tris_p1, base);

    // No triangle has been found yet -> create new one
    if ( !t1 )
    {
      Vertex&   v_new = add_vertex( p1 );
      Triangle& t_new = add_triangle( b1, b2, v_new );

      // Algorithm fails if new vertex or new triangle is invalid
      if ( remove_if_invalid(v_new, t_new) )
        return false;

      // Update the advancing front with new vertex
      update_front( base, v_new, t_new );

      // Keep track of the new triangle
      t1 = &t_new;
    }

    // ****** Create second triangle *******
    Vertex& d1 = t1->v3();
    Vertex& d2 = t1->v2();

    Edge* diag = front_.get_edge(d1, d2);

    // --> created triangle t1 is not located on the front
    if ( !diag )
      return true;

    TriVector new_tris_p2 {};

    // Find all vertices in vicinity of p2 
    VertexVector vertex_candidates_p2 
      = find_local_vertices(p2, r2, wide_search );

    // Create potential triangles with all found vertices
    check_vertex_candidates(vertex_candidates_p2, *diag, new_tris_p2);

    // If potential triangles have been found, choose the best one
    // --> this function also updates the advancing front
    Triangle* t2 = choose_best_triangle(new_tris_p2, *diag);

    // No triangle has been found yet -> create new one
    if ( !t2 )
    {
      Vertex&   v_new = add_vertex( p2 );
      Triangle& t_new = add_triangle(d1, d2, v_new);

      // Algorithm fails if new vertex or new triangle is invalid
      if ( remove_if_invalid(v_new, t_new) )
        return true;

      // Update the advancing front with new vertex
      update_front( *diag, v_new, t_new );

      // Keep track of the new triangle
      t2 = &t_new;
    }

    // ****** Merge the newly created triangles  *******
    // Gather vertices
    Vertex& q1 = t1->v1();
    Vertex& q2 = t1->v2();
    Vertex& q3 = t2->v3();
    Vertex& q4 = t2->v1();

    // Remove internal edge
    Edge* e_rem = intr_edges_.get_edge(q2, q4);
    remove_interior_edge( *e_rem );

    // Remove old triangular elements
    remove_triangle( *t1 );
    remove_triangle( *t2 );

    // Create new quadrilateral element
    Quad& q_new = add_quad( q1, q2, q3, q4 );
    q_new.is_active( true );

    return true;

  } // Mesh::advance_front_quad()

  /*------------------------------------------------------------------
  | Let the advancing front create a new triangle 
  ------------------------------------------------------------------*/
  bool advance_front_triangle(Edge& base, bool wide_search=false)
  {
    // Constants
    const double f1 = CONSTANTS.base_height_factor();
    const double f2 = CONSTANTS.mesh_range_factor();

    // Obtain the position of a potential new vertex and the radius 
    // to search around it for potential existing vertices
    const double height = domain_->size_function( base.xy() );
    const Vec2d  v_xy   = base.xy() + f1 * height * base.normal();
    const double r      = domain_->size_function( v_xy );

    // Find all vertices in the vicinity of the current base edge
    VertexVector vertex_candidates 
      = find_local_vertices(v_xy, f2 * r, wide_search);

    // Create potential triangles with all found vertices
    TriVector new_triangles {};
    check_vertex_candidates(vertex_candidates, base, new_triangles);

    // If potential triangles have been found, choose the best one
    if ( choose_best_triangle(new_triangles, base) )
      return true;

    // Check if a potential triangle can be created with the base edge
    // and a newly created vertex
    Vertex& b1 = base.v1();
    Vertex& b2 = base.v2();

    Vertex&   v_new = get_base_vertex( base );
    Triangle& t_new = add_triangle(b1, b2, v_new);
    
    // Algorithm fails if new vertex or new triangle is invalid
    if ( remove_if_invalid(v_new, t_new) )
      return false;
    
    // Update the advancing front with new vertex
    update_front( base, v_new, t_new );

    return true;

  } // Mesh::advance_front_triangle()


  /*------------------------------------------------------------------
  | This function is part of the advancing front loop.
  | For a given location <xy> and a respective search range <dist>,
  | all vertices that are located in this vicinity are put into a 
  | vector and the sorted ascendingly towards <xy>.
  | If the flag <wide_search> is activated, the search range is 
  | enlarged by a preset factor, in order to finde more vertex 
  | candidates.
  ------------------------------------------------------------------*/
  VertexVector find_local_vertices(const Vec2d& xy, double dist, 
                                   bool wide_search=false )
  {
    if (wide_search)
      dist *= CONSTANTS.wide_search_factor();

    // Get vertices in vicinity of xy  
    VertexVector vertex_candidates = verts_.get_items(xy, dist);

    // Sort vertices in ascending order towards xy
    std::sort( vertex_candidates.begin(), vertex_candidates.end(), 
    [xy] ( const Vertex* a, const Vertex* b )
    {
      return ( (a->xy()-xy).length_squared() 
             < (b->xy()-xy).length_squared() );
    });

    return std::move( vertex_candidates ); 

  } // find_local_vertices() 

  /*------------------------------------------------------------------
  | This function is part of the advancing front loop.
  | We loop over a given set of <vertex_candidates> and check
  | if we can create a possible triangle with the current base edge
  | (<b1>,<b2>) and the given vertices.
  | If it is possible, new triangles are created and pushed back to 
  | the vector <new_triangles>.
  ------------------------------------------------------------------*/
  void check_vertex_candidates(const VertexVector& vertex_candidates,
                               Edge& base, TriVector& new_triangles)
  {
    for ( Vertex* v : vertex_candidates )
    {
      // Skip vertices that are not located on the advancing front
      if ( !v->on_front() )
        continue;

      // Skip vertices that are colinear to the current base edge
      if ( orientation( base.v1().xy(), base.v2().xy(), v->xy() )
          == Orientation::CL )
        continue;

      // Create new potential triangle 
      Triangle& t_new = add_triangle( base.v1(), base.v2(), *v );

      // Check if new potential triangle is valid
      if ( !remove_if_invalid(t_new) )
        new_triangles.push_back( &t_new );

    }

  } // check_vertex_candidates()

  /*------------------------------------------------------------------
  | This function is part of the advancing front loop.
  | We sort a given vector of <new_triangles> in descending order
  | according to the triangle quality.
  | Finally, the advancing front is updated with the triangle of best
  | quality and all other triangles are removed.
  ------------------------------------------------------------------*/
  Triangle* choose_best_triangle(TriVector& new_triangles,
                                 Edge&      base)
  {
    if ( new_triangles.size() < 1 )
      return nullptr;

    DEBUG_LOG(
      "VALID TRIANGLES IN NEIGHBORHOOD: " 
       << new_triangles.size()
    );

    std::sort( new_triangles.begin(), new_triangles.end(),
    [this] ( Triangle* t1, Triangle* t2 )
    {
      const double h1 = domain_->size_function( t1->xy() );
      const double h2 = domain_->size_function( t2->xy() );
      const double q1 = t1->quality(h1);
      const double q2 = t2->quality(h2);

      return ( q1 > q2 );
    });

    Triangle* new_tri = new_triangles[0];
    Vertex&   v_adj   = new_tri->v3();

    update_front( base, v_adj, *new_tri );

    for (int i = 1; i < new_triangles.size(); i++)
      tris_.remove( *new_triangles[i] );

    return new_tri;

  } // choose_best_triangle()

  /*------------------------------------------------------------------
  | Check if a triangle is valid. If yes, return true - 
  | else return false.
  ------------------------------------------------------------------*/
  bool triangle_is_valid(const Triangle& tri)
  {
    const double rho   = domain_->size_function( tri.xy() );
    const double range = 2.0 * rho;

    DEBUG_LOG("CHECK NEW TRIANGLE: " << tri);

    if ( !tri.is_valid() )
      return false;

    if ( tri.intersects_front( front_, range ) )
    {
      DEBUG_LOG("  > FRONT INTERSECTION");
      return false;
    }

    if ( tri.intersects_vertex( verts_, range ) )
    {
      DEBUG_LOG("  > VERTEX INTERSECTION");
      return false;
    }

    if ( tri.intersects_triangle( tris_, range ) )
    {
      DEBUG_LOG("  > TRIANGLE INTERSECTION");
      return false;
    }

    if ( tri.intersects_quad( quads_, range ) )
    {
      DEBUG_LOG("  > QUAD INTERSECTION");
      return false;
    }

    if ( tri.quality(rho) < CONSTANTS.min_cell_quality() )
    {
      DEBUG_LOG("  > BAD TRIANGLE QUALITY");
      return false;
    }

    if ( tri.max_angle() > CONSTANTS.max_cell_angle() )
    {
      DEBUG_LOG("  > BAD MAXIMUM ANGLE");
      return false;
    }


    DEBUG_LOG("  > VALID");

    return true;

  } // Mesh::triangle_is_valid()

  /*------------------------------------------------------------------
  | Check if a vertex is valid. If yes, return true - 
  | else return false.
  ------------------------------------------------------------------*/
  bool vertex_is_valid(const Vertex& v)
  {
    const double rho   = domain_->size_function( v.xy() );
    const double range = 2.0 * rho;

    DEBUG_LOG("CHECK NEW VERTEX: " << v);

    if ( !domain_->is_inside( v ) )
    {
      DEBUG_LOG("  > OUTSIDE DOMAIN");
      return false;
    }

    if ( v.intersects_facet(tris_, range) )
    {
      DEBUG_LOG("  > TRIANGLE INTERSECTION");
      return false;
    }

    if ( v.intersects_facet(quads_, range) )
    {
      DEBUG_LOG("  > QUAD INTERSECTION");
      return false;
    }

    DEBUG_LOG("  > VALID");

    return true;

  } // Mesh::vertex_is_valid()

  /*------------------------------------------------------------------
  | Check if mesh entities are invalid. If yes, remove them and 
  | return true. Else, simply return false.
  ------------------------------------------------------------------*/
  bool remove_if_invalid(Vertex& v)
  {
    if ( vertex_is_valid(v) )
      return false;
    remove_vertex(v);
    return true;
  } 

  bool remove_if_invalid(Triangle& t)
  {
    if ( triangle_is_valid(t) )
      return false;
    remove_triangle(t);
    return true;
  } 

  bool remove_if_invalid(Vertex& v, Triangle& t)
  {
    if ( !vertex_is_valid(v) || !triangle_is_valid(t) )
    {
      remove_triangle(t);
      remove_vertex(v);
      return true;
    }
    return false;
  } 

  bool remove_if_invalid(Vertex& v, Triangle& t1, Triangle& t2)
  {
    if (  !vertex_is_valid(v) 
       || !triangle_is_valid(t1) || !triangle_is_valid(t2) )
    {
      remove_triangle(t1);
      remove_triangle(t2);
      remove_vertex(v);
      return true;
    }
    return false;
  }
=======
  , quads_      { qtree_scale, qtree_items, qtree_depth }
  , tris_       { qtree_scale, qtree_items, qtree_depth }
  { }
>>>>>>> fa0899f5faedbc3de2d30dba4c8c9fc7b7288940

  /*------------------------------------------------------------------
  | Move constructor
  ------------------------------------------------------------------*/
<<<<<<< HEAD
  Vertex& get_base_vertex(const Edge& base)
  {
    // Half of the factor h for height of equlateral triangle
    // h := sqrt(3) / 2  -   h_fac := h / 2
    constexpr double h_fac = 0.4330127019; 
    const double v_fac = h_fac * CONSTANTS.base_vertex_factor();

    // Obtain size function value at the centroid of an equlateral
    // triangle, created from the current base edge
    Vec2d c = base.xy() + base.normal() * base.length() * v_fac;
    const double rho = domain_->size_function(c);

    // Coordinate of new vertex 
    Vec2d xy = base.xy() + base.normal() * rho;

    return add_vertex( xy );

  } // get_base_vertex() 
=======
  Mesh(Mesh&& mesh)
  : mesh_id_    { mesh.mesh_id_ }
  , elem_color_ { mesh.elem_color_ }
  , mesh_area_  { mesh.mesh_area_ }
  , verts_      { std::move( mesh.verts_      )}
  , quads_      { std::move( mesh.quads_      )}
  , tris_       { std::move( mesh.tris_       )}
  , intr_edges_ { std::move( mesh.intr_edges_ )}
  , bdry_edges_ { std::move( mesh.bdry_edges_ )}
  {}
>>>>>>> fa0899f5faedbc3de2d30dba4c8c9fc7b7288940

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
    EdgeVector front_edges {};

    for ( const auto& e_ptr : intr_edges_ )
      if ( NullFacet::is_null( e_ptr->facet_l() ) || 
           NullFacet::is_null( e_ptr->facet_r() )  )
        front_edges.push_back( e_ptr.get() );

    for ( const auto& e_ptr : bdry_edges_ )
      if ( NullFacet::is_null( e_ptr->facet_l() ) )
        front_edges.push_back( e_ptr.get() );

    return std::move(front_edges);
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

<<<<<<< HEAD
  /*------------------------------------------------------------------
  | Export the mesh to a text file
  ------------------------------------------------------------------*/
  void write_to_txt_file(const std::string& path) const
  {
    std::ofstream outfile;

    std::string file_name = path;

    if(file_name.substr(file_name.find_last_of(".") + 1) != "txt") 
      file_name += ".txt";

    outfile.open( file_name );

    outfile << (*this);

    outfile.close();

  } // Mesh::write_to_txt_file()

  /*------------------------------------------------------------------
  | Export the mesh to a vtu file
  ------------------------------------------------------------------*/
  void write_to_vtu_file(const std::string& path) const
  {
    std::string file_name = path;

    if(file_name.substr(file_name.find_last_of(".") + 1) != "vtu") 
      file_name += ".vtu";

    // Create data structure for VTU format
    std::vector<double> points {};
    std::vector<size_t> connectivity {};
    std::vector<size_t> offsets {};
    std::vector<size_t> types {};
    std::vector<double> size_function {};
    std::vector<int>    element_color {};
    std::vector<double> edge_length {};
    std::vector<double> max_angle {};
    std::vector<double> cell_quality {};

    size_t i_offset = 0;

    for ( const auto& v_ptr : verts_ )
    {
      points.push_back( v_ptr->xy().x );
      points.push_back( v_ptr->xy().y );
      points.push_back( 0.0 );

      size_function.push_back( domain_->size_function(v_ptr->xy()) ); 
    }

    for ( const auto& q_ptr : quads_ )
    {
      connectivity.push_back( q_ptr->v1().index() );
      connectivity.push_back( q_ptr->v2().index() );
      connectivity.push_back( q_ptr->v3().index() );
      connectivity.push_back( q_ptr->v4().index() );

      i_offset += 4;
      offsets.push_back( i_offset );

      /// Type == 9 -> VTK_QUAD
      types.push_back( 9 );

      element_color.push_back( q_ptr->color() );

      edge_length.push_back( q_ptr->max_edge_length() );
      max_angle.push_back( q_ptr->max_angle() * 180. / M_PI );

      const double s = domain_->size_function(q_ptr->xy());
      cell_quality.push_back( q_ptr->quality(s) );
    }

    for ( const auto& t_ptr : tris_ )
    {
      connectivity.push_back( t_ptr->v1().index() );
      connectivity.push_back( t_ptr->v2().index() );
      connectivity.push_back( t_ptr->v3().index() );

      i_offset += 3;
      offsets.push_back( i_offset );

      /// Type == 5 -> VTK_TRIANGLE
      types.push_back( 5 );

      element_color.push_back( t_ptr->color() );

      edge_length.push_back( t_ptr->max_edge_length() );
      max_angle.push_back( t_ptr->max_angle() * 180. / M_PI );

      const double s = domain_->size_function(t_ptr->xy());
      cell_quality.push_back( t_ptr->quality(s) );
    }


    VtuWriter writer { points, connectivity, offsets, types };

    writer.add_point_data( size_function, "size_function", 1 );
    writer.add_cell_data( element_color, "element_color", 1 );
    writer.add_cell_data( edge_length, "edge_length", 1 );
    writer.add_cell_data( max_angle, "max_angle", 1 );
    writer.add_cell_data( cell_quality, "cell_quality", 1 );

    writer.write( file_name );

  } // Mesh::write_to_vtu_file()

=======
>>>>>>> fa0899f5faedbc3de2d30dba4c8c9fc7b7288940

private:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  int        mesh_id_        { 0 };
  int        elem_color_     { DEFAULT_ELEMENT_COLOR };
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
