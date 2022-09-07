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

#include "Vec2.h"
#include "utils.h"
#include "VtkIO.h"
#include "ProgressBar.h"

#include "Vertex.h"
#include "Triangle.h"
#include "Quad.h"
#include "Front.h"
#include "Boundary.h"
#include "Domain.h"
#include "QuadLayer.h"

/*********************************************************************
* ToDo:
* - Add boolean to indicate if mesh elements have been generated yet
* - Allow addition of neighbor meshes only when no mesh elements
*   are generated yet
* - Make ASSERT similar to CHECK
* - User may only provide positive mesh color values
*********************************************************************/

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* The mesh 
*********************************************************************/
class Mesh
{
  using EdgePair       = std::pair<const Edge*,const Edge*>;
  using EdgePairVector = std::vector<EdgePair>;
  using EdgeVector     = std::vector<Edge*>;
  using VertexVector   = std::vector<Vertex*>;
  using TriVector      = std::vector<Triangle*>;
  using QuadVector     = std::vector<Quad*>;
  using DoubleVector   = std::vector<double>;
  using BoolVector     = std::vector<bool>;
  using IntVector      = std::vector<int>;
  using Vec2dVector    = std::vector<Vec2d>;
  using MeshVector     = std::vector<Mesh*>;
  using BdryEdgeConn   = std::vector<std::vector<EdgeVector>>;
  using NbrMeshConn    = std::vector<std::vector<EdgeVector>>;
  using FrontData      = std::pair<EdgeVector,BoolVector>;

public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Mesh(Domain&   domain,
       int       mesh_id=CONSTANTS.default_mesh_id(),
       int       element_color=CONSTANTS.default_element_color(),
       double    qtree_scale=ContainerQuadTreeScale,
       size_t    qtree_items=ContainerQuadTreeItems, 
       size_t    qtree_depth=ContainerQuadTreeDepth)
  : domain_     { &domain }
  , mesh_id_    { ABS(mesh_id) }
  , elem_color_ { ABS(element_color) }
  , verts_      { qtree_scale, qtree_items, qtree_depth }
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

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  Vertex& get_base_vertex(const Edge& base)
  {
    // Half of the factor h for height of equlateral triangle
    // h := sqrt(3) / 2  -   h_fac := h / 2
    constexpr double h_fac = 0.4330127019; 
    const double v_fac = CONSTANTS.base_vertex_factor();

    // Obtain size function value at the centroid of an equlateral
    // triangle, created from the current base edge
    Vec2d c = base.xy() + base.normal() * base.length() * h_fac;
    const double rho = domain_->size_function(c);

    // Coordinate of new vertex 
    Vec2d xy = base.xy() + base.normal() * rho * v_fac;

    return add_vertex( xy );

  } // get_base_vertex() 

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  void update_front( Edge& base, Vertex& v_new, Triangle& t_new )
  {
    // Get advancing front edges adjacent to vertex
    // -> First two vertices of new triangle tri are always
    //    the base edge vertices
    Edge* e1 = front_.get_edge(v_new, base.v1());
    Edge* e2 = front_.get_edge(v_new, base.v2());

    // *** Both edges are connected to vertex ***
    //     -> No new edge must be created
    //     -> Base vertex v1 no longer on the advancing front
    //     -> Base vertex v2 no longer on the advancing front
    //     -> vertex no longer on the advancing front
    if ( e1 && e2 )
    {
      ASSERT( v_new.on_front(), "Grid structure corrupted." );

      base.v1().on_front( false );
      base.v2().on_front( false );
      v_new.on_front( false );

      if ( e1->is_interior() )
        intr_edges_.add_edge(e1->v1(), e1->v2());
      if ( e2->is_interior() )
        intr_edges_.add_edge(e2->v1(), e2->v2());

      front_.remove( *e1 );
      front_.remove( *e2 );
    }
    // *** First edge is connected to vertex ***
    //     -> New edge between second base vertex and vertex
    //     -> Base vertex v1 no longer on the advancing front
    //     -> vertex is part of the advancing front
    else if ( e1 && !e2 )
    {
      ASSERT( v_new.on_front(), "Grid structure corrupted." );

      base.v1().on_front( false );

      if ( e1->is_interior() )
        intr_edges_.add_edge(e1->v1(), e1->v2());

      front_.remove( *e1 );
      front_.add_edge(v_new, base.v2());
    }
    // *** Second edge is connected to vertex ***
    //     -> New edge between first base vertex and vertex
    //     -> Base vertex v2 no longer on the advancing front
    //     -> vertex is part of the advancing front
    else if ( !e1 && e2 )
    {
      ASSERT( v_new.on_front(), "Grid structure corrupted." );

      base.v2().on_front( false );

      if ( e2->is_interior() )
        intr_edges_.add_edge(e2->v1(), e2->v2());

      front_.remove( *e2 );
      front_.add_edge(base.v1(), v_new);
    }
    // *** Both edges are not connected to vertex ***
    //     -> Create two new edges
    //     -> Vertex now part of the advancing front
    else
    {
      v_new.on_front( true );
      front_.add_edge(base.v1(), v_new);
      front_.add_edge(v_new, base.v2());
    }

    // If current base is not at the boundary, add it to the 
    // interior edge list
    if ( base.is_interior() )
      intr_edges_.add_edge(base.v1(), base.v2());

    // Remove base edge
    front_.remove( base );

    // Mark new triangle as active
    t_new.is_active( true );

    // Add element area to the total mesh area
    mesh_area_ += t_new.area();

  } // update_front() 

  /*------------------------------------------------------------------
  | Every vertex gets assigned its neighboring vertices and these 
  | are then sorted by means of ascending angles
  | This function requires, that all the interior edges and 
  | boundary edges of the mesh have been generated 
  | -> element adjacency of the edges is not required here
  ------------------------------------------------------------------*/
  void setup_vertex_connectivity()
  {
    // Remove all current vertex-to-vertex connectivities
    for ( auto& v_ptr : verts_ )
      v_ptr->vertices().clear();

    // Get vertex-to-vertex connectivity from interior edges
    for ( const auto& e : intr_edges_ )
    {
      Vertex* v1 = &(e->v1());
      Vertex* v2 = &(e->v2());

      v1->vertices().push_back( v2 );
      v2->vertices().push_back( v1 );
    }

    // Get vertex-to-vertex connectivity from boundary edges
    for ( const auto& e : bdry_edges_ )
    {
      Vertex* v1 = &(e->v1());
      Vertex* v2 = &(e->v2());

      v1->vertices().push_back( v2 );
      v2->vertices().push_back( v1 );
    }

    // Sort local vertex connectivities by ascending angle
    for ( auto& v_ptr : verts_ )
    {
      const Vec2d xy = v_ptr->xy();

      std::sort( v_ptr->vertices().begin(), v_ptr->vertices().end(),
      [xy] ( Vertex* v1, Vertex* v2 )
      {
        const Vec2d dxy1 = v1->xy() - xy;
        const Vec2d dxy2 = v2->xy() - xy;
        const double a1 = std::atan2(dxy1.y, dxy1.x);
        const double a2 = std::atan2(dxy2.y, dxy2.x);

        return ( a1 < a2 );
      });
    }

  } // Mesh::setup_vertex_connectivity()

  /*------------------------------------------------------------------
  | Initialize the connectivity between facets and facets, as well  
  | as between edges and facets
  ------------------------------------------------------------------*/
  void setup_facet_connectivity()
  {
    Facet* f1 = nullptr;
    Facet* f2 = nullptr;
    
    // Setup connectivity for interor edges
    for ( const auto& e_ptr : intr_edges_ )
    {
      const Vertex& v1 = e_ptr->v1();
      const Vertex& v2 = e_ptr->v2();

      int idx1 {-1};
      int idx2 {-1};

      for ( auto f : v1.facets() )
      {
        int idx = f->get_edge_index(v1, v2);

        if ( idx < 0 ) continue;

        if ( idx1 < 0 )
        {
          idx1 = idx;
          f1 = f;
        }
        else if ( idx2 < 0 )
        {
          idx2 = idx;
          f2 = f;
          break;
        }
      }

      // Setup connectivity between facets
      f1->neighbor( idx1, f2 );
      f2->neighbor( idx2, f1 );

      // Setup connectivity between internal edge and facets
      if ( is_left( v1.xy(), v2.xy(), f1->xy() ) )
      {
        e_ptr->facet_l( f1 );
        e_ptr->facet_r( f2 );
      }
      else
      {
        e_ptr->facet_l( f2 );
        e_ptr->facet_r( f1 );
      }
    }

    // Setup connectivity for boundary edges
    for ( const auto& e_ptr : bdry_edges_ )
    {
      const Vertex& v1 = e_ptr->v1();
      const Vertex& v2 = e_ptr->v2();

      for ( auto f : v1.facets() )
      {
        int idx = f->get_edge_index(v1, v2);

        if ( idx < 0 ) continue;

        f1 = f;
        break;
      }

      // Setup connectivity between internal edge and facets
      if ( is_left( v1.xy(), v2.xy(), f1->xy() ) )
        e_ptr->facet_l( f1 );

    }
    
  } // Mesh::setup_facet_connectivity

  /*------------------------------------------------------------------
  | Clean up quad elements:
  | It may be, that some adjacent quads share two internal edges
  | These elements will be merged in the next step.
  |
  | -> This function requires a finalized mesh and that the functions
  |    setup_vertex_connectivity() and setup_facet_connectivity()
  |    have been called before.
  | 
  |     v3                       vp   
  |       x---------------------x
  |       | \                   |
  |       |   \ e2     q_nbr    |
  |       |     \               |
  |       |       \  v2         |
  |       |         x           |
  |       |           \         |
  |       |   q_cur     \  e1   |
  |       |               \     |
  |       |                 \   |
  |       |                   \ |
  |       x---------------------x
  |     v4                       v1
  ------------------------------------------------------------------*/
  void clean_double_quad_edges()
  {
    // Initialize all quad markers
    for ( auto& q_cur : quads_ )
      q_cur->marker( false );

    std::vector<std::pair<Quad*,Quad*>>     quads_to_remove;
    std::vector<std::pair<Edge*,Edge*>>     edges_to_remove;
    std::vector<Vertex*>                    verts_to_remove;
    std::vector<std::pair<Vertex*,Vertex*>> opposing_vertices;

    for ( auto& q_cur : quads_ )
    {
      if ( q_cur->marker() )
        continue;

      // Loop over all edges of the current quad
      for ( int i_vert = 0; i_vert < 4; ++i_vert )
      {
        // Vertices of the current quad
        Vertex& v1 = q_cur->vertex( i_vert );
        Vertex& v2 = q_cur->vertex( MOD(i_vert + 1, 4) );
        Vertex& v3 = q_cur->vertex( MOD(i_vert + 2, 4) );
        Vertex& v4 = q_cur->vertex( MOD(i_vert + 3, 4) );

        // Get the neighboring facets of the current successive
        // quad edges (v1,v2) and (v2,v3)
        int idx_1 = q_cur->get_edge_index(v1, v2);
        int idx_2 = q_cur->get_edge_index(v2, v3);

        Facet* nbr_1 = q_cur->neighbor(idx_1);
        Facet* nbr_2 = q_cur->neighbor(idx_2);

        // Proceed, if no neighbors are found (nullptr) or if 
        // neighbors of the adjacent edges differ, or if the neighbor
        // has already been added
        if ( !nbr_1 || !nbr_2 || nbr_1 != nbr_2 || nbr_1->marker() )
          continue;

        // In this stage, we address only quad / quad connections
        if ( nbr_1->n_vertices() < 4 )
          continue;

        // Now we can cast the facet to a quad
        Quad* q_nbr = static_cast<Quad*>(nbr_1);

        // Get the internal edges adjacent to both current quads
        Edge* e1 = intr_edges_.get_edge(v1, v2);
        Edge* e2 = intr_edges_.get_edge(v2, v3);

        ASSERT( (e1 != e2), "INVALID DATA STRUCTURE");
        ASSERT( ( e1->v1() == e2->v1() || e1->v1() == e2->v2() 
               || e1->v2() == e2->v1() || e1->v2() == e2->v2() ),
            "WRONG EDGES FOUND.");

        // Mark the current quads, such that these won't get chosen
        // in upcoming loops
        q_cur->marker( true );
        q_nbr->marker( true );

        // Add elements to the removal vectors
        quads_to_remove.push_back( {q_cur.get(), q_nbr} ); 
        edges_to_remove.push_back( {e1, e2} );
        verts_to_remove.push_back( &v2 );

        
        // We still need the vertex of the neighboring quad, that 
        // is located on the opposite of the current edge segments
        // --> Use internal edge definition of quads
        int idx_op = q_nbr->get_edge_index(v2,v3);
        Vertex& v_op = q_nbr->vertex( idx_op );

        opposing_vertices.push_back( {&v4, &v_op} );

        ASSERT( (v_op != v1), "BAD DATA STRUCTURE");
        ASSERT( (v_op != v2), "BAD DATA STRUCTURE");
        ASSERT( (v_op != v3), "BAD DATA STRUCTURE");

        // At this point we can break the inner loop over 
        // the quad edges
        break;
      }
    }

    // Merge the marked edges and created triangles from the 
    // adjacent marked quads
    for (size_t i = 0; i < quads_to_remove.size(); ++i)
    {
      Edge* e1 = edges_to_remove[i].first;
      Edge* e2 = edges_to_remove[i].second;

      Vertex* v2 = verts_to_remove[i];

      // Get correct vertex order of edge segments
      Vertex& v1 = (e1->v2() == *v2) 
                 ? e1->v1() : e1->v2();
      Vertex& v3 = (e2->v1() == *v2) 
                 ? e2->v2() : e2->v1();

      Vertex* o1 = opposing_vertices[i].first;
      Vertex* o2 = opposing_vertices[i].second;

      // Create new quad 
      Quad& q_new = add_quad( *o1, v1, *o2, v3 );
      q_new.is_active(true);

    }

    // Removal of old quads
    for (size_t i = 0; i < quads_to_remove.size(); ++i)
    {
      Quad* q1 = quads_to_remove[i].first;
      Quad* q2 = quads_to_remove[i].second;

      remove_quad( *q1 );
      remove_quad( *q2 );
    }

    // Removal of old interior edges
    for (size_t i = 0; i < edges_to_remove.size(); ++i)
    {
      Edge* e1 = edges_to_remove[i].first;
      Edge* e2 = edges_to_remove[i].second;

      remove_interior_edge( *e1 );
      remove_interior_edge( *e2 );
    }

    // Removal of old vertices
    for (size_t i = 0; i < verts_to_remove.size(); ++i)
    {
      Vertex* v = verts_to_remove[i];
      remove_vertex( *v );
    }

    // Re-initialize facet-to-facet connectivity
    setup_facet_connectivity();

  } // Mesh::clean_double_quad_edges()

  /*------------------------------------------------------------------
  | Clean up triangle elements:
  | It may be, that some triangles share two internal edges with a 
  | single quad element. 
  | These triangles and quads will be removed in this step and then
  | replaced by a single triangle.
  |
  | -> This function requires a finalized mesh and that the functions
  |    setup_vertex_connectivity() and setup_facet_connectivity()
  |    have been called before.
  | -> This function should be called after the function 
  |    clean_double_quad_edges(), because such triangles might be 
  |    genereated during the latter function.
  |
  |        v3
  |       x
  |       | \
  |       |\  \
  |       | \   \
  |       |  \    \
  |       |   \     \
  |       |    \      \
  |       |     \ t_nbr \
  |       |      \        \
  |       |       x--.      \
  |       |      v2   --.     \
  |       |              --.    \
  |       |   q_cur         --.   \
  |       |                    ---  \
  |       x--------------------------x
  |     v4                          v1
  |
  ------------------------------------------------------------------*/
  void clean_double_triangle_edges()
  {
    // Initialize all quad markers to false
    for ( auto& q_cur : quads_ )
      q_cur->marker(false);

    std::vector<std::pair<Quad*,Triangle*>> elements_to_remove;
    std::vector<std::pair<Edge*,Edge*>>     edges_to_remove;
    std::vector<Vertex*>                    verts_to_remove;

    for ( auto& q_cur : quads_ )
    {
      if ( q_cur->marker() )
        continue;

      // Loop over all edges of the current quad
      for ( int i_vert = 0; i_vert < 4; ++i_vert )
      {
        // Vertices of the current quad
        Vertex& v1 = q_cur->vertex( i_vert );
        Vertex& v2 = q_cur->vertex( MOD(i_vert + 1, 4) );
        Vertex& v3 = q_cur->vertex( MOD(i_vert + 2, 4) );

        // Get the neighboring facets of the current successive
        // quad edges (v1,v2) and (v2,v3)
        int idx_1 = q_cur->get_edge_index(v1, v2);
        int idx_2 = q_cur->get_edge_index(v2, v3);

        Facet* nbr_1 = q_cur->neighbor(idx_1);
        Facet* nbr_2 = q_cur->neighbor(idx_2);

        // Proceed, if no neighbors are found (nullptr) or if 
        // neighbors of the adjacent edges differ, or if the neighbor
        // has already been added
        if ( !nbr_1 || !nbr_2 || nbr_1 != nbr_2 || nbr_1->marker() )
          continue;

        // In this stage, we address only quad / triangle connections
        if ( nbr_1->n_vertices() > 3 )
          continue;

        // Now we can cast the facet to a triangle
        Triangle* t_nbr = static_cast<Triangle*>(nbr_1);

        // Get the internal edges adjacent to both current elements
        Edge* e1 = intr_edges_.get_edge(v1, v2);
        Edge* e2 = intr_edges_.get_edge(v2, v3);

        ASSERT( (e1 != e2), "INVALID DATA STRUCTURE");
        ASSERT( ( e1->v1() == e2->v1() || e1->v1() == e2->v2() 
               || e1->v2() == e2->v1() || e1->v2() == e2->v2() ),
            "WRONG EDGES FOUND.");

        // Mark the current quads, such that the won't get chosen
        // in upcoming loops
        q_cur->marker( true );
        t_nbr->marker( true );

        // Add elements to the removal vectors
        elements_to_remove.push_back( {q_cur.get(), t_nbr} ); 
        edges_to_remove.push_back( {e1, e2} );
        verts_to_remove.push_back( &v2 );

        // At this point we can break the inner loop over 
        // the quad edges
        break;
      }
    }

    // Create a new triangle from the remaining vertices
    for ( size_t i = 0; i < elements_to_remove.size(); ++i )
    {
      Quad*   q = elements_to_remove[i].first;
      Vertex* v = verts_to_remove[i];

      int id_v = q->get_vertex_index( *v );

      ASSERT( (id_v > -1), "BAD DATA STRUCTURE" );

      // Get the remaining three vertices of the quad 
      Vertex& v1 = q->vertex( MOD(id_v+1, 4) );
      Vertex& v2 = q->vertex( MOD(id_v+2, 4) );
      Vertex& v3 = q->vertex( MOD(id_v+3, 4) );

      // Create new triangle 
      Triangle& t_new = add_triangle( v1, v2, v3 );
      t_new.is_active(true);
    }

    // Removal of old elements
    for ( size_t i = 0; i < elements_to_remove.size(); ++i )
    {
      Quad*     q = elements_to_remove[i].first;
      Triangle* t = elements_to_remove[i].second;
      
      remove_quad( *q );
      remove_triangle( *t );
    }

    // Removal of old interiord edges
    for (size_t i = 0; i < edges_to_remove.size(); ++i)
    {
      Edge* e1 = edges_to_remove[i].first;
      Edge* e2 = edges_to_remove[i].second;

      remove_interior_edge( *e1 );
      remove_interior_edge( *e2 );
    }

    // Removal of old vertices
    for (size_t i = 0; i < verts_to_remove.size(); ++i)
    {
      Vertex* v = verts_to_remove[i];
      remove_vertex( *v );
    }

    // Re-initialize facet-to-facet connectivity
    setup_facet_connectivity();

  } // Mesh::clean_double_triangle_edges()

  /*------------------------------------------------------------------
  | Algorithm to create a single quad layer for a connected list of
  | advancing front edges that start with v_start and end with v_end
  ------------------------------------------------------------------*/
  bool add_quad_layer(Vertex*& v_start_in, Vertex*& v_end_in, 
                      double height)
  {
    // Find closest vertices to given input vertices
    Vertex* v_start     = nullptr;
    Vertex* v_end       = nullptr;
    double d2_start_min = 1.0E+10;
    double d2_end_min   = 1.0E+10;

    for ( const auto& v : verts_ )
    {
      double d2_start = (v_start_in->xy() - v->xy()).length_squared();
      double d2_end   = (v_end_in->xy() - v->xy()).length_squared();

      if (d2_start < d2_start_min)
      {
        v_start = v.get();
        d2_start_min = d2_start;
      }

      if (d2_end < d2_end_min)
      {
        v_end = v.get();
        d2_end_min = d2_end;
      }
    }

    if (!v_start || !v_end)
    {
      LOG(ERROR) << 
      "Failed to create quad layer for mesh " << mesh_id_ << ", "
      "due to invalid provided starting or ending vertices.";
      return false;
    }

    // Get advancing front edges adjacent to input vertices
    Edge* e_start = front_.get_edge(*v_start, 1); 
    Edge* e_end   = front_.get_edge(*v_end, 2); 

    if ( !e_start || !e_end )
    {
      LOG(ERROR) << 
      "During the generation of a quad layer for mesh " << mesh_id_ << 
      " it was not possible to locate an advancing front edge that " <<
      "is adjacent to the given input vertex, that was provided for " <<
      "the layer generation. Thus, the process is aborted.";
      return false;
    }

    // Check if given front segments can be traversed and if closed
    bool is_closed = (v_start == v_end);

    if ( !front_.is_traversable(*e_start, *e_end) )
    {
      LOG(ERROR) << 
      "During the generation of a quad layer for mesh " << mesh_id_ << 
      " it was not possible to traverse the advancing front with "
      "the given input vertex. Thus, the process is aborted.";
      return false;
    }

    // For closed quad layers, try not to start at sharp angle edges
    if ( is_closed )
    {
      const Vec2d& v1 = e_end->v1().xy();
      const Vec2d& v2 = e_end->v2().xy();
      const Vec2d& v3 = e_start->v2().xy();

      const double ang = angle(v1-v2, v3-v2);

      Edge *e_next = e_start->get_next_edge();

      if ( e_next && ang <= CONSTANTS.quad_layer_angle() )
      {
        e_end = e_start;
        e_start = e_next; 
      }
    }

    // Triangulate front edges with critical angles
    //prepare_quad_layer_front(e_start, e_end, height);

    // Create the quad layer structure, which keeps track of the target
    // vertex coordinates, that are projected from the base vertex 
    // coordinates
    QuadLayer quad_layer { e_start, e_end, is_closed, height };
    quad_layer.smooth_heights( *domain_ );
    quad_layer.setup_vertex_projection( verts_, front_, bdry_edges_ );

    // For each base edge in the quad layer, try to create a quad
    // element with its given projected coordinates
    create_quad_layer_elements( quad_layer );

    // Triangulate the quad layer based edges, where the generation
    // of quads did not succeed
    finish_quad_layer( quad_layer );

    // Remove deleted entities
    clear_waste();

    // Set pointers to new start and ending vertices
    int i = 0;
    int n = quad_layer.n_bases();

    do 
    {
      v_start_in = quad_layer.p1()[i];

      if ( is_closed )
        v_end_in = v_start_in;
      else
        v_end_in = quad_layer.p2()[MOD(i-1,n)]; 

      ++i;

    } while ( !(v_start_in->on_front()) && !(v_end_in->on_front()) );

    return true;

  } // add_quad_layer()

  /*------------------------------------------------------------------
  | For each QuadProjection, create a triangle with its base 
  | vertices (b1,b2) and a vertex p1, which is either located in 
  | the vicinity of the base edge or which is otherwise generated at 
  | the projected coordinate of the base vertex b1
  |   
  |           p1            p2
  |          x-------------x-------------
  |          | \           | \          |
  |          |   \         |   \        |
  |          |     \       |     \      |
  |          |       \     |       \    |
  |          |         \   |         \  |
  |          |    base   \ |           \|
  | ---------x-------------x------------x-------
  |           b1            b2
  |   
  ------------------------------------------------------------------*/
  void create_quad_layer_elements(QuadLayer& quad_layer)
  {
    auto& b1        = quad_layer.b1();
    auto& b2        = quad_layer.b2();

    auto& p1        = quad_layer.p1();
    auto& p2        = quad_layer.p2();

    auto& p1_xy     = quad_layer.p1_xy();
    auto& p2_xy     = quad_layer.p2_xy();

    auto& heights   = quad_layer.heights();
    auto& bases     = quad_layer.bases();

    int  n_bases   = quad_layer.n_bases();

    for ( int i = 0; i < n_bases; ++i )
    {
      // Search radius for vertices in the vicinity of the 
      // projected coordinates
      const double r = CONSTANTS.quad_layer_range() * heights[i];

      // Create first triangle (b1,b2,p1)
      Edge* base = bases[i];

      if (!base->in_container())
        continue;

      Triangle* t1 = add_quad_layer_triangle(base, p1_xy[i], r);

      if ( t1 ) 
        p1[i] = &(t1->v3());
      else
        continue;

      // Create second triangle (p1,b2,p2)
      base = front_.get_edge( *p1[i], *b2[i] );

      if ( !base ) 
        continue;

      Triangle* t2 = add_quad_layer_triangle(base, p2_xy[i], r);

      if ( t2 )
        p2[i] = &(t2->v3());
      else
        continue;

      // Merge both triangles t1 & t2 to a quad
      // --> First remove the interior edge between these triangles
      Edge* e_rem = intr_edges_.get_edge( *b2[i], *p1[i] );

      if ( e_rem ) 
        remove_interior_edge( *e_rem );
      else
        continue;

      // Remove old triangular elements
      remove_triangle( *t1 );
      remove_triangle( *t2 );
      t1 = nullptr;
      t2 = nullptr;

      // Create new quadrilateral element
      Quad& q_new = add_quad( *b1[i], *b2[i], *p2[i], *p1[i] );
      q_new.is_active( true );

    }

  } // Mesh::create_quad_layer_elements() 
  
  /*------------------------------------------------------------------
  | This is a helper function for the generation of a triangle 
  | during the quad layer generation
  ------------------------------------------------------------------*/
  Triangle* add_quad_layer_triangle(Edge* base, 
                                    const Vec2d& xy, double r)
  {
    TriVector new_tris {};

    Triangle* tri = nullptr;

    Vertex& v1 = base->v1();
    Vertex& v2 = base->v2();

    // Look for vertices in the vicinity of the projected coordinate 
    // xy, which could be used to construct a new triangle
    VertexVector vertex_candidates = find_local_vertices(xy, r);

    // Create potential triangles with all found vertices
    check_vertex_candidates(vertex_candidates, *base, new_tris);

    // If potential triangles have been found, use the best one
    tri = choose_best_triangle(new_tris, *base);

    // If no proper triangle has been found, create a new vertex
    // at the projected position xy and create a new triangle 
    // with this new vertex
    if ( !tri )
    {
      Vertex& v_new = add_vertex( xy );
      tri = &( add_triangle(v1, v2, v_new) );

      // If the created entities are invalid, clean up
      if ( remove_if_invalid(v_new, *tri) )
      {
        tri = nullptr;
      }
      // Otherwise, update the advancing front with the new triangle
      else
      {
        update_front( *base, v_new, *tri );
        v_new.is_fixed( true );
      }
    }

    return tri;

  } // Mesh::add_quad_layer_triangle()

  /*------------------------------------------------------------------
  | In some cases, gaps might be formed during the previous quad layer
  | generation steps. In this function, these gaps are closed with 
  | triangular elements.
  |
  |              p1[i]
  |      v      x 
  |     x       :
  |             :
  |  p2[i-1]    :  
  |   x.........x-------------x
  |             | b1[i]        b2[i]
  |             |           
  |             |
  |             |
  |             x
  |               
  ------------------------------------------------------------------*/
  void finish_quad_layer(QuadLayer& quad_layer)
  {
    auto& b1        = quad_layer.b1();

    auto& p1        = quad_layer.p1();
    auto& p2        = quad_layer.p2();

    int  n_bases   = quad_layer.n_bases();

    for ( int i = 1; i < n_bases; ++i )
    {
      if ( !p1[i] || !p2[i-1] || p1[i] == p2[i-1] )
        continue;

      Vertex& a = *p2[i-1];
      Vertex& b = *b1[i];
      Vertex& c = *p1[i];

      const Vec2d l1 = a.xy()-b.xy();
      const Vec2d l2 = c.xy()-b.xy();
      const double alpha = angle(l1,l2);

      // Don't add a new vertex and instead only connect (a,b,c)
      if ( alpha <= CONSTANTS.quad_layer_angle() )
      {
        Triangle* t_new = &( add_triangle(a, b, c) );

        if ( !remove_if_invalid(*t_new) )
        {
          Edge* base = front_.get_edge( b, c );
          update_front( *base, a, *t_new );
        }
      }
      // Create new vertex and then generate two triangles
      else
      {
        const Vec2d v_xy = b.xy() + l1 + l2;

        Vertex& v_new = add_vertex( v_xy );

        Triangle* t1_new = &( add_triangle(a, b, v_new) );
        Triangle* t2_new = &( add_triangle(b, c, v_new) );

        if ( !remove_if_invalid(v_new, *t1_new, *t2_new) )
        {
          Edge* base = nullptr;

          base = front_.get_edge( a, b );
          update_front( *base, v_new, *t1_new );

          base = front_.get_edge( b, c );
          update_front( *base, v_new, *t2_new );

          v_new.is_fixed( true );
        }
      }
    }

  } // Mesh::finish_quad_layer()

  /*------------------------------------------------------------------
  | Triangulate front edges with critical angles
  | 
  |                      p2
  |                      x
  |              v       |
  |             o        |
  |               .      | e_next 
  |                 .    |
  |                   .  |
  |      p1             .|
  |     x----------------x
  |           e_cur       c
  | 
  ------------------------------------------------------------------*/
  void prepare_quad_layer_front(Edge*& e_start, Edge*& e_end, 
                                const double height)
  {
    Edge* e_cur  = e_start;

    Edge* e_last = e_end->get_next_edge();
    Edge* e_prev = e_start->get_prev_edge();

    int edge_count = 0;
    int n_edges = static_cast<int>(front_.size());

    do 
    {
      Edge* e_next = e_cur->get_next_edge();

      ASSERT( e_next, "INVALID DATA STRUCTURE" );
      if ( !e_next )
        break;

      Vertex& p1 = e_cur->v1();
      Vertex&  c = e_cur->v2();
      Vertex& p2 = e_next->v2();

      const Vec2d d1 = p1.xy() - c.xy();
      const Vec2d d2 = p2.xy() - c.xy();

      const double alpha = angle(d1,d2);
      const double delta = ( p2.xy() - p1.xy() ).length();

      const double l = 0.5 * ( e_cur->length() + e_next->length() );
      const double h = MIN( l, height );

      if (   !( is_lefton(p1.xy(), p2.xy(), c.xy()) ) 
          && delta <= h * CONSTANTS.quad_layer_factor() 
          && alpha <= CONSTANTS.quad_layer_angle() )
      {
        Edge* e_buf = e_next->get_next_edge();

        if ( e_next == e_last )
          e_last = e_buf;

        if ( delta <= h )
        {
          const Vec2d v_xy = 0.5 * (p1.xy() + p2.xy());
          Vertex& v_new = add_vertex( v_xy );

          Triangle& t1_new = add_triangle(p1, c, v_new);
          Triangle& t2_new = add_triangle(c, p2, v_new);

          if ( !remove_if_invalid(v_new, t1_new, t2_new) )
          {
            update_front( *e_cur, v_new, t1_new );
            update_front( *e_next, v_new, t2_new );

            v_new.is_fixed( true );
          }

        }
        else
        {
          Triangle& t_new = add_triangle(p1, c, p2);

          if ( !remove_if_invalid(t_new) )
            update_front( *e_cur, p2, t_new );
        }

        e_cur = e_buf;
        ++edge_count;
      }
      else
      {
        e_cur = e_next;
        ++edge_count;
      }

    } while(  ( e_cur )
           && ( edge_count < n_edges )
           && ( e_cur != e_last ) );

    e_start = e_prev->get_next_edge();
    e_end   = e_last->get_prev_edge();

  } // Mesh::prepare_quad_layer_front()

  /*------------------------------------------------------------------
  | Collect all boundary / neighbor-mesh edges that will be used 
  | for the creation of the advancing front
  | Additionally, return also an array of booleans, which indicate
  | the orientation of the given edges.
  | Edges which result from neighbor mesh are oriented in the 
  | opposite direction. This can be used to identify these edges
  | from normal boundary edges (e.g. for the advancing front 
  | refinement).
  ------------------------------------------------------------------*/
  FrontInitData collect_front_edges()
  {
    FrontInitData front_data {};

    for ( const auto& boundary : *domain_ )
    {
      EdgeVector edges {};
      IntVector  markers {};
      BoolVector is_oriented {};

      for ( const auto& e : boundary->edges() )
      {
        EdgeVector nbr_edges = get_neighbor_mesh_edges(*e);

        if ( nbr_edges.size() > 0 )
        {
          for ( Edge* nbr_e : nbr_edges )
          {
            edges.push_back( nbr_e );
            is_oriented.push_back( false );
            markers.push_back( nbr_e->marker() );
          }
        }
        else
        {
          edges.push_back( e.get() ) ;
          is_oriented.push_back( true );
          markers.push_back( e->marker() );
        }
      }

      front_data.edges.push_back( edges );
      front_data.is_oriented.push_back( is_oriented );
      front_data.markers.push_back( markers );

    }

    return std::move(front_data); 

  } // Mesh::collect_front_edges()

  /*------------------------------------------------------------------
  | For a given edge <e>, search all boundary edges of a 
  | neighboring mesh, that are contained within <e>.
  | This function is used upon the initialization of the advancing
  | front.
  ------------------------------------------------------------------*/
  EdgeVector get_neighbor_mesh_edges(const Edge& e)
  {
    EdgeVector   nbr_edges {};

    const Vec2d& c = e.xy();
    double       r = CONSTANTS.edge_search_factor() * e.length();

    const Vec2d& v = e.v1().xy();
    const Vec2d& w = e.v2().xy();

    for ( Mesh* m : neighbor_meshes_ )
    {
      ASSERT( m, "Neighbor mesh data structure is invalid." );

      // Get edges in the vicinity of the given edge
      EdgeVector edges_in_vicinity = m->get_bdry_edges(c, r);

      // Get all edges, that are actually located on the segment
      // of the current domain boundary edge 
      for ( Edge* e_vicinity : edges_in_vicinity )
      {
        const Vec2d& p = e_vicinity->v1().xy();
        const Vec2d& q = e_vicinity->v2().xy();

        if ( in_on_segment(v,w,p) && in_on_segment(v,w,q) )
          nbr_edges.push_back( e_vicinity );
      }

      // In case edges have been found, sort them in 
      // ascending direction to the stating vertex
      std::sort(nbr_edges.begin(), nbr_edges.end(), 
      [v](Edge* e1, Edge* e2)
      {
        double delta_1 = (e1->xy() - v).length_squared();
        double delta_2 = (e2->xy() - v).length_squared();
        return (delta_1 < delta_2);
      });

      // If edges have been located, terminate the search
      // --> First come, first serve
      if ( nbr_edges.size() > 0 )
        break;
    }

    return std::move(nbr_edges);
  }



  /*------------------------------------------------------------------
  | Get edge 
  ------------------------------------------------------------------*/
  Edge* get_edge(const Vertex& v1, const Vertex& v2, bool dir=false)
  const 
  {
    Edge* found = nullptr;

    found = intr_edges_.get_edge(v1, v2, dir);

    if ( found != nullptr )
      return found;

    found = bdry_edges_.get_edge(v1, v2, dir);

    if ( found != nullptr )
      return found;

    found = front_.get_edge(v1, v2, dir);

    if ( found != nullptr )
      return found;

    return found;

  } // Mesh::get_edge()

  /*------------------------------------------------------------------
  | This function removes an entity and makes sure, that the 
  | removal succeeded
  ------------------------------------------------------------------*/
  inline void remove_vertex(Vertex& v)
  {
    bool removed = verts_.remove( v );
    ASSERT( removed, "Failed to remove vertex.");
    (void) removed;
  }

  inline void remove_triangle(Triangle& t)
  {
    bool removed = tris_.remove( t );
    ASSERT( removed, "Failed to remove triangle.");
    (void) removed;
  }

  inline void remove_quad(Quad& q)
  {
    bool removed = quads_.remove( q );
    ASSERT( removed, "Failed to remove quad.");
    (void) removed;
  }

  inline void remove_interior_edge(Edge& e)
  {
    bool removed = intr_edges_.remove( e );
    ASSERT( removed, "Failed to remove interior edge.");
    (void) removed;
  }

  inline void remove_boundary_edge(Edge& e)
  {
    bool removed = bdry_edges_.remove( e );
    ASSERT( removed, "Failed to remove interior edge.");
    (void) removed;
  }

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


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Domain*    domain_      {nullptr};
  int        mesh_id_     { 0 };
  int        elem_color_  { CONSTANTS.default_element_color() };

  bool       mesh_initialized_  { false };
  bool       mesh_completed_    { false };

  Vertices   verts_;
  Triangles  tris_;
  Quads      quads_;
  Front      front_;

  EdgeList   intr_edges_ { Orientation::NONE };
  EdgeList   bdry_edges_ { Orientation::NONE };

  double     mesh_area_ { 0.0 };

protected:
  MeshVector neighbor_meshes_ {};
  MeshVector merged_meshes_   {};


}; // Mesh


/*********************************************************************
* Print out the mesh to std::cout
*********************************************************************/
inline std::ostream& operator<<(std::ostream& os, const Mesh& mesh)
{
  os << "MESH " << mesh.id() << "\n";

  os << "VERTICES " << mesh.vertices().size() << "\n";
  for ( const auto& v_ptr : mesh.vertices() )
  {
    os << std::setprecision(5) << std::fixed 
              << v_ptr->xy().x << "," 
              << v_ptr->xy().y << "\n";
  }

  os << "INTERIOREDGES " << mesh.interior_edges().size() << "\n";
  for ( const auto& e_ptr : mesh.interior_edges() )
  {
    auto v1_index = e_ptr->v1().index();
    auto v2_index = e_ptr->v2().index();

    auto fl_index = ( e_ptr->facet_l() != nullptr ) 
                  ?   e_ptr->facet_l()->index()
                  :   -1;
    auto fr_index = ( e_ptr->facet_r() != nullptr ) 
                  ?   e_ptr->facet_r()->index()
                  :   -1;

    os << std::setprecision(0) << std::fixed 
      << std::setw(4) << v1_index << "," 
      << std::setw(4) << v2_index << ","
      << std::setw(4) << fl_index << ","
      << std::setw(4) << fr_index << "\n";
  }


  // During the output of all boundary edges, we track the 
  // occurence of interface edges to other meshes
  size_t n_interface_edges = 0;

  os << "BOUNDARYEDGES " << mesh.boundary_edges().size() << "\n";
  for ( const auto& e_ptr : mesh.boundary_edges() )
  {
    if (e_ptr->twin_edge() != nullptr)
      ++n_interface_edges;

    auto v1_index = e_ptr->v1().index();
    auto v2_index = e_ptr->v2().index();

    auto fl_index = ( e_ptr->facet_l() != nullptr ) 
                  ?   e_ptr->facet_l()->index()
                  :   -1;
    auto marker   = e_ptr->marker();

    os << std::setprecision(0) << std::fixed 
      << std::setw(4) << v1_index << "," 
      << std::setw(4) << v2_index << ","
      << std::setw(4) << fl_index << ","
      << std::setw(4) << marker << "\n";
  }

  os << "INTERFACEEDGES " << n_interface_edges << "\n";
  for ( const auto& e_ptr : mesh.boundary_edges() )
  {
    if (e_ptr->twin_edge() == nullptr)
      continue;

    auto v1_index = e_ptr->v1().index();
    auto v2_index = e_ptr->v2().index();

    auto fl_index = ( e_ptr->facet_l() != nullptr ) 
                  ?   e_ptr->facet_l()->index()
                  :   -1;

    auto e_twin = e_ptr->twin_edge();
    auto f_twin = e_twin->facet_l();
    auto fl_twin_index = ( f_twin != nullptr )
                       ?   f_twin->index()
                       :   -1;
    auto fl_twin_color = ( f_twin != nullptr )
                         ?   f_twin->color()
                         :   -1;

    os << std::setprecision(0) << std::fixed 
      << std::setw(4) << v1_index << "," 
      << std::setw(4) << v2_index << ","
      << std::setw(4) << fl_index << ","
      << std::setw(4) << fl_twin_index << ","
      << std::setw(4) << fl_twin_color << "\n";
  }

  os << "FRONT " << mesh.front().size() << "\n";
  for ( const auto& e_ptr : mesh.front() )
  {
    auto v1_index = e_ptr->v1().index();
    auto v2_index = e_ptr->v2().index();
    auto marker   = e_ptr->marker();

    os << std::setprecision(0) << std::fixed 
      << std::setw(4) << v1_index << "," 
      << std::setw(4) << v2_index << ","
      << std::setw(4) << marker << "\n";
  }

  os << "QUADS " << mesh.quads().size() << "\n";
  for ( const auto& q_ptr : mesh.quads() )
  {
    auto v1_index = q_ptr->v1().index();
    auto v2_index = q_ptr->v2().index();
    auto v3_index = q_ptr->v3().index();
    auto v4_index = q_ptr->v4().index();
    auto color    = q_ptr->color();

    os << std::setprecision(0) << std::fixed
      << std::setw(4) << v1_index << ","
      << std::setw(4) << v2_index << ","
      << std::setw(4) << v3_index << ","
      << std::setw(4) << v4_index << ","
      << std::setw(4) << color << "\n";
  }

  os << "TRIANGLES " << mesh.triangles().size() << "\n";
  for ( const auto& t_ptr : mesh.triangles() )
  {
    auto v1_index = t_ptr->v1().index();
    auto v2_index = t_ptr->v2().index();
    auto v3_index = t_ptr->v3().index();
    auto color    = t_ptr->color();

    os << std::setprecision(0) << std::fixed
      << std::setw(4) << v1_index << ","
      << std::setw(4) << v2_index << ","
      << std::setw(4) << v3_index << ","
      << std::setw(4) << color << "\n";
  }

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
       << mesh.domain().size_function( v_ptr->xy() ) << "\n";
  }

  return os;
} 

} // namespace TQAlgorithm
} // namespace TQMesh
