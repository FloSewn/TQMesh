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

using namespace CppUtils;

/*********************************************************************
* 
*********************************************************************/
class MeshingStrategy 
{
public:
  using VertexVector   = std::vector<Vertex*>;
  using TriVector      = std::vector<Triangle*>;

  /*------------------------------------------------------------------
  | Constructor / Destructor 
  ------------------------------------------------------------------*/
  MeshingStrategy(Mesh& mesh, const Domain& domain)
  : mesh_ { mesh }
  , domain_ { domain }
  , front_update_ {mesh, domain, front_} 
  {}

  virtual ~MeshingStrategy() {}

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  Mesh& mesh() { return mesh_; }
  bool show_progress() const { return show_progress_; }

  /*------------------------------------------------------------------
  | Triangulate a given initialized mesh structure
  ------------------------------------------------------------------*/
  virtual bool generate_elements() = 0;

protected:

  /*------------------------------------------------------------------
  | Initialize the advancing front structure
  ------------------------------------------------------------------*/
  Edge* init_advancing_front(bool sort_edges=true)
  {
    front_.init_front(mesh_);
    Edge* base = front_.set_base_first();
    ASSERT( base, "MeshingStrategy::generate_elements(): "
      "Invalid advancing front structure.");

    if (sort_edges)
      front_.sort_edges( false );

    return base;
  }

  /*------------------------------------------------------------------
  | Remove invalid mesh edges (which stem from previous
  | meshing approaches)
  ------------------------------------------------------------------*/
  void remove_invalid_mesh_edges()
  {
    auto invalid_mesh_edges = mesh_.get_front_edges();
    for ( auto& e_ptr : invalid_mesh_edges )
      mesh_.remove_interior_edge( *e_ptr );
  }

  /*------------------------------------------------------------------
  | Prepare the mesh for the output - this consists mainly of the 
  | addition of remaining advancing front edges, such that the 
  | mesh process can be continued in a subsequent step
  ------------------------------------------------------------------*/
  void add_remaining_front_edges_to_mesh()
  {
    // Add remaining front edges to the mesh 
    for ( auto& e_ptr : front_.edges() )
    {
      Vertex& v1 = e_ptr->v1();
      Vertex& v2 = e_ptr->v2();

      v1.remove_property( VertexProperty::on_front );
      v2.remove_property( VertexProperty::on_front );

      if ( mesh_.interior_edges().get_edge(v1, v2) )
        continue;

      if ( mesh_.boundary_edges().get_edge(v1, v2) )
        continue;

      mesh_.add_interior_edge(e_ptr->v1(), e_ptr->v2());
    }
  }

  /*------------------------------------------------------------------
  | Update the progress bar
  ------------------------------------------------------------------*/
  void update_progress_bar()
  {
    if ( !show_progress_ )
      return;
    double state = std::ceil(100.0 * mesh_.area() / domain_.area());
    progress_bar_.update( static_cast<int>(state) );
    progress_bar_.show( LOG_PROPERTIES.get_ostream(INFO) );
  }

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Mesh&         mesh_;
  const Domain& domain_;
  FrontUpdate   front_update_;

  Front         front_ {};
  ProgressBar   progress_bar_ {};
  bool          show_progress_ { false };

}; // MeshingStrategy

} // namespace TQMesh
