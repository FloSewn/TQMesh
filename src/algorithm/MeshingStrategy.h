/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "ProgressBar.h"

#include "Vertex.h"
#include "Edge.h"
#include "Front.h"
#include "Domain.h"
#include "Mesh.h"
#include "FrontUpdate.h"

namespace TQMesh {
namespace TQAlgorithm {

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
      const Vertex& v1 = e_ptr->v1();
      const Vertex& v2 = e_ptr->v2();

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

}; // MeshingStrategy

} // namespace TQAlgorithm
} // namespace TQMesh
