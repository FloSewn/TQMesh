/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/

#include <iostream>
#include <cassert>

#include <TQMeshConfig.h>

#include "run_examples.h"

#include "Vec2.h"
#include "Log.h"

#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Mesh.h"
#include "Smoother.h"

using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
*
*********************************************************************/
void run_example_6()
{
  /*------------------------------------------------------------------
  | Define the size function
  ------------------------------------------------------------------*/
  UserSizeFunction f_outer = [](const Vec2d& p) 
  { 
    return 0.5;
  };

  Domain outer_domain   { f_outer };

  /*------------------------------------------------------------------
  | Exterior boundary
  ------------------------------------------------------------------*/
  Boundary&  b_outer_ext = outer_domain.add_exterior_boundary();

  Vertex& v0_out = outer_domain.add_vertex(  0.0,  0.0, 1.0, 1.0 );
  Vertex& v1_out = outer_domain.add_vertex(  5.0,  0.0, 1.0, 1.0 );
  Vertex& v2_out = outer_domain.add_vertex(  5.0,  5.0, 1.0, 1.0 );
  Vertex& v3_out = outer_domain.add_vertex(  0.0,  5.0, 1.0, 1.0 );

  b_outer_ext.add_edge( v0_out, v1_out, 1 );
  b_outer_ext.add_edge( v1_out, v2_out, 1 );
  b_outer_ext.add_edge( v2_out, v3_out, 1 );
  b_outer_ext.add_edge( v3_out, v0_out, 1 );

  /*------------------------------------------------------------------
  | Interior boundary
  ------------------------------------------------------------------*/
  Boundary&  b_outer_int = outer_domain.add_interior_boundary();

  Vertex& v4_out = outer_domain.add_vertex(  1.5,  1.5, 1.0, 1.0 );
  Vertex& v5_out = outer_domain.add_vertex(  3.5,  1.5, 1.0, 1.0 );
  Vertex& v6_out = outer_domain.add_vertex(  3.5,  3.5, 1.0, 1.0 );
  Vertex& v7_out = outer_domain.add_vertex(  1.5,  3.5, 1.0, 1.0 );

  b_outer_int.add_edge( v4_out, v7_out, 2 );
  b_outer_int.add_edge( v7_out, v6_out, 2 );
  b_outer_int.add_edge( v6_out, v5_out, 2 );
  b_outer_int.add_edge( v5_out, v4_out, 2 );

  /*------------------------------------------------------------------
  | Initialize the mesh
  ------------------------------------------------------------------*/
  int outer_mesh_id = 1;
  int outer_mesh_color = 1;
  Mesh outer_mesh { outer_domain, outer_mesh_id, outer_mesh_color };
  outer_mesh.init_advancing_front();

  /*------------------------------------------------------------------
  |
  ------------------------------------------------------------------*/
  outer_mesh.create_quad_layers(v0_out, v0_out, 2, 0.05, 1.5);

  /*------------------------------------------------------------------
  |
  ------------------------------------------------------------------*/
  outer_mesh.triangulate();


  /*------------------------------------------------------------------
  |
  ------------------------------------------------------------------*/
  Smoother smoother {};
  smoother.smooth(outer_domain, outer_mesh, 4);


  /*------------------------------------------------------------------
  | Define the size function
  ------------------------------------------------------------------*/
  UserSizeFunction f_inner = [](const Vec2d& p) 
  { 
    return 0.02 + 0.5 * ( pow(p.x-2.5, 2) + pow(p.y-2.5, 2) );
  };

  Domain inner_domain   { f_inner };


  /*------------------------------------------------------------------
  | Interior boundary
  ------------------------------------------------------------------*/
  Boundary&  b_inner = inner_domain.add_exterior_boundary();

  Vertex& v4_in = inner_domain.add_vertex(  1.5,  1.5, 1.0, 1.0 );
  Vertex& v5_in = inner_domain.add_vertex(  3.5,  1.5, 1.0, 1.0 );
  Vertex& v6_in = inner_domain.add_vertex(  3.5,  3.5, 1.0, 1.0 );
  Vertex& v7_in = inner_domain.add_vertex(  1.5,  3.5, 1.0, 1.0 );

  b_inner.add_edge( v4_in, v5_in, 3 );
  b_inner.add_edge( v5_in, v6_in, 3 );
  b_inner.add_edge( v6_in, v7_in, 3 );
  b_inner.add_edge( v7_in, v4_in, 3 );


  /*------------------------------------------------------------------
  | Initialize the mesh
  ------------------------------------------------------------------*/
  int inner_mesh_id = 2;
  int inner_mesh_color = 2;
  Mesh inner_mesh { inner_domain, inner_mesh_id, inner_mesh_color };

  inner_mesh.add_neighbor_mesh( outer_mesh );

  inner_mesh.init_advancing_front();

  /*------------------------------------------------------------------
  |
  ------------------------------------------------------------------*/
  inner_mesh.triangulate();

  /*------------------------------------------------------------------
  |
  ------------------------------------------------------------------*/
  inner_mesh.merge_neighbor_mesh( outer_mesh );

  /*------------------------------------------------------------------
  | Smooth the mesh for four iterations
  ------------------------------------------------------------------*/
  smoother.smooth(inner_domain, inner_mesh, 4);

  /*------------------------------------------------------------------
  | Finally, the mesh is exportet to a file in VTU format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/aux/example_data/Example_6" };

  inner_mesh.write_to_file( file_name, ExportType::txt );

} // run_example_6()
