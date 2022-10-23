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
* This example shows how to generate and merge several meshes with
* different colors and size function.
* 
*   x-----------------------x
*   |  Outer                |
*   |  mesh                 |
*   |       x-------x       |
*   |       | Inner |       |
*   |       | mesh  |       |
*   |       |       |       |
*   |       x-------x       |
*   |                       |
*   |                       |
*   x-----------------------x
*
*********************************************************************/
void run_example_5()
{
  /*------------------------------------------------------------------
  | Define the size function of the outer mesh
  ------------------------------------------------------------------*/
  UserSizeFunction f_outer = [](const Vec2d& p) 
  { 
    return 0.5;
  };

  Domain outer_domain   { f_outer };

  /*------------------------------------------------------------------
  | Exterior boundary of the outer mesh
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
  | Interior boundary of the outer mesh
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
  | Initialize the outer mesh
  ------------------------------------------------------------------*/
  int outer_mesh_id = 1;
  int outer_mesh_color = 1;
  Mesh outer_mesh { outer_domain, outer_mesh_id, outer_mesh_color };
  outer_mesh.init_advancing_front();

  /*------------------------------------------------------------------
  | Create two quad layers for the outer mesh
  ------------------------------------------------------------------*/
  outer_mesh.create_quad_layers(v0_out, v0_out, 2, 0.05, 1.5);

  /*------------------------------------------------------------------
  | Generate the remaining elements of the outer mesh
  ------------------------------------------------------------------*/
  outer_mesh.triangulate();

  /*------------------------------------------------------------------
  | Smooth the elements of the outer mesh
  ------------------------------------------------------------------*/
  Smoother smoother {};
  smoother.smooth(outer_domain, outer_mesh, 4);



  /*------------------------------------------------------------------
  | Define the size function of the inner mesh
  ------------------------------------------------------------------*/
  UserSizeFunction f_inner = [](const Vec2d& p) 
  { 
    return 0.02 + 0.5 * ( pow(p.x-2.5, 2) + pow(p.y-2.5, 2) );
  };

  Domain inner_domain   { f_inner };


  /*------------------------------------------------------------------
  | Exterior boundary of the inner mesh
  | This boundary overlaps with the interior boundary of the outer 
  | mesh. It is thus important, that it features the same number of 
  | edges and vertex coordinates - but it must be defined in 
  | the opposite direction.
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
  | Initialize the inner mesh
  ------------------------------------------------------------------*/
  int inner_mesh_id = 2;
  int inner_mesh_color = 2;
  Mesh inner_mesh { inner_domain, inner_mesh_id, inner_mesh_color };

  /*------------------------------------------------------------------
  | Before the initialization of the advancing front, we pass 
  | the outer mesh as a "neighbor" to the inner mesh  
  | In this way, the inner mesh's exterior boundary edges will be 
  | conforming to the outer mesh's interior boundary edges
  ------------------------------------------------------------------*/
  inner_mesh.add_neighbor_mesh( outer_mesh );

  inner_mesh.init_advancing_front();

  /*------------------------------------------------------------------
  | Generate the inner mesh
  ------------------------------------------------------------------*/
  inner_mesh.triangulate();

  /*------------------------------------------------------------------
  | Finally, merge both meshes. In this way, the outer mesh's 
  | elements will be added to the inner mesh, whereas the outer 
  | mesh won't change.
  ------------------------------------------------------------------*/
  inner_mesh.merge_neighbor_mesh( outer_mesh );

  /*------------------------------------------------------------------
  | Smooth the final mesh for four iterations
  ------------------------------------------------------------------*/
  smoother.smooth(inner_domain, inner_mesh, 4);

  /*------------------------------------------------------------------
  | Finally, the mesh is exportet to a file in VTU format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/example_data/Example_5" };
  LOG(INFO) << "Writing mesh output to: " << file_name << ".txt";

  inner_mesh.write_to_file( file_name, ExportType::txt );

} // run_example_5()
