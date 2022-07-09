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
* This example covers the generation of a simple triangular mesh
*********************************************************************/
void run_example_3()
{
  // Define a variable size function
  UserSizeFunction f1 = [](const Vec2d& p) { return 0.5; };
  UserSizeFunction f2 = [](const Vec2d& p) { return 0.5; };
  UserSizeFunction f3 = [](const Vec2d& p) { return 0.5; };

  // First domain = major domain 
  Domain domain_1 { f1, 20.0 };
  Domain domain_2 { f2, 20.0 };
  Domain domain_3 { f3, 20.0 };

  Boundary&  b_1 = domain_1.add_exterior_boundary();
  Boundary&  b_2 = domain_2.add_exterior_boundary();
  Boundary&  b_3 = domain_3.add_exterior_boundary();


  // Build exterior boundary of domain 1
  Vertex& v1_1 = domain_1.add_vertex(  0.0,  0.0 );
  Vertex& v2_1 = domain_1.add_vertex(  1.0,  0.0 );
  Vertex& v3_1 = domain_1.add_vertex(  1.0,  0.5 );
  Vertex& v4_1 = domain_1.add_vertex(  1.0,  1.0 );
  Vertex& v5_1 = domain_1.add_vertex(  0.0,  1.0 );

  b_1.add_edge( v1_1, v2_1, 1 );
  b_1.add_edge( v2_1, v3_1, 1 );
  b_1.add_edge( v3_1, v4_1, 1 );
  b_1.add_edge( v4_1, v5_1, 1 );
  b_1.add_edge( v5_1, v1_1, 1 );


  // Build exterior boundary of domain 2
  Vertex& v1_2 = domain_2.add_vertex(  1.0,  0.0 );
  Vertex& v2_2 = domain_2.add_vertex(  1.5,  0.0 );
  Vertex& v3_2 = domain_2.add_vertex(  1.5,  0.5 );
  Vertex& v4_2 = domain_2.add_vertex(  1.0,  0.5 );

  b_2.add_edge( v1_2, v2_2, 2 );
  b_2.add_edge( v2_2, v3_2, 2 );
  b_2.add_edge( v3_2, v4_2, 2 );
  b_2.add_edge( v4_2, v1_2, 2 );


  // Build exterior boundary of domain 3
  Vertex& v1_3 = domain_3.add_vertex(  1.0,  0.5 );
  Vertex& v2_3 = domain_3.add_vertex(  1.5,  0.5 );
  Vertex& v3_3 = domain_3.add_vertex(  1.5,  1.0 );
  Vertex& v4_3 = domain_3.add_vertex(  1.0,  1.0 );

  b_3.add_edge( v1_3, v2_3, 3 );
  b_3.add_edge( v2_3, v3_3, 3 );
  b_3.add_edge( v3_3, v4_3, 3 );
  b_3.add_edge( v4_3, v1_3, 3 );

  // Create meshes
  Mesh mesh_1 { domain_1, 0, 0, 5.0 };
  mesh_1.init_advancing_front();
  mesh_1.triangulate();
  //mesh_1.refine_to_quads();

  Mesh mesh_2 { domain_2, 1, 1, 5.0 };
  mesh_2.add_neighbor_mesh( mesh_1 );
  mesh_2.init_advancing_front();
  mesh_2.triangulate();

  Mesh mesh_3 { domain_3, 2, 2, 5.0 };
  mesh_3.add_neighbor_mesh( mesh_1 );
  mesh_3.add_neighbor_mesh( mesh_2 );
  mesh_3.init_advancing_front();
  mesh_3.triangulate();

  LOG(INFO) << "MERGE MESH: ";
  mesh_1.merge_neighbor_mesh( mesh_2 );
  mesh_1.merge_neighbor_mesh( mesh_3 );
  mesh_1.refine_to_quads();

  /*------------------------------------------------------------------
  | Finally, the mesh is exportet to a file in VTU format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/aux/example_data/Example_3" };

  mesh_1.write_to_file( file_name, ExportType::txt );

} // run_example_3()
