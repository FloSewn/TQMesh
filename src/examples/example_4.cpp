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
void run_example_4()
{
  // Define a variable size function
  UserSizeFunction f1 = [](const Vec2d& p) { return 1.2; };
  UserSizeFunction f2 = [](const Vec2d& p) { return 1.2; };

  // First domain = major domain 
  Domain domain_1 { f1, 20.0 };

  Boundary&  b_ext_1 = domain_1.add_exterior_boundary();
  Boundary&  b_int_1 = domain_1.add_interior_boundary();

  // Build exterior boundary of domain 1
  Vertex& v1 = domain_1.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain_1.add_vertex(  8.0,  0.0 );
  Vertex& v3 = domain_1.add_vertex(  8.0,  8.0 );
  Vertex& v4 = domain_1.add_vertex(  0.0,  8.0 );

  b_ext_1.add_edge( v1, v2, 1 );
  b_ext_1.add_edge( v2, v3, 1 );
  b_ext_1.add_edge( v3, v4, 1 );
  b_ext_1.add_edge( v4, v1, 1 );

  // Build interior boundary of domain 1
  // --> This will be the exterior boundary of domain 2
  Vertex& v5 = domain_1.add_vertex(  2.0,  2.0 );
  Vertex& v6 = domain_1.add_vertex(  2.0,  6.0 );
  Vertex& v7 = domain_1.add_vertex(  6.0,  6.0 );
  Vertex& v8 = domain_1.add_vertex(  6.0,  2.0 );

  b_int_1.add_edge( v5, v6, 2 ); // These markers will not be needed 
  b_int_1.add_edge( v6, v7, 2 ); // after everything is meshed
  b_int_1.add_edge( v7, v8, 2 );
  b_int_1.add_edge( v8, v5, 2 );

  // Create the mesh for domain 1
  Mesh mesh_1 { domain_1, 0, 0, 50.0 };
  mesh_1.init_advancing_front();

  mesh_1.triangulate();
  mesh_1.merge_triangles_to_quads();
  //mesh_1.refine_to_quads();


  // Second domain = sub-domain of domain_1
  Domain domain_2 { f2, 20.0 };
  Boundary&  b_ext_2 = domain_2.add_exterior_boundary();

  // Build exterior boundary of domain 2
  Vertex& v9  = domain_2.add_vertex(  2.0,  2.0 );
  Vertex& v10 = domain_2.add_vertex(  6.0,  2.0 );
  Vertex& v11 = domain_2.add_vertex(  6.0,  6.0 );
  Vertex& v12 = domain_2.add_vertex(  2.0,  6.0 );

  b_ext_2.add_edge(  v9, v10, 2 );
  b_ext_2.add_edge( v10, v11, 2 );
  b_ext_2.add_edge( v11, v12, 2 );
  b_ext_2.add_edge( v12,  v9, 2 );
  
  // Create the mesh for domain 2
  Mesh mesh_2 { domain_2, 1, 1, 50.0 };
  mesh_2.add_neighbor_mesh( mesh_1 );
  mesh_2.init_advancing_front();
  mesh_2.triangulate();

  LOG(INFO) << "MERGE MESH: ";
  mesh_1.merge_neighbor_mesh( mesh_2 );
  //mesh_1.refine_to_quads();

  /*------------------------------------------------------------------
  | Finally, the mesh is exportet to a file in VTU format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/aux/example_data/Example_4" };

  mesh_1.write_to_file( file_name, ExportType::txt );

} // run_example_4()
