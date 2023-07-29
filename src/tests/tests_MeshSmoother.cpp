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

#include "tests.h"

#include "VecND.h"
#include "Testing.h"

#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Mesh.h"
#include "MeshSmoother.h"

namespace MeshSmootherTests 
{
using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* Test MeshSmoother::smooth() for a pure triangle mesh
*********************************************************************
void tri_mesh()
{
  // Log debug messages to specified output-file
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/test_data/MeshSmootherTests.tri_mesh.log" };
  LOG_PROPERTIES.set_info_ostream( TO_FILE, file_name );
  LOG_PROPERTIES.set_debug_ostream( TO_FILE, file_name );

  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 1.5; 
  };

  Domain domain   { f, 50.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0, 0.1, 2.5 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0, 1.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0, 1.0 );
  Vertex& v4 = domain.add_vertex(  0.0,  5.0, 1.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 2 );
  b_ext.add_edge( v3, v4, 3 );
  b_ext.add_edge( v4, v1, 4 );

  // Create the mesh
  Mesh mesh { domain, 0, 0, 50.0 };
  MeshSmoother smoother {};

  mesh.init_advancing_front();
  mesh.triangulate();
  mesh.refine_to_quads();
  smoother.smooth(domain, mesh, 6, 0.5, 0.75, 0.95);

  // Export mesh
  file_name = source_dir + "/auxiliary/test_data/MeshSmootherTests.tri_mesh.txt";

  mesh.write_to_file( file_name, ExportType::txt );

} // tri_mesh()  */

} // namespace MeshSmootherTests


/*********************************************************************
* Run tests for: MeshSmoother.h
*********************************************************************/
void run_tests_MeshSmoother()
{
  //MeshSmootherTests::tri_mesh();

  // Reset debug logging ostream
  CppUtils::LOG_PROPERTIES.set_info_ostream( CppUtils::TO_COUT );
  CppUtils::LOG_PROPERTIES.set_debug_ostream( CppUtils::TO_COUT );

} // run_tests_MeshSmoother()
