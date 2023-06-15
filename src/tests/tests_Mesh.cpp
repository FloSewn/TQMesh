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
#include "Timer.h"
#include "Container.h"

#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Mesh.h"

namespace MeshTests 
{
using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* Test mesh initialization
*********************************************************************/
void initialization()
{
  // Log debug messages to specified output-file
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/test_data/MeshTests.initialization.log" };
  LOG_PROPERTIES.set_info_ostream( TO_FILE, file_name );
  LOG_PROPERTIES.set_debug_ostream( TO_FILE, file_name );

  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { return 1.0 + 0.15*sqrt(p.x*p.y); };

  Domain           domain   { f, 10.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();
  Boundary&  b_int = domain.add_interior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0 );
  Vertex& v4 = domain.add_vertex(  0.0,  5.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v1, 1 );

  // Build interior boundary
  Vertex& v5 = domain.add_vertex(  2.5,  2.0, 0.2);
  Vertex& v6 = domain.add_vertex(  2.0,  3.5 );
  Vertex& v7 = domain.add_vertex(  3.0,  2.5 );
  Vertex& v8 = domain.add_vertex(  3.0,  2.0 );

  b_int.add_edge( v5, v6, 2 );
  b_int.add_edge( v6, v7, 2 );
  b_int.add_edge( v7, v8, 2 );
  b_int.add_edge( v8, v5, 2 );

  // Create the mesh
  Mesh mesh { domain, 0, 0, 10.0 };

  mesh.init_advancing_front();

  CHECK( EQ(mesh.front().area(), domain.area()) );

} // initialization()

/*********************************************************************
* Test Mesh::triangulate()
*********************************************************************/
void triangulate()
{
  // Log debug messages to specified output-file
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/test_data/MeshTests.triangulate.log" };
  LOG_PROPERTIES.set_info_ostream( TO_FILE, file_name );
  LOG_PROPERTIES.set_debug_ostream( TO_FILE, file_name );

  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { return 0.5; };

  Domain domain   { f, 20.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0, 0.9 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0, 0.9 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0, 0.9 );
  Vertex& v4 = domain.add_vertex( 10.0,  5.0, 0.9 );
  Vertex& v5 = domain.add_vertex( 10.0, 10.0, 0.9 );
  Vertex& v6 = domain.add_vertex(  0.0,  5.0, 0.9 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v5, 1 );
  b_ext.add_edge( v5, v6, 1 );
  b_ext.add_edge( v6, v1, 1 );

  // Create the mesh
  Mesh mesh { domain, 0, 0, 50.0 };
  mesh.init_advancing_front();

  mesh.triangulate();

  // Assertions
  CHECK( EQ(mesh.area(), domain.area()) );


  // Export mesh
  file_name =  source_dir + "/auxiliary/test_data/MeshTests.triangulate.txt";

  mesh.write_to_file( file_name, ExportType::txt );

} // triangulate() */

/*********************************************************************
* Test Mesh::pave()
*********************************************************************/
void pave()
{
  // Log debug messages to specified output-file
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/test_data/MeshTests.pave.log" };
  LOG_PROPERTIES.set_info_ostream( TO_FILE, file_name );
  LOG_PROPERTIES.set_debug_ostream( TO_FILE, file_name );

  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { return 0.5; };

  Domain domain   { f, 20.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0, 0.9 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0, 0.9 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0, 0.9 );
  Vertex& v4 = domain.add_vertex( 10.0,  5.0, 0.9 );
  Vertex& v5 = domain.add_vertex( 10.0, 10.0, 0.9 );
  Vertex& v6 = domain.add_vertex(  0.0,  5.0, 0.9 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v5, 1 );
  b_ext.add_edge( v5, v6, 1 );
  b_ext.add_edge( v6, v1, 1 );

  // Create the mesh
  Mesh mesh { domain, 0, 0, 50.0 };

  mesh.init_advancing_front();
  mesh.pave();

  // Assertions
  CHECK( EQ(mesh.area(), domain.area()) );

  // Export the mesh
  file_name = source_dir + "/auxiliary/test_data/MeshTests.pave.txt";

  mesh.write_to_file( file_name, ExportType::txt );

} // pave() */

} // namespace MeshTests


/*********************************************************************
* Run tests for: Mesh.h
*********************************************************************/
void run_tests_Mesh()
{
  MeshTests::initialization();
  MeshTests::triangulate();
  MeshTests::pave();

  // Reset debug logging ostream
  CppUtils::LOG_PROPERTIES.set_info_ostream( CppUtils::TO_COUT );
  CppUtils::LOG_PROPERTIES.set_debug_ostream( CppUtils::TO_COUT );

} // run_tests_Mesh()
