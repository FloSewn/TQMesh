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
#include "Cleanup.h"

namespace MeshTests 
{
using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* Test mesh initialization
*********************************************************************/
void initialization()
{
  int mesh_id = 0;
  int element_color = 0;

  Mesh mesh { mesh_id, element_color };

  Vertex& v1 = mesh.add_vertex({0.0, 0.0});
  Vertex& v2 = mesh.add_vertex({1.0, 0.0});
  Vertex& v3 = mesh.add_vertex({1.0, 1.0});

  Triangle& t1 = mesh.add_triangle(v1, v2, v3);

  mesh.add_boundary_edge(v1, v2, 0);
  mesh.add_boundary_edge(v2, v3, 1);
  mesh.add_interior_edge(v3, v1);

  Cleanup::assign_mesh_indices(mesh);
  Cleanup::setup_vertex_connectivity(mesh);

  LOG(INFO) << mesh;

  CHECK( v1.index() == 0 );
  CHECK( v2.index() == 1 );
  CHECK( v3.index() == 2 );
  CHECK( t1.index() == 0 );

} // initialization()

/*********************************************************************
* Test cleanup of double quad edges
*********************************************************************/
void clear_double_quad_edges()
{
  int mesh_id = 0;
  int element_color = 0;

  Mesh mesh { mesh_id, element_color };

  Vertex& v1 = mesh.add_vertex({3.0, 0.0});
  Vertex& v2 = mesh.add_vertex({6.0, 0.0});
  Vertex& v3 = mesh.add_vertex({0.0, 3.0});
  Vertex& v4 = mesh.add_vertex({3.0, 3.0});
  Vertex& v5 = mesh.add_vertex({6.0, 3.0});
  Vertex& v6 = mesh.add_vertex({4.0, 4.0});
  Vertex& v7 = mesh.add_vertex({3.0, 6.0});
  Vertex& v8 = mesh.add_vertex({6.0, 6.0});

  Triangle& t1 = mesh.add_triangle(v3, v4, v7);
  Quad& q1 = mesh.add_quad(v1, v2, v5, v4);
  Quad& q2 = mesh.add_quad(v4, v5, v8, v6);
  Quad& q3 = mesh.add_quad(v4, v6, v8, v7);

  mesh.add_interior_edge(v4, v5);
  Edge& intr_edge_4_6 = mesh.add_interior_edge(v4, v6);
  mesh.add_interior_edge(v8, v6);
  mesh.add_interior_edge(v7, v4);

  mesh.add_boundary_edge(v1, v2, 1);
  mesh.add_boundary_edge(v2, v5, 2);
  mesh.add_boundary_edge(v5, v8, 2);
  mesh.add_boundary_edge(v8, v7, 3);
  mesh.add_boundary_edge(v7, v3, 3);
  mesh.add_boundary_edge(v3, v4, 4);
  mesh.add_boundary_edge(v4, v1, 4);

  CHECK( mesh.n_elements() == 4 );
  CHECK( mesh.n_quads() == 3 );
  CHECK( mesh.n_triangles() == 1 );
  CHECK( v4.facets().size() == 4 );
  CHECK( v6.facets().size() == 2 );

  CHECK( Cleanup::check_mesh_validity(mesh) );

  Cleanup::assign_mesh_indices(mesh);
  Cleanup::setup_vertex_connectivity(mesh);
  Cleanup::setup_facet_connectivity(mesh);

  CHECK( q1.index() == 0 );
  CHECK( q2.index() == 1 );
  CHECK( q3.index() == 2 );
  CHECK( t1.index() == 3 );

  CHECK( intr_edge_4_6.facet_l() == &q3 );
  CHECK( intr_edge_4_6.facet_r() == &q2 );

  CHECK( v4.vertices().size() == 5 );


  Cleanup::clear_double_quad_edges(mesh, false);

  CHECK( mesh.n_elements() == 3 );
  CHECK( mesh.n_quads() == 2 );
  CHECK( mesh.n_triangles() == 1 );

  LOG_PROPERTIES.set_info_header( "" );
  LOG(INFO) << mesh;
  LOG_PROPERTIES.set_info_header( "  " );

} // clear_double_quad_edges()


/*********************************************************************
* Test cleanup of double quad edges
*********************************************************************/
void clear_double_triangle_edges()
{
  int mesh_id = 0;
  int element_color = 0;

  Mesh mesh { mesh_id, element_color };

  Vertex& v1 = mesh.add_vertex({3.0, 0.0});
  Vertex& v2 = mesh.add_vertex({6.0, 0.0});
  Vertex& v3 = mesh.add_vertex({0.0, 3.0});
  Vertex& v4 = mesh.add_vertex({3.0, 3.0});
  Vertex& v5 = mesh.add_vertex({6.0, 3.0});
  Vertex& v6 = mesh.add_vertex({4.0, 4.0});
  Vertex& v7 = mesh.add_vertex({3.0, 6.0});

  Triangle& t1 = mesh.add_triangle(v3, v4, v7);
  Triangle& t2 = mesh.add_triangle(v4, v5, v6);

  Quad& q1 = mesh.add_quad(v1, v2, v5, v4);
  Quad& q2 = mesh.add_quad(v4, v6, v5, v7);

  mesh.add_interior_edge(v4, v5);
  mesh.add_interior_edge(v4, v6);
  mesh.add_interior_edge(v5, v6);
  mesh.add_interior_edge(v7, v4);

  mesh.add_boundary_edge(v1, v2, 1);
  mesh.add_boundary_edge(v2, v5, 2);
  mesh.add_boundary_edge(v5, v7, 2);
  mesh.add_boundary_edge(v7, v3, 3);
  mesh.add_boundary_edge(v3, v4, 4);
  mesh.add_boundary_edge(v4, v1, 4);

  CHECK( mesh.n_elements() == 4 );
  CHECK( mesh.n_quads() == 2 );
  CHECK( mesh.n_triangles() == 2 );
  CHECK( v4.facets().size() == 4 );
  CHECK( v6.facets().size() == 2 );

  CHECK( Cleanup::check_mesh_validity(mesh) );

  Cleanup::assign_mesh_indices(mesh);
  Cleanup::setup_vertex_connectivity(mesh);
  Cleanup::setup_facet_connectivity(mesh);

  Cleanup::clear_double_triangle_edges(mesh, false);

  LOG_PROPERTIES.set_info_header( "" );
  LOG(INFO) << mesh;
  LOG_PROPERTIES.set_info_header( "  " );


} // clear_double_triangle_edges()

/*********************************************************************
* Test mesh initialization
*********************************************************************
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

} // initialization() */

/*********************************************************************
* Test Mesh::triangulate()
*********************************************************************
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
*********************************************************************
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
  adjust_logging_output_stream("MeshTests.initialization.log");
  MeshTests::initialization();

  adjust_logging_output_stream("MeshTests.clear_double_quad_edges.log");
  MeshTests::clear_double_quad_edges();

  adjust_logging_output_stream("MeshTests.clear_double_triangle_edges.log");
  MeshTests::clear_double_triangle_edges();

  //MeshTests::triangulate();
  //MeshTests::pave();

  // Reset debug logging ostream
  adjust_logging_output_stream("COUT");

} // run_tests_Mesh()
