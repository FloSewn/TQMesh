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
#include "MeshInitializer.h"
#include "FrontTriangulation.h"
#include "FrontQuadLayering.h"
//#include "FrontPaving.h"
#include "Smoother.h"

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

  LOG(INFO) << "\n" << mesh;

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

  mesh.add_triangle(v3, v4, v7);
  mesh.add_triangle(v4, v5, v6);

  mesh.add_quad(v1, v2, v5, v4);
  mesh.add_quad(v4, v6, v5, v7);

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

  LOG(INFO) << "\n" << mesh;

} // clear_double_triangle_edges()

/*********************************************************************
* Test Mesh::triangulate()
*********************************************************************/
void triangulate()
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { return 3.; };

  double quadtree_scale = 20.0;
  Domain domain   { f, quadtree_scale };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0 );
  Vertex& v4 = domain.add_vertex( 10.0,  5.0 );
  Vertex& v5 = domain.add_vertex( 10.0, 10.0 );
  Vertex& v6 = domain.add_vertex(  0.0,  5.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 2 );
  b_ext.add_edge( v3, v4, 2 );
  b_ext.add_edge( v4, v5, 3 );
  b_ext.add_edge( v5, v6, 4 );
  b_ext.add_edge( v6, v1, 4 );

  // Create the mesh
  MeshInitializer initializer {};

  Mesh mesh = initializer.create_empty_mesh(domain);
  initializer.prepare_mesh(mesh, domain);

  FrontTriangulation triangulation {mesh, domain};

  //CHECK( triangulation.generate_elements() );
    
  triangulation.n_elements(5);
  CHECK( triangulation.generate_elements() );

  triangulation.n_elements(4);
  CHECK( triangulation.generate_elements() );

  triangulation.n_elements(3);
  CHECK( triangulation.generate_elements() );

  // Assertions
  CHECK( ABS(mesh.area() - domain.area()) < 1.0E-07 );

  // Export mesh
  Cleanup::assign_size_function_to_vertices(mesh, domain);
  Cleanup::assign_mesh_indices(mesh);
  Cleanup::setup_vertex_connectivity(mesh);
  Cleanup::setup_facet_connectivity(mesh);
  LOG(DEBUG) << "\n" << mesh;

  // Check mesh stats
  CHECK( mesh.n_elements() == 12 );
  CHECK( mesh.n_vertices() == 14 );
  CHECK( mesh.n_interior_edges() == 11 );
  CHECK( mesh.n_boundary_edges() == 14 );

  // Check vertex-to-facet-connectivity  
  CHECK( mesh.vertices()[0].facets().size() == 4 );
  CHECK( mesh.vertices()[1].facets().size() == 1 );
  CHECK( mesh.vertices()[2].facets().size() == 3 );
  CHECK( mesh.vertices()[3].facets().size() == 1 );
  CHECK( mesh.vertices()[4].facets().size() == 4 );
  CHECK( mesh.vertices()[5].facets().size() == 3 );
  CHECK( mesh.vertices()[6].facets().size() == 4 );
  CHECK( mesh.vertices()[7].facets().size() == 1 );
  CHECK( mesh.vertices()[8].facets().size() == 3 );
  CHECK( mesh.vertices()[9].facets().size() == 1 );
  CHECK( mesh.vertices()[10].facets().size() == 3 );
  CHECK( mesh.vertices()[11].facets().size() == 3 );
  CHECK( mesh.vertices()[12].facets().size() == 4 );
  CHECK( mesh.vertices()[13].facets().size() == 1 );

  // Check edge-to-facet connectivity (boundary edges)
  CHECK( mesh.get_edge(mesh.vertices()[0], mesh.vertices()[1])->facet_l()->index() == 0 );
  CHECK( mesh.get_edge(mesh.vertices()[1], mesh.vertices()[2])->facet_l()->index() == 0 );
  CHECK( mesh.get_edge(mesh.vertices()[2], mesh.vertices()[3])->facet_l()->index() == 1 );
  CHECK( mesh.get_edge(mesh.vertices()[3], mesh.vertices()[4])->facet_l()->index() == 1 );
  CHECK( mesh.get_edge(mesh.vertices()[4], mesh.vertices()[5])->facet_l()->index() == 10 );
  CHECK( mesh.get_edge(mesh.vertices()[5], mesh.vertices()[6])->facet_l()->index() == 2 );
  CHECK( mesh.get_edge(mesh.vertices()[6], mesh.vertices()[7])->facet_l()->index() == 3 );
  CHECK( mesh.get_edge(mesh.vertices()[7], mesh.vertices()[8])->facet_l()->index() == 3 );
  CHECK( mesh.get_edge(mesh.vertices()[8], mesh.vertices()[9])->facet_l()->index() == 11 );
  CHECK( mesh.get_edge(mesh.vertices()[9], mesh.vertices()[10])->facet_l()->index() == 11 );
  CHECK( mesh.get_edge(mesh.vertices()[10], mesh.vertices()[11])->facet_l()->index() == 7 );
  CHECK( mesh.get_edge(mesh.vertices()[11], mesh.vertices()[12])->facet_l()->index() == 9 );
  CHECK( mesh.get_edge(mesh.vertices()[12], mesh.vertices()[13])->facet_l()->index() == 4 );
  CHECK( mesh.get_edge(mesh.vertices()[13], mesh.vertices()[0])->facet_l()->index() == 4 );

  // Check edge-to-facet connectivity (interior edges)
  CHECK( mesh.get_edge(mesh.vertices()[0], mesh.vertices()[2])->facet_l()->index() == 5 );
  CHECK( mesh.get_edge(mesh.vertices()[0], mesh.vertices()[2])->facet_r()->index() == 0 );

  CHECK( mesh.get_edge(mesh.vertices()[2], mesh.vertices()[4])->facet_l()->index() == 5 );
  CHECK( mesh.get_edge(mesh.vertices()[2], mesh.vertices()[4])->facet_r()->index() == 1 );

  CHECK( mesh.get_edge(mesh.vertices()[0], mesh.vertices()[4])->facet_l()->index() == 6 );
  CHECK( mesh.get_edge(mesh.vertices()[0], mesh.vertices()[4])->facet_r()->index() == 5 );

  CHECK( mesh.get_edge(mesh.vertices()[0], mesh.vertices()[12])->facet_l()->index() == 6 );
  CHECK( mesh.get_edge(mesh.vertices()[0], mesh.vertices()[12])->facet_r()->index() == 4 );

  CHECK( mesh.get_edge(mesh.vertices()[4], mesh.vertices()[12])->facet_l()->index() == 10 );
  CHECK( mesh.get_edge(mesh.vertices()[4], mesh.vertices()[12])->facet_r()->index() == 6 );

  CHECK( mesh.get_edge(mesh.vertices()[5], mesh.vertices()[12])->facet_l()->index() == 10 );
  CHECK( mesh.get_edge(mesh.vertices()[5], mesh.vertices()[12])->facet_r()->index() == 9 );

  CHECK( mesh.get_edge(mesh.vertices()[5], mesh.vertices()[11])->facet_l()->index() == 9 );
  CHECK( mesh.get_edge(mesh.vertices()[5], mesh.vertices()[11])->facet_r()->index() == 2 );

  CHECK( mesh.get_edge(mesh.vertices()[6], mesh.vertices()[11])->facet_l()->index() == 7 );
  CHECK( mesh.get_edge(mesh.vertices()[6], mesh.vertices()[11])->facet_r()->index() == 2 );

  CHECK( mesh.get_edge(mesh.vertices()[6], mesh.vertices()[10])->facet_l()->index() == 8 );
  CHECK( mesh.get_edge(mesh.vertices()[6], mesh.vertices()[10])->facet_r()->index() == 7 );

  CHECK( mesh.get_edge(mesh.vertices()[8], mesh.vertices()[10])->facet_l()->index() == 11 );
  CHECK( mesh.get_edge(mesh.vertices()[8], mesh.vertices()[10])->facet_r()->index() == 8 );

  CHECK( mesh.get_edge(mesh.vertices()[8], mesh.vertices()[6])->facet_l()->index() == 8 );
  CHECK( mesh.get_edge(mesh.vertices()[8], mesh.vertices()[6])->facet_r()->index() == 3 );

} // triangulate() 


/*********************************************************************
* Test Mesh::quadlayer()
*********************************************************************/
void quad_layer()
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { return 0.2; };

  double quadtree_scale = 20.0;
  Domain domain   { f, quadtree_scale };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0 );
  Vertex& v4 = domain.add_vertex( 10.0,  5.0 );
  Vertex& v5 = domain.add_vertex( 10.0, 10.0 );
  Vertex& v6 = domain.add_vertex(  0.0,  5.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 2 );
  b_ext.add_edge( v3, v4, 2 );
  b_ext.add_edge( v4, v5, 3 );
  b_ext.add_edge( v5, v6, 4 );
  b_ext.add_edge( v6, v1, 4 );

  // Create the mesh
  MeshInitializer initializer {};

  Mesh mesh = initializer.create_empty_mesh(domain);
  initializer.prepare_mesh(mesh, domain);

  FrontQuadLayering quadlayering {mesh, domain};
  quadlayering.n_layers( 4 );
  quadlayering.first_height( 0.35 );
  quadlayering.growth_rate( 1.0 );
  quadlayering.starting_position( 0.0, 2.5 );
  quadlayering.ending_position( 7.5, 5.0 );

  CHECK( quadlayering.generate_elements() );

  /*
  FrontTriangulation triangulation {mesh, domain};
  triangulation.n_elements(0);
  CHECK( triangulation.generate_elements() );

  Smoother smoother {};
  smoother.smooth(domain, mesh, 2);
  */

  // Export mesh
  Cleanup::assign_size_function_to_vertices(mesh, domain);
  Cleanup::assign_mesh_indices(mesh);
  Cleanup::setup_vertex_connectivity(mesh);
  Cleanup::setup_facet_connectivity(mesh);
  LOG(DEBUG) << "\n" << mesh;

} // quad_layer()

/*********************************************************************
* Test Mesh::pave()
*********************************************************************
void pave()
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { return 0.5; };

  double quadtree_scale = 20.0;
  Domain domain   { f, quadtree_scale };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0 );
  Vertex& v4 = domain.add_vertex( 10.0,  5.0 );
  Vertex& v5 = domain.add_vertex( 10.0, 10.0 );
  Vertex& v6 = domain.add_vertex(  0.0,  5.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 2 );
  b_ext.add_edge( v3, v4, 2 );
  b_ext.add_edge( v4, v5, 3 );
  b_ext.add_edge( v5, v6, 4 );
  b_ext.add_edge( v6, v1, 4 );

  // Create the mesh
  MeshInitializer initializer {};

  Mesh mesh = initializer.create_empty_mesh(domain);
  initializer.prepare_mesh(mesh, domain);

  FrontPaving paving {mesh, domain};

  CHECK( paving.generate_elements(0) );
  //CHECK( paving.generate_elements(3) );
  //CHECK( paving.generate_elements(3) );

  // Assertions
  //CHECK( ABS(mesh.area() - domain.area()) < 1.0E-07 );

  // Export mesh
  Cleanup::assign_size_function_to_vertices(mesh, domain);
  Cleanup::assign_mesh_indices(mesh);
  Cleanup::setup_vertex_connectivity(mesh);
  Cleanup::setup_facet_connectivity(mesh);
  LOG(DEBUG) << "\n" << mesh;

} // pave()  */

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

  adjust_logging_output_stream("MeshTests.triangulate.log");
  MeshTests::triangulate();

  adjust_logging_output_stream("MeshTests.quad_layer.log");
  MeshTests::quad_layer();

  //adjust_logging_output_stream("MeshTests.pave.log");
  //MeshTests::pave();

  // Reset debug logging ostream
  adjust_logging_output_stream("COUT");

} // run_tests_Mesh()
