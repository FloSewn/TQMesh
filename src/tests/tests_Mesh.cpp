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
#include "TestBuilder.h"

#include "VecND.h"
#include "Testing.h"
#include "Timer.h"
#include "Container.h"

#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Mesh.h"
#include "Cleanup.h"
#include "MeshBuilder.h"
#include "FrontTriangulation.h"
#include "FrontQuadLayering.h"
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
  MeshBuilder mesh_builder {};

  Mesh mesh = mesh_builder.create_empty_mesh(domain);
  CHECK( mesh_builder.prepare_mesh(mesh, domain) );

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
  MeshBuilder mesh_builder {};

  Mesh mesh = mesh_builder.create_empty_mesh(domain);
  CHECK( mesh_builder.prepare_mesh(mesh, domain) );

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
  MeshBuilder mesh_builder {};

  Mesh mesh = mesh_builder.create_empty_mesh(domain);
  CHECK( mesh_builder.prepare_mesh(mesh, domain) );

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
* Standard tests for Mesh::triangulate()
*********************************************************************/
void triangulate_standard_tests(const std::string& test_name)
{
  bool success = true;

  UserSizeFunction f = [](const Vec2d& p) { return 0.5; };
   
  TestBuilder test_builder { test_name, f};
  Domain& domain = test_builder.domain();

  // Create the mesh
  MeshBuilder mesh_builder {};

  Mesh mesh = mesh_builder.create_empty_mesh( domain );
  success &= mesh_builder.prepare_mesh(mesh, domain);
  CHECK( success );

  if ( !success )
    return;

  FrontTriangulation triangulation {mesh, domain};

  success &= triangulation.generate_elements();
  CHECK( success );

  success &= EQ(mesh.area(), domain.area());
  CHECK( success );

  if ( !success )
    LOG(ERROR) << "triangulate_" << test_name;

  Cleanup::merge_degenerate_triangles(mesh);

  // Export mesh
  Cleanup::assign_size_function_to_vertices(mesh, domain);
  Cleanup::assign_mesh_indices(mesh);

  if ( success )
  {
    Cleanup::setup_vertex_connectivity(mesh);
    Cleanup::setup_facet_connectivity(mesh);
  }

  LOG(DEBUG) << "\n" << mesh;

} // triangulate_standard_tests()



} // namespace MeshTests


/*********************************************************************
* Run tests for: Mesh.h
*********************************************************************/
void run_tests_Mesh()
{
  adjust_logging_output_stream("MeshTests.initialization.log");
  MeshTests::initialization();

  adjust_logging_output_stream("MeshTests.triangulate.log");
  MeshTests::triangulate();
  
  //adjust_logging_output_stream("MeshTests.pave.log");
  //MeshTests::pave();

  adjust_logging_output_stream("MeshTests.quad_layer.log");
  MeshTests::quad_layer();

  std::vector<std::string> standard_tests {
    "UnitSquare", "UnitCircle", "RefinedTriangle", 
    "FixedVertices", "TriangleSquareCircle",
    //"LakeSuperior",
  };

  for ( auto test_name : standard_tests )
  {
    adjust_logging_output_stream("MeshTests.triangulate_" + test_name + ".log");
    MeshTests::triangulate_standard_tests(test_name);
  }

  // Reset debug logging ostream
  adjust_logging_output_stream("COUT");

} // run_tests_Mesh()
