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
#include "MeshBuilder.h"
#include "RefinementStrategy.h"
#include "MeshCleanup.h"
#include "TriangulationStrategy.h"
#include "QuadLayerStrategy.h"
#include "SmoothingStrategy.h"

namespace MeshTests 
{
using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* Test mesh initialization
*********************************************************************/
void initialization()
{
  // Define dummy domain & mesh builder
  UserSizeFunction f = [](const Vec2d& p) { return 1.0; };
  Domain domain { f, 5.0 };
  MeshBuilder mesh_builder {};
  Mesh mesh = mesh_builder.create_empty_mesh(domain, 0, 0, 5.0);

  Vertex& v1 = mesh.add_vertex({0.0, 0.0});
  Vertex& v2 = mesh.add_vertex({1.0, 0.0});
  Vertex& v3 = mesh.add_vertex({1.0, 1.0});

  Triangle& t1 = mesh.add_triangle(v1, v2, v3);

  mesh.add_boundary_edge(v1, v2, 0);
  mesh.add_boundary_edge(v2, v3, 1);
  mesh.add_interior_edge(v3, v1);

  MeshCleanup::assign_mesh_indices(mesh);

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

  TestBuilder test_builder { "NormalStepAndSharpEdge", f};
  Domain& domain = test_builder.domain();

  // Create the mesh
  MeshBuilder mesh_builder {};

  Mesh mesh = mesh_builder.create_empty_mesh(domain);
  CHECK( mesh_builder.prepare_mesh(mesh, domain) );

  TriangulationStrategy triangulation {mesh, domain};

  //CHECK( triangulation.generate_elements() );
    
  triangulation.n_elements(5);
  CHECK( triangulation.generate_elements() );

  triangulation.n_elements(4);
  CHECK( triangulation.generate_elements() );

  triangulation.n_elements(3);
  CHECK( triangulation.generate_elements() );

  // Assertions
  CHECK( EQ(mesh.area(), domain.area(), 1E-08) );

  // Export mesh
  MeshCleanup::assign_size_function_to_vertices(mesh, domain);
  MeshCleanup::assign_mesh_indices(mesh);
  MeshCleanup::setup_facet_connectivity(mesh);
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

  TestBuilder test_builder { "SharpStepAndSharpEdge", f};
  Domain& domain = test_builder.domain();

  // Create the mesh
  MeshBuilder mesh_builder {};

  Mesh mesh = mesh_builder.create_empty_mesh(domain);
  CHECK( mesh_builder.prepare_mesh(mesh, domain) );

  QuadLayerStrategy quad_layer {mesh, domain};
  quad_layer.n_layers( 4 );
  quad_layer.first_height( 0.15 );
  quad_layer.growth_rate( 1.0 );
  quad_layer.starting_position( 0.0, 2.5 );
  quad_layer.ending_position( 7.5, 5.0 );

  CHECK( quad_layer.generate_elements() );

  TriangulationStrategy triangulation {mesh, domain};
  triangulation.n_elements(0);
  CHECK( triangulation.generate_elements() );

  MixedSmoothingStrategy smoother {mesh, domain};
  smoother.smooth(2);

  // Export mesh
  MeshCleanup::assign_size_function_to_vertices(mesh, domain);
  MeshCleanup::assign_mesh_indices(mesh);
  MeshCleanup::setup_facet_connectivity(mesh);
  LOG(DEBUG) << "\n" << mesh;

} // quad_layer()

/*********************************************************************
* Test exhaustive_search_triangulation()
*********************************************************************/
void exhaustive_search_triangulation()
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { return 0.5; };

  TestBuilder test_builder { "SharpStepAndSharpEdge", f};
  Domain& domain = test_builder.domain();

  // Create the mesh
  MeshBuilder mesh_builder {};
  
  Mesh mesh = mesh_builder.create_empty_mesh(domain);
  CHECK( mesh_builder.prepare_mesh(mesh, domain) );

  TriangulationStrategy triangulation {mesh, domain};
  triangulation.n_elements();
  CHECK( triangulation.generate_elements_exhaustive() );

  // Export mesh
  MeshCleanup::assign_size_function_to_vertices(mesh, domain);
  MeshCleanup::assign_mesh_indices(mesh);
  MeshCleanup::setup_facet_connectivity(mesh);
  LOG(DEBUG) << "\n" << mesh;

} // exhaustive_search_triangulation()


/*********************************************************************
* Test refinement to quads
*********************************************************************/
void refine_to_quads()
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { return 0.3; };

  double quadtree_scale = 50.0;
  Domain domain   { f, quadtree_scale };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0 );
  Vertex& v4 = domain.add_vertex( 10.0,  5.0, 0.1, 0.1 );
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

  // Create quad layers
  QuadLayerStrategy quad_layer {mesh, domain};
  quad_layer.n_layers( 3 );
  quad_layer.first_height( 0.35 );
  quad_layer.growth_rate( 1.0 );
  quad_layer.starting_position( 0.0, 2.5 );
  quad_layer.ending_position( 7.5, 5.0 );

  CHECK( quad_layer.generate_elements() );

  // Refinement
  QuadRefinementStrategy refinement{mesh, domain};
  CHECK( refinement.refine() );

  // Create triangulation
  TriangulationStrategy triangulation {mesh, domain};
  triangulation.n_elements(0);

  CHECK( triangulation.generate_elements() );

  MeshCleanup::clear_double_quad_edges(mesh);
  MeshCleanup::clear_double_triangle_edges(mesh);
  MeshCleanup::merge_degenerate_triangles(mesh);

  // Refinement
  CHECK( refinement.refine() );

  // Smooth grid
  MixedSmoothingStrategy smoother {mesh, domain};
  smoother.smooth(2);

  // Export mesh
  MeshCleanup::assign_size_function_to_vertices(mesh, domain);
  MeshCleanup::assign_mesh_indices(mesh);
  MeshCleanup::setup_facet_connectivity(mesh);
  LOG(DEBUG) << "\n" << mesh;

} // refine_to_quads()

/*********************************************************************
* Test triangle to quad merge operation
*********************************************************************/
void merge_triangles_to_quads()
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) { return 0.3; };

  TestBuilder test_builder { "TriangleSquareCircle", f};
  Domain& domain = test_builder.domain();

  // Create the mesh
  MeshBuilder mesh_builder {};

  Mesh mesh = mesh_builder.create_empty_mesh( domain );
  CHECK( mesh_builder.prepare_mesh(mesh, domain) );

  TriangulationStrategy triangulation {mesh, domain};

  CHECK( triangulation.generate_elements() );
  
  MeshCleanup::merge_triangles_to_quads(mesh);

  // Smooth grid
  MixedSmoothingStrategy smoother {mesh, domain};
  smoother.smooth(2);

  // Export mesh
  MeshCleanup::assign_size_function_to_vertices(mesh, domain);
  MeshCleanup::assign_mesh_indices(mesh);
  MeshCleanup::setup_facet_connectivity(mesh);
  LOG(DEBUG) << "\n" << mesh;

} // merge_triangles_to_quads()


/*********************************************************************
* Test very small refinement
*********************************************************************/
void small_refinement()
{
  UserSizeFunction f = [](const Vec2d& p) { return 0.5; };
   
  TestBuilder test_builder { "UnitCircle", f};
  Domain& domain = test_builder.domain();

  domain.add_fixed_vertex(0.0, 0.0, 0.0005, 1.0);

  // Create the mesh
  MeshBuilder mesh_builder {};

  Mesh mesh = mesh_builder.create_empty_mesh( domain );
  CHECK( mesh_builder.prepare_mesh(mesh, domain) );

  TriangulationStrategy triangulation {mesh, domain};

  CHECK( triangulation.generate_elements() );
  CHECK( EQ(mesh.area(), domain.area(), 1E-07) );

  // Export mesh
  MeshCleanup::assign_size_function_to_vertices(mesh, domain);
  MeshCleanup::assign_mesh_indices(mesh);
  MeshCleanup::setup_facet_connectivity(mesh);

  LOG(DEBUG) << "\n" << mesh;

} // small_refinement()

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
  MeshCleanup::assign_size_function_to_vertices(mesh, domain);
  MeshCleanup::assign_mesh_indices(mesh);
  MeshCleanup::setup_facet_connectivity(mesh);
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

  TriangulationStrategy triangulation {mesh, domain};

  success &= triangulation.generate_elements();
  CHECK( success );

  MixedSmoothingStrategy smoother {mesh, domain};
  smoother.smooth(2);

  success &= EQ(mesh.area(), domain.area(), 1E+16);
  CHECK( success );

  if ( !success )
    LOG(ERROR) << "triangulate_" << test_name;

  MeshCleanup::merge_degenerate_triangles(mesh);

  // Export mesh
  MeshCleanup::assign_size_function_to_vertices(mesh, domain);
  MeshCleanup::assign_mesh_indices(mesh);

  if ( success )
  {
    MeshCleanup::setup_facet_connectivity(mesh);
  }

  LOG(DEBUG) << "\n" << mesh;

} // triangulate_standard_tests()

/*********************************************************************
* Test CSV import
*********************************************************************/
void csv_import()
{
  // Define dummy domain & mesh builder
  UserSizeFunction f = [](const Vec2d& p) { return 1.0; };
  Domain domain { f };

  std::string b_ext_file { TQMESH_SOURCE_DIR };
  b_ext_file += "/auxiliary/test_data/ExteriorBoundary.csv";

  Boundary& b_ext = domain.add_exterior_boundary();
  b_ext.set_shape_from_csv(b_ext_file);


  std::string b_int_1_file { TQMESH_SOURCE_DIR };
  b_int_1_file += "/auxiliary/test_data/InteriorBoundary_1.csv";

  Boundary& b_int_1 = domain.add_interior_boundary();
  b_int_1.set_shape_from_csv(b_int_1_file);

  std::string b_int_2_file { TQMESH_SOURCE_DIR };
  b_int_2_file += "/auxiliary/test_data/InteriorBoundary_2.csv";

  Boundary& b_int_2 = domain.add_interior_boundary();
  b_int_2.set_shape_from_csv(b_int_2_file);
    
    
  // Create the mesh
  MeshBuilder mesh_builder {};
  Mesh mesh = mesh_builder.create_empty_mesh( domain );

  CHECK( mesh_builder.prepare_mesh(mesh, domain) );

  TriangulationStrategy triangulation {mesh, domain};

  CHECK( triangulation.generate_elements() );
  CHECK( EQ(mesh.area(), domain.area(), 1E-07) );

  // Smooth grid
  MixedSmoothingStrategy smoother {mesh, domain};
  smoother.smooth(2);

  // Export mesh
  MeshCleanup::assign_size_function_to_vertices(mesh, domain);
  MeshCleanup::assign_mesh_indices(mesh);
  MeshCleanup::setup_facet_connectivity(mesh);

  LOG(DEBUG) << "\n" << mesh;

} // csv_import()

/*********************************************************************
* Test CSV import
*********************************************************************/
void bad_csv_import()
{
  // Define dummy domain & mesh builder
  UserSizeFunction f_outer = [](const Vec2d& p) { return 0.5; };
  UserSizeFunction f_inner = [](const Vec2d& p) { return 0.5; };
  Domain domain_outer { f_outer };
  Domain domain_inner { f_inner };

  std::string bdry_file { TQMESH_SOURCE_DIR };
  bdry_file += "/auxiliary/test_data/BadCSVInput.csv";

  domain_outer.add_exterior_boundary().set_shape_circle(1, {0.77, 0.09}, 0.14, 60);
  domain_outer.add_interior_boundary().set_shape_from_csv(bdry_file);

  domain_inner.add_exterior_boundary().set_shape_from_csv(bdry_file);

  // Create the outer mesh
  MeshBuilder mesh_builder {};

  Mesh mesh_outer = mesh_builder.create_empty_mesh( domain_outer, 1 );

  CHECK( mesh_builder.prepare_mesh(mesh_outer, domain_outer) );


  TriangulationStrategy triangulation_outer {mesh_outer, domain_outer};
  CHECK( triangulation_outer.generate_elements() );
  CHECK( EQ(mesh_outer.area(), domain_outer.area(), 1E-07) );


  // Create the inner mesh
  Mesh mesh_inner = mesh_builder.create_empty_mesh( domain_inner, 2 );

  CHECK( mesh_builder.prepare_mesh(mesh_inner, domain_inner) );

  TriangulationStrategy triangulation_inner {mesh_inner, domain_inner};
  CHECK( triangulation_inner.generate_elements() );
  CHECK( EQ(mesh_inner.area(), domain_inner.area(), 1E-07) );


  // Smooth grid
  MixedSmoothingStrategy smoother_outer {mesh_outer, domain_outer};
  smoother_outer.smooth(2);

  // Export mesh
  MeshCleanup::assign_size_function_to_vertices(mesh_outer, domain_outer);
  MeshCleanup::assign_mesh_indices(mesh_outer);
  MeshCleanup::setup_facet_connectivity(mesh_outer);

  LOG(DEBUG) << "\n" << mesh_outer;
  LOG(DEBUG) << "\n" << mesh_inner;


} // bad_csv_import()


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

  adjust_logging_output_stream("MeshTests.quad_layer.log");
  MeshTests::quad_layer();

  adjust_logging_output_stream("MeshTests.small_refinement.log");
  MeshTests::small_refinement();

  adjust_logging_output_stream("MeshTests.exhaustive_search_triangulation.log");
  MeshTests::exhaustive_search_triangulation();

  adjust_logging_output_stream("MeshTests.merge_triangles_to_quads.log");
  MeshTests::merge_triangles_to_quads();

  std::vector<std::string> standard_tests {
    "UnitSquare", "UnitCircle", "RefinedTriangle", 
    "FixedVertices", 
    "TriangleSquareCircle",
    "NormalStepAndSharpEdge", "SharpStepAndSharpEdge",
    //"LakeSuperior",
  };

  for ( auto test_name : standard_tests )
  {
    adjust_logging_output_stream("MeshTests.triangulate_" + test_name + ".log");
    MeshTests::triangulate_standard_tests(test_name);
  }

  adjust_logging_output_stream("MeshTests.csv_import.log");
  MeshTests::csv_import();
   
  adjust_logging_output_stream("MeshTests.bad_csv_import.log");
  MeshTests::bad_csv_import();

  // Reset debug logging ostream
  adjust_logging_output_stream("COUT");

} // run_tests_Mesh()
