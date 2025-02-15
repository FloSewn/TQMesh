/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <TQMeshConfig.h>
#include "STLHeaders.h"
#include "CppUtils.h"

#include "tests.h"
#include "TestBuilder.h"

#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Mesh.h"
#include "MeshGenerator.h"
#include "MeshMerger.h"
#include "RefinementStrategy.h"
#include "MeshCleanup.h"
#include "MeshChecker.h"
#include "SmoothingStrategy.h"
#include "EntityChecks.h"

namespace MeshGeneratorTests 
{
using namespace CppUtils;
using namespace TQMesh;

/*********************************************************************
* Test mesh initialization for two meshes of different colors and 
* size functions
*
*               x-----------------------x
*               |  Outer                |
*               |  mesh                 |
*               |       x-------x       |
*               |       | Inner |       |
*               |       | mesh  |       |
*               |       |       |       |
*               |       x-------x       |
*               |                       |
*               |                       |
*               x-----------------------x
*
*
*********************************************************************/
void initialization()
{
  // Define a variable size function
  UserSizeFunction f_outer = [](const Vec2d& p) 
  { return 0.5; };

  TQMeshSetup::get_instance().set_quadtree_scale( 10.0 );

  // Define the outer domain
  Domain outer_domain { f_outer };

  // Define vertices of the exterior boundary
  Vertex& v1 = outer_domain.add_vertex(  0.0,  0.0 );
  Vertex& v2 = outer_domain.add_vertex(  5.0,  0.0 );
  Vertex& v3 = outer_domain.add_vertex(  5.0,  5.0 );
  Vertex& v4 = outer_domain.add_vertex(  0.0,  5.0 );
 
  // Define vertices of the interior boundary
  Vertex& v5 = outer_domain.add_vertex(  2.5,  2.0, 0.2);
  Vertex& v6 = outer_domain.add_vertex(  2.0,  3.5 );
  Vertex& v7 = outer_domain.add_vertex(  3.0,  2.5 );
  Vertex& v8 = outer_domain.add_vertex(  3.0,  2.0 );

  // Build exterior boundary
  int exterior_edge_color = 1;
  Boundary& outer_exterior_bdry = outer_domain.add_exterior_boundary();
  outer_exterior_bdry.add_edge( v1, v2, exterior_edge_color );
  outer_exterior_bdry.add_edge( v2, v3, exterior_edge_color );
  outer_exterior_bdry.add_edge( v3, v4, exterior_edge_color );
  outer_exterior_bdry.add_edge( v4, v1, exterior_edge_color );

  // Build interior boundary
  int interior_edge_color = 2;
  Boundary& outer_interior_bdry = outer_domain.add_interior_boundary();
  outer_interior_bdry.add_edge( v5, v6, interior_edge_color );
  outer_interior_bdry.add_edge( v6, v7, interior_edge_color );
  outer_interior_bdry.add_edge( v7, v8, interior_edge_color );
  outer_interior_bdry.add_edge( v8, v5, interior_edge_color );

  // Setup the generator
  MeshGenerator generator {};
  generator.new_mesh( outer_domain );
    
    
} // initialization()

/*********************************************************************
* Test mesh initialization for two meshes of different colors and 
* size functions
*********************************************************************/
void mesh_initializer()
{
  // Define size functions
  //UserSizeFunction f_1 = [](const Vec2d& p) { return 0.3 + 0.4*p.x; };
  UserSizeFunction f_1 = [](const Vec2d& p) { return 2.5; };
  //UserSizeFunction f_2 = [](const Vec2d& p) { return 2.5 - MIN((p.x-5.0) * 0.3, 2.0); };
  UserSizeFunction f_2 = [](const Vec2d& p) { return 2.5; };

  TQMeshSetup::get_instance().set_quadtree_scale( 20.0 );

  // Define domains
  Domain domain_1 { f_1 };
  Domain domain_2 { f_2 };

  // Define boundary vertices 
  Vertex& v1_1 = domain_1.add_vertex(  0.0,  0.0 );
  Vertex& v2_1 = domain_1.add_vertex(  5.0,  0.0 );
  Vertex& v3_1 = domain_1.add_vertex(  5.0,  5.0 );
  Vertex& v4_1 = domain_1.add_vertex(  0.0,  5.0 );

  Vertex& v1_2 = domain_2.add_vertex(  5.0,  0.0 );
  Vertex& v2_2 = domain_2.add_vertex( 10.0,  0.0 );
  Vertex& v3_2 = domain_2.add_vertex( 10.0,  5.0 );
  Vertex& v4_2 = domain_2.add_vertex(  5.0,  5.0 );

  // Build domain boundaries
  int edge_color = 1;
  Boundary& bdry_1 = domain_1.add_exterior_boundary();
  Boundary& bdry_2 = domain_2.add_exterior_boundary();
  
  bdry_1.add_edge( v1_1, v2_1, edge_color );
  bdry_1.add_edge( v2_1, v3_1, edge_color );
  bdry_1.add_edge( v3_1, v4_1, edge_color );
  bdry_1.add_edge( v4_1, v1_1, edge_color );

  bdry_2.add_edge( v1_2, v2_2, edge_color );
  bdry_2.add_edge( v2_2, v3_2, edge_color );
  bdry_2.add_edge( v3_2, v4_2, edge_color );
  bdry_2.add_edge( v4_2, v1_2, edge_color );

  // Generate mesh 1
  MeshGenerator generator {};
  Mesh& mesh_1 = generator.new_mesh( domain_1, 1, 1);

  CHECK(
    generator.quad_layer_generation(mesh_1)
      .n_layers( 1 )
      .first_height( 0.20 )
      .growth_rate( 1.0 )
      .starting_position( 0.0, 0.0 )
      .ending_position( 0.0, 0.0 )
      .generate_elements()
  );

  CHECK( generator.triangulation(mesh_1).generate_elements() );
  CHECK( mesh_1.n_quads() == 8 );
  CHECK( EntityChecks::check_mesh_validity( mesh_1 ) );



  // Generate mesh 2
  Mesh& mesh_2 = generator.new_mesh( domain_2, 2, 2);

  generator.quad_layer_generation(mesh_2)
    .n_layers( 1 )
    .first_height( 0.20 )
    .growth_rate( 1.0 )
    .starting_position( 0.0, 0.0 )
    .ending_position( 0.0, 0.0 );
  CHECK( generator.quad_layer_generation(mesh_2).generate_elements() );

  CHECK( generator.triangulation(mesh_2).generate_elements() );
  CHECK( mesh_2.n_quads() == 8 );
  CHECK( EntityChecks::check_mesh_validity( mesh_2 ) );

  // Refinement of mesh 1 must fail, since it is connected to 
  // mesh 2
  CHECK( !generator.quad_refinement(mesh_1).refine() );

  // Merge both meshes
  CHECK( generator.merge_meshes( mesh_1, mesh_2 ) );
  CHECK( EntityChecks::check_mesh_validity( mesh_1 ) );
  CHECK( !generator.is_valid( mesh_2 ) );

  //MeshCleanup::merge_triangles_to_quads(mesh_1);
  //MeshCleanup::merge_degenerate_triangles(mesh_1);
    
  CHECK( generator.tri2quad_modification(mesh_1).modify() );

  // Refinement of mesh 1 must work now, since it has been
  // merged with mesh 2
  CHECK( generator.quad_refinement(mesh_1).refine() );
  CHECK( generator.quad_refinement(mesh_1).refine() );
  CHECK( generator.quad_refinement(mesh_1).refine() );

  // Smoothing of mesh 1
  //CHECK( generator.mixed_smoothing(mesh_1).smooth(2) );

  // Check mesh validity
  CHECK( EntityChecks::check_mesh_validity( mesh_1 ) );

  // Write mesh to vtu file
  std::string source_directory { TQMESH_SOURCE_DIR };
  std::string filepath { source_directory 
    + "/auxiliary/test_data/" 
    + "MeshGeneratorTests.mesh_initializer" 
  };

  CHECK( generator.write_mesh( mesh_1, filepath, MeshExportType::VTU ) );
  CHECK( !generator.write_mesh( mesh_2, filepath, MeshExportType::VTU ) );

} // mesh_initializer()

/*********************************************************************
* Test multiple neighboring meshes
*********************************************************************/
void multiple_neighbors()
{
  UserSizeFunction f_c  = [](const Vec2d& p) { return 2.5; };
  UserSizeFunction f_n  = [](const Vec2d& p) { return 2.5; };
  UserSizeFunction f_ne = [](const Vec2d& p) { return 2.5; };
  UserSizeFunction f_e  = [](const Vec2d& p) { return 2.5; };
  UserSizeFunction f_se = [](const Vec2d& p) { return 2.5; };

  TQMeshSetup::get_instance().set_quadtree_scale( 25.0 );

  Domain domain_c  { f_c  };
  Domain domain_n  { f_n  };
  Domain domain_ne { f_ne };
  Domain domain_e  { f_e  };
  Domain domain_se { f_se };

  // Center mesh
  Vertex& v1_c = domain_c.add_vertex( -2.5, -2.5 );
  Vertex& v2_c = domain_c.add_vertex(  2.5, -2.5 );
  Vertex& v3_c = domain_c.add_vertex(  2.5,  2.5 );
  Vertex& v4_c = domain_c.add_vertex( -2.5,  2.5 );

  // North mesh
  Vertex& v1_n = domain_n.add_vertex( -2.5,  2.5 );
  Vertex& v2_n = domain_n.add_vertex(  2.5,  2.5 );
  Vertex& v3_n = domain_n.add_vertex(  2.5,  5.0 );
  Vertex& v4_n = domain_n.add_vertex( -2.5,  5.0 );

  // North east mesh
  Vertex& v1_ne = domain_ne.add_vertex(  2.5,  2.5 );
  Vertex& v2_ne = domain_ne.add_vertex(  5.0,  2.5 );
  Vertex& v3_ne = domain_ne.add_vertex(  5.0,  5.0 );
  Vertex& v4_ne = domain_ne.add_vertex(  2.5,  5.0 );

  // East mesh
  Vertex& v1_e = domain_e.add_vertex(  2.5, -2.5 );
  Vertex& v2_e = domain_e.add_vertex(  5.0, -2.5 );
  Vertex& v3_e = domain_e.add_vertex(  5.0,  2.5 );
  Vertex& v4_e = domain_e.add_vertex(  2.5,  2.5 );

  // South east mesh
  Vertex& v1_se = domain_se.add_vertex(  2.5, -5.0 );
  Vertex& v2_se = domain_se.add_vertex(  5.0, -5.0 );
  Vertex& v3_se = domain_se.add_vertex(  5.0, -2.5 );
  Vertex& v4_se = domain_se.add_vertex(  2.5, -2.5 );


  Boundary& bdry_c = domain_c.add_exterior_boundary();
  bdry_c.add_edge( v1_c, v2_c, 1 );
  bdry_c.add_edge( v2_c, v3_c, 1 );
  bdry_c.add_edge( v3_c, v4_c, 1 );
  bdry_c.add_edge( v4_c, v1_c, 1 );

  Boundary& bdry_n = domain_n.add_exterior_boundary();
  bdry_n.add_edge( v1_n, v2_n, 2 );
  bdry_n.add_edge( v2_n, v3_n, 2 );
  bdry_n.add_edge( v3_n, v4_n, 2 );
  bdry_n.add_edge( v4_n, v1_n, 2 );

  Boundary& bdry_ne = domain_ne.add_exterior_boundary();
  bdry_ne.add_edge( v1_ne, v2_ne, 3 );
  bdry_ne.add_edge( v2_ne, v3_ne, 3 );
  bdry_ne.add_edge( v3_ne, v4_ne, 3 );
  bdry_ne.add_edge( v4_ne, v1_ne, 3 );

  Boundary& bdry_e = domain_e.add_exterior_boundary();
  bdry_e.add_edge( v1_e, v2_e, 4 );
  bdry_e.add_edge( v2_e, v3_e, 4 );
  bdry_e.add_edge( v3_e, v4_e, 4 );
  bdry_e.add_edge( v4_e, v1_e, 4 );

  Boundary& bdry_se = domain_se.add_exterior_boundary();
  bdry_se.add_edge( v1_se, v2_se, 5 );
  bdry_se.add_edge( v2_se, v3_se, 5 );
  bdry_se.add_edge( v3_se, v4_se, 5 );
  bdry_se.add_edge( v4_se, v1_se, 5 );




  // Setup the generator
  MeshGenerator generator {};
  Mesh& mesh_c  = generator.new_mesh( domain_c,  1, 1);
  Mesh& mesh_n  = generator.new_mesh( domain_n,  2, 2);
  Mesh& mesh_ne = generator.new_mesh( domain_ne, 3, 3);
  Mesh& mesh_e  = generator.new_mesh( domain_e,  4, 4);
  Mesh& mesh_se = generator.new_mesh( domain_se, 5, 5);


  // Generate meshes
  CHECK( generator.triangulation(mesh_c).generate_elements() );
  CHECK( generator.triangulation(mesh_n).generate_elements() );
  CHECK( generator.triangulation(mesh_ne).generate_elements() );
  CHECK( generator.triangulation(mesh_e).generate_elements() );
  CHECK( generator.triangulation(mesh_se).generate_elements() );

  CHECK( generator.merge_meshes( mesh_c, mesh_n) );
  CHECK( generator.merge_meshes( mesh_c, mesh_ne) );
  CHECK( generator.merge_meshes( mesh_c, mesh_e) );
  CHECK( generator.merge_meshes( mesh_c, mesh_se) );

  CHECK( generator.quad_refinement(mesh_c).refine() );
  CHECK( generator.quad_refinement(mesh_c).refine() );
  CHECK( generator.quad_refinement(mesh_c).refine() );


  MeshCleanup::assign_size_function_to_vertices(mesh_c, domain_c);
  MeshCleanup::assign_mesh_indices(mesh_c);
  MeshCleanup::setup_facet_connectivity(mesh_c);

  LOG(INFO) << "\n" << mesh_c;


} // multiple_neighbors()

/*********************************************************************
* Test quad layer generation near global mesh size
*********************************************************************/
void quad_layer_near_mesh_size()
{
  UserSizeFunction f = [](const Vec2d& p) { return 0.1; };

  std::vector<Vec2d> exterior_vertex_coordinates { 
    { 0.0, 0.0 },
    { 1.0, 0.0 },
    { 1.0, 2.0 },
    { 0.0, 2.0 } 
  };

  std::vector<int> exterior_edge_colors ( 4, 1 );

  TQMeshSetup::get_instance().set_quadtree_scale( 8.0 );

  Domain domain { f };
  Boundary& b_ext = domain.add_exterior_boundary();
  b_ext.set_shape_from_coordinates( exterior_vertex_coordinates, 
                                    exterior_edge_colors );

  MeshGenerator generator {};
  Mesh& mesh = generator.new_mesh( domain );

  generator.quad_layer_generation(mesh)
    .n_layers(40)
    .first_height(0.005)
    .growth_rate(1.1)
    .starting_position( {0.0,0.0})
    .ending_position( {1.0,0.0})
    .generate_elements();

  generator.triangulation(mesh).generate_elements();
  generator.mixed_smoothing(mesh).smooth(2);

  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string filename 
  { source_dir + "/auxiliary/test_data/MeshGeneratorTests.quad_layer_near_mesh_size" };

  generator.write_mesh(mesh, filename, MeshExportType::VTU);

} // quad_layer_near_mesh_size()

/*********************************************************************
* This is the same domain as in the example of the square-in-channel.
* However, we use some slightly different parameters, which (in the 
* current state of the code) lead to a failed mesh attempt.
* -> The meshing process did not fail but a missing call to 
* "MeshCleanup::setup_facet_connectivity()" after the quad refinement
* led to an invalid internal mesh structure that resulted in the 
* occurence of advancing front edges in the final mesh.
*********************************************************************/
void quad_refinement()
{
  UserSizeFunction f = [](const Vec2d& p) { return 0.1; };

  Domain domain { f };

  Boundary&  b_ext = domain.add_exterior_boundary();

  Vertex& v0 = domain.add_vertex(  0.0,  0.0 ); // Vertex coordinates
  Vertex& v1 = domain.add_vertex(  4.0,  0.0 ); // ...
  Vertex& v2 = domain.add_vertex(  4.0,  1.0 );
  Vertex& v3 = domain.add_vertex(  0.0,  1.0 );

  b_ext.add_edge( v0, v1, 2 ); // 2: color for bottom edge
  b_ext.add_edge( v1, v2, 3 ); // 3: color for right edge
  b_ext.add_edge( v2, v3, 2 ); // 2: color for top edge
  b_ext.add_edge( v3, v0, 1 ); // 1: color for left edge


  Boundary&  b_int = domain.add_interior_boundary();


  // Apply a refinement at the interior boundary vertices to a local
  // mesh size of 0.05 - within a range of 0.2:
  Vertex& v4 = domain.add_vertex( 0.35, 0.35, 0.05, 0.2 );
  Vertex& v5 = domain.add_vertex( 0.35, 0.65, 0.05, 0.2 );
  Vertex& v6 = domain.add_vertex( 0.65, 0.65, 0.05, 0.2 );
  Vertex& v7 = domain.add_vertex( 0.65, 0.35, 0.05, 0.2 );

  b_int.add_edge( v4, v5, 4 ); // Use color 4 for all interior edges
  b_int.add_edge( v5, v6, 4 ); // ...
  b_int.add_edge( v6, v7, 4 );
  b_int.add_edge( v7, v4, 4 );


  MeshGenerator generator {};
  Mesh& mesh = generator.new_mesh( domain );


  auto& quad_layer = generator.quad_layer_generation(mesh);

  quad_layer.n_layers(3);                // Number of quad layers
  quad_layer.first_height(0.01);         // First layer height
  quad_layer.growth_rate(2.0);           // Growth rate between layers
  quad_layer.starting_position(v0.xy()); // Start and ending
  quad_layer.ending_position(v1.xy());   // coordinates of the layer
  quad_layer.generate_elements();

  // ... or apply everything directly
  generator.quad_layer_generation(mesh)
    .n_layers(3)
    .first_height(0.01)
    .growth_rate(2.0)
    .starting_position(v2.xy())
    .ending_position(v3.xy())
    .generate_elements();

  generator.quad_layer_generation(mesh)
    .n_layers(3)
    .first_height(0.01)
    .growth_rate(1.3)
    .starting_position(v4.xy())
    .ending_position(v4.xy())
    .generate_elements();

  generator.triangulation(mesh).generate_elements();
  CHECK(mesh.get_front_edges().size() == 0);

  generator.tri2quad_modification(mesh).modify();
  CHECK(mesh.get_front_edges().size() == 0);

  generator.quad_refinement(mesh).refine();

  // This call is required to update the facet connectivity, 
  // otherwise the mesh will return some "nonexistent" advancing 
  // front edges in "get_front_edges()"
  MeshCleanup::setup_facet_connectivity(mesh);
  CHECK(mesh.get_front_edges().size() == 0);

  generator.mixed_smoothing(mesh)
    .epsilon(0.7)               // Change the smoothing strength
    .quad_layer_smoothing(true) // Enable smoothing to quad layers
    .smooth(3);                 // Smooth for three iterations
  CHECK(mesh.get_front_edges().size() == 0);


  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string filename 
  { source_dir + "/auxiliary/test_data/MeshGeneratorTests.quad_refinement" };

  generator.write_mesh(mesh, filename, MeshExportType::TXT);

} // quad_refinement()

/*********************************************************************
* Test fixed interior edges 
*********************************************************************/
void fixed_interior_edges()
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { return 0.5; };

  TQMeshSetup::get_instance().set_quadtree_scale( 10.0 );

  // Define the domain
  Domain domain { f };

  // Define vertices of the exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0 );
  Vertex& v4 = domain.add_vertex(  0.0,  5.0 );

  // Build exterior boundary
  int exterior_edge_color = 1;
  Boundary& exterior_bdry = domain.add_exterior_boundary();
  exterior_bdry.add_edge( v1, v2, exterior_edge_color );
  exterior_bdry.add_edge( v2, v3, exterior_edge_color );
  exterior_bdry.add_edge( v3, v4, exterior_edge_color );
  exterior_bdry.add_edge( v4, v1, exterior_edge_color );

  // Add fixed vertices
  Vertex& v1_f = domain.add_fixed_vertex(2.5, 1.5, 0.05, 1.5);
  Vertex& v2_f = domain.add_fixed_vertex(2.5, 3.5, 0.05, 1.5);

  // Add fixed edges 
  domain.add_fixed_edge( v1_f, v2_f );

  // Setup the generator
  MeshGenerator generator {};
  Mesh& mesh = generator.new_mesh( domain );
  generator.triangulation(mesh).generate_elements();

  // Check mesh 
  MeshChecker mesh_checker {mesh, domain};
  CHECK( mesh_checker.check_completeness() );

  // Export mesh
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string filename 
  { source_dir + "/auxiliary/test_data/MeshGeneratorTests.fixed_interior_edges" };

  generator.write_mesh(mesh, filename, MeshExportType::TXT);

} // fixed_interior_edges() */


/*********************************************************************
* Test fixed edges 
*********************************************************************/
void fixed_edges()
{
  UserSizeFunction f = [](const Vec2d& p) { return 0.3; };

  Domain domain { f };

  // Vertices
  Vertex& v0 = domain.add_vertex(0.0,  0.0 );
  Vertex& v1 = domain.add_vertex(3.0,  0.0 );
  Vertex& v2 = domain.add_vertex(6.0,  0.0 );
  Vertex& v3 = domain.add_vertex(7.0,  5.0 );
  Vertex& v4 = domain.add_vertex(0.0,  6.0 );

  Boundary&  b_ext = domain.add_exterior_boundary();
  b_ext.add_edge( v0, v1, 1 );
  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 2 );
  b_ext.add_edge( v3, v4, 3 );
  b_ext.add_edge( v4, v0, 4 );

  Vertex& v5 = domain.add_vertex(2.0,  2.0 );
  Vertex& v6 = domain.add_vertex(4.0,  4.0 );
  Vertex& v7 = domain.add_vertex(4.0,  2.0 );

  Boundary&  b_int = domain.add_interior_boundary();
  b_int.add_edge( v5, v6, 5 );
  b_int.add_edge( v6, v7, 5 );
  b_int.add_edge( v7, v5, 5 );

  // Fixed vertices
  Vertex& v8 = domain.add_fixed_vertex(1.5, 4,  0.05, 1.0);

  // Define fixed edges
  domain.add_fixed_edge( v0, v5 );
  domain.add_fixed_edge( v1, v6 );
  domain.add_fixed_edge( v5, v8 );
  domain.add_fixed_edge( v6, v8 );

  // Setup the generator
  MeshGenerator generator {};
  Mesh& mesh = generator.new_mesh( domain );
  generator.triangulation(mesh).generate_elements();

  generator.mixed_smoothing(mesh)
    .epsilon(0.7)
    .smooth(5);

  // Check mesh 
  MeshChecker mesh_checker { mesh, domain };
  CHECK( mesh_checker.check_completeness() );

  // Export mesh
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string filename 
  { source_dir + "/auxiliary/test_data/MeshGeneratorTests.fixed_edges" };

  generator.write_mesh(mesh, filename, MeshExportType::TXT);

} // fixed_edges()


} // namespace MeshGeneratorTests

/*********************************************************************
* Run tests for: MeshGenerator.h
*********************************************************************/
void run_tests_MeshGenerator()
{
  //adjust_logging_output_stream("MeshGeneratorTests.initialization.log");
  //MeshGeneratorTests::initialization();

  //adjust_logging_output_stream("MeshGeneratorTests.mesh_initializer.log");
  //MeshGeneratorTests::mesh_initializer();

  //adjust_logging_output_stream("MeshGeneratorTests.multiple_neighbors.log");
  //MeshGeneratorTests::multiple_neighbors();

  //adjust_logging_output_stream("MeshGeneratorTests.quad_layer_near_mesh_size.log");
  //MeshGeneratorTests::quad_layer_near_mesh_size();

  //adjust_logging_output_stream("MeshGeneratorTests.quad_refinement.log");
  //MeshGeneratorTests::quad_refinement();

  //adjust_logging_output_stream("MeshGeneratorTests.fixed_interior_edges.log");
  //MeshGeneratorTests::fixed_interior_edges();

  adjust_logging_output_stream("MeshGeneratorTests.fixed_edges.log");
  MeshGeneratorTests::fixed_edges();

  // Reset debug logging ostream
  adjust_logging_output_stream("COUT");

} // run_tests_Mesh()
