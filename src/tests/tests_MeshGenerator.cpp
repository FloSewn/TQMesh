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
#include "MeshGenerator.h"
#include "MeshMerger.h"
#include "RefinementStrategy.h"
#include "MeshCleanup.h"
#include "SmoothingStrategy.h"
#include "EntityChecks.h"

namespace MeshGeneratorTests 
{
using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

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

  double quadtree_scale = 10.0;

  // Define the outer domain
  Domain outer_domain { f_outer, quadtree_scale };

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
  int exterior_edge_marker = 1;
  Boundary& outer_exterior_bdry = outer_domain.add_exterior_boundary();
  outer_exterior_bdry.add_edge( v1, v2, exterior_edge_marker );
  outer_exterior_bdry.add_edge( v2, v3, exterior_edge_marker );
  outer_exterior_bdry.add_edge( v3, v4, exterior_edge_marker );
  outer_exterior_bdry.add_edge( v4, v1, exterior_edge_marker );

  // Build interior boundary
  int interior_edge_marker = 2;
  Boundary& outer_interior_bdry = outer_domain.add_interior_boundary();
  outer_interior_bdry.add_edge( v5, v6, interior_edge_marker );
  outer_interior_bdry.add_edge( v6, v7, interior_edge_marker );
  outer_interior_bdry.add_edge( v7, v8, interior_edge_marker );
  outer_interior_bdry.add_edge( v8, v5, interior_edge_marker );

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

  // Define domains
  double quadtree_scale = 20.0;
  Domain domain_1 { f_1, quadtree_scale };
  Domain domain_2 { f_2, quadtree_scale };

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
  int edge_marker = 1;
  Boundary& bdry_1 = domain_1.add_exterior_boundary();
  Boundary& bdry_2 = domain_2.add_exterior_boundary();
  
  bdry_1.add_edge( v1_1, v2_1, edge_marker );
  bdry_1.add_edge( v2_1, v3_1, edge_marker );
  bdry_1.add_edge( v3_1, v4_1, edge_marker );
  bdry_1.add_edge( v4_1, v1_1, edge_marker );

  bdry_2.add_edge( v1_2, v2_2, edge_marker );
  bdry_2.add_edge( v2_2, v3_2, edge_marker );
  bdry_2.add_edge( v3_2, v4_2, edge_marker );
  bdry_2.add_edge( v4_2, v1_2, edge_marker );

  // Generate mesh 1
  MeshGenerator generator {};
  Mesh& mesh_1 = generator.new_mesh( domain_1, 1, 1);

  generator.set_algorithm(mesh_1, MeshingAlgorithm::QuadLayer);
  generator.quad_layer().n_layers( 1 );
  generator.quad_layer().first_height( 0.20 );
  generator.quad_layer().growth_rate( 1.0 );
  generator.quad_layer().starting_position( 0.0, 0.0 );
  generator.quad_layer().ending_position( 0.0, 0.0 );
  CHECK( generator.generate_mesh_elements() );

  generator.set_algorithm(mesh_1, MeshingAlgorithm::Triangulation);
  CHECK( generator.generate_mesh_elements() );
  CHECK( mesh_1.n_quads() == 8 );

  // Generate mesh 2
  Mesh& mesh_2 = generator.new_mesh( domain_2, 2, 2);

  generator.set_algorithm(mesh_2, MeshingAlgorithm::QuadLayer);
  generator.quad_layer().n_layers( 1 );
  generator.quad_layer().first_height( 0.20 );
  generator.quad_layer().growth_rate( 1.0 );
  generator.quad_layer().starting_position( 0.0, 0.0 );
  generator.quad_layer().ending_position( 0.0, 0.0 );
  CHECK( generator.generate_mesh_elements() );

  generator.set_algorithm(mesh_2, MeshingAlgorithm::Triangulation);
  CHECK( generator.generate_mesh_elements() );
  CHECK( mesh_2.n_quads() == 8 );

  // Refinement of mesh 1 must fail, since it is connected to 
  // mesh 2
  CHECK( !RefinementStrategy::refine_to_quads( mesh_1 ) );

  // Merge both meshes
  MeshMerger merger { mesh_1, mesh_2 };
  CHECK( merger.merge() );

  CHECK( EntityChecks::check_mesh_validity( mesh_1 ) );

  //MeshCleanup::merge_triangles_to_quads(mesh_1);
  //MeshCleanup::merge_degenerate_triangles(mesh_1);
    
  // Refinement of mesh 1 must work now, since it has been
  // merged with mesh 2
  CHECK( RefinementStrategy::refine_to_quads( mesh_1 ) );
  CHECK( RefinementStrategy::refine_to_quads( mesh_1 ) );
  CHECK( RefinementStrategy::refine_to_quads( mesh_1 ) );

  // Smoothing of mesh 1
  generator.set_algorithm(mesh_1, SmoothingAlgorithm::Mixed);
  CHECK( generator.smooth_mesh_elements(2) );

  CHECK( EntityChecks::check_mesh_validity( mesh_1 ) );

  // Prepare for output
  MeshCleanup::assign_size_function_to_vertices(mesh_1, domain_1);
  MeshCleanup::assign_mesh_indices(mesh_1);
  MeshCleanup::setup_facet_connectivity(mesh_1);

  LOG(INFO) << "\n" << mesh_1;

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

  Domain domain_c  { f_c,  25.0 };
  Domain domain_n  { f_n,  25.0 };
  Domain domain_ne { f_ne, 25.0 };
  Domain domain_e  { f_e,  25.0 };
  Domain domain_se { f_se, 25.0 };

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
  generator.set_algorithm(mesh_c, MeshingAlgorithm::Triangulation);
  CHECK( generator.generate_mesh_elements() );

  generator.set_algorithm(mesh_n, MeshingAlgorithm::Triangulation);
  CHECK( generator.generate_mesh_elements() );

  generator.set_algorithm(mesh_ne, MeshingAlgorithm::Triangulation);
  CHECK( generator.generate_mesh_elements() );

  generator.set_algorithm(mesh_e, MeshingAlgorithm::Triangulation);
  CHECK( generator.generate_mesh_elements() );

  generator.set_algorithm(mesh_se, MeshingAlgorithm::Triangulation);
  CHECK( generator.generate_mesh_elements() );


  MeshMerger merger_1 { mesh_c, mesh_n };
  CHECK( merger_1.merge() );

  MeshMerger merger_2 { mesh_c, mesh_ne };
  CHECK( merger_2.merge() );

  MeshMerger merger_3 { mesh_c, mesh_e };
  CHECK( merger_3.merge() );

  MeshMerger merger_4 { mesh_c, mesh_se };
  CHECK( merger_4.merge() );



  MeshCleanup::assign_size_function_to_vertices(mesh_c, domain_c);
  MeshCleanup::assign_mesh_indices(mesh_c);
  MeshCleanup::setup_facet_connectivity(mesh_c);

  LOG(INFO) << "\n" << mesh_c;


} // multiple_neighbors()

} // namespace MeshGeneratorTests

/*********************************************************************
* Run tests for: MeshGenerator.h
*********************************************************************/
void run_tests_MeshGenerator()
{
  adjust_logging_output_stream("MeshGeneratorTests.initialization.log");
  MeshGeneratorTests::initialization();

  adjust_logging_output_stream("MeshGeneratorTests.mesh_initializer.log");
  MeshGeneratorTests::mesh_initializer();

  adjust_logging_output_stream("MeshGeneratorTests.multiple_neighbors.log");
  MeshGeneratorTests::multiple_neighbors();

  // Reset debug logging ostream
  adjust_logging_output_stream("COUT");

} // run_tests_Mesh()
