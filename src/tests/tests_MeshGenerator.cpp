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
#include "MeshInitializer.h"
#include "Cleanup.h"
#include "FrontInitializer.h"

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
  //MeshGenerator mesh_generator { outer_domain };
    
    
} // initialization()

/*********************************************************************
* Test mesh initialization for two meshes of different colors and 
* size functions
*********************************************************************/
void mesh_initializer()
{
  // Define size functions
  UserSizeFunction f_1 = [](const Vec2d& p) { return 1.0; };
  UserSizeFunction f_2 = [](const Vec2d& p) { return 0.5; };

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

  MeshInitializer initializer {};

  Mesh mesh_2 = initializer.create_empty_mesh(domain_2);
  initializer.prepare_mesh(mesh_2, domain_2);
  initializer.add_mesh_and_domain(mesh_2, domain_2);
  Cleanup::assign_mesh_indices(mesh_2);

  Mesh mesh_1 = initializer.create_empty_mesh(domain_1);
  initializer.prepare_mesh(mesh_1, domain_1);
  initializer.add_mesh_and_domain(mesh_1, domain_1);
  Cleanup::assign_mesh_indices(mesh_1);

  LOG(INFO) << "\n" << mesh_1;


} // mesh_initializer()

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

  // Reset debug logging ostream
  adjust_logging_output_stream("COUT");

} // run_tests_Mesh()
