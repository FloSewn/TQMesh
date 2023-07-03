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
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) { return 1.0; };

  // Define the outer domain
  double quadtree_scale = 10.0;
  Domain domain { f, quadtree_scale };

  // Define boundary vertices 
  Vertex& v1 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0 );
  Vertex& v4 = domain.add_vertex(  0.0,  5.0 );

  // Build domain boundary
  int exterior_edge_marker = 1;
  Boundary& exterior_bdry = domain.add_exterior_boundary();
  exterior_bdry.add_edge( v1, v2, exterior_edge_marker );
  exterior_bdry.add_edge( v2, v3, exterior_edge_marker );
  exterior_bdry.add_edge( v3, v4, exterior_edge_marker );
  exterior_bdry.add_edge( v4, v1, exterior_edge_marker );


  MeshInitializer initializer {};

  Mesh mesh = initializer.create_empty_mesh(domain);

  initializer.prepare_mesh(mesh, domain);

  Cleanup::assign_mesh_indices(mesh);

  LOG(INFO) << "\n" << mesh;

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
