/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <iostream>
#include <cassert>

#include "TQMesh.h"
#include "run_examples.h"

using namespace CppUtils;
using namespace TQMesh;

/*********************************************************************
* This example covers the generation of a mesh for a square within 
* a channel, a test case that is often considered as benchmark for 
* computational fluid dynamics solvers.
* The square and the channel boundaries will be discretized with 
* dedicated quad layers and we will utilize a method to obtain an
* all-quad mesh.
*********************************************************************/
bool square_in_channel()
{
  /*------------------------------------------------------------------
  | Define the size function and the domain structure
  ------------------------------------------------------------------*/
  UserSizeFunction f = [](const Vec2d& p) { return 0.1; };

  Domain domain { f };

  /*------------------------------------------------------------------
  | We will build the following mesh boundaries:
  |
  |       v3                                             v2     
  |      *----------------------------------------------* 
  |      |      v5    v6                                |
  |      |       *---*                                  |
  |      |       |   |                                  |
  |      |       *---*                                  |
  |      |      v4    v7                                |
  |      *----------------------------------------------*
  |       v0                                             v1
  | 
  | In contrast to example 1, we can also generate boundaries
  | successively through the addition of vertices and edge segments.
  ------------------------------------------------------------------*/
  Boundary&  b_ext = domain.add_exterior_boundary();

  Vertex& v0 = domain.add_vertex(  0.0,  0.0 ); // Vertex coordinates
  Vertex& v1 = domain.add_vertex(  4.0,  0.0 ); // ...
  Vertex& v2 = domain.add_vertex(  4.0,  1.0 );
  Vertex& v3 = domain.add_vertex(  0.0,  1.0 );

  b_ext.add_edge( v0, v1, 2 ); // 2: Color for bottom edge
  b_ext.add_edge( v1, v2, 3 ); // 3: Color for right edge
  b_ext.add_edge( v2, v3, 2 ); // 2: Color for top edge
  b_ext.add_edge( v3, v0, 1 ); // 1: Color for left edge

  /*------------------------------------------------------------------
  | However, compared to providing the edges in terms of closed
  | polygonal chains (as in example 1), this construction requires 
  | us to maintain a correct orientation of all defined edges.
  | We must ensure that:
  |
  | - Exterior boundaries are defined in counter-clockwise orientation
  | 
  | - Interior boundaries are defined in clockwise orientation
  | 
  | We can inspect this difference in terms of the interior boundary
  | defined next: 
  ------------------------------------------------------------------*/
  Boundary&  b_int = domain.add_interior_boundary();

  // Apply a refinement at the interior boundary vertices to a local
  // mesh size of 0.05 - within a range of 0.2:
  Vertex& v4 = domain.add_vertex( 0.35, 0.35, 0.03, 0.25 );
  Vertex& v5 = domain.add_vertex( 0.35, 0.65, 0.03, 0.25 );
  Vertex& v6 = domain.add_vertex( 0.65, 0.65, 0.03, 0.25 );
  Vertex& v7 = domain.add_vertex( 0.65, 0.35, 0.03, 0.25 );

  b_int.add_edge( v4, v5, 4 ); // Use color 4 for all interior edges
  b_int.add_edge( v5, v6, 4 ); // ...
  b_int.add_edge( v6, v7, 4 );
  b_int.add_edge( v7, v4, 4 );

  /*------------------------------------------------------------------
  | Now we are ready to initialize the mesh structure
  ------------------------------------------------------------------*/
  MeshGenerator generator {};
  Mesh& mesh = generator.new_mesh( domain );

  /*------------------------------------------------------------------
  | Next we will create several quad layers at the following domain 
  | boundary edges:
  | 1) Edge (v0,v1)
  | 2) Edge (v2,v3)
  | 3) Edges (v4,v5), (v5,v6), (v6,v7), (v7,v4)
  |    -> by providing vertex v5 twice to "create_quad_layers()",
  |       all edge segments that are connected to v5 in a traversable
  |       group will be used for the quad layer generation
  ------------------------------------------------------------------*/
  // We can either define and generate quad layers in two successive
  // steps:
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

  /*------------------------------------------------------------------
  | Finally, we will create the mesh by first triangulating it.
  | Then we turn most triangles into quadrlilaterals using the 
  | "tri2quad" mesh modification.
  ------------------------------------------------------------------*/
  generator.triangulation(mesh).generate_elements();
  generator.tri2quad_modification(mesh).modify();

  /*------------------------------------------------------------------
  | In order to obtain a mesh that only consists of quad elements,
  | we will refine the mesh
  ------------------------------------------------------------------*/
  generator.quad_refinement(mesh).refine();

  /*------------------------------------------------------------------
  | To improve the mesh quality, we will smooth the mesh for 
  | a few iterations. For this example, we will allow elements within  
  | quad layers to be smoothed as well.
  ------------------------------------------------------------------*/
  generator.mixed_smoothing(mesh)
    .epsilon(0.7)               // Change the smoothing strength 
    .quad_layer_smoothing(true) // Enable smoothing to quad layers
    .smooth(3);                 // Smooth for three iterations

  /*------------------------------------------------------------------
  | Check if the meshing generation process succeeded
  ------------------------------------------------------------------*/
  MeshChecker checker { mesh, domain };
  if ( !checker.check_completeness() )
  {
    LOG(ERROR) << "Mesh generation failed";
    return false;
  }

  /*------------------------------------------------------------------
  | Finally, we export the mesh to a file in VTU / TXT format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/example_data/square_in_channel" };

  LOG(INFO) << "Writing mesh output to: " << file_name << ".vtu";
  generator.write_mesh(mesh, file_name, MeshExportType::VTU );

  LOG(INFO) << "Writing mesh output to: " << file_name << ".txt";
  generator.write_mesh(mesh, file_name, MeshExportType::TXT );

  return true;

} // square_in_channel()
