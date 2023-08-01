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

#include "run_examples.h"

#include "VecND.h"
#include "Log.h"

#include "MeshGenerator.h"

using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* This example covers the generation of a simple mixed 
* triangle / quad mesh which features quad layers at specified 
* boundaries.
*********************************************************************/
void run_example_2()
{
  /*------------------------------------------------------------------
  | Define the size function
  ------------------------------------------------------------------*/
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 0.1;
  };

  Domain domain   { f };

  /*------------------------------------------------------------------
  | Build the mesh domain boundaries
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
  ------------------------------------------------------------------*/
  Boundary&  b_ext = domain.add_exterior_boundary();
  Boundary&  b_int = domain.add_interior_boundary();

  // Exterior boundary
  Vertex& v0 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v1 = domain.add_vertex(  4.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  4.0,  1.0 );
  Vertex& v3 = domain.add_vertex(  0.0,  1.0 );

  Vertex& v4 = domain.add_vertex( 0.35, 0.35, 0.05, 0.2 );
  Vertex& v5 = domain.add_vertex( 0.35, 0.65, 0.05, 0.2 );
  Vertex& v6 = domain.add_vertex( 0.65, 0.65, 0.05, 0.2 );
  Vertex& v7 = domain.add_vertex( 0.65, 0.35, 0.05, 0.2 );

  b_ext.add_edge( v0, v1, 2 );
  b_ext.add_edge( v1, v2, 3 );
  b_ext.add_edge( v2, v3, 2 );
  b_ext.add_edge( v3, v0, 1 );

  // Interior boundary
  b_int.add_edge( v4, v5, 4 );
  b_int.add_edge( v5, v6, 4 );
  b_int.add_edge( v6, v7, 4 );
  b_int.add_edge( v7, v4, 4 );

  /*------------------------------------------------------------------
  | Initialize the mesh
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
  // You can either use this format to generate the quad layers
  auto& quad_layer = generator.quad_layer_generation(mesh);
  quad_layer.n_layers(3);                // The number of quad layers
  quad_layer.first_height(0.01);         // The first layer height
  quad_layer.growth_rate(2.0);           // The growth rate between layers
  quad_layer.starting_position(v0.xy()); // Start and ending 
  quad_layer.ending_position(v1.xy());   // coordinates of the layer 
  quad_layer.generate_elements();       

  // ... or this format
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
  | Smooth the mesh for four iterations
  ------------------------------------------------------------------*/
  generator.mixed_smoothing(mesh)
    .epsilon(0.9)  // This parameter controls the smoothing strength 
    .smooth(8);    // We apply 8 smoothing iterations

  /*------------------------------------------------------------------
  | Finally, the mesh is exportet to a file in TXT format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/example_data/Example_2" };
  LOG(INFO) << "Writing mesh output to: " << file_name << ".vtu";

  generator.write_mesh(mesh, file_name, MeshExportType::VTU );

} // run_example_2()
