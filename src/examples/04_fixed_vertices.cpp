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
* This example covers the mesh generation with fixed interior 
* vertices, as well as the usage of different mesh element colors
*********************************************************************/
void fixed_vertices()
{
  /*------------------------------------------------------------------
  | Define the size function and the domain structure
  ------------------------------------------------------------------*/
  UserSizeFunction f = [](const Vec2d& p) { return 0.2; };

  Domain domain   { f };

  /*------------------------------------------------------------------
  | Define the exterior boundary of the domain
  |
  |   v7                      v6
  |  x-----------------------x
  |  |                       |
  |  |                       |
  |  |     v2        v3      |
  |  |       x-------x       |
  |  |   x   |       |       |
  |  |       |       |       |
  |  |       |       |       |
  |  |       |       x-------x
  |  |   x   |       v4       v5
  |  |       |              
  |  |       |              
  |  |       |              
  |  x-------x              
  |  v0       v1
  ------------------------------------------------------------------*/
  std::vector<Vec2d> vertex_coordinates { 
    { 0.0,  0.0 },
    { 1.0,  0.0 },
    { 1.0,  2.0 },
    { 2.0,  2.0 },
    { 2.0,  1.0 },
    { 3.0,  1.0 },
    { 3.0,  3.0 },
    { 0.0,  3.0 },
  };

  std::vector<int> edge_markers ( vertex_coordinates.size(), 1 );

  Boundary& boundary = domain.add_exterior_boundary();
  boundary.set_shape_from_coordinates(vertex_coordinates, edge_markers);

  /*------------------------------------------------------------------
  | Here, we add two vertices in the interior of the domain,
  | which will be considered during the meshing process.
  | We reduce their scaling-factors, in order to refine the mesh
  | locally.
  ------------------------------------------------------------------*/
  domain.add_fixed_vertex(0.5, 1.0, 0.005, 0.5);
  domain.add_fixed_vertex(0.5, 2.0, 0.005, 0.5);

  /*------------------------------------------------------------------
  | Initialize the mesh. We also provide the index of the mesh,
  | as well as a color index that will be passed to all mesh elements
  ------------------------------------------------------------------*/
  int mesh_index = 0;
  int mesh_color = 0;
  MeshGenerator generator {};
  Mesh& mesh = generator.new_mesh( domain, mesh_index, mesh_color );

  /*------------------------------------------------------------------
  | Next, we add three quad layers. This will create several quad-
  | elements at the boundary edges from vertex v0 to v7.
  | All these quad elements will be given the current mesh color 
  | index, which is set to zero.
  ------------------------------------------------------------------*/
  generator.quad_layer_generation(mesh)
    .n_layers(3)
    .first_height(0.02)
    .growth_rate(1.5)
    .starting_position({0.0,0.0})
    .ending_position({0.0,3.0})
    .generate_elements();

  /*------------------------------------------------------------------
  | Now generate the remaining mesh elements. But before this,
  | we change the mesh's element color to one, such that all 
  | remaining elements will be given this color index
  ------------------------------------------------------------------*/
  mesh.element_color(1);
  generator.triangulation(mesh).generate_elements();

  /*------------------------------------------------------------------
  | Smooth the mesh for four iterations
  ------------------------------------------------------------------*/
  generator.mixed_smoothing(mesh).smooth(2); 

  /*------------------------------------------------------------------
  | Finally, the mesh is exportet to a file in TXT format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/example_data/fixed_vertices" };

  LOG(INFO) << "Writing mesh output to: " << file_name << ".vtu";
  generator.write_mesh(mesh, file_name, MeshExportType::VTU);

  LOG(INFO) << "Writing mesh output to: " << file_name << ".txt";
  generator.write_mesh(mesh, file_name, MeshExportType::TXT);


} // fixed_vertices()
