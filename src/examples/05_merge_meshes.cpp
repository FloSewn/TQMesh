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
* This example shows how to generate and merge several meshes with
* different colors and size functions.
* 
*   x-----------------------x
*   |  Outer                |
*   |  mesh                 |
*   |       x-------x       |
*   |       | Inner |       |
*   |       | mesh  |       |
*   |       |       |       |
*   |       x-------x       |
*   |                       |
*   |                       |
*   x-----------------------x
*
*********************************************************************/
void merge_meshes()
{
  /*------------------------------------------------------------------
  | Define the size function and domain of the outer mesh
  ------------------------------------------------------------------*/
  UserSizeFunction f_outer = [](const Vec2d& p) { return 0.35; };

  Domain outer_domain { f_outer };

  /*------------------------------------------------------------------
  | Exterior boundary of the outer mesh
  ------------------------------------------------------------------*/
  std::vector<Vec2d> outer_coords {
    { 0.0, 0.0 }, { 5.0, 0.0 }, { 5.0, 5.0 }, { 0.0, 5.0 },
  };

  std::vector<int> outer_markers ( outer_coords.size(), 1 );

  Boundary& bdry_outer_ext = outer_domain.add_exterior_boundary();
  bdry_outer_ext.set_shape_from_coordinates(outer_coords, outer_markers);

  /*------------------------------------------------------------------
  | Interior boundary of the outer mesh
  ------------------------------------------------------------------*/
  std::vector<Vec2d> inner_coords {
    { 1.5, 1.5 }, { 3.5, 1.5 }, { 3.5, 3.5 }, { 1.5, 3.5 },
  };

  std::vector<int> inner_markers ( inner_coords.size(), 2 );

  Boundary& bdry_outer_int = outer_domain.add_interior_boundary();
  bdry_outer_int.set_shape_from_coordinates(inner_coords, inner_markers);

  /*------------------------------------------------------------------
  | Initialize the outer mesh
  ------------------------------------------------------------------*/
  MeshGenerator generator {};

  int outer_mesh_id = 1;
  int outer_mesh_color = 1;
  Mesh& outer_mesh = generator.new_mesh(outer_domain, outer_mesh_id, 
                                        outer_mesh_color);

  /*------------------------------------------------------------------
  | Create two quad layers for the outer mesh
  ------------------------------------------------------------------*/
  generator.quad_layer_generation(outer_mesh)
    .n_layers(2)
    .first_height(0.05)
    .growth_rate(1.5)
    .starting_position({0.0,0.0})
    .ending_position({0.0,0.0})
    .generate_elements();

  generator.quad_layer_generation(outer_mesh)
    .n_layers(2)
    .first_height(0.05)
    .growth_rate(1.5)
    .starting_position({1.5,3.5})
    .ending_position({1.5,3.5})
    .generate_elements(); 

  /*------------------------------------------------------------------
  | Generate the remaining elements of the outer mesh
  ------------------------------------------------------------------*/
  generator.triangulation(outer_mesh).generate_elements();

  /*------------------------------------------------------------------
  | Smooth the elements of the outer mesh
  ------------------------------------------------------------------*/
  generator.mixed_smoothing(outer_mesh).smooth(2);   

  /*------------------------------------------------------------------
  | Define the size function of the inner mesh
  ------------------------------------------------------------------*/
  UserSizeFunction f_inner = [](const Vec2d& p) 
  { 
    return 0.02 + 0.2 * ( pow(p.x-2.5, 2) + pow(p.y-2.5, 2) );
  };

  Domain inner_domain   { f_inner };

  /*------------------------------------------------------------------
  | Exterior boundary of the inner mesh
  | This boundary overlaps with the interior boundary of the outer 
  | mesh. It is thus important, that it features the same number of 
  | edges and vertex coordinates.
  ------------------------------------------------------------------*/
  Boundary& bdry_inner = inner_domain.add_exterior_boundary();
  bdry_inner.set_shape_from_coordinates(inner_coords, inner_markers);

  /*------------------------------------------------------------------
  | Initialize the inner mesh
  ------------------------------------------------------------------*/
  int inner_mesh_id = 2;
  int inner_mesh_color = 2;
  Mesh& inner_mesh = generator.new_mesh(inner_domain, inner_mesh_id, 
                                        inner_mesh_color);

  /*------------------------------------------------------------------
  | Generate the inner mesh
  ------------------------------------------------------------------*/
  generator.triangulation(inner_mesh).generate_elements();

  /*------------------------------------------------------------------
  | Finally, merge both meshes. In this way, the outer mesh's 
  | elements will be added to the inner mesh. Be aware that the 
  | second argument of "merge_meshes()" (in this case the outer 
  | mesh) will be destroyed upon this function call - so you 
  | can not use it anymore.
  ------------------------------------------------------------------*/
  generator.merge_meshes(inner_mesh, outer_mesh);

  /*------------------------------------------------------------------
  | Smooth the final mesh for four iterations
  ------------------------------------------------------------------*/
  generator.mixed_smoothing(inner_mesh).smooth(2);   

  /*------------------------------------------------------------------
  | Finally, the mesh is exportet to a file in TXT format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/example_data/merge_meshes" };

  LOG(INFO) << "Writing mesh output to: " << file_name << ".vtu";
  generator.write_mesh(inner_mesh, file_name, MeshExportType::VTU);

  LOG(INFO) << "Writing mesh output to: " << file_name << ".txt";
  generator.write_mesh(inner_mesh, file_name, MeshExportType::TXT);

} // merge_meshes()
