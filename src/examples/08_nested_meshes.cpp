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
* This example shows how to generate and merge several nested meshes 
* with varying colors and size functions.
*
*     x--------------------------------x
*     |             Mesh 1             |
*     |   x------------------------x   |
*     |   |         Mesh 2         |   |
*     |   |   x----------------x   |   |
*     |   |   |     Mesh 3     |   |   |
*     |   |   |                |   |   |
*     |   |   |     x----x     |   |   |
*     |   |   |     |Mesh|     |   |   |
*     |   |   |     | 4  |     |   |   |
*     |   |   |     x----x     |   |   |
*     |   |   |                |   |   |
*     |   |   |                |   |   |
*     |   |   x----------------x   |   |
*     |   |                        |   |
*     |   x------------------------x   |
*     |                                |
*     x--------------------------------x
*
*********************************************************************/
void nested_meshes()
{
  MeshGenerator generator {};

  /*------------------------------------------------------------------
  | Define the size function and domain of mesh 1
  ------------------------------------------------------------------*/
  UserSizeFunction f_1 = [](const Vec2d& p) { return 5.0; };
  UserSizeFunction f_2 = [](const Vec2d& p) { return 2.0; };
  UserSizeFunction f_3 = [](const Vec2d& p) { return 2.0; };
  UserSizeFunction f_4 = [](const Vec2d& p) { return 1.0; };

  Domain domain_1 { f_1 };
  Domain domain_2 { f_2 };
  Domain domain_3 { f_3 };
  Domain domain_4 { f_4 };

  /*------------------------------------------------------------------
  | Exterior boundary of mesh 1
  ------------------------------------------------------------------*/
  std::vector<Vec2d> ext_coords_1 {
    { 0.0, 0.0 }, { 100.0, 0.0 }, { 100.0, 100.0 }, { 0.0, 100.0 },
  };
  std::vector<Vec2d> ext_coords_2 {
    { 10.0, 10.0 }, { 90.0, 10.0 }, { 90.0, 90.0 }, { 10.0, 90.0 },
  };
  std::vector<Vec2d> ext_coords_3 {
    { 20.0, 20.0 }, { 80.0, 20.0 }, { 80.0, 80.0 }, { 20.0, 80.0 },
  };
  std::vector<Vec2d> ext_coords_4 {
    { 40.0, 40.0 }, { 60.0, 40.0 }, { 60.0, 60.0 }, { 40.0, 60.0 },
  };

  std::vector<int> ext_markers_1 ( ext_coords_1.size(), 1 );
  std::vector<int> ext_markers_2 ( ext_coords_1.size(), 2 );
  std::vector<int> ext_markers_3 ( ext_coords_1.size(), 3 );
  std::vector<int> ext_markers_4 ( ext_coords_1.size(), 4 );

  Boundary& bdry_ext_1 = domain_1.add_exterior_boundary();
  Boundary& bdry_ext_2 = domain_2.add_exterior_boundary();
  Boundary& bdry_ext_3 = domain_3.add_exterior_boundary();
  Boundary& bdry_ext_4 = domain_4.add_exterior_boundary();

  bdry_ext_1.set_shape_from_coordinates(ext_coords_1, ext_markers_1);
  bdry_ext_2.set_shape_from_coordinates(ext_coords_2, ext_markers_2);
  bdry_ext_3.set_shape_from_coordinates(ext_coords_3, ext_markers_3);
  bdry_ext_4.set_shape_from_coordinates(ext_coords_4, ext_markers_4);

  /*------------------------------------------------------------------
  | Interior boundary of mesh 1
  ------------------------------------------------------------------*/
  std::vector<Vec2d> int_coords_1 {
    { 10.0, 10.0 }, { 90.0, 10.0 }, { 90.0, 90.0 }, { 10.0, 90.0 },
  };
  std::vector<Vec2d> int_coords_2 {
    { 20.0, 20.0 }, { 80.0, 20.0 }, { 80.0, 80.0 }, { 20.0, 80.0 },
  };
  std::vector<Vec2d> int_coords_3 {
    { 40.0, 40.0 }, { 60.0, 40.0 }, { 60.0, 60.0 }, { 40.0, 60.0 },
  };

  std::vector<int> int_markers_1 ( int_coords_1.size(), 2 );
  std::vector<int> int_markers_2 ( int_coords_1.size(), 3 );
  std::vector<int> int_markers_3 ( int_coords_1.size(), 4 );

  Boundary& bdry_int_1 = domain_1.add_interior_boundary();
  Boundary& bdry_int_2 = domain_2.add_interior_boundary();
  Boundary& bdry_int_3 = domain_3.add_interior_boundary();

  bdry_int_1.set_shape_from_coordinates(int_coords_1, int_markers_1);
  bdry_int_2.set_shape_from_coordinates(int_coords_2, int_markers_2);
  bdry_int_3.set_shape_from_coordinates(int_coords_3, int_markers_3);

  /*------------------------------------------------------------------
  | Create meshes
  ------------------------------------------------------------------*/
  int mesh_id_1    = 0;
  int mesh_id_2    = 1;
  int mesh_id_3    = 2;
  int mesh_id_4    = 3;

  int mesh_color_1 = 1;
  int mesh_color_2 = 1;
  int mesh_color_3 = 0;
  int mesh_color_4 = 0;

  Mesh& mesh_1 = generator.new_mesh(domain_1, mesh_id_1,  mesh_color_1);
  generator.triangulation(mesh_1).generate_elements();
  generator.mixed_smoothing(mesh_1).smooth(2);   

  Mesh& mesh_2 = generator.new_mesh(domain_2, mesh_id_2,  mesh_color_2);
  generator.triangulation(mesh_2).generate_elements();
  generator.mixed_smoothing(mesh_2).smooth(2);   

  Mesh& mesh_3 = generator.new_mesh(domain_3, mesh_id_3,  mesh_color_3);
  generator.triangulation(mesh_3).generate_elements();
  generator.mixed_smoothing(mesh_3).smooth(2);   

  Mesh& mesh_4 = generator.new_mesh(domain_4, mesh_id_4,  mesh_color_4);
  generator.triangulation(mesh_4).generate_elements();
  generator.mixed_smoothing(mesh_4).smooth(2);   

  /*------------------------------------------------------------------
  | Merge meshes
  ------------------------------------------------------------------*/
  generator.merge_meshes(mesh_1, mesh_2);
  generator.merge_meshes(mesh_1, mesh_3);
  generator.merge_meshes(mesh_1, mesh_4);

  /*------------------------------------------------------------------
  | Finally, the mesh is exportet to a file in TXT format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/example_data/nested_meshes" };

  LOG(INFO) << "Writing mesh output to: " << file_name << ".vtu";
  generator.write_mesh(mesh_1, file_name, MeshExportType::VTU);

  LOG(INFO) << "Writing mesh output to: " << file_name << ".txt";
  generator.write_mesh(mesh_1, file_name, MeshExportType::TXT);

} /* nested_meshes() */
