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
* This example covers the boundary definition via CSV files
*********************************************************************/
bool airfoil_from_csv()
{
  std::string csv_file { TQMESH_SOURCE_DIR };
  csv_file += "/auxiliary/test_data/Airfoil.csv";

  /*------------------------------------------------------------------
  | Define the size function
  ------------------------------------------------------------------*/
  UserSizeFunction f_outer = [](const Vec2d& p) { return 0.004; };

  Domain outer_domain { f_outer };
  outer_domain.add_exterior_boundary().set_shape_circle(1, {0.77, 0.09}, 0.14, 60);
  outer_domain.add_interior_boundary().set_shape_from_csv(csv_file);

  /*------------------------------------------------------------------
  | We use some additional fixed vertices to locally refine the mesh
  ------------------------------------------------------------------*/
  outer_domain.add_fixed_vertex(0.69095, 0.09625, 0.0020, 0.02);
  outer_domain.add_fixed_vertex(0.85985, 0.07582, 0.0008, 0.005);

  /*------------------------------------------------------------------
  | Initialize the outer mesh
  ------------------------------------------------------------------*/
  MeshGenerator generator {};

  int outer_mesh_id = 1;
  int outer_mesh_color = 1;
  Mesh& outer_mesh 
    = generator.new_mesh(outer_domain, outer_mesh_id, outer_mesh_color);

  /*------------------------------------------------------------------
  | Create quad layers 
  ------------------------------------------------------------------*/
  generator.quad_layer_generation(outer_mesh)
    .n_layers(10)
    .first_height(0.0004)
    .growth_rate(1.10)
    .starting_position({0.69132,0.09754})
    .ending_position({0.69132,0.09754})
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
  | Check if the meshing generation process succeeded
  ------------------------------------------------------------------*/
  MeshChecker outer_checker { outer_mesh, outer_domain };
  if ( !outer_checker.check_completeness() )
  {
    LOG(ERROR) << "Mesh generation failed";
    return false;
  }

  /*------------------------------------------------------------------
  | Initialize the inner mesh
  ------------------------------------------------------------------*/
  UserSizeFunction f_inner = [](const Vec2d& p) { return 0.003; };

  Domain inner_domain { f_inner };
  inner_domain.add_exterior_boundary().set_shape_from_csv(csv_file);

  /*------------------------------------------------------------------
  | Initialize the inner mesh
  ------------------------------------------------------------------*/
  int inner_mesh_id = 2;
  int inner_mesh_color = 2;
  Mesh& inner_mesh 
    = generator.new_mesh(inner_domain, inner_mesh_id, inner_mesh_color);

  /*------------------------------------------------------------------
  | Generate the inner mesh
  ------------------------------------------------------------------*/
  generator.triangulation(inner_mesh).generate_elements();

  /*------------------------------------------------------------------
  | Check if the meshing generation process succeeded
  ------------------------------------------------------------------*/
  MeshChecker inner_checker { inner_mesh, inner_domain };
  if ( !inner_checker.check_completeness() )
  {
    LOG(ERROR) << "Mesh generation failed";
    return false;
  }

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
  { source_dir + "/auxiliary/example_data/airfoil" };

  LOG(INFO) << "Writing mesh output to: " << file_name << ".vtu";
  generator.write_mesh(inner_mesh, file_name, MeshExportType::VTU);

  LOG(INFO) << "Writing mesh output to: " << file_name << ".txt";
  generator.write_mesh(inner_mesh, file_name, MeshExportType::TXT);

  return true;

} // airfoil_from_csv()
