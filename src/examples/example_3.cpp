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

#include "Vec2.h"
#include "Log.h"

#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Mesh.h"
#include "Smoother.h"

using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* 
*********************************************************************/
void run_example_3()
{
  /*------------------------------------------------------------------
  | Define the size function
  ------------------------------------------------------------------*/
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 0.4;
  };

  Domain domain   { f };

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  Boundary& b_ext = domain.add_exterior_boundary();
  b_ext.set_shape_rectangle( domain.vertices(), 1, {1.0,1.0}, 8.0, 8.0 );

  Boundary& b_shape_1 = domain.add_interior_boundary();
  b_shape_1.set_shape_circle( domain.vertices(), 2, {0.0,1.0}, 1.25, 30 );

  Boundary& b_shape_2 = domain.add_interior_boundary();
  b_shape_2.set_shape_square( domain.vertices(), 3, {3.0,2.5}, 1.75);

  Boundary& b_shape_3 = domain.add_interior_boundary();
  b_shape_3.set_shape_triangle( domain.vertices(), 4, {3.0,-0.5}, 1.75 );

  /*------------------------------------------------------------------
  | Initialize the mesh
  ------------------------------------------------------------------*/
  Mesh mesh { domain };
  mesh.init_advancing_front();

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  Vertex& v5 = domain.vertices()[4];
  mesh.create_quad_layers(v5, v5, 3, 0.05, 1.3);

  Vertex& v35 = domain.vertices()[34];
  Vertex& v37 = domain.vertices()[36];
  mesh.create_quad_layers(v37, v35, 1, 0.25, 1.0);

  Vertex& v39 = domain.vertices()[38];
  mesh.create_quad_layers(v39, v39, 1, 0.30, 1.0);

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  mesh.triangulate();

  /*------------------------------------------------------------------
  | Smooth the mesh for four iterations
  ------------------------------------------------------------------*/
  Smoother smoother {};
  smoother.smooth(domain, mesh, 4);

  /*------------------------------------------------------------------
  | Finally, the mesh is exportet to a file in VTU format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/aux/example_data/Example_3" };

  mesh.write_to_file( file_name, ExportType::txt );

} // run_example_3()
