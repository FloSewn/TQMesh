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

#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Mesh.h"
#include "Smoother.h"

using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* This example covers the mesh generation with boundary shapes
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
  | The exterior boundary will be generated by using a rectangular
  | shape. It is defined by its center point xy = (1,1), as well as 
  | its width = 8 and height = 8. 
  | A boundary marker with value = 1 will be applied to all boundary
  | edges.
  | For this rectangular boundary, four vertices will be added to
  | the domain.
  ------------------------------------------------------------------*/
  Boundary& b_ext = domain.add_exterior_boundary();
  b_ext.set_shape_rectangle( domain.vertices(), 1, {1.0,1.0}, 8.0, 8.0 );

  /*------------------------------------------------------------------
  | Next we will place a circular interior boundary insides the 
  | domain. Its center point is set to xy = (0,1) and its radius is 
  | set to 1.25.
  | The circular boundary will be created from 30 edge segments, all
  | of which will get the boundary marker with value = 2.
  | Accordingly, 30 new vertices will be added to the domain.
  ------------------------------------------------------------------*/
  Boundary& b_shape_1 = domain.add_interior_boundary();
  b_shape_1.set_shape_circle( domain.vertices(), 2, {0.0,1.0}, 1.25, 30 );

  /*------------------------------------------------------------------
  | Finally we will add two more interior boundaries: A suqare 
  | and a triangle. These will get the markers 3 and respectively, 
  | and both will feature the edge length of 1.75.
  ------------------------------------------------------------------*/
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
  | Here we will create some quad layers. Due to the boundary 
  | generation via shapes, we did not explicitly define any vertices
  | that could be passed as starting and ending vertices for the 
  | quad layer generation algorithm.
  | However, since we know how many vertices were created during the
  | generation of the given boundary shapes, we can use this 
  | information to get the respecive vertices.
  ------------------------------------------------------------------*/
  // The following domain vertex is located on the circular 
  // interior boundary
  Vertex& v4 = domain.vertices()[4];
  mesh.create_quad_layers(v4, v4, 3, 0.05, 1.3);

  // These next two vertices are located on the interior quad boundary
  Vertex& v34 = domain.vertices()[34];
  Vertex& v36 = domain.vertices()[36];
  mesh.create_quad_layers(v36, v34, 1, 0.25, 1.0);

  // This last vertices is located on the interior triangle boundary
  Vertex& v38 = domain.vertices()[38];
  mesh.create_quad_layers(v38, v38, 1, 0.30, 1.0);

  /*------------------------------------------------------------------
  | Now generate the mesh elements
  ------------------------------------------------------------------*/
  mesh.triangulate();

  /*------------------------------------------------------------------
  | Smooth the mesh for four iterations
  ------------------------------------------------------------------*/
  Smoother smoother {};
  smoother.smooth(domain, mesh, 4);

  /*------------------------------------------------------------------
  | Finally, the mesh is exportet to a file in TXT format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/example_data/Example_3" };
  LOG(INFO) << "Writing mesh output to: " << file_name << ".txt";

  mesh.write_to_file( file_name, ExportType::txt );

} // run_example_3()
