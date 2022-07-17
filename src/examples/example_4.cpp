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
* This example covers the mesh generation with fixed interior 
* vertices, as well as the usage of different mesh element colors
*********************************************************************/
void run_example_4()
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
  Boundary& b_ext = domain.add_exterior_boundary();

  Vertex& v0 = domain.add_vertex(  0.0,  0.0,  1.0,  1.0 );
  Vertex& v1 = domain.add_vertex(  1.0,  0.0,  1.0,  1.0 );
  Vertex& v2 = domain.add_vertex(  1.0,  2.0,  0.5,  1.5 );
  Vertex& v3 = domain.add_vertex(  2.0,  2.0,  0.5,  1.5 );
  Vertex& v4 = domain.add_vertex(  2.0,  1.0,  0.5,  1.2 );
  Vertex& v5 = domain.add_vertex(  3.0,  1.0,  0.5,  1.2 );
  Vertex& v6 = domain.add_vertex(  3.0,  3.0,  1.0,  1.0 );
  Vertex& v7 = domain.add_vertex(  0.0,  3.0,  1.0,  1.0 );

  b_ext.add_edge( v0, v1, 1 );
  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v5, 1 );
  b_ext.add_edge( v5, v6, 1 );
  b_ext.add_edge( v6, v7, 1 );
  b_ext.add_edge( v7, v0, 1 );

  /*------------------------------------------------------------------
  | Here, we add two vertices in the interior of the domain,
  | which will be considered during the meshing process.
  | We reduce their scaling-factors, in order to refine the mesh
  | locally.
  ------------------------------------------------------------------*/
  domain.add_fixed_vertex(0.5, 1.0, 0.1, 0.5);
  domain.add_fixed_vertex(0.5, 2.0, 0.1, 0.5);

  /*------------------------------------------------------------------
  | Initialize the mesh. We also provide the index of the mesh,
  | as well as a color index that will be passed to all mesh elements
  ------------------------------------------------------------------*/
  int mesh_index = 0;
  int element_color = 0;
  Mesh mesh { domain, mesh_index, element_color };
  mesh.init_advancing_front();

  /*------------------------------------------------------------------
  | Next, we add three quad layers. This will create several quad-
  | elements at the boundary edges from vertex v0 to v7.
  | All these quad elements will be given the current mesh color 
  | index, which is set to zero.
  ------------------------------------------------------------------*/
  mesh.create_quad_layers(v0, v7, 3, 0.02, 1.5);

  /*------------------------------------------------------------------
  | Now generate the remaining mesh elements. But before this,
  | we change the mesh's element color to one, such that all 
  | remaining elements will be given this color index
  ------------------------------------------------------------------*/
  mesh.element_color(1);
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
  { source_dir + "/aux/example_data/Example_4" };
  LOG(INFO) << "Writing mesh output to: " << file_name << ".txt";

  mesh.write_to_file( file_name, ExportType::txt );

} // run_example_4()
