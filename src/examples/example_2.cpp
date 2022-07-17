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
  Vertex& v0 = domain.add_vertex(  0.0,  0.0, 1.0, 1.0 );
  Vertex& v1 = domain.add_vertex(  4.0,  0.0, 1.0, 1.0 );
  Vertex& v2 = domain.add_vertex(  4.0,  1.0, 1.0, 1.0 );
  Vertex& v3 = domain.add_vertex(  0.0,  1.0, 1.0, 1.0 );

  Vertex& v4 = domain.add_vertex( 0.35, 0.35, 0.8, 1.2 );
  Vertex& v5 = domain.add_vertex( 0.35, 0.65, 0.8, 1.2 );
  Vertex& v6 = domain.add_vertex( 0.65, 0.65, 0.8, 1.2 );
  Vertex& v7 = domain.add_vertex( 0.65, 0.35, 0.8, 1.2 );

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
  Mesh mesh { domain };
  mesh.init_advancing_front();

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
  mesh.create_quad_layers(v0, v1, 3, 0.01, 2.0);
  mesh.create_quad_layers(v2, v3, 3, 0.01, 2.0);
  mesh.create_quad_layers(v4, v4, 3, 0.01, 1.3);

  /*------------------------------------------------------------------
  | Finally, we will create the mesh with the "paving()" method.
  | It will create a mixed mesh that consists mainly of quads and 
  | maybe some triangles 
  ------------------------------------------------------------------*/
  mesh.pave();

  /*------------------------------------------------------------------
  | In order to obtain a mesh that only consists of quad elements,
  | we will refine the mesh
  ------------------------------------------------------------------*/
  mesh.refine_to_quads();

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
  { source_dir + "/aux/example_data/Example_2" };
  LOG(INFO) << "Writing mesh output to: " << file_name << ".txt";

  mesh.write_to_file( file_name, ExportType::txt );

} // run_example_2()
