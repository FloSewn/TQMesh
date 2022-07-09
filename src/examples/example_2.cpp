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
* This example covers the generation of a simple triangular mesh
*********************************************************************/
void run_example_2()
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 0.1;
  };

  Domain domain   { f };

  Boundary&  b_ext = domain.add_exterior_boundary();
  Boundary&  b_int = domain.add_interior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0, 1.0, 1.0 );
  Vertex& v2 = domain.add_vertex(  4.0,  0.0, 1.0, 1.0 );
  Vertex& v3 = domain.add_vertex(  4.0,  1.0, 1.0, 1.0 );
  Vertex& v4 = domain.add_vertex(  0.0,  1.0, 1.0, 1.0 );

  Vertex& v5 = domain.add_vertex( 0.35, 0.35, 0.8, 1.2 );
  Vertex& v6 = domain.add_vertex( 0.35, 0.65, 0.8, 1.2 );
  Vertex& v7 = domain.add_vertex( 0.65, 0.65, 0.8, 1.2 );
  Vertex& v8 = domain.add_vertex( 0.65, 0.35, 0.8, 1.2 );

  b_ext.add_edge( v1, v2, 2 );
  b_ext.add_edge( v2, v3, 3 );
  b_ext.add_edge( v3, v4, 2 );
  b_ext.add_edge( v4, v1, 1 );

  // Build interior boundary
  b_int.add_edge( v5, v6, 4 );
  b_int.add_edge( v6, v7, 4 );
  b_int.add_edge( v7, v8, 4 );
  b_int.add_edge( v8, v5, 4 );

  // Create the mesh
  Mesh mesh { domain, 0, 0, 10.0 };
  mesh.init_advancing_front();

  mesh.create_quad_layers(v1, v2, 5, 0.005, 1.25);
  mesh.create_quad_layers(v3, v4, 5, 0.005, 1.25);
  mesh.create_quad_layers(v5, v5, 6, 0.005, 1.10);

  mesh.pave();

  /*------------------------------------------------------------------
  | Here we use a smoother, in order to improve the mesh quality.
  | The smoothing is applied for four iterations. 
  ------------------------------------------------------------------*/
  Smoother smoother {};
  smoother.smooth(domain, mesh, 4);

  /*------------------------------------------------------------------
  | Finally, the mesh is exportet to a file in VTU format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/aux/example_data/Example_2" };

  mesh.write_to_file( file_name, ExportType::vtu );


} // run_example_2()
