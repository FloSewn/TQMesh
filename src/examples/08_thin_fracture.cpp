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
using namespace TQMesh;
/*********************************************************************
* 
*********************************************************************/
void thin_fracture()
{
  UserSizeFunction f = [](const Vec2d& p) { return 4.0; };

  Domain domain { f };

  // Vertices
  Vertex& v0 = domain.add_vertex( 55.00, -140.0000 ); 
  Vertex& v1 = domain.add_vertex( 90.00, -140.0000 ); 
  Vertex& v2 = domain.add_vertex( 90.00, -110.0000 );
  Vertex& v3 = domain.add_vertex( 55.00, -110.0000 );

  Vertex& v4 = domain.add_vertex( 77.03125, -134.0878, 0.1, 0.035 ); 
  Vertex& v5 = domain.add_vertex( 77.04852, -134.0778, 0.1, 0.035 ); 
  Vertex& v6 = domain.add_vertex( 68.49935, -119.4510, 0.1, 0.035 );
  Vertex& v7 = domain.add_vertex( 68.48208, -119.4611, 0.1, 0.035 );

  Boundary&  b_ext = domain.add_exterior_boundary();
  b_ext.add_edge( v0, v1, 1 ); // 1: Marker for bottom edge
  b_ext.add_edge( v1, v2, 2 ); // 2: Marker for right edge
  b_ext.add_edge( v2, v3, 3 ); // 3: Marker for top edge
  b_ext.add_edge( v3, v0, 4 ); // 4: Marker for left edge

  Boundary&  b_int = domain.add_interior_boundary();
  b_int.add_edge( v4, v7, 5 ); // Use marker 5 for all interior edges
  b_int.add_edge( v7, v6, 5 ); // ...
  b_int.add_edge( v6, v5, 5 );
  b_int.add_edge( v5, v4, 5 );


  MeshGenerator generator {};
  Mesh& mesh = generator.new_mesh( domain );

  // ... or apply everything directly   
  generator.quad_layer_generation(mesh)
    .n_layers(2)
    .first_height(0.1)
    .growth_rate(1.3)
    .starting_position(v4.xy())
    .ending_position(v4.xy())
    .generate_elements();


  generator.triangulation(mesh).generate_elements();
  generator.tri2quad_modification(mesh).modify();
  generator.quad_refinement(mesh).refine();

  generator.mixed_smoothing(mesh)
    .epsilon(0.7)                
    .smooth(3);  

  /*------------------------------------------------------------------
  | Finally, we export the mesh to a file in VTU / TXT format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/example_data/thin_fracture" };

  LOG(INFO) << "Writing mesh output to: " << file_name << ".vtu";
  generator.write_mesh(mesh, file_name, MeshExportType::VTU );

  LOG(INFO) << "Writing mesh output to: " << file_name << ".txt";
  generator.write_mesh(mesh, file_name, MeshExportType::TXT );


} // thin_fracture()
