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
* This example covers the usage of fixed interior edges
*********************************************************************/
bool fixed_edges()
{
  UserSizeFunction f = [](const Vec2d& p) { return 0.35; };

  Domain domain { f };

  // Vertices
  Vertex& v0 = domain.add_vertex( 1.0,  0.0, 0.05, 0.1 ); 
  Vertex& v1 = domain.add_vertex( 6.0,  0.5, 0.05, 0.1 ); 
  Vertex& v2 = domain.add_vertex( 4.0,  5.0, 0.05, 0.1 );
  Vertex& v3 = domain.add_vertex(-1.0,  4.5, 0.05, 0.1 );

  Boundary&  b_ext = domain.add_exterior_boundary();
  b_ext.add_edge( v0, v1, 1 ); // 1: Color for bottom edge
  b_ext.add_edge( v1, v2, 2 ); // 2: Color for right edge
  b_ext.add_edge( v2, v3, 3 ); // 3: Color for top edge
  b_ext.add_edge( v3, v0, 4 ); // 4: Color for left edge

  // Fixed vertices
  Vertex& v4 = domain.add_fixed_vertex(2.5, 2.5, 0.05, 1.0);
  Vertex& v5 = domain.add_fixed_vertex(1.5, 1.5, 0.05, 1.0);
  Vertex& v6 = domain.add_fixed_vertex(3.5, 1.5, 0.05, 1.0);
  Vertex& v7 = domain.add_fixed_vertex(3.5, 3.5, 0.05, 1.0);
  Vertex& v8 = domain.add_fixed_vertex(1.5, 3.5, 0.05, 1.0);

  // Define fixed edges 
  domain.add_fixed_edge( v4, v5 );
  domain.add_fixed_edge( v4, v6 );
  domain.add_fixed_edge( v4, v7 );
  domain.add_fixed_edge( v4, v8 );

  domain.add_fixed_edge( v5, v6 );
  domain.add_fixed_edge( v6, v7 );
  domain.add_fixed_edge( v7, v8 );
  domain.add_fixed_edge( v8, v5 );

  // It is possible to generate fixed edges from normal vertices
  // and fixed vertices
  domain.add_fixed_edge( v0, v5 );
  domain.add_fixed_edge( v0, v6 );
  domain.add_fixed_edge( v1, v6 );
  domain.add_fixed_edge( v1, v7 );
  domain.add_fixed_edge( v2, v7 );
  domain.add_fixed_edge( v2, v8 );
  domain.add_fixed_edge( v3, v8 );
  domain.add_fixed_edge( v3, v5 );

  /*------------------------------------------------------------------
  | Mesh generation 
  ------------------------------------------------------------------*/
  MeshGenerator generator {};
  Mesh& mesh = generator.new_mesh( domain );

  generator.triangulation(mesh).generate_elements();

  //generator.quad_refinement(mesh).refine();

  generator.mixed_smoothing(mesh)
    .epsilon(0.7) 
    .smooth(5);  

  /*------------------------------------------------------------------
  | Check if the meshing generation process succeeded
  ------------------------------------------------------------------*/
  MeshChecker checker { mesh, domain };
  if ( !checker.check_completeness() )
  {
    LOG(ERROR) << "Mesh generation failed";
    return false;
  }

  /*------------------------------------------------------------------
  | Finally, we export the mesh to a file in VTU / TXT format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/example_data/fixed_edges" };

  LOG(INFO) << "Writing mesh output to: " << file_name << ".vtu";
  generator.write_mesh(mesh, file_name, MeshExportType::VTU );

  LOG(INFO) << "Writing mesh output to: " << file_name << ".txt";
  generator.write_mesh(mesh, file_name, MeshExportType::TXT );


  return true;

} // fixed_edges()
