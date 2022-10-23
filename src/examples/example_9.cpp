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
void run_example_9()
{
  /*------------------------------------------------------------------
  |
  |       v3                                             v2     
  |      *----------------------------------------------* 
  |      |                              ^               |
  |      |       .---.                  |               |
  |      |      D|   |                  |H              |
  |      |       '---' - -              |               |
  |      |               |h             v               |
  |      *----------------------------------------------*
  |       v0                                             v1
  |
  ------------------------------------------------------------------*/
  int marker_inlet  = 1;
  int marker_walls  = 2;
  int marker_outlet = 3;
  int marker_circle = 4;

  const double H = 0.41; // [m]
  const double D = 0.10; // [m]
  const double h = 0.15; // [m] 
  const double L = 2.20; // [m]
  const double R = 0.5 * D;

  const double refinement = 2.00;
  const double growth = 1.5;

  const double elem_size = refinement * H / 25.;
  const double layer_height = 0.35 * elem_size;
  const int n_circ = static_cast<int>( M_PI / asin(0.5*elem_size/D) );

  int    n_layers = 1;
  double layer_fac = (1.0 - pow(growth, n_layers)) / (1.0 - growth);
  double first_height = layer_height / layer_fac;

  while ( first_height / elem_size > 0.2 )
  {
    ++n_layers;
    layer_fac = (1.0 - pow(growth, n_layers)) / (1.0 - growth);
    first_height = layer_height / layer_fac;
  }

  /*------------------------------------------------------------------
  | Define the size function
  ------------------------------------------------------------------*/
  UserSizeFunction f = [elem_size,R,h](const Vec2d& p) 
  { 
    return elem_size;
  };

  Domain domain   { f };


  Boundary& b_ext = domain.add_exterior_boundary();

  // Exterior boundary
  Vertex& v0 = domain.add_vertex(  0.0,  0.0, 1.0, H/2 );
  Vertex& v1 = domain.add_vertex(    L,  0.0, 1.0, H/2 );
  Vertex& v2 = domain.add_vertex(    L,    H, 1.0, H/2 );
  Vertex& v3 = domain.add_vertex(  0.0,    H, 1.0, H/2 );

  b_ext.add_edge( v0, v1, marker_walls );
  b_ext.add_edge( v1, v2, marker_outlet );
  b_ext.add_edge( v2, v3, marker_walls );
  b_ext.add_edge( v3, v0, marker_inlet );

  // Interior boundary - Circle
  Boundary& b_circ = domain.add_interior_boundary();
  b_circ.set_shape_circle( domain.vertices(), 
                           marker_circle, 
                           {h+R, h+R}, R, n_circ,
                           1.0, 2.5*layer_height);

  /*------------------------------------------------------------------
  | Initialize the mesh
  ------------------------------------------------------------------*/
  Mesh mesh { domain };
  mesh.init_advancing_front();

  /*------------------------------------------------------------------
  | Quad layers
  ------------------------------------------------------------------*/
  Vertex& v_circ = domain.vertices()[4];
  mesh.create_quad_layers(v_circ, v_circ, 
                          n_layers, first_height, growth);
  mesh.create_quad_layers(v0, v1, 
                          n_layers, first_height, growth);
  mesh.create_quad_layers(v2, v3,
                          n_layers, first_height, growth);

  /*------------------------------------------------------------------
  | Use fixed vertex to maintain triangle size around circle
  ------------------------------------------------------------------*/
  //domain.add_fixed_vertex(h+R, h+R, 0.5, 0.5*H);
  CONSTANTS.base_vertex_factor(1.0);

  /*------------------------------------------------------------------
  | Triangulate
  ------------------------------------------------------------------*/
  mesh.triangulate();

  /*------------------------------------------------------------------
  | Smooth the mesh for four iterations
  ------------------------------------------------------------------*/
  Smoother smoother {};
  smoother.smooth(domain, mesh, 4, 0.5);

  /*------------------------------------------------------------------
  | Finally, the mesh is exportet to a file in VTU format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/example_data/Example_9" };
  LOG(INFO) << "Writing mesh output to: " << file_name << ".vtu";

  mesh.write_to_file( file_name, ExportType::vtu );

} // run_example_9()
