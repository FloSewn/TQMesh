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
* This example involves the generation of the TQMesh banner
*********************************************************************/
void run_example_6()
{
  /*------------------------------------------------------------------
  | Define the size function
  ------------------------------------------------------------------*/
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    const double x0 = 6.0; 
    const double y0 = 3.0; 
    const double Rx = 15.0;
    const double Ry = 8.0;
    const double R  = 9.5;
    const double h  = 0.20;
    const double alpha = 0.4 * M_PI;

    double rho = 0.55;
    const double rho_min = 0.025;

    for (int i=0; i<5; i++)
    {
      const double idbl = (double)i;
      const double xi = x0 + Rx * cos(0.5*M_PI+idbl*alpha);
      const double yi = y0 + Ry * sin(0.5*M_PI+idbl*alpha);
      const double dx = p.x - xi;
      const double dy = p.y - yi;
      const double ri = sqrt(dx*dx +dy*dy);
      const double rho_i = h * fabs(ri - R) + rho_min;

      rho = MIN(rho, rho_i);
    }

    return rho;
  };

  Domain domain   { f };

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  Boundary&  b_ext = domain.add_exterior_boundary();
  Boundary&  b_T   = domain.add_interior_boundary();
  Boundary&  b_Q   = domain.add_interior_boundary();
  Boundary&  b_M   = domain.add_interior_boundary();

  /*------------------------------------------------------------------
  | Vertices for exterior boundary 
  ------------------------------------------------------------------*/
  Vertex& b0  = domain.add_vertex( -22.0, -10.00 );
  Vertex& b1  = domain.add_vertex(  22.0, -10.00 );
  Vertex& b2  = domain.add_vertex(  22.0,  10.00 );
  Vertex& b3  = domain.add_vertex( -22.0,  10.00 );

  b_ext.add_edge( b0, b1, 1 );
  b_ext.add_edge( b1, b2, 1 );
  b_ext.add_edge( b2, b3, 1 );
  b_ext.add_edge( b3, b0, 1 );

  /*------------------------------------------------------------------
  | Vertices for letter "T"
  ------------------------------------------------------------------*/
  Vertex& T0  = domain.add_vertex( -14.0,  -6.00 );
  Vertex& T1  = domain.add_vertex( -14.0,   3.00 );
  Vertex& T2  = domain.add_vertex( -18.0,   3.00 );
  Vertex& T3  = domain.add_vertex( -18.0,   6.00 );
  Vertex& T4  = domain.add_vertex(  -7.0,   6.00 );
  Vertex& T5  = domain.add_vertex(  -7.0,   3.00 );
  Vertex& T6  = domain.add_vertex( -11.0,   3.00 );
  Vertex& T7  = domain.add_vertex( -11.0,  -6.00 );

  b_T.add_edge( T0, T1, 2 );
  b_T.add_edge( T1, T2, 2 );
  b_T.add_edge( T2, T3, 2 );
  b_T.add_edge( T3, T4, 2 );
  b_T.add_edge( T4, T5, 2 );
  b_T.add_edge( T5, T6, 2 );
  b_T.add_edge( T6, T7, 2 );
  b_T.add_edge( T7, T0, 2 );

  /*------------------------------------------------------------------
  | Vertices for letter "Q"
  ------------------------------------------------------------------*/
  Vertex& Q0  = domain.add_vertex(  -3.0,  -6.00  );
  Vertex& Q1  = domain.add_vertex(  -5.0,  -5.00  );
  Vertex& Q2  = domain.add_vertex(  -5.0,   4.00  );
  Vertex& Q3  = domain.add_vertex(  -3.0,   6.00  );
  Vertex& Q4  = domain.add_vertex(   3.0,   6.00  );
  Vertex& Q5  = domain.add_vertex(   5.0,   4.00  );
  Vertex& Q6  = domain.add_vertex(   5.0,  -4.00  );
  Vertex& Q7  = domain.add_vertex(   2.0,   0.00  );
  Vertex& Q8  = domain.add_vertex(   2.0,   3.00  );
  Vertex& Q9  = domain.add_vertex(  -2.0,   3.00, 0.3, 0.4 );
  Vertex& Q10 = domain.add_vertex(  -2.0,  -3.00  );
  Vertex& Q11 = domain.add_vertex(   2.0,  -3.00  );
  Vertex& Q12 = domain.add_vertex(   4.0,  -6.00  );

  b_Q.add_edge( Q0,  Q1, 3 );
  b_Q.add_edge( Q1,  Q2, 3 );
  b_Q.add_edge( Q2,  Q3, 3 );
  b_Q.add_edge( Q3,  Q4, 3 );
  b_Q.add_edge( Q4,  Q5, 3 );
  b_Q.add_edge( Q5,  Q6, 3 );
  b_Q.add_edge( Q6,  Q7, 3 );
  b_Q.add_edge( Q7,  Q8, 3 );
  b_Q.add_edge( Q8,  Q9, 3 );
  b_Q.add_edge( Q9,  Q10, 3 );
  b_Q.add_edge( Q10, Q11, 3 );
  b_Q.add_edge( Q11, Q12, 3 );
  b_Q.add_edge( Q12, Q0, 3 );

  /*------------------------------------------------------------------
  | Vertices for letter "M"
  ------------------------------------------------------------------*/
  Vertex& M0  = domain.add_vertex(   7.0,  -6.00 );
  Vertex& M1  = domain.add_vertex(   7.0,   6.00 );
  Vertex& M2  = domain.add_vertex(  10.0,   6.00 );
  Vertex& M3  = domain.add_vertex(  12.0,   3.00 );
  Vertex& M4  = domain.add_vertex(  14.0,   3.00 );
  Vertex& M5  = domain.add_vertex(  16.0,   6.00 );
  Vertex& M6  = domain.add_vertex(  19.0,   6.00 );
  Vertex& M7  = domain.add_vertex(  19.0,  -6.00 );
  Vertex& M8  = domain.add_vertex(  16.0,  -6.00 );
  Vertex& M9  = domain.add_vertex(  16.0,   1.00 );
  Vertex& M10 = domain.add_vertex(  14.0,  -1.00 );
  Vertex& M11 = domain.add_vertex(  12.0,  -1.00 );
  Vertex& M12 = domain.add_vertex(  10.0,   1.00 );
  Vertex& M13 = domain.add_vertex(  10.0,  -6.00 );

  b_M.add_edge( M0,   M1, 4 );
  b_M.add_edge( M1,   M2, 4 );
  b_M.add_edge( M2,   M3, 4 );
  b_M.add_edge( M3,   M4, 4 );
  b_M.add_edge( M4,   M5, 4 );
  b_M.add_edge( M5,   M6, 4 );
  b_M.add_edge( M6,   M7, 4 );
  b_M.add_edge( M7,   M8, 4 );
  b_M.add_edge( M8,   M9, 4 );
  b_M.add_edge( M9,   M10, 4 );
  b_M.add_edge( M10,  M11, 4 );
  b_M.add_edge( M11,  M12, 4 );
  b_M.add_edge( M12,  M13, 4 );
  b_M.add_edge( M13,  M0, 4 );

  /*------------------------------------------------------------------
  | Initialize the mesh
  ------------------------------------------------------------------*/
  MeshGenerator generator {};
  Mesh& mesh = generator.new_mesh( domain );

  /*------------------------------------------------------------------
  | Generate some quad layers 
  ------------------------------------------------------------------*/
  generator.quad_layer_generation(mesh)
    .n_layers(3)
    .first_height(0.1)
    .growth_rate(1.5)
    .starting_position(b0.xy())
    .ending_position(b0.xy())
    .generate_elements();

  generator.quad_layer_generation(mesh)
    .n_layers(3)
    .first_height(0.1)
    .growth_rate(1.3)
    .starting_position(Q0.xy())
    .ending_position(Q0.xy())
    .generate_elements();

  generator.quad_layer_generation(mesh)
    .n_layers(3)
    .first_height(0.1)
    .growth_rate(1.3)
    .starting_position(M0.xy())
    .ending_position(M0.xy())
    .generate_elements();

  generator.quad_layer_generation(mesh)
    .n_layers(3)
    .first_height(0.1)
    .growth_rate(1.3)
    .starting_position(T0.xy())
    .ending_position(T0.xy())
    .generate_elements();

  /*------------------------------------------------------------------
  | Now generate the mesh elements
  ------------------------------------------------------------------*/
  generator.triangulation(mesh).generate_elements();

  /*------------------------------------------------------------------
  | Smooth the mesh for four iterations
  ------------------------------------------------------------------*/
  //generator.mixed_smoothing(mesh).smooth(4);
  generator.mixed_smoothing(mesh)
    .epsilon(0.7) // The smoothing strength 
    .quad_layer_smoothing(true) // Enable smoothing to quad layers
    .smooth(3); // The number of smoothing iterations

  /*------------------------------------------------------------------
  | Finally, the mesh is exportet to a file in TXT format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/example_data/Example_6" };

  LOG(INFO) << "Writing mesh output to: " << file_name << ".txt";
  generator.write_mesh(mesh, file_name, MeshExportType::TXT);

  LOG(INFO) << "Writing mesh output to: " << file_name << ".vtu";
  generator.write_mesh(mesh, file_name, MeshExportType::VTU);

} // run_example_6()
