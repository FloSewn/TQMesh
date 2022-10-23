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
void run_example_8()
{
  /*------------------------------------------------------------------
  | First, we define the size function. This function describes
  | the size of the mesh elements with respect to their location 
  | in the domain. In this case, we use a constant size of 0.35
  ------------------------------------------------------------------*/
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 0.25;
  };

  /*------------------------------------------------------------------
  | Next, we need to define the domain. It requires the size function
  | as argument.
  ------------------------------------------------------------------*/
  Domain domain { f };

  /*------------------------------------------------------------------
  | Exterior boundary
  ------------------------------------------------------------------*/
  const double r = 1.0;

  Boundary&  b_ext = domain.add_exterior_boundary();

  std::vector<Vertex*> vertices {};

  Vertex& v0 = domain.add_vertex( -2.0*r,  0.0*r );
  vertices.push_back( &v0 );

  const double delta_ang = M_PI / 32;
  for (double ang = 0.0; ang <= M_PI; ang += delta_ang)
  {
    const double x = r * cos(M_PI-ang);
    const double y = r * sin( ang);

    Vertex& vi = domain.add_vertex( x, y, 1.0, 0.5 );
    vertices.push_back( &vi );
  }

  Vertex& v1 = domain.add_vertex(  2.0*r,  0.0*r );
  vertices.push_back( &v1 );

  Vertex& v2 = domain.add_vertex(  2.0*r,  2.0*r );
  vertices.push_back( &v2 );

  Vertex& v3 = domain.add_vertex( -2.0*r,  2.0*r );
  vertices.push_back( &v3 );


  int n = static_cast<int>( vertices.size() );
  for (size_t i = 0; i < vertices.size(); ++i)
  {
    int j = static_cast<int>( i );
    int k = MOD(j+1, n);
    Vertex* v_start = vertices[j];
    Vertex* v_end = vertices[k];

    LOG(INFO) << "EDGE: " << v_start->xy() << " -> " << v_end->xy();

    b_ext.add_edge( *v_start, *v_end, 1 );
  }

  /*------------------------------------------------------------------
  | Generate mesh
  ------------------------------------------------------------------*/
  Mesh mesh { domain };
  mesh.init_advancing_front();

  mesh.create_quad_layers(*vertices[1], *vertices[n-4], 3, 0.01, 1.6);

  mesh.triangulate();

  /*------------------------------------------------------------------
  | Smoother 
  ------------------------------------------------------------------*/
  Smoother smoother {};
  smoother.smooth(domain, mesh, 4);

  /*------------------------------------------------------------------
  | Export
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/example_data/Example_8" };
  LOG(INFO) << "Writing mesh output to: " << file_name << ".txt";

  mesh.write_to_file( file_name, ExportType::txt );

} // run_example_8()
