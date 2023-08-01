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
* This example covers the generation of a simple triangular mesh
*********************************************************************/
void run_example_1()
{
  /*------------------------------------------------------------------
  | First, we define the size function. This function describes
  | the size of the mesh elements with respect to their location 
  | in the domain. In this case, we use a constant size of 0.35
  ------------------------------------------------------------------*/
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 0.35;
  };

  /*------------------------------------------------------------------
  | Next, we need to define the domain. It requires the size function
  | as argument.
  ------------------------------------------------------------------*/
  Domain domain { f, 20.0 };

  /*------------------------------------------------------------------
  | Now we define the exteriour boundary of the domain. The boundary 
  | itself is created from four connected edges, which in turn are 
  | defined by the four vertices v0, v1, v2 and v3.
  | Our final exterior boundary will look like this:
  |
  |                   v3                 v2
  |                     *---------------*
  |                     |               |
  |                     |               |
  |                     |               |
  |                     |               |
  |                     |               |
  |                     |               |
  |                     |               |
  |                     *---------------*
  |                    v0                v1
  |
  | Notice the orientation of the edges (v0,v1), (v1,v2), (v2,v3) 
  | and (v3,v0). 
  | Exterior boundary edges must always be defined in counter-
  | clockwise direction.
  ------------------------------------------------------------------*/
  Boundary& b_ext = domain.add_exterior_boundary();

  Vertex& v0 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v1 = domain.add_vertex(  5.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  5.0,  5.0 );
  Vertex& v3 = domain.add_vertex(  0.0,  5.0 );

  b_ext.add_edge( v0, v1, 1 );
  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v0, 1 );

  /*------------------------------------------------------------------
  | In this step, we will create an interior boundary, which has the 
  | shape of a triangle.
  |
  |                   v3                 v2
  |                     *---------------*
  |                     |           v5  |
  |                     |          *    |
  |                     |         /|    |
  |                     |       /  |    |
  |                     |     /    |    |
  |                     |    *-----*    |
  |                     |   v4      v6  |
  |                     *---------------*
  |                    v0                v1
  |
  | This boundary is made up from the edges (v5,v6), (v6,v7), (v7,v5).
  | Interior boundary edges must always be defined in clockwise 
  | direction.
  | At vertex v4, we will refine the mesh locally by adjusting the 
  | size function locally to a value 0.05. This is applied within a 
  | distance of 0.2 to the vertex.
  ------------------------------------------------------------------*/
  Boundary& b_int = domain.add_interior_boundary();

  double size  = 0.05;
  double range = 0.2;

  Vertex& v4 = domain.add_vertex(  1.5,  1.5, size, range );
  Vertex& v5 = domain.add_vertex(  1.5,  3.5 );
  Vertex& v6 = domain.add_vertex(  3.5,  3.5 );

  b_int.add_edge( v4, v5, 2 );
  b_int.add_edge( v5, v6, 2 );
  b_int.add_edge( v6, v4, 2 );

  /*------------------------------------------------------------------
  | In the following lines, the mesh will be initialized  
  | triangulated.
  ------------------------------------------------------------------*/
  MeshGenerator generator {};
  Mesh& mesh = generator.new_mesh( domain );

  /*------------------------------------------------------------------
  | Here we generate the actual mesh elements using a triangulation
  ------------------------------------------------------------------*/
  generator.triangulation(mesh).generate_elements();

  /*------------------------------------------------------------------
  | Here we use a smoother, in order to improve the mesh quality.
  | The smoothing is applied for four iterations. 
  ------------------------------------------------------------------*/
  generator.mixed_smoothing(mesh).smooth(2);

  /*------------------------------------------------------------------
  | Finally, the mesh is exportet to a file in TXT format.
  ------------------------------------------------------------------*/
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string filename 
  { source_dir + "/auxiliary/example_data/Example_1.vtu" };
  LOG(INFO) << "Writing mesh output to: " << filename;

  generator.write_mesh(mesh, filename, MeshExportType::VTU);

} // run_example_1()
