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
  Domain domain { f };

  /*------------------------------------------------------------------
  | Now we define the exteriour boundary of the domain. The boundary 
  | itself is created from four connected edges, which in turn are 
  | defined by the four vertices v1, v2, v3 and v4.
  | Our final exterior boundary will look like this:
  |
  |                   v4                 v3
  |                     *<--------------*
  |                     |               ^
  |                     |               |
  |                     |               |
  |                     |               |
  |                     |               |
  |                     |               |
  |                     v               |
  |                     *-------------->*
  |                    v1                v2
  |
  | Notice the orientation of the edges (v1,v2), (v2,v3), (v3,v4) 
  | and (v4,v1). 
  | Exterior boundary edges must always be defined in counter-
  | clockwise direction.
  ------------------------------------------------------------------*/
  Boundary&  b_ext = domain.add_exterior_boundary();

  Vertex& v1 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0 );
  Vertex& v4 = domain.add_vertex(  0.0,  5.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v1, 1 );

  /*------------------------------------------------------------------
  | In this step, we will create an interior boundary, which has the 
  | shape of a triangle.
  |
  |                   v4                 v3
  |                     *<--------------*
  |                     |           v6  ^
  |                     |          *    |
  |                     |         /|    |
  |                     |       /  |    |
  |                     |     /    |    |
  |                     |    *-----*    |
  |                     v   v5      v7  |
  |                     *-------------->*
  |                    v1                v2
  |
  | This boundary is made up from the edges (v5,v6), (v6,v7), (v7,v5).
  | Interior boundary edges must always be defined in clockwise 
  | direction.
  | At vertex v5, we will refine the mesh locally by adjusting its 
  | sizing factor to 0.1. We also adjust the range to 1.5, in which 
  | this lower sizing will be applied. 
  ------------------------------------------------------------------*/
  Boundary&  b_int = domain.add_interior_boundary();

  double sizing = 0.1;
  double range  = 1.5;

  Vertex& v5 = domain.add_vertex(  1.5,  1.5, sizing, range );
  Vertex& v6 = domain.add_vertex(  1.5,  3.5 );
  Vertex& v7 = domain.add_vertex(  3.5,  3.5 );

  b_int.add_edge( v5, v6, 2 );
  b_int.add_edge( v6, v7, 2 );
  b_int.add_edge( v7, v5, 2 );

  /*------------------------------------------------------------------
  | In the following lines, the mesh will be initialized, as well 
  | as its advancing front structure (which is nescessary for the 
  | mesh generation)
  | After that, we simply triangulate the domain.
  ------------------------------------------------------------------*/
  Mesh mesh { domain };
  mesh.init_advancing_front();
  mesh.triangulate();

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
  { source_dir + "/aux/example_data/Example_1" };

  mesh.write_to_file( file_name, ExportType::txt );

} // run_example_1()
