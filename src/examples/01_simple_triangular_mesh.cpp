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
* This example covers the generation of a simple triangular mesh
*********************************************************************/
bool simple_triangular_mesh()
{

  /*------------------------------------------------------------------
  | First, we define the size function. This function describes
  | the size of the mesh elements with respect to their location 
  | in the domain. In this case, we use a constant size of 0.35
  ------------------------------------------------------------------*/
  UserSizeFunction f = [](const Vec2d& p) { return 0.35; };

  /*------------------------------------------------------------------
  | Next, we define the domain. It requires the size function
  | as argument. Before the domain gets instantinated, we set the 
  | scale of the underlying quadtree data structure to a value of 20.0 
  | This scale should generally be set slightly larger than the 
  | extent of the mesh that is going to be generated. It also has
  | to be done always prior to the generation of the domain. 
  ------------------------------------------------------------------*/
  TQMeshSetup::get_instance().set_quadtree_scale( 20.0 );
  Domain domain { f };

  /*------------------------------------------------------------------
  | Now we define the exteriour boundary of the domain. To do this, 
  | we provide it with a sequence of 2D vertex coordinates, which 
  | define the boundary in terms of a closed polygonal chain.
  |
  | Additionally, each boundary edge is associated to an integer, 
  | which we refer to as "colors". 
  | We provide them throgh an additional sequence of integers, that
  | correspond to the colors of the edges in the given polygonal 
  | chain.
  |
  | Our final exterior boundary will look like this:
  |
  |           (0.0, 5.0)                 (5.0, 5.0)
  |                     o---------------o
  |                     |      [2]      |
  |                     |               |
  |                     |               |
  |                     | [2]       [1] | 
  |                     |               |
  |                     |               |
  |                     |      [1]      |
  |                     o---------------o
  |           (0.0, 0.0)                 (5.0, 0.0)
  | 
  | The tuples denote boundary vertex coordinates and the numbers in 
  | brackets denote the corresponding edge colors.
  ------------------------------------------------------------------*/
  Boundary& b_ext = domain.add_exterior_boundary();

  std::vector<Vec2d> exterior_vertex_coordinates { 
    { 0.0, 0.0 },
    { 5.0, 0.0 },
    { 5.0, 5.0 },
    { 0.0, 5.0 } 
  };

  std::vector<int> exterior_edge_colors { 1, 1, 2, 2 };

  b_ext.set_shape_from_coordinates( exterior_vertex_coordinates, 
                                    exterior_edge_colors );

  /*------------------------------------------------------------------
  | In the next step, we will define an interior boundary of  
  | triangular shape: 
  |                     o---------------o
  |                     |               |
  |                     |          o    |
  |                     |         /|    |
  |                     |       /  |    |
  |                     |     /    |    |
  |                     |    o-----o    |
  |                     |               |
  |                     o---------------o
  |
  | This boundary is made up from the vertex sequence: 
  |       [(1.5,1.5), (1.5,3.5), (3.5,3.5)] 
  | All edges will obtain the color number "3". 
  | 
  | At vertex (1.5,1.5), we will refine the mesh locally. To do 
  | this, we provide an additional sequence of refinement properties,
  | where each entry corresponds to the local mesh size and a local
  | range scale in which the local mesh size will be of influence.
  | For vertices, where we do not want to apply a local mesh size,
  | we simply use a default values tuple of (0.0, 0.0)
  ------------------------------------------------------------------*/
  Boundary& b_int = domain.add_interior_boundary();

  std::vector<Vec2d> interior_vertex_coordinates { 
    { 1.5, 1.5 },
    { 1.5, 3.5 },
    { 3.5, 3.5 } 
  };

  std::vector<int> interior_edge_colors { 3, 3, 3 };

  std::vector<Vec2d> interior_vertex_properties { 
    { 0.05, 0.2 },
    { 0.00, 0.0 },
    { 0.00, 0.0 },
  };

  b_int.set_shape_from_coordinates( interior_vertex_coordinates, 
                                    interior_edge_colors,
                                    interior_vertex_properties );

  /*------------------------------------------------------------------
  | Now we are ready to initialize a new mesh structure
  ------------------------------------------------------------------*/
  MeshGenerator generator {};
  Mesh& mesh = generator.new_mesh( domain );

  /*------------------------------------------------------------------
  | We generate the actual mesh elements using a triangulation
  | algorithm
  ------------------------------------------------------------------*/
  generator.triangulation(mesh).generate_elements();

  /*------------------------------------------------------------------
  | Then we use a smoothing operation, in order to improve the 
  | mesh quality. The smoothing is applied for two iterations. 
  ------------------------------------------------------------------*/
  generator.mixed_smoothing(mesh).smooth(2);

  /*------------------------------------------------------------------
  | We can use a "MeshChecker" to verify if the triangulation 
  | succeeded
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
  std::string filename 
  { source_dir + "/auxiliary/example_data/simple_triangular_mesh" };

  LOG(INFO) << "Writing mesh output to: " << filename << ".vtu";
  generator.write_mesh(mesh, filename, MeshExportType::VTU);

  LOG(INFO) << "Writing mesh output to: " << filename << ".txt";
  generator.write_mesh(mesh, filename, MeshExportType::TXT);

  return true;

} // simple_triangular_mesh()
