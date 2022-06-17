#include <iostream>
#include <vector>
#include <list>
#include <cassert>
#include <iomanip>   
#include <cmath>
#include <cstdlib>
#include <memory>       

#include "run_tests.h"

#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Mesh.h"
#include "Smoother.h"

namespace SmootherTests
{

using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* Test Smoother::smooth() for a pure triangle mesh
*********************************************************************/
void Test_Smoother_TriMesh(bool export_mesh)
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 1.5; 
  };

  Domain domain   { f, 50.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0, 0.1, 2.5 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0, 1.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0, 1.0 );
  Vertex& v4 = domain.add_vertex(  0.0,  5.0, 1.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 2 );
  b_ext.add_edge( v3, v4, 3 );
  b_ext.add_edge( v4, v1, 4 );

  // Create the mesh
  Mesh mesh { domain, 0, 50.0 };
  Smoother smoother {};

  mesh.pave();
  mesh.refine_to_quads();
  mesh.refine_to_quads();
  mesh.refine_to_quads();
  smoother.smooth(domain, mesh, 6, 0.5, 0.75, 0.95);

  // Export the mesh
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();

    std::cout << mesh;
    //domain.export_size_function({0.0,0.0}, {5.0,5.0}, 100, 100);
  }

  DBG_MSG("Tests for Smoother::smooth() succeeded");

} // Test_Smoother_TriMesh() */

} // namespace SmootherTests


/*********************************************************************
* 
*********************************************************************/
void run_smoother_tests()
{
  MSG("\n#===== Smoother tests =====");

  //SmootherTests::Test_Smoother_TriMesh(false);
}

