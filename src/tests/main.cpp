#include <iostream>
#include <cstdlib>

#include "run_tests.h"


/*********************************************************************
* The main function
*********************************************************************/
int main()
{
  run_geometry_tests();
  run_qtree_tests(false);
  run_container_tests(false);
  run_vertex_tests();
  run_edgelist_tests();
  run_boundary_tests();
  run_sizefunction_tests();
  run_front_tests();
  run_triangle_tests();
  run_quad_tests();
  run_mesh_tests(false);

  return EXIT_SUCCESS;
}
