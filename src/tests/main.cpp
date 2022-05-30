#include <iostream>
#include <cstdlib>

#include "run_tests.h"


/*********************************************************************
* The main function
*********************************************************************/
int main()
{
  run_vertex_tests();
  run_edgelist_tests();
  run_boundary_tests();
  run_sizefunction_tests();
  run_front_tests();
  run_triangle_tests();
  run_quad_tests();
  run_mesh_tests(false);
  run_smoother_tests();

  return EXIT_SUCCESS;
}
