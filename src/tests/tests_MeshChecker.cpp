/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/

#include "tests.h"
#include "TestBuilder.h"
#include "TQMesh.h"

namespace CheckerTests 
{
using namespace CppUtils;
using namespace TQMesh;

/*********************************************************************
* 
*********************************************************************/
void final_check()
{
  bool success = true;

  UserSizeFunction f = [](const Vec2d& p) { return 0.5; };

  TestBuilder test_builder { "UnitSquare", f};

  Domain& domain = test_builder.domain();

  // Create the mesh
  MeshBuilder mesh_builder {};

  Mesh mesh = mesh_builder.create_empty_mesh( domain );
  success &= mesh_builder.prepare_mesh(mesh, domain);
  CHECK( success );

  TriangulationStrategy triangulation {mesh, domain};

  success &= triangulation.generate_elements();
  CHECK( success );

  MixedSmoothingStrategy smoother {mesh, domain};
  smoother.smooth(2);

  MeshChecker mesh_checker {mesh, domain};

  CHECK( mesh_checker.check_completeness() );

  LOG(DEBUG) << "\n" << mesh;

} // final_check()


} // namespace CheckerTests

/*********************************************************************
* Run tests for: MeshChecker.h
*********************************************************************/
void run_tests_MeshChecker()
{
  adjust_logging_output_stream("CheckerTests.final_check.log");
  CheckerTests::final_check();

  // Reset debug logging ostream
  adjust_logging_output_stream("COUT");

} // run_tests_MeshChecker()
