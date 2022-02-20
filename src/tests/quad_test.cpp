#include <iostream>
#include <vector>
#include <list>
#include <cassert>
#include <iomanip>   
#include <cmath>
#include <cstdlib>
#include <memory>       

#include "run_tests.h"

#include "utils.h"
#include "Vec2.h"
#include "Timer.h"
#include "Container.h"

#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Front.h"
#include "Quad.h"

namespace QuadTests
{

using namespace TQMesh::TQUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* Test Quad constructor
*********************************************************************/
void Test_Quad_initialization()
{
  Vertices vertices { };
  Quads    quads { };

  Vertex& v1 = vertices.push_back(  0.0,  0.0 );
  Vertex& v2 = vertices.push_back(  2.0,  0.0 );
  Vertex& v3 = vertices.push_back(  2.0,  2.0 );
  Vertex& v4 = vertices.push_back(  0.0,  2.0 );

  Quad& q1 = quads.push_back( v1, v2, v3, v4 );

  ASSERT( quads.size() == 1,
      "Quad initialization failed.");

  ASSERT( (q1.v1() == v1 && q1.vertex(0) == v1),
      "Quad initialization failed.");
  ASSERT( (q1.v2() == v2 && q1.vertex(1) == v2),
      "Quad initialization failed.");
  ASSERT( (q1.v3() == v3 && q1.vertex(2) == v3),
      "Quad initialization failed.");
  ASSERT( (q1.v4() == v4 && q1.vertex(3) == v4),
      "Quad initialization failed.");

  ASSERT( EQ( q1.area(), 4.0),
      "Quad initialization failed.");

  ASSERT( EQ( q1.max_angle(), 0.5*M_PI),
      "Quad initialization failed.");
  ASSERT( EQ( q1.min_angle(), 0.5*M_PI),
      "Quad initialization failed.");

  (void) v1,v2,v3;
  (void) q1;

  DBG_MSG("Tests for Quad initialization succeeded");

} // Test_Quad_initialization()

} // namespace QuadTests

/*********************************************************************
* 
*********************************************************************/
void run_quad_tests()
{
  MSG("\n#===== Quad tests =====");

  QuadTests::Test_Quad_initialization();

} // run_front_tests()
