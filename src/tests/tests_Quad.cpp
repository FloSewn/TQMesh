/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/

#include <iostream>
#include <cassert>

#include "tests.h"

#include "Testing.h"
#include "Vec2.h"
#include "Timer.h"
#include "Container.h"

#include "utils.h"
#include "Vec2.h"
#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Front.h"
#include "Quad.h"

namespace QuadTests 
{
using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* Test initialization
*********************************************************************/
void initialization()
{
  Vertices vertices { };
  Quads    quads { };

  Vertex& v1 = vertices.push_back(  0.0,  0.0 );
  Vertex& v2 = vertices.push_back(  2.0,  0.0 );
  Vertex& v3 = vertices.push_back(  2.0,  2.0 );
  Vertex& v4 = vertices.push_back(  0.0,  2.0 );

  Quad& q1 = quads.push_back( v1, v2, v3, v4 );

  CHECK( quads.size() == 1 );

  CHECK( (q1.v1() == v1 && q1.vertex(0) == v1) );
  CHECK( (q1.v2() == v2 && q1.vertex(1) == v2) );
  CHECK( (q1.v3() == v3 && q1.vertex(2) == v3) );
  CHECK( (q1.v4() == v4 && q1.vertex(3) == v4) );

  CHECK( EQ( q1.area(), 4.0) );

  CHECK( EQ( q1.max_angle(), 0.5*M_PI) );
  CHECK( EQ( q1.min_angle(), 0.5*M_PI) );

  (void) v1,v2,v3;
  (void) q1;

} // initialization()

} // namespace QuadTests


/*********************************************************************
* Run tests for: Quad.h
*********************************************************************/
void run_tests_Quad()
{
  QuadTests::initialization();

} // run_tests_Quad()
