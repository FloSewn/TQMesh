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

#include "Vec2.h"
#include "Vertex.h"
#include "Testing.h"

namespace VertexTests 
{
using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* Test 
*********************************************************************/
void constructor()
{
  Vertex v1 { 1.0, 1.0 };

  CHECK( true );

} // constructor()


} // namespace VertexTests


/*********************************************************************
* Run tests for: Vertex.h
*********************************************************************/
void run_tests_Vertex()
{
  VertexTests::constructor();

} // run_tests_Vertex()
