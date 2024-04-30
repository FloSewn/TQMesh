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

#include "VecND.h"
#include "MathUtility.h"
#include "Testing.h"

#include "Vertex.h"
#include "Triangle.h"
#include "Quad.h"

namespace VertexTests 
{
using namespace CppUtils;
using namespace TQMesh;

/*********************************************************************
* Test vertex constructor
*********************************************************************/
void constructor()
{
  Vertex v { 1.0, 1.0, 0.5, 0.1 };

  v.index(1);

  CHECK( EQ(v.mesh_size(), 0.5) );
  CHECK( EQ(v.size_range(), 0.1) );
  CHECK( v.index() == 1 );
  CHECK( !v.on_front() );
  CHECK( !v.on_boundary() );
  CHECK( !v.is_fixed() );

  CHECK( true );

} // constructor()


/*********************************************************************
* Test generation of quads / triangles
*********************************************************************/
void facets()
{
  Vertex v1 { 0.0, 0.0 };
  Vertex v2 { 1.0, 0.0 };
  Vertex v3 { 1.0, 1.0 };
  Vertex v4 { 0.0, 1.0 };

  Triangle t { v1, v2, v3 };
  Quad q { v1, v2, v3, v4 };

  CHECK( v1.is_adjacent(t) );
  CHECK( v2.is_adjacent(t) );
  CHECK( v3.is_adjacent(t) );
  CHECK(!v4.is_adjacent(t) );

  CHECK( v1.is_adjacent(q) );
  CHECK( v2.is_adjacent(q) );
  CHECK( v3.is_adjacent(q) );
  CHECK( v4.is_adjacent(q) );

} // facets()

/*********************************************************************
* Test vertex properties
*********************************************************************/
void properties()
{
  Vertex v { 0.0, 0.0 };
  CHECK( v.has_no_property() );

  v.add_property(VertexProperty::on_front);
  CHECK( !v.has_no_property() );
  CHECK( v.has_property(VertexProperty::on_front) );
  CHECK( !v.has_property(VertexProperty::on_boundary) );
  CHECK( !v.has_property(VertexProperty::is_fixed) );
  CHECK( !v.has_property(VertexProperty::in_quad_layer) );

  v.add_property(VertexProperty::is_fixed);
  CHECK( !v.has_no_property() );
  CHECK( v.has_property(VertexProperty::on_front) );
  CHECK( !v.has_property(VertexProperty::on_boundary) );
  CHECK( v.has_property(VertexProperty::is_fixed) );
  CHECK( !v.has_property(VertexProperty::in_quad_layer) );

  v.remove_property(VertexProperty::is_fixed);
  CHECK( !v.has_no_property() );
  CHECK( v.has_property(VertexProperty::on_front) );
  CHECK( !v.has_property(VertexProperty::on_boundary) );
  CHECK( !v.has_property(VertexProperty::is_fixed) );
  CHECK( !v.has_property(VertexProperty::in_quad_layer) );

  v.add_property(VertexProperty::on_boundary);
  v.add_property(VertexProperty::in_quad_layer);
  v.add_property(VertexProperty::is_fixed);
  v.remove_property(VertexProperty::on_front);

  CHECK( !v.has_no_property() );
  CHECK( v.has_property(VertexProperty::on_boundary) );
  CHECK( v.has_property(VertexProperty::in_quad_layer) );
  CHECK( v.has_property(VertexProperty::is_fixed) );
  CHECK( !v.has_property(VertexProperty::on_front) );


  Vertex b { 1.0, 1.0 };
  b.set_property( v.properties() );
  CHECK( !b.has_no_property() );
  CHECK( b.has_property(VertexProperty::on_boundary) );
  CHECK( b.has_property(VertexProperty::in_quad_layer) );
  CHECK( b.has_property(VertexProperty::is_fixed) );
  CHECK( !b.has_property(VertexProperty::on_front) );

  b.remove_property(VertexProperty::on_boundary);
  b.remove_property(VertexProperty::is_fixed);
  b.remove_property(VertexProperty::in_quad_layer);
  CHECK( b.has_no_property() );


  b.remove_property(VertexProperty::is_fixed);
  CHECK( b.has_no_property() );
  CHECK( !b.has_property(VertexProperty::on_front) );
  CHECK( !b.has_property(VertexProperty::on_boundary) );
  CHECK( !b.has_property(VertexProperty::is_fixed) );
  CHECK( !b.has_property(VertexProperty::in_quad_layer) );

  b.add_property( v.properties() );

  CHECK( !b.has_no_property() );
  CHECK( b.has_property(VertexProperty::on_boundary) );
  CHECK( b.has_property(VertexProperty::in_quad_layer) );
  CHECK( b.has_property(VertexProperty::is_fixed) );
  CHECK( !b.has_property(VertexProperty::on_front) );

} // properties()

} // namespace VertexTests


/*********************************************************************
* Run tests for: Vertex.h
*********************************************************************/
void run_tests_Vertex()
{
  VertexTests::constructor();
  VertexTests::facets();
  VertexTests::properties();

} // run_tests_Vertex()
