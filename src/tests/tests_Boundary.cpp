/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/

#include "tests.h"
#include "TQMesh.h"

namespace BoundaryTests 
{
using namespace CppUtils;
using namespace TQMesh;

/*********************************************************************
* Test Domain:is_inside()
*
*   x-------x-------x
*   |               |
*   |     x---x     |
*   x     |   |     x
*   |     x---x     |
*   |               |
*   x-------x-------x
*
*********************************************************************/
void is_inside()
{
  Domain     domain {};

  Boundary& b_ext = domain.add_exterior_boundary();
  Boundary& b_int = domain.add_interior_boundary();

  // Built exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  8.0,  0.0 );
  Vertex& v3 = domain.add_vertex( 16.0,  0.0 );
  Vertex& v4 = domain.add_vertex( 16.0,  6.0 );
  Vertex& v5 = domain.add_vertex( 16.0, 12.0 );
  Vertex& v6 = domain.add_vertex(  8.0, 12.0 );
  Vertex& v7 = domain.add_vertex(  0.0,  6.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v5, 1 );
  b_ext.add_edge( v5, v6, 1 );
  b_ext.add_edge( v6, v7, 1 );
  b_ext.add_edge( v7, v1, 1 );

  auto extent = domain.extent();
  CHECK( EQ(extent.first[0],   0.0 ) ); // x-minimum 
  CHECK( EQ(extent.first[1],  16.0 ) ); // x-maximum
  CHECK( EQ(extent.second[0],  0.0 ) ); // y-minimum 
  CHECK( EQ(extent.second[1], 12.0 ) ); // y-maximum

  // Built interior boundary
  Vertex& v8  = domain.add_vertex(  6.0,  4.0 );
  Vertex& v9  = domain.add_vertex(  6.0,  8.0 );
  Vertex& v10 = domain.add_vertex( 10.0,  8.0 );
  Vertex& v11 = domain.add_vertex( 10.0,  4.0 );

  b_int.add_edge( v8,  v9,  2 );
  b_int.add_edge( v9,  v10, 2 );
  b_int.add_edge( v10, v11, 2 );
  b_int.add_edge( v11, v8,  2 );


  // Test if vertices are inside the domain
  Vertex& v_in    = domain.add_vertex( 3.0, 2.0 );
  Vertex& v_out_1 = domain.add_vertex( 8.0, 6.0 );
  Vertex& v_out_2 = domain.add_vertex(-8.0, 6.0 );

  extent = domain.extent();
  CHECK( EQ(extent.first[0],  -8.0 ) ); // x-minimum 
  CHECK( EQ(extent.first[1],  16.0 ) ); // x-maximum
  CHECK( EQ(extent.second[0],  0.0 ) ); // y-minimum 
  CHECK( EQ(extent.second[1], 12.0 ) ); // y-maximum

  CHECK( domain.is_inside( v_in ) );

  CHECK( !(domain.is_inside( v_out_1 )) ); 

  CHECK( !(domain.is_inside( v_out_2 )) );

  (void) v_in, v_out_1, v_out_2;

} // is_inside()

/*********************************************************************
* Test interior / exterior boundary creation 
*
*   x---x---x
*   |       |
*   x---x---x
*
*********************************************************************/
void interior_exterior()
{
  Container<Vertex> vertices { };

  Vertex& v1 = vertices.push_back( 1.0, 1.0 );
  Vertex& v2 = vertices.push_back( 2.0, 1.0 );
  Vertex& v3 = vertices.push_back( 3.0, 1.0 );
  Vertex& v4 = vertices.push_back( 3.0, 2.0 );
  Vertex& v5 = vertices.push_back( 2.0, 2.0 );
  Vertex& v6 = vertices.push_back( 1.0, 2.0 );

  // Exterior boundaries are defined CCW
  Boundary extr_bdry { vertices, BdryType::EXTERIOR };

  extr_bdry.add_edge(v1,v2,1);
  extr_bdry.add_edge(v2,v3,1);
  extr_bdry.add_edge(v3,v4,2);
  extr_bdry.add_edge(v4,v5,3);
  extr_bdry.add_edge(v5,v6,3);
  extr_bdry.add_edge(v6,v1,4);

  // Check area
  CHECK( EQ(extr_bdry.area(),2.0) );

  // interior boundaries are defined CW
  Boundary intr_bdry { vertices, BdryType::INTERIOR };

  intr_bdry.add_edge(v1,v6,1);
  intr_bdry.add_edge(v6,v5,1);
  intr_bdry.add_edge(v5,v4,2);
  intr_bdry.add_edge(v4,v3,3);
  intr_bdry.add_edge(v3,v2,3);
  intr_bdry.add_edge(v2,v1,4);

  // Check area
  CHECK( EQ(intr_bdry.area(),-2.0) );

} // interior_exterior()

/*********************************************************************
* Test clearing of all boundary edges
*
*   x---x---x
*   |       |
*   x---x---x
*
*********************************************************************/
void clear_edges()
{
  Container<Vertex> vertices { };

  Vertex& v1 = vertices.push_back( 1.0, 1.0 );
  Vertex& v2 = vertices.push_back( 2.0, 1.0 );
  Vertex& v3 = vertices.push_back( 3.0, 1.0 );
  Vertex& v4 = vertices.push_back( 3.0, 2.0 );
  Vertex& v5 = vertices.push_back( 2.0, 2.0 );
  Vertex& v6 = vertices.push_back( 1.0, 2.0 );

  // Exterior boundaries are defined CCW
  Boundary extr_bdry { vertices, BdryType::EXTERIOR };

  extr_bdry.add_edge(v1,v2,1);
  extr_bdry.add_edge(v2,v3,1);
  extr_bdry.add_edge(v3,v4,2);
  extr_bdry.add_edge(v4,v5,3);
  extr_bdry.add_edge(v5,v6,3);
  extr_bdry.add_edge(v6,v1,4);

  extr_bdry.clear_edges();

  // CHECK that all edges are cleared
  CHECK( (extr_bdry.size() == 0) );

  // Assert that no vertex is connected to any edge
  for ( auto& v : vertices )
  {
    CHECK( (v->edges().size() == 0) );
    (void) v;
  }

} // clear_edges()


/*********************************************************************
* Test boundary shapes
*
*   x---x---x
*   |       |
*   x---x---x
*
*********************************************************************/
void shapes()
{
  Container<Vertex> vertices { };

  // Define exterior boundary with rectangular shape 
  Boundary extr_bdry { vertices, BdryType::EXTERIOR };

  // Define interior boundary with rectangular shape 
  Boundary intr_bdry { vertices, BdryType::INTERIOR };

} // shapes()


} // namespace BoundaryTests


/*********************************************************************
* Run tests for: Boundary.h
*********************************************************************/
void run_tests_Boundary()
{
  BoundaryTests::interior_exterior();
  BoundaryTests::is_inside();
  BoundaryTests::clear_edges();
  BoundaryTests::shapes();

} // run_tests_Boundary()
