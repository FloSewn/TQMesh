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
#include "Boundary.h"
#include "Domain.h"


namespace BoundaryTests
{

using namespace TQMesh::TQUtils;
using namespace TQMesh::TQAlgorithm;

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
void Test_Domain_is_inside()
{
  Vertices   vertices {};
  Domain     domain { vertices };

  Boundary& b_ext = domain.add_boundary( BdryType::EXTERIOR );
  Boundary& b_int = domain.add_boundary( BdryType::INTERIOR );

  // Built exterior boundary
  Vertex& v1 = vertices.push_back(  0.0,  0.0 );
  Vertex& v2 = vertices.push_back(  8.0,  0.0 );
  Vertex& v3 = vertices.push_back( 16.0,  0.0 );
  Vertex& v4 = vertices.push_back( 16.0,  6.0 );
  Vertex& v5 = vertices.push_back( 16.0, 12.0 );
  Vertex& v6 = vertices.push_back(  8.0, 12.0 );
  Vertex& v7 = vertices.push_back(  0.0,  6.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v5, 1 );
  b_ext.add_edge( v5, v6, 1 );
  b_ext.add_edge( v6, v7, 1 );
  b_ext.add_edge( v7, v1, 1 );


  // Built interior boundary
  Vertex& v8  = vertices.push_back(  6.0,  4.0 );
  Vertex& v9  = vertices.push_back(  6.0,  8.0 );
  Vertex& v10 = vertices.push_back( 10.0,  8.0 );
  Vertex& v11 = vertices.push_back( 10.0,  4.0 );

  b_int.add_edge( v8,  v9,  2 );
  b_int.add_edge( v9,  v10, 2 );
  b_int.add_edge( v10, v11, 2 );
  b_int.add_edge( v11, v8,  2 );


  // Test if vertices are inside the domain
  Vertex& v_in    = vertices.push_back( 3.0, 2.0 );
  Vertex& v_out_1 = vertices.push_back( 8.0, 6.0 );
  Vertex& v_out_2 = vertices.push_back(-8.0, 6.0 );

  ASSERT( domain.is_inside( v_in ), 
          "Domain::is_inside() failed." );

  ASSERT( !(domain.is_inside( v_out_1 )), 
          "Domain::is_inside() failed." );

  ASSERT( !(domain.is_inside( v_out_2 )), 
          "Domain::is_inside() failed." );

  (void) v_in, v_out_1, v_out_2;


  DBG_MSG("Tests for Domain::is_inside() succeeded");
  
} // Test_Domain_is_inside()

/*********************************************************************
* Test interior / exterior boundary creation 
*
*   x---x---x
*   |       |
*   x---x---x
*
*********************************************************************/
void Test_Boundary_interior_exterior()
{
  Container<Vertex> vertices { };

  Vertex& v1 = vertices.push_back( 1.0, 1.0 );
  Vertex& v2 = vertices.push_back( 2.0, 1.0 );
  Vertex& v3 = vertices.push_back( 3.0, 1.0 );
  Vertex& v4 = vertices.push_back( 3.0, 2.0 );
  Vertex& v5 = vertices.push_back( 2.0, 2.0 );
  Vertex& v6 = vertices.push_back( 1.0, 2.0 );

  // Exterior boundaries are defined CCW
  Boundary extr_bdry { BdryType::EXTERIOR };

  extr_bdry.add_edge(v1,v2,1);
  extr_bdry.add_edge(v2,v3,1);
  extr_bdry.add_edge(v3,v4,2);
  extr_bdry.add_edge(v4,v5,3);
  extr_bdry.add_edge(v5,v6,3);
  extr_bdry.add_edge(v6,v1,4);

  // Check area
  ASSERT( EQ(extr_bdry.area(),2.0) , 
          "Edge::compute_area() failed." );

  // interior boundaries are defined CW
  Boundary intr_bdry { BdryType::INTERIOR };

  intr_bdry.add_edge(v1,v6,1);
  intr_bdry.add_edge(v6,v5,1);
  intr_bdry.add_edge(v5,v4,2);
  intr_bdry.add_edge(v4,v3,3);
  intr_bdry.add_edge(v3,v2,3);
  intr_bdry.add_edge(v2,v1,4);

  // Check area
  ASSERT( EQ(intr_bdry.area(),-2.0) , 
          "Edge::compute_area() failed." );

  DBG_MSG("Tests for interior / exterior boundaries succeeded");

} // Test_Boundary_interior_exterior()

/*********************************************************************
* Test clearing of all boundary edges
*
*   x---x---x
*   |       |
*   x---x---x
*
*********************************************************************/
void Test_Boundary_clear_edges()
{
  Container<Vertex> vertices { };

  Vertex& v1 = vertices.push_back( 1.0, 1.0 );
  Vertex& v2 = vertices.push_back( 2.0, 1.0 );
  Vertex& v3 = vertices.push_back( 3.0, 1.0 );
  Vertex& v4 = vertices.push_back( 3.0, 2.0 );
  Vertex& v5 = vertices.push_back( 2.0, 2.0 );
  Vertex& v6 = vertices.push_back( 1.0, 2.0 );

  // Exterior boundaries are defined CCW
  Boundary extr_bdry { BdryType::EXTERIOR };

  extr_bdry.add_edge(v1,v2,1);
  extr_bdry.add_edge(v2,v3,1);
  extr_bdry.add_edge(v3,v4,2);
  extr_bdry.add_edge(v4,v5,3);
  extr_bdry.add_edge(v5,v6,3);
  extr_bdry.add_edge(v6,v1,4);

  extr_bdry.clear_edges();

  // Assert that all edges are cleared
  ASSERT( (extr_bdry.size() == 0) , 
          "Boundary::clear_edges() failed." );

  // Assert that no vertex is connected to any edge
  for ( auto& v : vertices )
  {
    ASSERT( (v->edges().size() == 0),
          "Boundary::clear_edges() failed." );
    (void) v;
  }

  DBG_MSG("Tests for Boundary::clear_edges() succeeded");

} // Test_Boundary_clear_edges()


/*********************************************************************
* Test boundary shapes
*
*   x---x---x
*   |       |
*   x---x---x
*
*********************************************************************/
void Test_Boundary_shapes()
{
  Container<Vertex> vertices { };

  // Define exterior boundary with rectangular shape 
  Boundary extr_bdry { BdryType::EXTERIOR };

  // Define interior boundary with rectangular shape 
  Boundary intr_bdry { BdryType::INTERIOR };

} // Test_Boundary_shapes()


} // namespace BoundaryTests

/*********************************************************************
* 
*********************************************************************/
void run_boundary_tests()
{
  MSG("\n#===== EdgeList tests =====");
  
  BoundaryTests::Test_Boundary_interior_exterior();
  BoundaryTests::Test_Domain_is_inside();
  BoundaryTests::Test_Boundary_clear_edges();
  BoundaryTests::Test_Boundary_shapes();

} // run_boundary_tests()
