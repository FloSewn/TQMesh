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
#include "Triangle.h"

namespace TriangleTests
{

using namespace TQMesh::TQUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* Test Triangle constructor
*********************************************************************/
void Test_Triangle_initialization()
{
  Vertices   vertices { };
  Triangles  triangles { };

  Vertex& v1 = vertices.push_back(  0.0,  0.0 );
  Vertex& v2 = vertices.push_back(  2.0,  0.0 );
  Vertex& v3 = vertices.push_back(  2.0,  2.0 );

  Triangle& t1 = triangles.push_back( v1, v2, v3 );

  ASSERT( (t1.get_edge_index(v1,v2) == 2),
      "Triangle::get_edge_index() failed.");
  ASSERT( (t1.get_edge_index(v2,v1) == 2),
      "Triangle::get_edge_index() failed.");

  ASSERT( (t1.get_edge_index(v2,v3) == 0),
      "Triangle::get_edge_index() failed.");
  ASSERT( (t1.get_edge_index(v3,v2) == 0),
      "Triangle::get_edge_index() failed.");

  ASSERT( (t1.get_edge_index(v3,v1) == 1),
      "Triangle::get_edge_index() failed.");
  ASSERT( (t1.get_edge_index(v1,v3) == 1),
      "Triangle::get_edge_index() failed.");

  ASSERT( triangles.size() == 1,
      "Triangle initialization failed.");

  ASSERT( (t1.v1() == v1 && t1.vertex(0) == v1),
      "Triangle initialization failed.");
  ASSERT( (t1.v2() == v2 && t1.vertex(1) == v2),
      "Triangle initialization failed.");
  ASSERT( (t1.v3() == v3 && t1.vertex(2) == v3),
      "Triangle initialization failed.");

  ASSERT( EQ( t1.area(), 2.0),
      "Triangle initialization failed.");

  ASSERT( EQ( t1.max_angle(), 0.5*M_PI),
      "Triangle initialization failed.");
  ASSERT( EQ( t1.min_angle(), 0.25*M_PI),
      "Triangle initialization failed.");

  double dr1 = (t1.v1().xy()-t1.circumcenter()).length();
  double dr2 = (t1.v1().xy()-t1.circumcenter()).length();
  double dr3 = (t1.v1().xy()-t1.circumcenter()).length();

  ASSERT( EQ( dr1, t1.circumradius() ), 
      "Triangle initialization failed.");
  ASSERT( EQ( dr2, t1.circumradius() ), 
      "Triangle initialization failed.");
  ASSERT( EQ( dr3, t1.circumradius() ), 
      "Triangle initialization failed.");

  (void) v1,v2,v3;
  (void) t1;
  (void) dr1, dr2, dr3;

  DBG_MSG("Tests for Triangle initialization succeeded");

} // Test_Triangle_initialization()

/*********************************************************************
* Test Triangle::intersects_vertex()
*********************************************************************/
void Test_Triangle_intersects_vertex()
{
  Vertices   vertices { };
  Triangles  triangles { };

  Vertex& v1 = vertices.push_back(  0.0,  0.0 );
  Vertex& v2 = vertices.push_back(  2.0,  0.0 );
  Vertex& v3 = vertices.push_back(  2.0,  2.0 );

  Vertex& v4 = vertices.push_back(  1.0,  0.0 );
  Vertex& v5 = vertices.push_back(  1.0,  0.5 );
  Vertex& v6 = vertices.push_back(  0.0,  2.0 );

  Triangle& t1 = triangles.push_back( v1, v2, v3 );

  ASSERT( t1.intersects_vertex(v4),
      "Triangle::intersects_vertex failed.");
  ASSERT( t1.intersects_vertex(v5),
      "Triangle::intersects_vertex failed.");
  ASSERT( !t1.intersects_vertex(v6),
      "Triangle::intersects_vertex failed.");
  ASSERT( !t1.intersects_vertex(v3),
      "Triangle::intersects_vertex failed.");


  (void) v1,v2,v3,v4,v5,v6;
  (void) t1;


  DBG_MSG("Tests for Triangle::intersects_vertex() succeeded");

} // Test_Triangle_intersects_vertex()

/*********************************************************************
* Test Triangle::intersects_domain()
*********************************************************************/
void Test_Triangle_intersects_domain()
{
  Vertices   vertices {};
  Triangles  triangles {};
  Domain     domain { vertices };

  Boundary& b_ext = domain.add_boundary( BdryType::EXTERIOR );
  Boundary& b_int = domain.add_boundary( BdryType::INTERIOR );

  // Built exterior boundary
  Vertex& v1 = vertices.push_back(  0.0,  0.0 );
  Vertex& v2 = vertices.push_back(  8.0,  0.0 );
  Vertex& v3 = vertices.push_back(  8.0,  4.0 );
  Vertex& v4 = vertices.push_back(  0.0,  4.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v1, 1 );

  // Built interior boundary
  Vertex& v5 = vertices.push_back(  6.0,  1.0 );
  Vertex& v6 = vertices.push_back(  6.0,  2.0 );
  Vertex& v7 = vertices.push_back(  7.0,  2.0 );
  Vertex& v8 = vertices.push_back(  7.0,  1.0 );

  b_int.add_edge( v5, v6, 2 );
  b_int.add_edge( v6, v7, 2 );
  b_int.add_edge( v7, v8, 2 );
  b_int.add_edge( v8, v5, 2 );


  // Built triangles
  Vertex& v_t11 = vertices.push_back(  1.0,  1.0 );
  Vertex& v_t12 = vertices.push_back(  3.0,  2.0 );
  Vertex& v_t13 = vertices.push_back(  1.0,  2.0 );
  Triangle& t1 = triangles.push_back( v_t11, v_t12, v_t13);

  Vertex& v_t21 = vertices.push_back(  2.0, -1.0 );
  Vertex& v_t22 = vertices.push_back(  4.0, -1.0 );
  Vertex& v_t23 = vertices.push_back(  3.0,  1.0 );
  Triangle& t2 = triangles.push_back( v_t21, v_t22, v_t23);

  Vertex& v_t31 = vertices.push_back(  5.0,  1.0 );
  Vertex& v_t32 = vertices.push_back(  7.0,  3.0 );
  Vertex& v_t33 = vertices.push_back(  5.0,  3.0 );
  Triangle& t3 = triangles.push_back( v_t31, v_t32, v_t33);

  Vertex& v_t41 = vertices.push_back(  3.0,  4.0 );
  Vertex& v_t42 = vertices.push_back(  4.0,  3.0 );
  Vertex& v_t43 = vertices.push_back(  4.0,  4.0 );
  Triangle& t4 = triangles.push_back( v_t41, v_t42, v_t43);


  ASSERT( !t1.intersects_domain( domain ), 
          "Triangle::intersects_domain() failed." );

  ASSERT( t2.intersects_domain( domain ), 
          "Triangle::intersects_domain() failed." );

  ASSERT( !t3.intersects_domain( domain ), 
          "Triangle::intersects_domain() failed." );

  ASSERT( !t4.intersects_domain( domain ), 
          "Triangle::intersects_domain() failed." );

  (void) t1,t2,t3,t4;


  DBG_MSG("Tests for Triangle::intersects_domain() succeeded");

} // Test_Triangle_intersects_domain()


/*********************************************************************
* Test Triangle::intersects_triangle()
*********************************************************************/
void Test_Triangle_intersects_triangle()
{
  Vertices   vertices { };
  Triangles  triangles { };

  Vertex& v1 = vertices.push_back(  0.0,  0.0 );
  Vertex& v2 = vertices.push_back(  2.0,  0.0 );
  Vertex& v3 = vertices.push_back(  2.0,  2.0 );

  Vertex& v4 = vertices.push_back(  0.0,  1.0 );
  Vertex& v5 = vertices.push_back(  2.0,  3.0 );
  Vertex& v6 = vertices.push_back(  0.0,  4.0 );

  Triangle& t1 = triangles.push_back( v1, v2, v3 );
  Triangle& t2 = triangles.push_back( v4, v5, v6 );

  t1.is_active(true);
  t2.is_active(true);

  ASSERT( !t1.intersects_triangle(triangles, 5.0),
      "Triangle::intersects_triangle failed.");

  Vertex& v7 = vertices.push_back( -1.0,  0.0 );
  Vertex& v8 = vertices.push_back(  1.0,  0.0 );
  Vertex& v9 = vertices.push_back( -2.0,  2.0 );

  Triangle& t3 = triangles.push_back( v7, v8, v9 );
  t3.is_active(true);

  ASSERT( t1.intersects_triangle(triangles, 5.0),
      "Triangle::intersects_triangle failed.");

  (void) v1,v2,v3,v4,v5,v6;
  (void) t1,t2,t3;


  DBG_MSG("Tests for Triangle::intersects_triangle() succeeded");

} // Test_Triangle_intersects_triangle()

} // namespace TriangleTests

/*********************************************************************
* 
*********************************************************************/
void run_triangle_tests()
{
  MSG("\n#===== Triangle tests =====");

  TriangleTests::Test_Triangle_initialization();
  TriangleTests::Test_Triangle_intersects_vertex();
  TriangleTests::Test_Triangle_intersects_domain();
  TriangleTests::Test_Triangle_intersects_triangle();


} // run_front_tests()
