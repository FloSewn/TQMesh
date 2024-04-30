/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/

#include "tests.h"
#include "TQMesh.h"

namespace TriangleTests 
{
using namespace CppUtils;
using namespace TQMesh;

/*********************************************************************
* Test initialization
*********************************************************************/
void initialization()
{
  Vertices   vertices { };
  Triangles  triangles { };

  Vertex& v1 = vertices.push_back(  0.0,  0.0 );
  Vertex& v2 = vertices.push_back(  2.0,  0.0 );
  Vertex& v3 = vertices.push_back(  2.0,  2.0 );

  Triangle& t1 = triangles.push_back( v1, v2, v3 );

  CHECK( (t1.get_edge_index(v1,v2) == 2) );
  CHECK( (t1.get_edge_index(v2,v1) == 2) );

  CHECK( (t1.get_edge_index(v2,v3) == 0) );
  CHECK( (t1.get_edge_index(v3,v2) == 0) );

  CHECK( (t1.get_edge_index(v3,v1) == 1) );
  CHECK( (t1.get_edge_index(v1,v3) == 1) );

  CHECK( triangles.size() == 1 );

  CHECK( (t1.v1() == v1 && t1.vertex(0) == v1) );
  CHECK( (t1.v2() == v2 && t1.vertex(1) == v2) );
  CHECK( (t1.v3() == v3 && t1.vertex(2) == v3) );

  CHECK( EQ( t1.area(), 2.0) );

  CHECK( EQ( t1.max_angle(), 0.5*M_PI) );
  CHECK( EQ( t1.min_angle(), 0.25*M_PI) );

  double dr1 = (t1.v1().xy()-t1.circumcenter()).norm();
  double dr2 = (t1.v1().xy()-t1.circumcenter()).norm();
  double dr3 = (t1.v1().xy()-t1.circumcenter()).norm();

  CHECK( EQ( dr1, t1.circumradius() ) ); 
  CHECK( EQ( dr2, t1.circumradius() ) ); 
  CHECK( EQ( dr3, t1.circumradius() ) ); 

  (void) v1,v2,v3;
  (void) t1;
  (void) dr1, dr2, dr3;

} // initialization()

/*********************************************************************
* Test vertex intersection
*********************************************************************/
void intersects_vertex()
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

  CHECK( t1.intersects_vertex(v4) );
  CHECK( t1.intersects_vertex(v5) );
  CHECK( !t1.intersects_vertex(v6) );
  CHECK( !t1.intersects_vertex(v3) );

  (void) v1,v2,v3,v4,v5,v6;
  (void) t1;

} // intersects_vertex()

/*********************************************************************
* Test domain intersection
*********************************************************************/
void intersects_domain()
{
  Triangles  triangles {};
  Domain     domain {};

  Boundary& b_ext = domain.add_exterior_boundary();
  Boundary& b_int = domain.add_interior_boundary();

  // Built exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  8.0,  0.0 );
  Vertex& v3 = domain.add_vertex(  8.0,  4.0 );
  Vertex& v4 = domain.add_vertex(  0.0,  4.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v1, 1 );

  // Built interior boundary
  Vertex& v5 = domain.add_vertex(  6.0,  1.0 );
  Vertex& v6 = domain.add_vertex(  6.0,  2.0 );
  Vertex& v7 = domain.add_vertex(  7.0,  2.0 );
  Vertex& v8 = domain.add_vertex(  7.0,  1.0 );

  b_int.add_edge( v5, v6, 2 );
  b_int.add_edge( v6, v7, 2 );
  b_int.add_edge( v7, v8, 2 );
  b_int.add_edge( v8, v5, 2 );

  // Built triangles
  Vertex& v_t11 = domain.add_vertex(  1.0,  1.0 );
  Vertex& v_t12 = domain.add_vertex(  3.0,  2.0 );
  Vertex& v_t13 = domain.add_vertex(  1.0,  2.0 );
  Triangle& t1 = triangles.push_back( v_t11, v_t12, v_t13);

  Vertex& v_t21 = domain.add_vertex(  2.0, -1.0 );
  Vertex& v_t22 = domain.add_vertex(  4.0, -1.0 );
  Vertex& v_t23 = domain.add_vertex(  3.0,  1.0 );
  Triangle& t2 = triangles.push_back( v_t21, v_t22, v_t23);

  Vertex& v_t31 = domain.add_vertex(  5.0,  1.0 );
  Vertex& v_t32 = domain.add_vertex(  7.0,  3.0 );
  Vertex& v_t33 = domain.add_vertex(  5.0,  3.0 );
  Triangle& t3 = triangles.push_back( v_t31, v_t32, v_t33);

  Vertex& v_t41 = domain.add_vertex(  3.0,  4.0 );
  Vertex& v_t42 = domain.add_vertex(  4.0,  3.0 );
  Vertex& v_t43 = domain.add_vertex(  4.0,  4.0 );
  Triangle& t4 = triangles.push_back( v_t41, v_t42, v_t43);


  CHECK( !t1.intersects_domain( domain ) );

  CHECK( t2.intersects_domain( domain ) );

  CHECK( !t3.intersects_domain( domain ) );

  CHECK( !t4.intersects_domain( domain ) );

  (void) t1,t2,t3,t4;

} // intersects_domain()

/*********************************************************************
* Test intersection with triangles
*********************************************************************/
void intersects_triangle()
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

  CHECK( !t1.intersects_triangle(triangles, 5.0) );

  Vertex& v7 = vertices.push_back( -1.0,  0.0 );
  Vertex& v8 = vertices.push_back(  1.0,  0.0 );
  Vertex& v9 = vertices.push_back( -2.0,  2.0 );

  Triangle& t3 = triangles.push_back( v7, v8, v9 );
  t3.is_active(true);

  CHECK( t1.intersects_triangle(triangles, 5.0) );

  (void) v1,v2,v3,v4,v5,v6;
  (void) t1,t2,t3;


} // intersects_triangle()

} // namespace TriangleTests


/*********************************************************************
* Run tests for: Triangle.h
*********************************************************************/
void run_tests_Triangle()
{
  TriangleTests::initialization();
  TriangleTests::intersects_vertex();
  TriangleTests::intersects_domain();
  TriangleTests::intersects_triangle();

} // run_tests_Triangle()
