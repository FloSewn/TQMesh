/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/

#include <iostream>
#include <fstream>
#include <cassert>

#include <TQMeshConfig.h>

#include "tests.h"
#include "TestInitializer.h"

#include "Testing.h"
#include "VecND.h"
#include "Timer.h"
#include "Container.h"

#include "utils.h"
#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Front.h"
#include "Mesh.h"
#include "FrontInitializer.h"

namespace FrontTests 
{
using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* 
*********************************************************************/
static inline void export_mesh_file(const Vertices& vertices,
                                    const Front& front,
                                    const Domain& domain)
{
  unsigned int v_index = 0;
  LOG_PROPERTIES.set_info_header("");
  LOG(INFO) << "MESH 0";
  LOG(INFO) << "VERTICES " << vertices.size();
  for ( const auto& v_ptr : vertices )
  {
    LOG(INFO) << std::setprecision(5) << std::fixed 
              << (*v_ptr).xy().x << "," << (*v_ptr).xy().y;
    (*v_ptr).index( v_index++ );
  }
  LOG(INFO) << "INTERIOREDGES 0";

  LOG(INFO) << "BOUNDARYEDGES " << front.size();
  for ( const auto& e : front )
    LOG(INFO) 
      << std::setprecision(0) << std::fixed 
      << std::setw(4) << e->v1().index() << "," 
      << std::setw(4) << e->v2().index() << ","
      << std::setw(4) << -1 << ","
      << std::setw(4) << e->marker();

  LOG(INFO) << "INTERFACEEDGES 0";
  LOG(INFO) << "FRONT 0";
  LOG(INFO) << "QUADS 0";
  LOG(INFO) << "TRIANGLES 0";
  LOG(INFO) << "QUADNEIGHBORS 0";
  LOG(INFO) << "TRIANGLENEIGHBORS 0";
  LOG(INFO) << "SIZEFUNCTION " << vertices.size();;
  for ( const auto& v_ptr : vertices )
    LOG(INFO) << std::setprecision(5) << std::fixed 
              << domain.size_function(v_ptr->xy());
  LOG_PROPERTIES.set_info_header("  ");

} // export_mesh_file()


/*********************************************************************
* Test the initialization of the front with a unit squre domain
*********************************************************************/
void test_UnitSquare()
{
  UserSizeFunction f = [](const Vec2d& p) { return 1.0; };

  TestInitializer initializer { "UnitSquare", f};

  // Collect data for the initialization of the advancing front
  FrontInitializer front_data { initializer.domain() };

  // Advancing front requires initialized vertex container
  Vertices vertices { 1.5 };
    
  // Create advancing front
  Front front { }; 
  front.init_front( initializer.domain(), front_data, vertices );
  
  CHECK( EQ(front.area(), initializer.domain().area()) );
  CHECK( EQ(front.area(), 1.0) );
  CHECK( front.edges().size() == 4 );
  for ( const auto& e_ptr : front.edges() )
    CHECK( EQ(e_ptr->length(), 1.0) );

  export_mesh_file(vertices, front, initializer.domain());

} // test_UnitSquare()

/*********************************************************************
* Test the initialization of the front with a unit squre domain
*********************************************************************/
void test_UnitCircle()
{
  UserSizeFunction f = [](const Vec2d& p) { return 1.0; };

  TestInitializer initializer { "UnitCircle", f };

  // Collect data for the initialization of the advancing front
  FrontInitializer front_data { initializer.domain() };

  // Advancing front requires initialized vertex container
  Vertices vertices { 2.5 };
    
  // Create advancing front
  Front front { }; 
  front.init_front( initializer.domain(), front_data, vertices );
  
  //CHECK( EQ(front.area(), initializer.domain().area()) );
  //CHECK( EQ(front.area(), 1.0) );
  //CHECK( front.edges().size() == 4 );
  //for ( const auto& e_ptr : front.edges() )
  //  CHECK( EQ(e_ptr->length(), 1.0) );

  export_mesh_file(vertices, front, initializer.domain());

} // test_UnitCircle()


/*********************************************************************
* Test the initialization of the front
*********************************************************************/
void initialization()
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { return 5.0; };
  //{ return 1.0 + 0.15*sqrt(p.x*p.x); };

  Domain domain { f, 10.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();
  Boundary&  b_int = domain.add_interior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0 );
  Vertex& v4 = domain.add_vertex(  0.0,  5.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v1, 1 );


  // Build interior boundary
  Vertex& v5 = domain.add_vertex(  2.5,  2.0 );
  Vertex& v6 = domain.add_vertex(  2.0,  3.5, 0.1 );
  Vertex& v7 = domain.add_vertex(  3.0,  2.5 );
  Vertex& v8 = domain.add_vertex(  3.0,  2.0 );

  b_int.add_edge( v5, v6, 2 );
  b_int.add_edge( v6, v7, 2 );
  b_int.add_edge( v7, v8, 2 );
  b_int.add_edge( v8, v5, 2 );

  // Collect data for the initialization of the advancing front
  FrontInitializer front_data { domain };

  // Advancing front requires initialized vertex container
  Vertices vertices { 10.0 };
    
  // Create advancing front
  Front front { }; 
  front.init_front( domain, front_data, vertices );
  
  CHECK( EQ(front.area(), domain.area()) );

  export_mesh_file(vertices, front, domain);

} // initialization()

/*********************************************************************
* Test Advacing Front edge sorting
*********************************************************************/
void sort_edges()
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { return 1.0 + 0.15*sqrt(p.x*p.y); };

  Domain           domain   { f, 10.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();
  Boundary&  b_int = domain.add_interior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0 );
  Vertex& v4 = domain.add_vertex(  0.0,  5.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v1, 1 );


  // Build interior boundary
  Vertex& v5 = domain.add_vertex(  2.5,  2.0, 0.2);
  Vertex& v6 = domain.add_vertex(  2.0,  3.5 );
  Vertex& v7 = domain.add_vertex(  3.0,  2.5 );
  Vertex& v8 = domain.add_vertex(  3.0,  2.0 );

  b_int.add_edge( v5, v6, 2 );
  b_int.add_edge( v6, v7, 2 );
  b_int.add_edge( v7, v8, 2 );
  b_int.add_edge( v8, v5, 2 );

  // Collect data for the initialization of the advancing front
  std::vector<Mesh*> dummy {};
  FrontInitializer front_data { domain, dummy };

  // Advancing front requires initialized vertex container
  Vertices vertices { 10.0 };

  // Create advancing front
  Front front { };
  front.init_front( domain, front_data, vertices );

  // Sort edges in ascending order
  front.sort_edges();

  // Check front edges to be arranged in ascending order
  double length = 0.0;
  for (const auto& e : front)
  {
    CHECK( (length <= (*e).length()) );
    length = (*e).length();
  }


  // Sort edges in descending order
  front.sort_edges(false);

  length = 1.0E+10;
  for (const auto& e : front)
  {
    CHECK( (length >= (*e).length()) );
    length = (*e).length();
  }

} // sort_edges()


/*********************************************************************
* Test the Advacing Front edge size
*********************************************************************/
void edge_size()
{
  const double edge_size = 0.23;

  // Define a variable size function
  UserSizeFunction f = [edge_size](const Vec2d& p) 
  { return edge_size; };

  Domain           domain   { f, 10.0 };

  // Build exterior boundary
  Boundary&  b_ext = domain.add_exterior_boundary();

  Vertex& v1 = domain.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain.add_vertex(  1.0,  0.0 );
  Vertex& v3 = domain.add_vertex(  1.0,  1.0 );
  Vertex& v4 = domain.add_vertex(  0.0,  1.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v1, 1 );


  // Collect data for the initialization of the advancing front
  std::vector<Mesh*> dummy {};
  FrontInitializer front_data { domain, dummy };

  // Advancing front requires initialized vertex container
  Vertices vertices {};

  // Create advancing front
  Front front { };
  front.init_front( domain, front_data, vertices );

  // Check if all edge lengths are more or less in accordance to
  // the global size parameter
  for ( const auto& e : front )
  {
    const double error = ABS(e->length() - edge_size) / edge_size;
    CHECK( error < 0.25 );
  }

} // edge_size()


} // namespace FrontTests


/*********************************************************************
* Run tests for: Front.h
*********************************************************************/
void run_tests_Front()
{
  adjust_logging_output_stream("FrontTests.test_UnitSquare.log");
  FrontTests::test_UnitSquare();

  adjust_logging_output_stream("FrontTests.test_UnitCirlce.log");
  FrontTests::test_UnitCircle();

  adjust_logging_output_stream("FrontTests.initialization.log");
  FrontTests::initialization();

  adjust_logging_output_stream("FrontTests.sort_edges.log");
  FrontTests::sort_edges();

  adjust_logging_output_stream("FrontTests.edge_size.log");
  FrontTests::edge_size();

  // Reset debug logging ostream
  adjust_logging_output_stream("COUT");

} // run_tests_Front()
