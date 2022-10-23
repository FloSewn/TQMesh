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

#include "Testing.h"
#include "Vec2.h"
#include "Timer.h"
#include "Container.h"

#include "utils.h"
#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Front.h"

namespace FrontTests 
{
using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* Helper function to create the data structure that is needed for 
* the initialization of the advancing front
*********************************************************************/
FrontInitData collect_front_data(const Domain& domain)
{
  FrontInitData front_data {};

  for ( const auto& boundary : domain )
  {
    std::vector<Edge*> edges {};
    std::vector<bool>  is_oriented {};
    std::vector<int>   markers {};

    for ( const auto& e : boundary->edges() )
    {
      edges.push_back( e.get() ) ;
      is_oriented.push_back( true );
      markers.push_back( e->marker() );
    }

    front_data.edges.push_back( edges );
    front_data.is_oriented.push_back( is_oriented );
    front_data.markers.push_back( markers );
  }

  return std::move( front_data );
}


/*********************************************************************
* Test the initialization of the front
*********************************************************************/
void initialization()
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { return 1.0 + 0.15*sqrt(p.x*p.y); };

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
  Vertex& v5 = domain.add_vertex(  2.5,  2.0, 0.2);
  Vertex& v6 = domain.add_vertex(  2.0,  3.5 );
  Vertex& v7 = domain.add_vertex(  3.0,  2.5 );
  Vertex& v8 = domain.add_vertex(  3.0,  2.0 );

  b_int.add_edge( v5, v6, 2 );
  b_int.add_edge( v6, v7, 2 );
  b_int.add_edge( v7, v8, 2 );
  b_int.add_edge( v8, v5, 2 );

  // Collect data for the initialization of the advancing front
  FrontInitData front_data = collect_front_data( domain );

  // Advancing front requires initialized vertex container
  Vertices vertices { 10.0 };

  for ( const auto& v_ptr : domain.vertices() )
    vertices.push_back( v_ptr->xy(), v_ptr->sizing(), v_ptr->range() );
    
  // Create advancing front
  Front front { }; 
  front.init_front_edges( domain, front_data, vertices );
  
  CHECK( EQ(front.area(), domain.area()) );


  // Export the advancing front data 
  std::string source_dir { TQMESH_SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/test_data/FrontTests.initialization.txt" };

  std::ofstream outfile;
  outfile.open( file_name );

  unsigned int v_index = 0;
  outfile << "VERTICES " << vertices.size() << std::endl;
  for ( const auto& v_ptr : vertices )
  {
    outfile << std::setprecision(5) << std::fixed 
              << (*v_ptr).xy().x << "," 
              << (*v_ptr).xy().y << std::endl;
    (*v_ptr).index( v_index++ );
  }

  outfile << "EDGES " << front.size() << "\n";
  for ( const auto& e : front )
    outfile 
      << std::setprecision(0) << std::fixed 
      << std::setw(4) << e->v1().index() << "," 
      << std::setw(4) << e->v2().index() << ","
      << std::setw(4) << e->marker() << "\n";

  domain.export_size_function( outfile, {0.0,0.0}, {5.0,5.0}, 100, 100);

  outfile << "QTREE-LEAFS " << vertices.quad_tree().n_leafs() << std::endl;
  outfile << vertices.quad_tree();

  outfile.close();

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
  FrontInitData front_data = collect_front_data( domain );

  // Advancing front requires initialized vertex container
  Vertices vertices { 10.0 };

  for ( const auto& v_ptr : domain.vertices() )
    vertices.push_back( v_ptr->xy(), v_ptr->sizing(), v_ptr->range() );

  // Create advancing front
  Front front { };
  front.init_front_edges( domain, front_data, vertices );

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



} // namespace FrontTests


/*********************************************************************
* Run tests for: Front.h
*********************************************************************/
void run_tests_Front()
{
  FrontTests::initialization();
  FrontTests::sort_edges();

} // run_tests_Front()
