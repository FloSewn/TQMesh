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

namespace FrontTests
{

using namespace TQMesh::TQUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* Test Front initialization
*********************************************************************/
void Test_Front_initialization(bool export_data)
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) { return 1.0 + 0.15*sqrt(p.x*p.y); };

  Vertices         vertices { 10.0 };
  Domain           domain   { vertices, f };

  Boundary&  b_ext = domain.add_boundary( BdryType::EXTERIOR );
  Boundary&  b_int = domain.add_boundary( BdryType::INTERIOR );

  // Build exterior boundary
  Vertex& v1 = vertices.push_back(  0.0,  0.0 );
  Vertex& v2 = vertices.push_back(  5.0,  0.0 );
  Vertex& v3 = vertices.push_back(  5.0,  5.0 );
  Vertex& v4 = vertices.push_back(  0.0,  5.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v1, 1 );


  // Build interior boundary
  Vertex& v5 = vertices.push_back(  2.5,  2.0, 0.2);
  Vertex& v6 = vertices.push_back(  2.0,  3.5 );
  Vertex& v7 = vertices.push_back(  3.0,  2.5 );
  Vertex& v8 = vertices.push_back(  3.0,  2.0 );

  b_int.add_edge( v5, v6, 2 );
  b_int.add_edge( v6, v7, 2 );
  b_int.add_edge( v7, v8, 2 );
  b_int.add_edge( v8, v5, 2 );

  // Create advancing front
  Front front { domain, vertices };
  
  ASSERT( EQ(front.area(), domain.area()),
      "Front initialization failed.");

  if (export_data)
  {
    unsigned int v_index = 0;
    std::cout << "VERTICES " << vertices.size() << std::endl;
    for ( const auto& v_ptr : vertices )
    {
      std::cout << std::setprecision(5) << std::fixed 
                << (*v_ptr).xy().x << "," 
                << (*v_ptr).xy().y << std::endl;
      (*v_ptr).index( v_index++ );
    }

    std::cout << "EDGES " << front.size() << "\n";
    for ( const auto& e : front )
      std::cout 
        << std::setprecision(0) << std::fixed 
        << std::setw(4) << e->v1().index() << "," 
        << std::setw(4) << e->v2().index() << ","
        << std::setw(4) << e->marker() << "\n";

    domain.export_size_function({0.0,0.0}, {5.0,5.0}, 100, 100);

    std::cout << "QTREE-LEAFS " << vertices.qtree().n_leafs() << std::endl;
    std::cout << vertices.qtree();
  }

  DBG_MSG("Tests for Front initialization succeeded");


} // Test_Front_initialization()


/*********************************************************************
* Test Advacing Front edge sorting
*********************************************************************/
void Test_Front_sort_edges()
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) { return 1.0 + 0.15*sqrt(p.x*p.y); };

  Vertices         vertices { 10.0 };
  Domain           domain   { vertices, f };

  Boundary&  b_ext = domain.add_boundary( BdryType::EXTERIOR );
  Boundary&  b_int = domain.add_boundary( BdryType::INTERIOR );

  // Build exterior boundary
  Vertex& v1 = vertices.push_back(  0.0,  0.0 );
  Vertex& v2 = vertices.push_back(  5.0,  0.0 );
  Vertex& v3 = vertices.push_back(  5.0,  5.0 );
  Vertex& v4 = vertices.push_back(  0.0,  5.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v1, 1 );


  // Build interior boundary
  Vertex& v5 = vertices.push_back(  2.5,  2.0, 0.2);
  Vertex& v6 = vertices.push_back(  2.0,  3.5 );
  Vertex& v7 = vertices.push_back(  3.0,  2.5 );
  Vertex& v8 = vertices.push_back(  3.0,  2.0 );

  b_int.add_edge( v5, v6, 2 );
  b_int.add_edge( v6, v7, 2 );
  b_int.add_edge( v7, v8, 2 );
  b_int.add_edge( v8, v5, 2 );

  // Create advancing front
  Front front { domain, vertices };

  // Sort edges in ascending order
  front.sort_edges();

  // Check front edges to be arranged in ascending order
  double length = 0.0;
  for (const auto& e : front)
  {
    ASSERT( (length <= (*e).length()),
        "Front::sort_edges() failed.");
    length = (*e).length();
  }


  // Sort edges in descending order
  front.sort_edges(false);

  length = 1.0E+10;
  for (const auto& e : front)
  {
    ASSERT( (length >= (*e).length()),
        "Front::sort_edges() failed.");
    length = (*e).length();
  }

  DBG_MSG("Tests for Front::sort_edges() succeeded");

} // Test_Front_sort_edges()


} // namespace FrontTests

/*********************************************************************
* 
*********************************************************************/
void run_front_tests()
{
  MSG("\n#===== Front tests =====");

  FrontTests::Test_Front_initialization(false);
  FrontTests::Test_Front_sort_edges();

} // run_front_tests()
