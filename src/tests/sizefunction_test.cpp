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

namespace SizeFunctionTests
{

using namespace TQMesh::TQUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* Test SizeFunction evaluation
*********************************************************************/
void Test_SizeFunction_evaluation(bool export_sizefun)
{


  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) { return 1.0 + 0.15*sqrt(p.x*p.y); };

  Domain           domain   { f };

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
  Vertex& v5 = domain.add_vertex(  2.5,  2.0, 0.1);
  Vertex& v6 = domain.add_vertex(  2.0,  3.5 );
  Vertex& v7 = domain.add_vertex(  3.0,  2.5 );
  Vertex& v8 = domain.add_vertex(  3.0,  2.0 );

  b_int.add_edge( v5, v6, 2 );
  b_int.add_edge( v6, v7, 2 );
  b_int.add_edge( v7, v8, 2 );
  b_int.add_edge( v8, v5, 2 );


  /*
  ASSERT( EQ( domain.size_function( {0.5,0.5} ), 1.0 ),
      "Domain size function evaluation failed." );
  */

  if (export_sizefun)
  {
    unsigned int v_index = 0;
    std::cout << "VERTICES " << domain.vertices().size() << std::endl;
    for ( const auto& v_ptr : domain.vertices() )
    {
      std::cout << std::setprecision(5) << std::fixed 
                << (*v_ptr).xy().x << "," 
                << (*v_ptr).xy().y << std::endl;
      (*v_ptr).index( v_index++ );
    }

    std::cout << "EDGES " << b_ext.size()  + b_int.size() << "\n";
    for ( const auto& e : b_ext )
      std::cout 
        << std::setprecision(0) << std::fixed 
        << std::setw(4) << e->v1().index() << "," 
        << std::setw(4) << e->v2().index() << ","
        << std::setw(4) << e->marker() << "\n";

    for ( const auto& e : b_int )
      std::cout 
        << std::setprecision(0) << std::fixed 
        << std::setw(4) << e->v1().index() << "," 
        << std::setw(4) << e->v2().index() << ","
        << std::setw(4) << e->marker() << "\n";

    domain.export_size_function({0.0,0.0}, {5.0,5.0}, 100, 100);
  }

  DBG_MSG("Tests for SizeFunction evaluation succeeded");


} // Test_SizeFunction_evaluation()


} // namespace SizeFunctionTests

/*********************************************************************
* 
*********************************************************************/
void run_sizefunction_tests()
{
  MSG("\n#===== SizeFunction tests =====");

  SizeFunctionTests::Test_SizeFunction_evaluation(false);

} // run_edge_list_tests()
