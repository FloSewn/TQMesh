/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <string>

#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Boundary.h"
#include "SizeFunction.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* 
*********************************************************************/
class TestInitializer
{
public:

  TestInitializer(const std::string& testcase, 
                  UserSizeFunction f)
  : domain_ { f }
  {
    if      ( testcase == "UnitSquare" )
      init_UnitSquare();
    else if ( testcase == "UnitCircle" )
      init_UnitCircle();
    else
      TERMINATE("TestInitializer(): No valid test case defined.");

  }

  /*------------------------------------------------------------------
  | Getter
  ------------------------------------------------------------------*/
  const Domain& domain() const { return domain_; }
  Domain& domain() { return domain_; }

private:

  /*------------------------------------------------------------------
  | UnitSquare
  ------------------------------------------------------------------*/
  void init_UnitSquare()
  {
    domain_.quad_tree_scale( 1.5 );

    Boundary& b_ext = domain_.add_exterior_boundary();

    Vertex& v1 = domain_.add_vertex( -0.5, -0.5 );
    Vertex& v2 = domain_.add_vertex(  0.5, -0.5 );
    Vertex& v3 = domain_.add_vertex(  0.5,  0.5 );
    Vertex& v4 = domain_.add_vertex( -0.5,  0.5 );

    b_ext.add_edge( v1, v2, 1 );
    b_ext.add_edge( v2, v3, 2 );
    b_ext.add_edge( v3, v4, 3 );
    b_ext.add_edge( v4, v1, 4 );
  } 

  /*------------------------------------------------------------------
  | UnitCircle
  ------------------------------------------------------------------*/
  void init_UnitCircle()
  {
    domain_.quad_tree_scale( 2.5 );

    int marker        = 1;
    Vec2d center      = { 0.0, 0.0 };
    double radius     = 1.0;
    size_t n_segments = 40; 
    double mesh_size  = 1.0;
    double size_range = 1.0;

    Boundary& b_ext = domain_.add_exterior_boundary();

    b_ext.set_shape_circle(domain_.vertices(), marker, 
                           center, radius, n_segments, 
                           mesh_size, size_range );
  }

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Domain domain_;

}; // TestInitializer

} // namespace TQAlgorithm
} // namespace TQMesh
