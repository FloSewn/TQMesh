/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <TQMeshConfig.h>
#include "STLHeaders.h"
#include "CppUtils.h"

#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Boundary.h"

namespace TQMesh {

using namespace CppUtils;

/*********************************************************************
* 
*********************************************************************/
class TestBuilder
{
public:

  TestBuilder(const std::string& testcase, 
              UserSizeFunction f)
  : domain_ { f }
  {
    if      ( testcase == "UnitSquare" )
      init_UnitSquare();
    else if ( testcase == "UnitCircle" )
      init_UnitCircle();
    else if ( testcase == "TriangleSquareCircle" )
      init_TriangleSquareCircle();
    else if ( testcase == "RefinedTriangle" )
      init_RefinedTriangle();
    else if ( testcase == "FixedVertices" )
      init_FixedVertices();
    else if ( testcase == "LakeSuperior" )
      init_LakeSuperior();
    else if ( testcase == "SharpStepAndSharpEdge" )
      init_SharpStepAndSharpEdge();
    else if ( testcase == "NormalStepAndSharpEdge" )
      init_NormalStepAndSharpEdge();
    else
      TERMINATE("TestBuilder(): No valid test case defined.");

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

  } // TestBuilder::init_UnitSquare()

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

    b_ext.set_shape_circle(marker, center, radius, n_segments, 
                           mesh_size, size_range );

  } // TestBuilder::init_UnitCircle()

  /*------------------------------------------------------------------
  | TriQuadCircle
  ------------------------------------------------------------------*/
  void init_TriangleSquareCircle()
  {
    domain_.quad_tree_scale( 20.0 );

    // Exterior boundary
    Boundary& b_ext = domain_.add_exterior_boundary();
    b_ext.set_shape_rectangle(1, {1.0, 1.0}, 10.0, 8.0);

    // Interior circle boundary
    Boundary& b_circle = domain_.add_interior_boundary();
    b_circle.set_shape_circle(2, {0.0,1.0}, 1.25, 40 );

    // Interior square boundary
    Boundary& b_square = domain_.add_interior_boundary();
    b_square.set_shape_square(3, {3.0,2.5}, 1.75);

    // Interior triangle boundary
    Boundary& b_triangle = domain_.add_interior_boundary();
    b_triangle.set_shape_triangle(4, {3.0,-0.5}, 1.75 );

  } // TestBuilder::init_TriangleSquareCircle()

  /*------------------------------------------------------------------
  | RefinedTriangle
  ------------------------------------------------------------------*/
  void init_RefinedTriangle()
  {
    domain_.quad_tree_scale( 10.0 );

    // Exterior boundary
    Boundary&  b_ext = domain_.add_exterior_boundary();
    b_ext.set_shape_square(1, {2.5,2.5}, 5.0);

    // Interior triangle boundary
    Boundary&  b_int = domain_.add_interior_boundary();

    Vertex& v1 = domain_.add_vertex(  1.5,  1.5, 0.025, 0.7 );
    Vertex& v2 = domain_.add_vertex(  1.5,  3.5 );
    Vertex& v3 = domain_.add_vertex(  3.5,  3.5 );

    b_int.add_edge( v1, v2, 2 );
    b_int.add_edge( v2, v3, 2 );
    b_int.add_edge( v3, v1, 2 );

  } // TestBuilder::init_RefinedTriangle()

  /*------------------------------------------------------------------
  | RefinedTriangle
  ------------------------------------------------------------------*/
  void init_FixedVertices()
  {
    domain_.quad_tree_scale( 10.0 );

    Boundary& b_ext = domain_.add_exterior_boundary();

    Vertex& v0 = domain_.add_vertex(  0.0,  0.0 );
    Vertex& v1 = domain_.add_vertex(  1.0,  0.0 );
    Vertex& v2 = domain_.add_vertex(  1.0,  2.0 );
    Vertex& v3 = domain_.add_vertex(  2.0,  2.0 );
    Vertex& v4 = domain_.add_vertex(  2.0,  1.0 );
    Vertex& v5 = domain_.add_vertex(  3.0,  1.0 );
    Vertex& v6 = domain_.add_vertex(  3.0,  3.0 );
    Vertex& v7 = domain_.add_vertex(  0.0,  3.0 );

    b_ext.add_edge( v0, v1, 1 );
    b_ext.add_edge( v1, v2, 1 );
    b_ext.add_edge( v2, v3, 1 );
    b_ext.add_edge( v3, v4, 1 );
    b_ext.add_edge( v4, v5, 1 );
    b_ext.add_edge( v5, v6, 1 );
    b_ext.add_edge( v6, v7, 1 );
    b_ext.add_edge( v7, v0, 1 );

    // Fixed vertices
    domain_.add_fixed_vertex(0.5, 1.0, 0.01, 0.7);
    domain_.add_fixed_vertex(0.5, 2.0, 0.01, 0.7);

  } // TestBuilder::init_FixedVertices()

  /*------------------------------------------------------------------
  | Lake Superior model
  ------------------------------------------------------------------*/
  void init_LakeSuperior()
  {
    domain_.quad_tree_scale( 20.0 );

    std::string source_directory { TQMESH_SOURCE_DIR };
    std::string filepath { source_directory + "/auxiliary/test_data/LakeSuperior.txt" };
    std::ifstream input_file( filepath );

    if ( !input_file.is_open() )
      TERMINATE("Error opening \"" + filepath + "\".");

    std::vector<Vec2d> coordinates;

    std::string line;

    while (std::getline(input_file, line)) 
    {
      std::stringstream ss(line);
      Vec2d coord;
      char comma;

      if (ss >> coord.x >> comma >> coord.y && comma == ',') 
        coordinates.push_back(coord);
      else
        LOG(ERROR) << "Error parsing line: " << line;
    }

    input_file.close();

    std::vector<Vertex*> v {};
    int n_vertices = static_cast<int>( coordinates.size() );

    // Initialize the boundary
    Boundary& b_ext = domain_.add_exterior_boundary();

    for ( int i = 0; i < n_vertices; ++i )
    {
      Vertex& v_new = domain_.add_vertex(coordinates[i].x, 
                                         coordinates[i].y);
      v.push_back( &v_new );
    }

    for ( int i = 0; i < n_vertices; ++i )
    {
      Vertex& v0 = *v[i];
      Vertex& v1 = *v[MOD(i+1,n_vertices)];
      b_ext.add_edge( v0, v1, 1 );
    }

  } // TestBuilder::init_LakeSuperior()


  /*------------------------------------------------------------------
  | A sharp step and a sharp edge
  ------------------------------------------------------------------*/
  void init_SharpStepAndSharpEdge()
  {
    domain_.quad_tree_scale( 30.0 );

    Boundary& b_ext = domain_.add_exterior_boundary();

    Vertex& v0 = domain_.add_vertex(  0.0,  0.0 );
    Vertex& v1 = domain_.add_vertex(  5.0,  0.0, 0.1, 0.1 );
    Vertex& v2 = domain_.add_vertex(  4.0,  5.0 );
    Vertex& v3 = domain_.add_vertex( 10.0,  3.0 );
    Vertex& v4 = domain_.add_vertex( 10.0, 10.0 );
    Vertex& v5 = domain_.add_vertex(  0.0,  5.0 );

    b_ext.add_edge( v0, v1, 1 );
    b_ext.add_edge( v1, v2, 2 );
    b_ext.add_edge( v2, v3, 2 );
    b_ext.add_edge( v3, v4, 3 );
    b_ext.add_edge( v4, v5, 4 );
    b_ext.add_edge( v5, v0, 4 );

  } // TestBuilder::init_SharpStepAndSharpEdge()

  /*------------------------------------------------------------------
  | A normal step and a sharp edge
  ------------------------------------------------------------------*/
  void init_NormalStepAndSharpEdge()
  {
    domain_.quad_tree_scale( 30.0 );

    Boundary& b_ext = domain_.add_exterior_boundary();

    Vertex& v0 = domain_.add_vertex(  0.0,  0.0 );
    Vertex& v1 = domain_.add_vertex(  5.0,  0.0 );
    Vertex& v2 = domain_.add_vertex(  5.0,  5.0 );
    Vertex& v3 = domain_.add_vertex( 10.0,  5.0 );
    Vertex& v4 = domain_.add_vertex( 10.0, 10.0 );
    Vertex& v5 = domain_.add_vertex(  0.0,  5.0 );

    b_ext.add_edge( v0, v1, 1 );
    b_ext.add_edge( v1, v2, 2 );
    b_ext.add_edge( v2, v3, 2 );
    b_ext.add_edge( v3, v4, 3 );
    b_ext.add_edge( v4, v5, 4 );
    b_ext.add_edge( v5, v0, 4 );

  } // TestBuilder::init_NormalStepAndSharpEdge()

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Domain domain_;

}; // TestBuilder

} // namespace TQMesh
