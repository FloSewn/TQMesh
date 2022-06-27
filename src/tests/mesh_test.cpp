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
#include "Mesh.h"

namespace MeshTests
{

using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/*********************************************************************
* Container for benchmark data
*********************************************************************/
struct BenchmarkContainer
{
  std::vector<double>   h;
  std::vector<double>   L;

  std::vector<int>      n_vertices;
  std::vector<int>      n_tris;
  std::vector<int>      n_quads;
  std::vector<int>      n_intr_edges;
  std::vector<int>      n_bdry_edges;

  std::vector<double>   quad_layer_time;
  std::vector<double>   meshing_time;
  std::vector<double>   smoothing_time;

}; // BenchmarkContainer


/*********************************************************************
* Test mesh initialization
*********************************************************************/
void Test_Mesh_initialization()
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) { return 1.0 + 0.15*sqrt(p.x*p.y); };

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

  // Create the mesh
  Mesh mesh { domain, 0, 0, 10.0 };

  ASSERT( EQ(mesh.front().area(), domain.area()),
      "Front initialization failed.");

  DBG_MSG("Tests for Mesh initialization succeeded");

} // Test_Mesh_initialization() */

/*********************************************************************
* Test Mesh::triangulate()
*********************************************************************/
void Test_Mesh_triangulate(bool export_mesh)
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 2.5; 
  };

  Domain domain   { f, 20.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0, 0.9 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0, 0.9 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0, 0.9 );
  Vertex& v4 = domain.add_vertex( 10.0,  5.0, 0.9 );
  Vertex& v5 = domain.add_vertex( 10.0, 10.0, 0.9 );
  Vertex& v6 = domain.add_vertex(  0.0,  5.0, 0.9 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v5, 1 );
  b_ext.add_edge( v5, v6, 1 );
  b_ext.add_edge( v6, v1, 1 );

  // Create the mesh
  Mesh mesh { domain, 0, 0, 50.0 };

  mesh.triangulate();

  // Assertions
  ASSERT( EQ(mesh.area(), domain.area()),
      "Mesh::mesh_area() failed.");

  // Export the mesh
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();

    std::cout << mesh;
    //domain.export_size_function({0.0,0.0}, {5.0,5.0}, 100, 100);
  }

  DBG_MSG("Tests for Mesh::triangulate() succeeded");

} // Test_Mesh_triangulate() */

/*********************************************************************
* Test Mesh::pave()
*********************************************************************/
void Test_Mesh_pave(bool export_mesh)
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 0.5; 
  };

  Domain domain   { f, 20.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0, 0.9 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0, 0.9 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0, 0.9 );
  Vertex& v4 = domain.add_vertex( 10.0,  5.0, 0.9 );
  Vertex& v5 = domain.add_vertex( 10.0, 10.0, 0.9 );
  Vertex& v6 = domain.add_vertex(  0.0,  5.0, 0.9 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v5, 1 );
  b_ext.add_edge( v5, v6, 1 );
  b_ext.add_edge( v6, v1, 1 );

  // Create the mesh
  Mesh mesh { domain, 0, 0, 50.0 };

  mesh.pave();

  // Assertions
  ASSERT( EQ(mesh.area(), domain.area()),
      "Mesh::mesh_area() failed.");

  // Export the mesh
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();

    std::cout << mesh;
    //domain.export_size_function({0.0,0.0}, {5.0,5.0}, 100, 100);
  }

  DBG_MSG("Tests for Mesh::pave() succeeded");

} // Test_Mesh_pave() */


/*********************************************************************
* Test Mesh::triangulate_quad_layer()
*********************************************************************/
void Test_Mesh_refine_to_quads(bool export_mesh)
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 2.5; 
  };

  Domain domain   { f, 50.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0, 1.0 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0, 1.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0, 1.0 );
  Vertex& v4 = domain.add_vertex(  0.0,  5.0, 1.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 2 );
  b_ext.add_edge( v3, v4, 3 );
  b_ext.add_edge( v4, v1, 4 );

  // Create the mesh
  Mesh mesh { domain, 0, 0, 50.0 };

  mesh.create_quad_layers(v1, v1, 1, 1.0, 1.5);

  mesh.pave();

  //MSG("============== REFINE ==============");
  mesh.refine_to_quads();
  mesh.refine_to_quads();
  mesh.refine_to_quads();

  // Export the mesh
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();

    std::cout << mesh;
    //domain.export_size_function({0.0,0.0}, {5.0,5.0}, 100, 100);
  }

  DBG_MSG("Tests for Mesh::triangulate_quad_layer() succeeded");

} // Test_Mesh_refine_to_quads() */


/*********************************************************************
* Test Mesh::advance_front()
*********************************************************************/
void Test_Mesh_advance_front_quad(bool export_mesh)
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 0.35;
  };

  Domain domain   { f, 50.0 };

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
  Vertex& v5 = domain.add_vertex(  1.5,  1.5, 0.6, 0.8 );
  Vertex& v6 = domain.add_vertex(  1.5,  3.5 );
  Vertex& v7 = domain.add_vertex(  3.5,  3.5 );

  b_int.add_edge( v5, v6, 2 );
  b_int.add_edge( v6, v7, 2 );
  b_int.add_edge( v7, v5, 2 );

  // Create the mesh
  Mesh mesh { domain, 0, 0, 50.0 };

  mesh.pave();
  //mesh.triangulate();
  mesh.merge_triangles_to_quads();

  // Export the mesh
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();

    std::cout << mesh;
    //domain.export_size_function({0.0,0.0}, {5.0,5.0}, 100, 100);
  }


  DBG_MSG("Tests for Mesh::advance_front_quad() succeeded");

} // Test_Mesh_advance_front_quad() */

/*********************************************************************
* Test Mesh::create_simple_hex_layers()
*********************************************************************/
void Test_Mesh_create_simple_hex_layers(bool export_mesh)
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    //return 1.5 - MIN(0.2*(p.x-3.0), 1.0); 
    //return 0.35; 
    return 0.75; 
  };

  Domain domain   { f, 150.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  2.0,  0.0, 1.0, 1.0 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0, 1.0, 1.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  6.0, 1.0, 1.0 );
  Vertex& v4 = domain.add_vertex(  8.0,  4.0, 1.0, 1.0 );
  Vertex& v5 = domain.add_vertex( 14.0, 10.0, 1.0, 4.0 );
  Vertex& v6 = domain.add_vertex(  0.0,  8.0, 1.0, 1.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 2 );
  b_ext.add_edge( v3, v4, 3 );
  b_ext.add_edge( v4, v5, 4 );
  b_ext.add_edge( v5, v6, 5 );
  b_ext.add_edge( v6, v1, 6 );

  domain.add_fixed_vertex(8.0, 8.0, 0.3, 2.0);

  // Create the mesh
  Mesh mesh { domain, 0, 0, 150.0 };

  // Create quad layers
  mesh.create_quad_layers(v1, v6, 3, 0.1, 1.2);

  mesh.triangulate();
  //mesh.pave();
  //mesh.merge_triangles_to_quads();
  mesh.refine_to_quads();

  // Export the mesh
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();

    std::cout << mesh;
    domain.export_size_function({0.0,0.0}, {14.0,10.0}, 100, 100);
  }

  DBG_MSG("Tests for Mesh::create_simple_hex_layers() succeeded");


} // Test_Mesh_create_simple_hex_layers() */

/*********************************************************************
* Test Mesh::add_quad_layer()
*********************************************************************/
void Test_Mesh_add_quad_layer(bool export_mesh)
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 1.0; 
  };

  Domain domain   { f, 50.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0, 1.0, 0.5 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0, 1.0, 0.5 );
  Vertex& v3 = domain.add_vertex(  5.0,  5.0, 1.0, 0.5 );
  Vertex& v4 = domain.add_vertex(  0.0,  5.0, 1.0, 0.5 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v1, 1 );

  // Create the mesh
  Mesh mesh { domain, 0, 0, 50.0 };

  mesh.create_quad_layers(v1, v1, 4, 0.1, 1.0);

  mesh.triangulate();
  //mesh.pave();
  //mesh.merge_triangles_to_quads();

  // Export the mesh
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();

    std::cout << mesh;
    //domain.export_size_function({0.0,0.0}, {5.0,5.0}, 100, 100);
  }

  DBG_MSG("Tests for Mesh::add_quad_layer() succeeded");


} // Test_Mesh_add_quad_layer() */


/*********************************************************************
* Test Mesh::add_quad_layer_step()
*********************************************************************/
void Test_Mesh_add_quad_layer_step(bool export_mesh)
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 1.5; 
  };

  Domain domain   { f, 50.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0, 1.0 );
  Vertex& v2 = domain.add_vertex(  5.0,  0.0, 1.0 );
  Vertex& v3 = domain.add_vertex(  5.0,  2.5, 1.0 );
  Vertex& v4 = domain.add_vertex(  7.5,  2.5, 1.0 );
  Vertex& v5 = domain.add_vertex(  7.5,  5.0, 1.0 );
  Vertex& v6 = domain.add_vertex(  0.0,  5.0, 1.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v5, 1 );
  b_ext.add_edge( v5, v6, 1 );
  b_ext.add_edge( v6, v1, 1 );

  // Create the mesh
  Mesh mesh { domain, 0, 0, 50.0 };

  mesh.triangulate();
  //mesh.pave();

  // Export the mesh
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();

    std::cout << mesh;
    //domain.export_size_function({0.0,0.0}, {5.0,5.0}, 100, 100);
  }

  DBG_MSG("Tests for Mesh::add_quad_layer() (step) succeeded");


} // Test_Mesh_add_quad_layer_step() */


/*********************************************************************
* Create a simple wedge mesh
*********************************************************************/
void Test_Mesh_wedge(bool export_mesh)
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return (p.x < 300) ? 1.5 : 3.0; 
  };

  Domain domain   { f, 2000.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v0  = domain.add_vertex(   0.0,   0.00,  1.0,  1.0 );
  Vertex& v1  = domain.add_vertex( 300.0,   0.00,  1.0,  1.0 );
  Vertex& v2  = domain.add_vertex( 600.0,   0.00,  1.0,  1.0 );
  Vertex& v3  = domain.add_vertex( 600.0,  25.40,  1.0,  1.0 );
  Vertex& v4  = domain.add_vertex( 300.0,  25.40,  1.0,  1.0 );
  Vertex& v5  = domain.add_vertex(   0.0,  25.40,  1.0,  1.0 );
  Vertex& v6  = domain.add_vertex(-150.0,  25.40,  1.0,  1.0 );
  Vertex& v7  = domain.add_vertex(-150.0,   6.60,  1.0,  1.0 );
  Vertex& v8  = domain.add_vertex(  -1.0,   6.60,  1.0,  1.0 );
  Vertex& v9  = domain.add_vertex(   0.0,   6.35,  1.0,  1.0 );
  Vertex& v10 = domain.add_vertex(-150.0,   6.35,  1.0,  1.0 );
  Vertex& v11 = domain.add_vertex(-150.0,   0.00,  1.0,  1.0 );

  b_ext.add_edge(  v0,  v1, 1 );
  b_ext.add_edge(  v1,  v2, 1 );
  b_ext.add_edge(  v2,  v3, 2 );
  b_ext.add_edge(  v3,  v4, 3 );
  b_ext.add_edge(  v4,  v5, 3 );
  b_ext.add_edge(  v5,  v6, 3 );
  b_ext.add_edge(  v6,  v7, 4 );
  b_ext.add_edge(  v7,  v8, 5 );
  b_ext.add_edge(  v8,  v9, 5 );
  b_ext.add_edge(  v9, v10, 5 );
  b_ext.add_edge( v10, v11, 6 );
  b_ext.add_edge( v11,  v0, 1 );

  // Create the mesh
  Mesh mesh { domain, 0, 0, 2000.0 };

  mesh.triangulate();

  // Export the mesh
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();

    std::cout << mesh;
    domain.export_size_function({0.0,0.0}, {14.0,10.0}, 100, 100);
  }

  DBG_MSG("Tests for Mesh::wedge() succeeded");


} // Test_Mesh_wedge() */


/*********************************************************************
* Create a the TQM banner 
*********************************************************************/
void Test_Mesh_banner(bool export_mesh)
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    const double x0 = 6.0; //9.0;
    const double y0 = 3.0; //3.0;
    const double Rx = 15.0;
    const double Ry = 8.0;
    const double R  = 9.5;
    const double h  = 0.20;
    const double alpha = 0.4 * M_PI;

    double rho = 0.55;
    const double rho_min = 0.125;

    for (int i=0; i<5; i++)
    {
      const double idbl = (double)i;
      const double xi = x0 + Rx * cos(0.5*M_PI+idbl*alpha);
      const double yi = y0 + Ry * sin(0.5*M_PI+idbl*alpha);
      const double dx = p.x - xi;
      const double dy = p.y - yi;
      const double ri = sqrt(dx*dx +dy*dy);
      const double rho_i = h * fabs(ri - R) + rho_min;

      rho = MIN(rho, rho_i);
    }

    return rho;
  };

  Domain domain   { f, 60.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();
  Boundary&  b_T   = domain.add_interior_boundary();
  Boundary&  b_Q   = domain.add_interior_boundary();
  Boundary&  b_M   = domain.add_interior_boundary();

  // Vertices for exterior boundary
  Vertex& b0  = domain.add_vertex( -22.0, -10.00,  0.5,  2.0 );
  Vertex& b1  = domain.add_vertex(  22.0, -10.00,  0.5,  2.0 );
  Vertex& b2  = domain.add_vertex(  22.0,  10.00,  0.5,  2.0 );
  Vertex& b3  = domain.add_vertex( -22.0,  10.00,  0.5,  2.0 );

  // Vertices for letter "T"
  Vertex& T0  = domain.add_vertex( -14.0,  -6.00,  1.0,  1.0 );
  Vertex& T1  = domain.add_vertex( -14.0,   3.00,  0.7,  2.0 );
  Vertex& T2  = domain.add_vertex( -18.0,   3.00,  1.0,  1.0 );
  Vertex& T3  = domain.add_vertex( -18.0,   6.00,  1.0,  1.0 );
  Vertex& T4  = domain.add_vertex(  -7.0,   6.00,  0.7,  2.0 );
  Vertex& T5  = domain.add_vertex(  -7.0,   3.00,  0.7,  2.0 );
  Vertex& T6  = domain.add_vertex( -11.0,   3.00,  1.0,  1.0 );
  Vertex& T7  = domain.add_vertex( -11.0,  -6.00,  1.0,  1.0 );

  // Vertices for letter "Q"
  Vertex& Q0  = domain.add_vertex(  -3.0,  -6.00,  1.0,  1.0 );
  Vertex& Q1  = domain.add_vertex(  -5.0,  -5.00,  1.0,  1.0 );
  Vertex& Q2  = domain.add_vertex(  -5.0,   4.00,  1.0,  1.0 );
  Vertex& Q3  = domain.add_vertex(  -3.0,   6.00,  1.0,  1.0 );
  Vertex& Q4  = domain.add_vertex(   3.0,   6.00,  1.0,  1.0 );
  Vertex& Q5  = domain.add_vertex(   5.0,   4.00,  1.0,  1.0 );
  Vertex& Q6  = domain.add_vertex(   5.0,  -4.00,  1.0,  1.0 );
  Vertex& Q7  = domain.add_vertex(   2.0,   0.00,  1.0,  1.0 );
  Vertex& Q8  = domain.add_vertex(   2.0,   3.00,  1.0,  1.0 );
  Vertex& Q9  = domain.add_vertex(  -2.0,   3.00,  1.0,  1.0 );
  Vertex& Q10 = domain.add_vertex(  -2.0,  -3.00,  1.0,  1.0 );
  Vertex& Q11 = domain.add_vertex(   2.0,  -3.00,  0.7,  2.0 );
  Vertex& Q12 = domain.add_vertex(   4.0,  -6.00,  1.0,  1.0 );

  // Vertices for letter "M"
  Vertex& M0  = domain.add_vertex(   7.0,  -6.00,  1.0,  1.0 );
  Vertex& M1  = domain.add_vertex(   7.0,   6.00,  1.0,  1.0 );
  Vertex& M2  = domain.add_vertex(  10.0,   6.00,  1.0,  1.0 );
  Vertex& M3  = domain.add_vertex(  12.0,   3.00,  1.0,  1.0 );
  Vertex& M4  = domain.add_vertex(  14.0,   3.00,  1.0,  1.0 );
  Vertex& M5  = domain.add_vertex(  16.0,   6.00,  1.0,  1.0 );
  Vertex& M6  = domain.add_vertex(  19.0,   6.00,  1.0,  1.0 );
  Vertex& M7  = domain.add_vertex(  19.0,  -6.00,  1.0,  1.0 );
  Vertex& M8  = domain.add_vertex(  16.0,  -6.00,  1.0,  1.0 );
  Vertex& M9  = domain.add_vertex(  16.0,   1.00,  0.7,  2.0 );
  Vertex& M10 = domain.add_vertex(  14.0,  -1.00,  1.0,  1.0 );
  Vertex& M11 = domain.add_vertex(  12.0,  -1.00,  1.0,  1.0 );
  Vertex& M12 = domain.add_vertex(  10.0,   1.00,  0.7,  2.0 );
  Vertex& M13 = domain.add_vertex(  10.0,  -6.00,  1.0,  1.0 );

  b_ext.add_edge( b0, b1, 1 );
  b_ext.add_edge( b1, b2, 1 );
  b_ext.add_edge( b2, b3, 1 );
  b_ext.add_edge( b3, b0, 1 );

  b_T.add_edge( T0, T1, 2 );
  b_T.add_edge( T1, T2, 2 );
  b_T.add_edge( T2, T3, 2 );
  b_T.add_edge( T3, T4, 2 );
  b_T.add_edge( T4, T5, 2 );
  b_T.add_edge( T5, T6, 2 );
  b_T.add_edge( T6, T7, 2 );
  b_T.add_edge( T7, T0, 2 );

  b_Q.add_edge( Q0,  Q1, 3 );
  b_Q.add_edge( Q1,  Q2, 3 );
  b_Q.add_edge( Q2,  Q3, 3 );
  b_Q.add_edge( Q3,  Q4, 3 );
  b_Q.add_edge( Q4,  Q5, 3 );
  b_Q.add_edge( Q5,  Q6, 3 );
  b_Q.add_edge( Q6,  Q7, 3 );
  b_Q.add_edge( Q7,  Q8, 3 );
  b_Q.add_edge( Q8,  Q9, 3 );
  b_Q.add_edge( Q9,  Q10, 3 );
  b_Q.add_edge( Q10, Q11, 3 );
  b_Q.add_edge( Q11, Q12, 3 );
  b_Q.add_edge( Q12, Q0, 3 );

  b_M.add_edge( M0,   M1, 4 );
  b_M.add_edge( M1,   M2, 4 );
  b_M.add_edge( M2,   M3, 4 );
  b_M.add_edge( M3,   M4, 4 );
  b_M.add_edge( M4,   M5, 4 );
  b_M.add_edge( M5,   M6, 4 );
  b_M.add_edge( M6,   M7, 4 );
  b_M.add_edge( M7,   M8, 4 );
  b_M.add_edge( M8,   M9, 4 );
  b_M.add_edge( M9,   M10, 4 );
  b_M.add_edge( M10,  M11, 4 );
  b_M.add_edge( M11,  M12, 4 );
  b_M.add_edge( M12,  M13, 4 );
  b_M.add_edge( M13,  M0, 4 );

  // Create the mesh
  Mesh mesh { domain, 0, 0, 60.0 };

  // Create quad layers
  mesh.create_quad_layers(b0, b0, 2, 0.5, 1.0);
  mesh.create_quad_layers(Q0, Q0, 2, 0.3, 1.0);
  mesh.create_quad_layers(M0, M0, 2, 0.3, 1.0);
  mesh.create_quad_layers(T0, T0, 2, 0.3, 1.0);

  // Create triangular mesh
  //mesh.pave();
  mesh.triangulate();
  //mesh.merge_triangles_to_quads();

  // Export the mesh
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();

    std::cout << mesh;
    domain.export_size_function({0.0,0.0}, {14.0,10.0}, 100, 100);
  }

  DBG_MSG("Tests for Mesh::banner() succeeded");


} // Test_Mesh_banner() */


/*********************************************************************
* Test: Create mesh for a vortex shedding test case
*********************************************************************/
void Test_Mesh_vortex_shedding(bool export_mesh)
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 0.1;
  };

  Domain domain   { f, 50.0 };

  Boundary&  b_ext = domain.add_exterior_boundary();
  Boundary&  b_int = domain.add_interior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0,  0.0, 1.0, 1.0 );
  Vertex& v2 = domain.add_vertex(  4.0,  0.0, 1.0, 1.0 );
  Vertex& v3 = domain.add_vertex(  4.0,  1.0, 1.0, 1.0 );
  Vertex& v4 = domain.add_vertex(  0.0,  1.0, 1.0, 1.0 );

  Vertex& v5 = domain.add_vertex( 0.35, 0.35, 0.8, 1.3 );
  Vertex& v6 = domain.add_vertex( 0.35, 0.65, 0.8, 1.3 );
  Vertex& v7 = domain.add_vertex( 0.65, 0.65, 0.8, 1.3 );
  Vertex& v8 = domain.add_vertex( 0.65, 0.35, 0.8, 1.3 );

  b_ext.add_edge( v1, v2, 2 );
  b_ext.add_edge( v2, v3, 3 );
  b_ext.add_edge( v3, v4, 2 );
  b_ext.add_edge( v4, v1, 1 );

  // Build interior boundary
  b_int.add_edge( v5, v6, 4 );
  b_int.add_edge( v6, v7, 4 );
  b_int.add_edge( v7, v8, 4 );
  b_int.add_edge( v8, v5, 4 );

  // Create the mesh
  Mesh mesh { domain, 0, 0, 10.0 };

  mesh.create_quad_layers(v1, v2, 2, 0.05, 1.0);
  mesh.create_quad_layers(v3, v4, 2, 0.05, 1.0);
  mesh.create_quad_layers(v5, v5, 2, 0.01, 1.0);

  //mesh.triangulate();
  mesh.pave();
  //mesh.merge_triangles_to_quads();
  //mesh.refine_to_quads();

  // Export the mesh
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();

    std::cout << mesh;
    domain.export_size_function({0.0,0.0}, {5.0,5.0}, 100, 100);
  }


  DBG_MSG("Tests for Mesh::advance_front_quad() succeeded");

} // Test_Mesh_vortex_shedding() */

/*********************************************************************
* Test: Create mesh for a vortex shedding test case
*********************************************************************/
void Test_Mesh_create_bdry_shape_mesh(bool export_mesh)
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 0.5;
  };

  Domain domain   { f, 50.0 };

  Boundary&  b_ext     = domain.add_exterior_boundary();
  Boundary&  b_shape_1 = domain.add_interior_boundary();
  Boundary&  b_shape_2 = domain.add_interior_boundary();
  Boundary&  b_shape_3 = domain.add_interior_boundary();

  b_ext.set_shape_rectangle( domain.vertices(), 1, {6.0,1.0}, 16.0, 8.0 );
  b_shape_1.set_shape_circle( domain.vertices(), 2, {1.0,1.0}, 0.5, 30 );
  b_shape_2.set_shape_square( domain.vertices(), 3, {2.0,1.5}, 0.75);
  b_shape_3.set_shape_triangle( domain.vertices(), 4, {2.0,0.5}, 0.75 );


  // Create the mesh
  Mesh mesh { domain, 0, 0, 40.0 };

  //mesh.create_quad_layers(v1, v2, 2, 0.05, 1.0);
  //mesh.create_quad_layers(v3, v4, 2, 0.05, 1.0);
  //mesh.create_quad_layers(v5, v5, 2, 0.01, 1.0);

  mesh.triangulate();
  //mesh.pave();
  //mesh.merge_triangles_to_quads();

  mesh.refine_to_quads();

  // Export the mesh
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();

    std::cout << mesh;
    domain.export_size_function({0.0,0.0}, {5.0,5.0}, 100, 100);
  }

  DBG_MSG("Tests for Mesh::triangulate() succeeded");

} // Test_Mesh_create_bdry_shape_mesh()


/*********************************************************************
* Test: Create meshing with multiple domains
*********************************************************************/
void Test_Mesh_multiple_domains(bool export_mesh)
{
  // Define a variable size function
  UserSizeFunction f1 = [](const Vec2d& p) { return 0.5; };
  UserSizeFunction f2 = [](const Vec2d& p) { return 0.2; };

  // First domain = major domain 
  Domain domain_1 { f1, 20.0 };

  Boundary&  b_ext_1 = domain_1.add_exterior_boundary();
  Boundary&  b_int_1 = domain_1.add_interior_boundary();

  // Build exterior boundary of domain 1
  Vertex& v1 = domain_1.add_vertex(  0.0,  0.0 );
  Vertex& v2 = domain_1.add_vertex(  8.0,  0.0 );
  Vertex& v3 = domain_1.add_vertex(  8.0,  8.0 );
  Vertex& v4 = domain_1.add_vertex(  0.0,  8.0 );

  b_ext_1.add_edge( v1, v2, 1 );
  b_ext_1.add_edge( v2, v3, 1 );
  b_ext_1.add_edge( v3, v4, 1 );
  b_ext_1.add_edge( v4, v1, 1 );

  // Build interior boundary of domain 1
  // --> This will be the exterior boundary of domain 2
  Vertex& v5 = domain_1.add_vertex(  2.0,  2.0 );
  Vertex& v6 = domain_1.add_vertex(  2.0,  6.0 );
  Vertex& v7 = domain_1.add_vertex(  6.0,  6.0 );
  Vertex& v8 = domain_1.add_vertex(  6.0,  2.0 );

  b_int_1.add_edge( v5, v6, 2 ); // These markers will not be needed 
  b_int_1.add_edge( v6, v7, 2 ); // after everything is meshed
  b_int_1.add_edge( v7, v8, 2 );
  b_int_1.add_edge( v8, v5, 2 );

  // Create the mesh for domain 1
  // --> Here the edges of domain 1 are discretized into further 
  //     segments, according to the size function of domain 1
  Mesh mesh_1 { domain_1, 0, 0, 50.0 };

  mesh_1.triangulate();
  mesh_1.merge_triangles_to_quads();

  mesh_1.refine_to_quads();

  // ToDo:
  // --> Every domain edge must know its sub-vertices of the mesh
  //     as well as the corresponding boundary edges of the mesh
  // --> This must happen during the advancing front initialization
  //
  // --> Structure for every domain edge:
  //     * vector of pointers to sub-edges
  //       (both arranged in direction of the edge)
  //     * Pointer to possible twin edge of another domain
  //  



  //=======================================================
  // Second domain = sub-domain of domain_1
  // --> domain_1 is given as constructor argument instead 
  //     of size function
  // --> it's size function is used
  // --> the remaining arguments can be given by user
  //
  // ToDo: 
  // --> Everytime an entity gets added, we must check
  //     in the background, if this entity is a copy
  //     of an entity in domain_1
  //     * Vertices: by coordinate
  //     * Edges: by start & ending vertex coordinate
  //     -> Use get_item() in case of large domains?
  //
  //
  Domain domain_2 { f2, 20.0 };
  Boundary&  b_ext_2 = domain_2.add_exterior_boundary();

  // Build exterior boundary of domain 2
  // --> Must somehow handle, that vertex scale properties 
  //     of Domain 1 and Domain 2 should be the same!
  Vertex& v9  = domain_2.add_vertex(  2.0,  2.0 );
  Vertex& v10 = domain_2.add_vertex(  6.0,  2.0 );
  Vertex& v11 = domain_2.add_vertex(  6.0,  6.0 );
  Vertex& v12 = domain_2.add_vertex(  2.0,  6.0 );

  b_ext_2.add_edge(  v9, v10, 2 );
  b_ext_2.add_edge( v10, v11, 2 );
  b_ext_2.add_edge( v11, v12, 2 );
  b_ext_2.add_edge( v12,  v9, 2 );
  
  // Create the mesh for domain 2
  // --> Here the edges of domain 2 are discretized into further 
  //     segments, by using the corresponding entities of domain 1
  //     or by placing new vertices according to its size function
  Mesh mesh_2 { domain_2, mesh_1, 1, 1, 50.0 };

  mesh_2.triangulate();
  mesh_2.refine_to_quads();



  // Merge meshes
  // --> combine both meshes to one 
  // --> Boundaries must keep their markers
  // --> Boundary edges at mesh interfaces must be turned into 
  //     interior edges
  // --> In case of entitiy duplicates, use always vertices or 
  //     edges of the first mesh
  // --> If vertices,  interior & boundary edges and elements
  //     are initialized, we can get the full connectivity by 
  //     the mesh's functions
  //Mesh mesh_3 { mesh_1, mesh_2 };


  // Export the mesh
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    //mesh_1.assign_mesh_indices();
    //mesh_2.assign_mesh_indices();

    //std::cout << mesh_1;
    //std::cout << mesh_2;
    //domain.export_size_function({0.0,0.0}, {5.0,5.0}, 100, 100);

    std::string export_prefix = 
      "/home/florian/datadisk/Code/C++-Code/TQMesh/build/mesh";

    mesh_1.write_to_file( export_prefix + "_1.txt", ExportType::txt );
    mesh_2.write_to_file( export_prefix + "_2.txt", ExportType::txt );
  }

  DBG_MSG("Tests for Mesh_multiple_domains() succeeded");

} // Test_Mesh_multiple_domains()


/*********************************************************************
* Test: Benchmark the meshing process
*********************************************************************/
void Test_Mesh_benchmark(double h, double L, 
                         BenchmarkContainer& bm_data,
                         bool export_mesh=false)
{
  // Define a variable size function
  UserSizeFunction f = [h](const Vec2d& p) 
  { 
    return h;
  }; 

  Timer            timer {};
  Domain           domain   { f, 50.0*L };

  Boundary&  b_ext = domain.add_exterior_boundary();
  Boundary&  b_int = domain.add_interior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  2.0*L,  0.0*L, 1.0, 1.0 );
  Vertex& v2 = domain.add_vertex(  5.0*L,  0.0*L, 1.0, 1.0 );
  Vertex& v3 = domain.add_vertex(  5.0*L,  6.0*L, 1.0, 1.0 );
  Vertex& v4 = domain.add_vertex(  8.0*L,  4.0*L, 1.0, 1.0 );
  Vertex& v5 = domain.add_vertex( 14.0*L, 10.0*L, 1.0, 4.0 );
  Vertex& v6 = domain.add_vertex(  0.0*L,  8.0*L, 1.0, 1.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 2 );
  b_ext.add_edge( v3, v4, 3 );
  b_ext.add_edge( v4, v5, 4 );
  b_ext.add_edge( v5, v6, 5 );
  b_ext.add_edge( v6, v1, 6 );

  // Build interior boundary
  Vertex& v7 = domain.add_vertex(  2.0*L,  6.0*L, 1.0, 1.0 );
  Vertex& v8 = domain.add_vertex(  3.0*L,  6.0*L, 1.0, 1.0 );
  Vertex& v9 = domain.add_vertex(  3.0*L,  3.0*L, 1.0, 1.0 );

  b_int.add_edge( v7, v8, 1 );
  b_int.add_edge( v8, v9, 1 );
  b_int.add_edge( v9, v7, 1 );

  // Add a fixed vertex
  domain.add_fixed_vertex(8.0*L, 8.0*L, 0.3, 2.0);

  /*------------------------------------------------------------------
  | Create the mesh & discretize the front
  ------------------------------------------------------------------*/
  timer.count();

  Mesh mesh { domain, 0, 0, 40.0*L };

  /*------------------------------------------------------------------
  | Create the quad layer
  ------------------------------------------------------------------*/
  timer.count();

  mesh.create_quad_layers(v1, v6, 3, 0.2*L*h, 1.2);

  /*------------------------------------------------------------------
  | Triangulate the mesh
  ------------------------------------------------------------------*/
  timer.count();

  mesh.triangulate();

  /*------------------------------------------------------------------
  | Smoothing
  ------------------------------------------------------------------*/
  timer.count();

  /*------------------------------------------------------------------
  | Finalize benchmark - output times to user
  ------------------------------------------------------------------*/
  timer.count();

  size_t n_vertices   = mesh.vertices().size();
  size_t n_tris       = mesh.triangles().size();
  size_t n_quads      = mesh.quads().size();
  size_t n_intr_edges = mesh.interior_edges().size();
  size_t n_bdry_edges = mesh.boundary_edges().size();

  MSG( "\n--------------------------------------------" );
  MSG( "h=" << h << ", L=" << L );
  MSG( "Number of vertices:             " << n_vertices );
  MSG( "Number of triangles:            " << n_tris );
  MSG( "Number of quads:                " << n_quads );
  MSG( "Number of interior edges:       " << n_intr_edges );
  MSG( "Number of boundary edges:       " << n_bdry_edges );


  double t_tot = timer.delta(0) + timer.delta(1) 
               + timer.delta(2) + timer.delta(3);

  MSG( "Time for mesh initialization:   " 
      << std::setprecision(5) << std::fixed << timer.delta(0) << "s" );
  MSG( "Time for quad layer generation: " 
      << std::setprecision(5) << std::fixed << timer.delta(1) << "s" );
  MSG( "Time for triangulation:         " 
      << std::setprecision(5) << std::fixed << timer.delta(2) << "s" );
  MSG( "Time for mesh smoothing:        " 
      << std::setprecision(5) << std::fixed << timer.delta(3) << "s" );
  MSG( "Total meshing time:             " 
      << std::setprecision(5) << std::fixed << t_tot << "s" );

  //------------------------------------------------------------------
  // Update benchmark container
  //------------------------------------------------------------------
  bm_data.h.push_back( h );
  bm_data.L.push_back( L );

  bm_data.n_vertices.push_back( n_vertices );
  bm_data.n_tris.push_back( n_tris );
  bm_data.n_quads.push_back( n_quads );
  bm_data.n_intr_edges.push_back( n_intr_edges );
  bm_data.n_bdry_edges.push_back( n_bdry_edges );

  bm_data.quad_layer_time.push_back( timer.delta(1) );
  bm_data.meshing_time.push_back( timer.delta(2) );
  bm_data.smoothing_time.push_back( timer.delta(3) );
               
  /*------------------------------------------------------------------
  | Export the mesh
  ------------------------------------------------------------------*/
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();
    std::cout << mesh;
    domain.export_size_function({0.0,0.0}, {14.0,10.0}, 100, 100);
  }


} // Test_Mesh_benchmark()


/*********************************************************************
* Test: Benchmark the meshing process
*********************************************************************/
void Test_Mesh_benchmark_TMesh(double h, double L, 
                               BenchmarkContainer& bm_data,
                               bool export_mesh=false)
{
  // Define a variable size function
  UserSizeFunction f = [h](const Vec2d& p) 
  { 
    return h;
  };

  Timer            timer {};
  Domain           domain   { f, 20.0*L };

  Boundary&  b_ext = domain.add_exterior_boundary();

  // Build exterior boundary
  Vertex& v1 = domain.add_vertex(  0.0*L,  0.0*L, 1.0, 1.0 );
  Vertex& v2 = domain.add_vertex(  6.0*L,  0.0*L, 1.0, 1.0 );
  Vertex& v3 = domain.add_vertex(  8.0*L,  3.0*L, 1.0, 1.0 );
  Vertex& v4 = domain.add_vertex(  8.0*L,  6.0*L, 1.0, 1.0 );
  Vertex& v5 = domain.add_vertex(  8.0*L,  9.0*L, 1.0, 1.0 );
  Vertex& v6 = domain.add_vertex(  3.0*L,  9.0*L, 1.0, 1.0 );
  Vertex& v7 = domain.add_vertex( -1.0*L,  9.0*L, 1.0, 1.0 );
  Vertex& v8 = domain.add_vertex(  1.0*L,  3.0*L, 1.0, 1.0 );
  Vertex& v9 = domain.add_vertex( -3.0*L,  6.0*L, 1.0, 1.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v5, 1 );
  b_ext.add_edge( v5, v6, 2 );
  b_ext.add_edge( v6, v7, 2 );
  b_ext.add_edge( v7, v8, 2 );
  b_ext.add_edge( v8, v9, 2 );
  b_ext.add_edge( v9, v1, 2 );

  //------------------------------------------------------------------
  // Create the mesh & discretize the front
  //------------------------------------------------------------------
  timer.count();

  Mesh mesh { domain, 0, 0, 20.0*L };

  //------------------------------------------------------------------
  // Create the quad layer
  //------------------------------------------------------------------
  timer.count();

  //------------------------------------------------------------------
  // Triangulate the mesh
  //------------------------------------------------------------------
  timer.count();

  mesh.triangulate();
  
  //------------------------------------------------------------------
  // Smoothing
  //------------------------------------------------------------------
  timer.count();

  //------------------------------------------------------------------
  // Finalize benchmark - output times to user
  //------------------------------------------------------------------
  timer.count();

  size_t n_vertices   = mesh.vertices().size();
  size_t n_tris       = mesh.triangles().size();
  size_t n_quads      = mesh.quads().size();
  size_t n_intr_edges = mesh.interior_edges().size();
  size_t n_bdry_edges = mesh.boundary_edges().size();

  MSG( "");
  MSG( "--------------------------------------------" );
  MSG( "h=" << h << ", L=" << L );
  MSG( "Number of vertices:             " << n_vertices );
  MSG( "Number of triangles:            " << n_tris );
  MSG( "Number of quads:                " << n_quads );
  MSG( "Number of interior edges:       " << n_intr_edges );
  MSG( "Number of boundary edges:       " << n_bdry_edges );


  double t_tot = timer.delta(0) + timer.delta(1) 
               + timer.delta(2) + timer.delta(3);

  MSG( "Time for mesh initialization:   " 
      << std::setprecision(5) << std::fixed << timer.delta(0) << "s" );
  MSG( "Time for quad layer generation: " 
      << std::setprecision(5) << std::fixed << timer.delta(1) << "s" );
  MSG( "Time for triangulation:         " 
      << std::setprecision(5) << std::fixed << timer.delta(2) << "s" );
  MSG( "Time for mesh smoothing:        " 
      << std::setprecision(5) << std::fixed << timer.delta(3) << "s" );
  MSG( "Total meshing time:             " 
      << std::setprecision(5) << std::fixed << t_tot << "s" );

  //------------------------------------------------------------------
  // Update benchmark container
  //------------------------------------------------------------------
  bm_data.h.push_back( h );
  bm_data.L.push_back( L );

  bm_data.n_vertices.push_back( n_vertices );
  bm_data.n_tris.push_back( n_tris );
  bm_data.n_quads.push_back( n_quads );
  bm_data.n_intr_edges.push_back( n_intr_edges );
  bm_data.n_bdry_edges.push_back( n_bdry_edges );

  bm_data.quad_layer_time.push_back( timer.delta(1) );
  bm_data.meshing_time.push_back( timer.delta(2) );
  bm_data.smoothing_time.push_back( timer.delta(3) );
               
  //------------------------------------------------------------------
  // Export the mesh
  //------------------------------------------------------------------
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();
    std::cout << mesh;
    domain.export_size_function({0.0,0.0}, {14.0,10.0}, 100, 100);
  }

} // Test_Mesh_benchmark_TMesh()  */

} // namespace MeshTests


/*********************************************************************
* 
*********************************************************************/
void run_mesh_tests(bool benchmark)
{
  MSG("\n#===== Mesh tests =====");

  MeshTests::Test_Mesh_initialization();
  //MeshTests::Test_Mesh_triangulate(true);
  //MeshTests::Test_Mesh_pave(false);
  //MeshTests::Test_Mesh_refine_to_quads(true);
  //MeshTests::Test_Mesh_advance_front_quad(true);
  //MeshTests::Test_Mesh_add_quad_layer_step(true);
  //MeshTests::Test_Mesh_wedge(true);
  //MeshTests::Test_Mesh_banner(true);
    
  //MeshTests::Test_Mesh_create_simple_hex_layers(true);
  //MeshTests::Test_Mesh_vortex_shedding(true);
    
  //MeshTests::Test_Mesh_create_bdry_shape_mesh(true);

  //MeshTests::Test_Mesh_add_quad_layer(true);
    
  MeshTests::Test_Mesh_multiple_domains(true);


  if ( benchmark )
  {
    MSG("\n#===== TQMesh benchmarks =====");

    MeshTests::BenchmarkContainer bm_data {};

    double h = 0.5;
    double L = 1.0;

    for (size_t i = 0; i < 6; i++)
    {
      MeshTests::Test_Mesh_benchmark(h, L, bm_data, false);
      //MeshTests::Test_Mesh_benchmark_TMesh(h, L, bm_data, false);
      L *= 2.0;
    }

    // Output data 
    for ( size_t i = 0; i < bm_data.h.size(); ++i )
    {
      std::cout << std::setprecision(5) << std::fixed 
        << bm_data.h[i] << ", " 
        << bm_data.L[i] << ", " 
        << bm_data.n_vertices[i] << ", " 
        << bm_data.n_tris[i] << ", " 
        << bm_data.n_quads[i] << ", " 
        << bm_data.n_intr_edges[i] << ", " 
        << bm_data.n_bdry_edges[i] << ", " 
        << bm_data.quad_layer_time[i] << ", " 
        << bm_data.meshing_time[i] << ", " 
        << bm_data.smoothing_time[i] << "\n";
    }

  }


} // run_mesh_tests()

