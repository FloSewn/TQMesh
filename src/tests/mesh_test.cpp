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

using namespace TQMesh::TQUtils;
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

  // Create the mesh
  Mesh mesh { domain, 10.0 };

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

  Vertices         vertices { 20.0 };
  Domain           domain   { vertices, f };

  Boundary&  b_ext = domain.add_boundary( BdryType::EXTERIOR );

  // Build exterior boundary
  Vertex& v1 = vertices.push_back(  0.0,  0.0, 0.9 );
  Vertex& v2 = vertices.push_back(  5.0,  0.0, 0.9 );
  Vertex& v3 = vertices.push_back(  5.0,  5.0, 0.9 );
  Vertex& v4 = vertices.push_back( 10.0,  5.0, 0.9 );
  Vertex& v5 = vertices.push_back( 10.0, 10.0, 0.9 );
  Vertex& v6 = vertices.push_back(  0.0,  5.0, 0.9 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v5, 1 );
  b_ext.add_edge( v5, v6, 1 );
  b_ext.add_edge( v6, v1, 1 );

  // Create the mesh
  Mesh mesh { domain, 50.0 };

  mesh.triangulate();
  mesh.smoothing(4, 0.9);

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

  Vertices         vertices { 20.0 };
  Domain           domain   { vertices, f };

  Boundary&  b_ext = domain.add_boundary( BdryType::EXTERIOR );

  // Build exterior boundary
  Vertex& v1 = vertices.push_back(  0.0,  0.0, 0.9 );
  Vertex& v2 = vertices.push_back(  5.0,  0.0, 0.9 );
  Vertex& v3 = vertices.push_back(  5.0,  5.0, 0.9 );
  Vertex& v4 = vertices.push_back( 10.0,  5.0, 0.9 );
  Vertex& v5 = vertices.push_back( 10.0, 10.0, 0.9 );
  Vertex& v6 = vertices.push_back(  0.0,  5.0, 0.9 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v5, 1 );
  b_ext.add_edge( v5, v6, 1 );
  b_ext.add_edge( v6, v1, 1 );

  // Create the mesh
  Mesh mesh { domain, 50.0 };

  mesh.pave();
  mesh.smoothing(4, 0.9);

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
void Test_Mesh_triangulate_quad_layer(bool export_mesh)
{
  // Define a variable size function
  UserSizeFunction f = [](const Vec2d& p) 
  { 
    return 2.5; 
  };

  Vertices         vertices { 20.0 };
  Domain           domain   { vertices, f };

  Boundary&  b_ext = domain.add_boundary( BdryType::EXTERIOR );

  // Build exterior boundary
  Vertex& v1 = vertices.push_back(  0.0,  0.0, 0.9 );
  Vertex& v2 = vertices.push_back(  5.0,  0.0, 0.9 );
  Vertex& v3 = vertices.push_back(  5.0,  5.0, 0.9 );
  Vertex& v4 = vertices.push_back( 10.0,  5.0, 0.9 );
  Vertex& v5 = vertices.push_back( 10.0, 10.0, 0.9 );
  Vertex& v6 = vertices.push_back(  0.0,  9.0, 0.9 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 2 );
  b_ext.add_edge( v3, v4, 3 );
  b_ext.add_edge( v4, v5, 4 );
  b_ext.add_edge( v5, v6, 5 );
  b_ext.add_edge( v6, v1, 6 );

  // Create the mesh
  Mesh mesh { domain, 50.0 };

  mesh.create_quad_layers(v1, v2, 1, 0.25, 1.5);

  mesh.triangulate();
  mesh.smoothing(4, 0.9);

  // Export the mesh
  if (export_mesh)
  {
    // Make sure that all vertex / element indices are assigned
    mesh.assign_mesh_indices();

    std::cout << mesh;
    //domain.export_size_function({0.0,0.0}, {5.0,5.0}, 100, 100);
  }

  DBG_MSG("Tests for Mesh::triangulate_quad_layer() succeeded");

} // Test_Mesh_triangulate_quad_layer() */


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

  Vertices         vertices { 50.0 };
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
  Vertex& v5 = vertices.push_back(  1.5,  1.5, 0.6, 0.8 );
  Vertex& v6 = vertices.push_back(  1.5,  3.5 );
  Vertex& v7 = vertices.push_back(  3.5,  3.5 );

  b_int.add_edge( v5, v6, 2 );
  b_int.add_edge( v6, v7, 2 );
  b_int.add_edge( v7, v5, 2 );

  // Create the mesh
  Mesh mesh { domain, 50.0 };

  mesh.pave();
  //mesh.triangulate();
  mesh.smoothing(4, 0.9);
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
    return 0.35; 
    //return 0.5; 
  };

  Vertices         vertices { 150.0 };
  Domain           domain   { vertices, f };

  Boundary&  b_ext = domain.add_boundary( BdryType::EXTERIOR );

  // Build exterior boundary
  Vertex& v1 = vertices.push_back(  2.0,  0.0, 1.0, 1.0 );
  Vertex& v2 = vertices.push_back(  5.0,  0.0, 1.0, 1.0 );
  Vertex& v3 = vertices.push_back(  5.0,  6.0, 1.0, 1.0 );
  Vertex& v4 = vertices.push_back(  8.0,  4.0, 1.0, 1.0 );
  Vertex& v5 = vertices.push_back( 14.0, 10.0, 1.0, 4.0 );
  Vertex& v6 = vertices.push_back(  0.0,  8.0, 1.0, 1.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 2 );
  b_ext.add_edge( v3, v4, 3 );
  b_ext.add_edge( v4, v5, 4 );
  b_ext.add_edge( v5, v6, 5 );
  b_ext.add_edge( v6, v1, 6 );

  domain.add_fixed_vertex(8.0, 8.0, 0.3, 2.0);

  // Create the mesh
  Mesh mesh { domain, 150.0 };

  // Create quad layers
  mesh.create_quad_layers(v1, v6, 3, 0.1, 1.2);

  mesh.triangulate();
  //mesh.pave();
  //mesh.merge_triangles_to_quads();
  mesh.smoothing(4, 0.9);

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

  Vertices         vertices { 50.0 };
  Domain           domain   { vertices, f };

  Boundary&  b_ext = domain.add_boundary( BdryType::EXTERIOR );

  // Build exterior boundary
  Vertex& v1 = vertices.push_back(  0.0,  0.0, 1.0, 0.5 );
  Vertex& v2 = vertices.push_back(  5.0,  0.0, 1.0, 0.5 );
  Vertex& v3 = vertices.push_back(  5.0,  5.0, 1.0, 0.5 );
  Vertex& v4 = vertices.push_back(  0.0,  5.0, 1.0, 0.5 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v1, 1 );

  // Create the mesh
  Mesh mesh { domain, 50.0 };

  mesh.create_quad_layers(v1, v1, 4, 0.1, 1.0);

  //mesh.triangulate();
  //mesh.pave();
  //mesh.merge_triangles_to_quads();
  mesh.smoothing(4, 0.9);

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

  Vertices         vertices { 50.0 };
  Domain           domain   { vertices, f };

  Boundary&  b_ext = domain.add_boundary( BdryType::EXTERIOR );

  // Build exterior boundary
  Vertex& v1 = vertices.push_back(  0.0,  0.0, 1.0 );
  Vertex& v2 = vertices.push_back(  5.0,  0.0, 1.0 );
  Vertex& v3 = vertices.push_back(  5.0,  2.5, 1.0 );
  Vertex& v4 = vertices.push_back(  7.5,  2.5, 1.0 );
  Vertex& v5 = vertices.push_back(  7.5,  5.0, 1.0 );
  Vertex& v6 = vertices.push_back(  0.0,  5.0, 1.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 1 );
  b_ext.add_edge( v3, v4, 1 );
  b_ext.add_edge( v4, v5, 1 );
  b_ext.add_edge( v5, v6, 1 );
  b_ext.add_edge( v6, v1, 1 );

  // Create the mesh
  Mesh mesh { domain, 50.0 };

  //mesh.triangulate();
  //mesh.pave();
  mesh.smoothing(4, 0.9);

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

  Vertices         vertices { 2000.0 };
  Domain           domain   { vertices, f };

  Boundary&  b_ext = domain.add_boundary( BdryType::EXTERIOR );

  // Build exterior boundary
  Vertex& v0  = vertices.push_back(   0.0,   0.00,  1.0,  1.0 );
  Vertex& v1  = vertices.push_back( 300.0,   0.00,  1.0,  1.0 );
  Vertex& v2  = vertices.push_back( 600.0,   0.00,  1.0,  1.0 );
  Vertex& v3  = vertices.push_back( 600.0,  25.40,  1.0,  1.0 );
  Vertex& v4  = vertices.push_back( 300.0,  25.40,  1.0,  1.0 );
  Vertex& v5  = vertices.push_back(   0.0,  25.40,  1.0,  1.0 );
  Vertex& v6  = vertices.push_back(-150.0,  25.40,  1.0,  1.0 );
  Vertex& v7  = vertices.push_back(-150.0,   6.60,  1.0,  1.0 );
  Vertex& v8  = vertices.push_back(  -1.0,   6.60,  1.0,  1.0 );
  Vertex& v9  = vertices.push_back(   0.0,   6.35,  1.0,  1.0 );
  Vertex& v10 = vertices.push_back(-150.0,   6.35,  1.0,  1.0 );
  Vertex& v11 = vertices.push_back(-150.0,   0.00,  1.0,  1.0 );

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
  Mesh mesh { domain, 2000.0 };

  mesh.triangulate();
  mesh.smoothing(8, 0.9);

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

  Vertices         vertices { 60.0 };
  Domain           domain   { vertices, f };

  Boundary&  b_ext = domain.add_boundary( BdryType::EXTERIOR );
  Boundary&  b_T   = domain.add_boundary( BdryType::INTERIOR );
  Boundary&  b_Q   = domain.add_boundary( BdryType::INTERIOR );
  Boundary&  b_M   = domain.add_boundary( BdryType::INTERIOR );

  // Vertices for exterior boundary
  Vertex& b0  = vertices.push_back( -22.0, -10.00,  0.5,  2.0 );
  Vertex& b1  = vertices.push_back(  22.0, -10.00,  0.5,  2.0 );
  Vertex& b2  = vertices.push_back(  22.0,  10.00,  0.5,  2.0 );
  Vertex& b3  = vertices.push_back( -22.0,  10.00,  0.5,  2.0 );

  // Vertices for letter "T"
  Vertex& T0  = vertices.push_back( -14.0,  -6.00,  1.0,  1.0 );
  Vertex& T1  = vertices.push_back( -14.0,   3.00,  0.7,  2.0 );
  Vertex& T2  = vertices.push_back( -18.0,   3.00,  1.0,  1.0 );
  Vertex& T3  = vertices.push_back( -18.0,   6.00,  1.0,  1.0 );
  Vertex& T4  = vertices.push_back(  -7.0,   6.00,  0.7,  2.0 );
  Vertex& T5  = vertices.push_back(  -7.0,   3.00,  0.7,  2.0 );
  Vertex& T6  = vertices.push_back( -11.0,   3.00,  1.0,  1.0 );
  Vertex& T7  = vertices.push_back( -11.0,  -6.00,  1.0,  1.0 );

  // Vertices for letter "Q"
  Vertex& Q0  = vertices.push_back(  -3.0,  -6.00,  1.0,  1.0 );
  Vertex& Q1  = vertices.push_back(  -5.0,  -5.00,  1.0,  1.0 );
  Vertex& Q2  = vertices.push_back(  -5.0,   4.00,  1.0,  1.0 );
  Vertex& Q3  = vertices.push_back(  -3.0,   6.00,  1.0,  1.0 );
  Vertex& Q4  = vertices.push_back(   3.0,   6.00,  1.0,  1.0 );
  Vertex& Q5  = vertices.push_back(   5.0,   4.00,  1.0,  1.0 );
  Vertex& Q6  = vertices.push_back(   5.0,  -4.00,  1.0,  1.0 );
  Vertex& Q7  = vertices.push_back(   2.0,   0.00,  1.0,  1.0 );
  Vertex& Q8  = vertices.push_back(   2.0,   3.00,  1.0,  1.0 );
  Vertex& Q9  = vertices.push_back(  -2.0,   3.00,  1.0,  1.0 );
  Vertex& Q10 = vertices.push_back(  -2.0,  -3.00,  1.0,  1.0 );
  Vertex& Q11 = vertices.push_back(   2.0,  -3.00,  0.7,  2.0 );
  Vertex& Q12 = vertices.push_back(   4.0,  -6.00,  1.0,  1.0 );

  // Vertices for letter "M"
  Vertex& M0  = vertices.push_back(   7.0,  -6.00,  1.0,  1.0 );
  Vertex& M1  = vertices.push_back(   7.0,   6.00,  1.0,  1.0 );
  Vertex& M2  = vertices.push_back(  10.0,   6.00,  1.0,  1.0 );
  Vertex& M3  = vertices.push_back(  12.0,   3.00,  1.0,  1.0 );
  Vertex& M4  = vertices.push_back(  14.0,   3.00,  1.0,  1.0 );
  Vertex& M5  = vertices.push_back(  16.0,   6.00,  1.0,  1.0 );
  Vertex& M6  = vertices.push_back(  19.0,   6.00,  1.0,  1.0 );
  Vertex& M7  = vertices.push_back(  19.0,  -6.00,  1.0,  1.0 );
  Vertex& M8  = vertices.push_back(  16.0,  -6.00,  1.0,  1.0 );
  Vertex& M9  = vertices.push_back(  16.0,   1.00,  0.7,  2.0 );
  Vertex& M10 = vertices.push_back(  14.0,  -1.00,  1.0,  1.0 );
  Vertex& M11 = vertices.push_back(  12.0,  -1.00,  1.0,  1.0 );
  Vertex& M12 = vertices.push_back(  10.0,   1.00,  0.7,  2.0 );
  Vertex& M13 = vertices.push_back(  10.0,  -6.00,  1.0,  1.0 );

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
  Mesh mesh { domain, 60.0 };

  // Create quad layers
  mesh.create_quad_layers(b0, b0, 2, 0.5, 1.0);
  mesh.create_quad_layers(Q0, Q0, 2, 0.3, 1.0);
  mesh.create_quad_layers(M0, M0, 2, 0.3, 1.0);
  mesh.create_quad_layers(T0, T0, 2, 0.3, 1.0);

  // Create triangular mesh
  //mesh.pave();
  mesh.triangulate();
  //mesh.merge_triangles_to_quads();
  mesh.smoothing(8, 0.9);

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

  Vertices         vertices { 50.0 };
  Domain           domain   { vertices, f };

  Boundary&  b_ext = domain.add_boundary( BdryType::EXTERIOR );
  Boundary&  b_int = domain.add_boundary( BdryType::INTERIOR );

  // Build exterior boundary
  Vertex& v1 = vertices.push_back(  0.0,  0.0, 1.0, 1.0 );
  Vertex& v2 = vertices.push_back(  4.0,  0.0, 1.0, 1.0 );
  Vertex& v3 = vertices.push_back(  4.0,  1.0, 1.0, 1.0 );
  Vertex& v4 = vertices.push_back(  0.0,  1.0, 1.0, 1.0 );

  Vertex& v5 = vertices.push_back( 0.35, 0.35, 0.8, 1.3 );
  Vertex& v6 = vertices.push_back( 0.35, 0.65, 0.8, 1.3 );
  Vertex& v7 = vertices.push_back( 0.65, 0.65, 0.8, 1.3 );
  Vertex& v8 = vertices.push_back( 0.65, 0.35, 0.8, 1.3 );

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
  Mesh mesh { domain, 10.0 };

  mesh.create_quad_layers(v1, v2, 2, 0.05, 1.0);
  mesh.create_quad_layers(v3, v4, 2, 0.05, 1.0);
  mesh.create_quad_layers(v5, v5, 2, 0.01, 1.0);

  mesh.triangulate();
  //mesh.pave();
  //mesh.merge_triangles_to_quads();
  mesh.smoothing(4, 0.9);

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
  Vertices         vertices { 50.0*L };
  Domain           domain   { vertices, f };

  Boundary&  b_ext = domain.add_boundary( BdryType::EXTERIOR );
  Boundary&  b_int = domain.add_boundary( BdryType::INTERIOR );

  // Build exterior boundary
  Vertex& v1 = vertices.push_back(  2.0*L,  0.0*L, 1.0, 1.0 );
  Vertex& v2 = vertices.push_back(  5.0*L,  0.0*L, 1.0, 1.0 );
  Vertex& v3 = vertices.push_back(  5.0*L,  6.0*L, 1.0, 1.0 );
  Vertex& v4 = vertices.push_back(  8.0*L,  4.0*L, 1.0, 1.0 );
  Vertex& v5 = vertices.push_back( 14.0*L, 10.0*L, 1.0, 4.0 );
  Vertex& v6 = vertices.push_back(  0.0*L,  8.0*L, 1.0, 1.0 );

  b_ext.add_edge( v1, v2, 1 );
  b_ext.add_edge( v2, v3, 2 );
  b_ext.add_edge( v3, v4, 3 );
  b_ext.add_edge( v4, v5, 4 );
  b_ext.add_edge( v5, v6, 5 );
  b_ext.add_edge( v6, v1, 6 );

  // Build interior boundary
  Vertex& v7 = vertices.push_back(  2.0*L,  6.0*L, 1.0, 1.0 );
  Vertex& v8 = vertices.push_back(  3.0*L,  6.0*L, 1.0, 1.0 );
  Vertex& v9 = vertices.push_back(  3.0*L,  3.0*L, 1.0, 1.0 );

  b_int.add_edge( v7, v8, 1 );
  b_int.add_edge( v8, v9, 1 );
  b_int.add_edge( v9, v7, 1 );

  // Add a fixed vertex
  domain.add_fixed_vertex(8.0*L, 8.0*L, 0.3, 2.0);

  /*------------------------------------------------------------------
  | Create the mesh & discretize the front
  ------------------------------------------------------------------*/
  timer.count();

  Mesh mesh { domain, 40.0*L };

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

  mesh.smoothing(4, 0.9);

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
  Vertices         vertices { 20.0*L };
  Domain           domain   { vertices, f };

  Boundary&  b_ext = domain.add_boundary( BdryType::EXTERIOR );

  // Build exterior boundary
  Vertex& v1 = vertices.push_back(  0.0*L,  0.0*L, 1.0, 1.0 );
  Vertex& v2 = vertices.push_back(  6.0*L,  0.0*L, 1.0, 1.0 );
  Vertex& v3 = vertices.push_back(  8.0*L,  3.0*L, 1.0, 1.0 );
  Vertex& v4 = vertices.push_back(  8.0*L,  6.0*L, 1.0, 1.0 );
  Vertex& v5 = vertices.push_back(  8.0*L,  9.0*L, 1.0, 1.0 );
  Vertex& v6 = vertices.push_back(  3.0*L,  9.0*L, 1.0, 1.0 );
  Vertex& v7 = vertices.push_back( -1.0*L,  9.0*L, 1.0, 1.0 );
  Vertex& v8 = vertices.push_back(  1.0*L,  3.0*L, 1.0, 1.0 );
  Vertex& v9 = vertices.push_back( -3.0*L,  6.0*L, 1.0, 1.0 );

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

  Mesh mesh { domain, 20.0*L };

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

  mesh.smoothing(4, 0.9);

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

  //MeshTests::Test_Mesh_initialization();
  //MeshTests::Test_Mesh_triangulate(false);
  //MeshTests::Test_Mesh_pave(false);
  //MeshTests::Test_Mesh_triangulate_quad_layer(false);
  //MeshTests::Test_Mesh_advance_front_quad(true);
  //MeshTests::Test_Mesh_add_quad_layer_step(true);
  //MeshTests::Test_Mesh_wedge(true);
  MeshTests::Test_Mesh_banner(true);
  //
  //MeshTests::Test_Mesh_create_simple_hex_layers(true);
  //MeshTests::Test_Mesh_vortex_shedding(true);
  //
  //MeshTests::Test_Mesh_add_quad_layer(true);

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

