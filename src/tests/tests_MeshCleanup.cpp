/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/

#include <iostream>
#include <cassert>

#include <TQMeshConfig.h>

#include "tests.h"
#include "TestBuilder.h"

#include "VecND.h"
#include "Testing.h"
#include "Timer.h"
#include "Container.h"
 
#include "Vertex.h"
#include "Edge.h"
#include "Domain.h"
#include "Mesh.h"
#include "MeshCleanup.h"
#include "EntityChecks.h"

namespace CleanupTests 
{
using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;


/*********************************************************************
* Test cleanup of double quad edges
*********************************************************************/
void clear_double_quad_edges()
{
  int mesh_id = 0;
  int element_color = 0;

  Mesh mesh { mesh_id, element_color };

  Vertex& v1 = mesh.add_vertex({3.0, 0.0});
  Vertex& v2 = mesh.add_vertex({6.0, 0.0});
  Vertex& v3 = mesh.add_vertex({0.0, 3.0});
  Vertex& v4 = mesh.add_vertex({3.0, 3.0});
  Vertex& v5 = mesh.add_vertex({6.0, 3.0});
  Vertex& v6 = mesh.add_vertex({4.0, 4.0});
  Vertex& v7 = mesh.add_vertex({3.0, 6.0});
  Vertex& v8 = mesh.add_vertex({6.0, 6.0});

  Triangle& t1 = mesh.add_triangle(v3, v4, v7);
  Quad& q1 = mesh.add_quad(v1, v2, v5, v4);
  Quad& q2 = mesh.add_quad(v4, v5, v8, v6);
  Quad& q3 = mesh.add_quad(v4, v6, v8, v7);

  mesh.add_interior_edge(v4, v5);
  Edge& intr_edge_4_6 = mesh.add_interior_edge(v4, v6);
  mesh.add_interior_edge(v8, v6);
  mesh.add_interior_edge(v7, v4);

  mesh.add_boundary_edge(v1, v2, 1);
  mesh.add_boundary_edge(v2, v5, 2);
  mesh.add_boundary_edge(v5, v8, 2);
  mesh.add_boundary_edge(v8, v7, 3);
  mesh.add_boundary_edge(v7, v3, 3);
  mesh.add_boundary_edge(v3, v4, 4);
  mesh.add_boundary_edge(v4, v1, 4);

  CHECK( mesh.n_elements() == 4 );
  CHECK( mesh.n_quads() == 3 );
  CHECK( mesh.n_triangles() == 1 );
  CHECK( v4.facets().size() == 4 );
  CHECK( v6.facets().size() == 2 );

  CHECK( EntityChecks::check_mesh_validity(mesh) );

  MeshCleanup::assign_mesh_indices(mesh);
  MeshCleanup::setup_facet_connectivity(mesh);

  CHECK( q1.index() == 0 );
  CHECK( q2.index() == 1 );
  CHECK( q3.index() == 2 );
  CHECK( t1.index() == 3 );

  CHECK( intr_edge_4_6.facet_l() == &q3 );
  CHECK( intr_edge_4_6.facet_r() == &q2 );

  MeshCleanup::clear_double_quad_edges(mesh, false);

  CHECK( mesh.n_elements() == 3 );
  CHECK( mesh.n_quads() == 2 );
  CHECK( mesh.n_triangles() == 1 );

  LOG(INFO) << "\n" << mesh;

} // clear_double_quad_edges()

/*********************************************************************
* Test cleanup of double quad edges
*********************************************************************/
void clear_double_triangle_edges()
{
  int mesh_id = 0;
  int element_color = 0;

  Mesh mesh { mesh_id, element_color };

  Vertex& v1 = mesh.add_vertex({3.0, 0.0});
  Vertex& v2 = mesh.add_vertex({6.0, 0.0});
  Vertex& v3 = mesh.add_vertex({0.0, 3.0});
  Vertex& v4 = mesh.add_vertex({3.0, 3.0});
  Vertex& v5 = mesh.add_vertex({6.0, 3.0});
  Vertex& v6 = mesh.add_vertex({4.0, 4.0});
  Vertex& v7 = mesh.add_vertex({3.0, 6.0});

  mesh.add_triangle(v3, v4, v7);
  mesh.add_triangle(v4, v5, v6);

  mesh.add_quad(v1, v2, v5, v4);
  mesh.add_quad(v4, v6, v5, v7);

  mesh.add_interior_edge(v4, v5);
  mesh.add_interior_edge(v4, v6);
  mesh.add_interior_edge(v5, v6);
  mesh.add_interior_edge(v7, v4);

  mesh.add_boundary_edge(v1, v2, 1);
  mesh.add_boundary_edge(v2, v5, 2);
  mesh.add_boundary_edge(v5, v7, 2);
  mesh.add_boundary_edge(v7, v3, 3);
  mesh.add_boundary_edge(v3, v4, 4);
  mesh.add_boundary_edge(v4, v1, 4);

  CHECK( mesh.n_elements() == 4 );
  CHECK( mesh.n_quads() == 2 );
  CHECK( mesh.n_triangles() == 2 );
  CHECK( v4.facets().size() == 4 );
  CHECK( v6.facets().size() == 2 );

  CHECK( EntityChecks::check_mesh_validity(mesh) );

  MeshCleanup::assign_mesh_indices(mesh);
  MeshCleanup::setup_facet_connectivity(mesh);

  MeshCleanup::clear_double_triangle_edges(mesh, false);

  LOG(INFO) << "\n" << mesh;

} // clear_double_triangle_edges()

/*********************************************************************
* Test merge of degenerate triangles
*********************************************************************/
void merge_degenerate_triangles()
{
  int mesh_id = 0;
  int element_color = 0;

  Mesh mesh { mesh_id, element_color };

  Vertex& v1 = mesh.add_vertex({0.0, 0.0});
  Vertex& v2 = mesh.add_vertex({6.0, 0.0});
  Vertex& v3 = mesh.add_vertex({3.0, 3.0});
  Vertex& v4 = mesh.add_vertex({3.0, 6.0});

  mesh.add_triangle(v1, v2, v3);
  mesh.add_triangle(v2, v4, v3);
  mesh.add_triangle(v4, v1, v3);

  mesh.add_interior_edge(v1, v3);
  mesh.add_interior_edge(v2, v3);
  mesh.add_interior_edge(v3, v4);

  mesh.add_boundary_edge(v1, v2, 1);
  mesh.add_boundary_edge(v2, v4, 2);
  mesh.add_boundary_edge(v4, v1, 3);

  CHECK( mesh.n_elements() == 3 );
  CHECK( mesh.n_triangles() == 3 );
  CHECK( v3.facets().size() == 3 );

  CHECK( EntityChecks::check_mesh_validity(mesh) );

  MeshCleanup::merge_degenerate_triangles(mesh);

  CHECK( mesh.n_elements() == 1 );
  CHECK( mesh.n_triangles() == 1 );

  LOG(INFO) << "\n" << mesh;

} // merge_degenerate_triangles()

} // namespace CleanupTests


/*********************************************************************
* Run tests for: MeshCleanup.h
*********************************************************************/
void run_tests_MeshCleanup()
{
  adjust_logging_output_stream("CleanupTests.clear_double_quad_edges.log");
  CleanupTests::clear_double_quad_edges();

  adjust_logging_output_stream("CleanupTests.clear_double_triangle_edges.log");
  CleanupTests::clear_double_triangle_edges();

  adjust_logging_output_stream("CleanupTests.merge_degenerate_triangles.log");
  CleanupTests::merge_degenerate_triangles();

  // Reset debug logging ostream
  adjust_logging_output_stream("COUT");

} // run_tests_MeshCleanup()
