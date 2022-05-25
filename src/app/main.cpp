/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <iostream>
#include <cstdlib>

#include "utils.h"

#include "Container.h"

#include "Vertex.h"
#include "Edge.h"
#include "Boundary.h"
#include "Domain.h"
#include "Mesh.h"
#include "Smoother.h"

#include "ParaReader.h"
#include "size_function.h"

using namespace TQMesh::TQUtils;
using namespace TQMesh::TQAlgorithm;

int main(int argc, char* argv[])
{
  /*------------------------------------------------------------------
  | Handle command line arguments
  ------------------------------------------------------------------*/
  if ( argc < 2 )
  {
    MSG("TQMesh <ParameterFile>");
    return EXIT_SUCCESS;
  }

  /*------------------------------------------------------------------
  | Initialize parameter reader
  ------------------------------------------------------------------*/
  CppUtils::ParaReader reader { argv[1] };

  reader.new_scalar_parameter<std::string>(
      "size_function", "Element size:");

  reader.new_scalar_parameter<std::string>(
      "algorithm", "Meshing algorithm:");

  reader.new_scalar_parameter<int>(
      "quad_refinements", "Number of quad refinements:");

  reader.new_list_parameter<double>(
      "quad_layers", "Add quad layers:", 5);

  reader.new_list_parameter<double>(
      "vertices", "Define vertices:", "End vertices", 4);

  reader.new_list_parameter<double>(
      "fixed_vertices", "Define fixed vertices:", "End fixed vertices", 4);

  reader.new_list_parameter<int>(
      "extr_bdry", "Define exterior boundary:", "End exterior boundary", 3);

  reader.new_list_parameter<int>(
      "intr_bdry", "Define interior boundary:", "End interior boundary", 3);

  reader.new_list_parameter<double>(
      "intr_bdry_rect", "Define interior rectangular boundary:", 5);
  
  reader.new_list_parameter<double>(
      "intr_bdry_circ", "Define interior circular boundary:", 5);

  /*------------------------------------------------------------------
  | Initial query of mandatory parameters
  ------------------------------------------------------------------*/
  reader.query<std::string>("size_function");
  reader.query<double>("vertices");
  reader.query<int>("extr_bdry");

  if ( !reader.found("size_function") )
  {
    MSG("[ERROR] Invalid mesh size definition.");
    return EXIT_FAILURE;
  }

  if ( !reader.found("vertices") )
  {
    MSG("[ERROR] Invalid definition of boundary vertices.");
    return EXIT_FAILURE;
  }

  if ( !reader.found("extr_bdry") )
  {
    MSG("[ERROR] Invalid definition of the exterior boundary.");
    return EXIT_FAILURE;
  }

  /*------------------------------------------------------------------
  | Initialize element size function
  ------------------------------------------------------------------*/
  UserSizeFunction size_fun = init_size_function( 
      reader.get_value<std::string>("size_function") 
  );

  /*------------------------------------------------------------------
  | Query boundary vertices and obtain domain extent
  ------------------------------------------------------------------*/
  Vec2d xy_min {  DBL_MAX,  DBL_MAX };
  Vec2d xy_max { -DBL_MAX, -DBL_MAX };

  auto para_vertices = reader.get_parameter<double>("vertices");

  for ( size_t i = 0; i < para_vertices.rows(); ++i )
  {
    const double x = para_vertices.get_value(0, i);
    const double y = para_vertices.get_value(1, i);

    xy_min.x = MIN(xy_min.x, x);
    xy_min.y = MIN(xy_min.y, y);

    xy_max.x = MAX(xy_max.x, x);
    xy_max.y = MAX(xy_max.y, y);
  }

  const double dx = xy_max.x - xy_min.x;
  const double dy = xy_max.y - xy_min.y;
  const double domain_size = ABS(10.0 * MAX(dx, dy));

  /*------------------------------------------------------------------
  | Initialize vertices and exterior boundary
  ------------------------------------------------------------------*/
  Domain domain { size_fun, domain_size };

  // init vertices
  for ( size_t i = 0; i < para_vertices.rows(); ++i )
  {
    double x = para_vertices.get_value(0, i);
    double y = para_vertices.get_value(1, i);
    double s = para_vertices.get_value(2, i);
    double r = para_vertices.get_value(3, i);

    domain.add_vertex(x,y,s,r);
  }

  /*------------------------------------------------------------------
  | Initialize exterior boundary definitions
  ------------------------------------------------------------------*/
  Vertices& vertices = domain.vertices();

  Boundary& b_ext = domain.add_exterior_boundary();

  auto para_extr_bdry = reader.get_parameter<int>("extr_bdry");

  for ( size_t i = 0; i < para_extr_bdry.rows(); ++i )
  {
    size_t i1 = para_extr_bdry.get_value(0, i);
    size_t i2 = para_extr_bdry.get_value(1, i);
    int    m  = para_extr_bdry.get_value(2, i);

    Vertex& v1 = vertices[i1];
    Vertex& v2 = vertices[i2];

    b_ext.add_edge( v1, v2, m );
  }


  /*------------------------------------------------------------------
  | Query and initialize interior boundary definitions
  | (optional parameter)
  ------------------------------------------------------------------*/
  while( reader.query<int>("intr_bdry") )
  {
    Boundary& b_int = domain.add_interior_boundary();

    auto para_intr_bdry = reader.get_parameter<int>("intr_bdry");

    for ( size_t i = 0; i < para_intr_bdry.rows(); ++i )
    {
      size_t i1 = para_intr_bdry.get_value(0, i);
      size_t i2 = para_intr_bdry.get_value(1, i);
      int     m = para_intr_bdry.get_value(2, i);

      Vertex& v1 = vertices[i1];
      Vertex& v2 = vertices[i2];

      b_int.add_edge( v1, v2, m );
    }
  }

  /*------------------------------------------------------------------
  | Query and initialize interior boundary rectangles
  | (optional parameter)
  ------------------------------------------------------------------*/
  while( reader.query<double>("intr_bdry_rect") )
  {
    Boundary& b_int = domain.add_interior_boundary();

    auto para_intr_bdry_rect 
      = reader.get_parameter<double>("intr_bdry_rect");

    int    m = static_cast<int>( para_intr_bdry_rect.get_value(0) );
    double x = para_intr_bdry_rect.get_value(1);
    double y = para_intr_bdry_rect.get_value(2);
    double w = para_intr_bdry_rect.get_value(3);
    double h = para_intr_bdry_rect.get_value(4);

    b_int.set_shape_rectangle( vertices, m, {x,y}, w, h );
  }

  /*------------------------------------------------------------------
  | Query and initialize interior boundary circles
  | (optional parameter)
  ------------------------------------------------------------------*/
  while( reader.query<double>("intr_bdry_circ") )
  {
    Boundary& b_int = domain.add_interior_boundary();

    auto para_intr_bdry_circ 
      = reader.get_parameter<double>("intr_bdry_circ");

    int    m = static_cast<int>( para_intr_bdry_circ.get_value(0) );
    double x = para_intr_bdry_circ.get_value(1);
    double y = para_intr_bdry_circ.get_value(2);
    double r = para_intr_bdry_circ.get_value(3);
    double n = static_cast<int>( para_intr_bdry_circ.get_value(4) );

    b_int.set_shape_circle( vertices, m, {x,y}, r, n );
  }

  /*------------------------------------------------------------------
  | Query and intialize fixed vertex definitions
  | (optional parameter)
  ------------------------------------------------------------------*/
  if ( reader.query<double>("fixed_vertices") )
  {
    auto para_fixed_vertices 
      = reader.get_parameter<double>("fixed_vertices");

    for ( size_t i = 0; i < para_fixed_vertices.rows(); ++i )
    {
      const double x = para_fixed_vertices.get_value(0, i);
      const double y = para_fixed_vertices.get_value(1, i);
      const double s = para_fixed_vertices.get_value(2, i);
      const double r = para_fixed_vertices.get_value(3, i);

      domain.add_fixed_vertex(x,y,s,r);
    }
  }

  /*------------------------------------------------------------------
  | Query quad layer definitions and store vertices for creation
  | of quad layers
  | (optional parameter)
  ------------------------------------------------------------------*/
  std::vector<std::pair<Vertex&,Vertex&>> quad_layer_vertices {};
  std::vector<int> quad_layer_numbers {};
  std::vector<double> quad_layer_heights {};
  std::vector<double> quad_layer_growth {};

  while( reader.query<double>("quad_layers") )
  {
    auto para_quad_layers 
      = reader.get_parameter<double>("quad_layers");

      int   i1 = static_cast<int>( para_quad_layers.get_value(0) );
      int   i2 = static_cast<int>( para_quad_layers.get_value(1) );
      int    n = static_cast<int>( para_quad_layers.get_value(2) );
      double h = para_quad_layers.get_value(3);
      double g = para_quad_layers.get_value(4);

      quad_layer_vertices.push_back( {vertices[i1], vertices[i2]} );
      quad_layer_numbers.push_back( n );
      quad_layer_heights.push_back( h );
      quad_layer_growth.push_back( g );
  }

  /*------------------------------------------------------------------
  | Query meshing algorithm
  | (optional parameter)
  ------------------------------------------------------------------*/
  std::string meshing_algorithm = "Triangulation";

  if ( reader.query<std::string>("algorithm") )
    meshing_algorithm = reader.get_value<std::string>("algorithm");

  /*------------------------------------------------------------------
  | Query mesh refinements
  | (optional parameter)
  ------------------------------------------------------------------*/
  int n_quad_refinements = 0;

  if ( reader.query<int>("quad_refinements") )
    n_quad_refinements = reader.get_value<int>("quad_refinements");



  /*------------------------------------------------------------------
  | Create the mesh
  ------------------------------------------------------------------*/
  Mesh mesh { domain, domain_size };

  /*------------------------------------------------------------------
  | Create quad layers
  ------------------------------------------------------------------*/
  for ( size_t i = 0; i < quad_layer_vertices.size(); ++i )
  {
    Vertex& v1 = quad_layer_vertices[i].first;
    Vertex& v2 = quad_layer_vertices[i].second;
    int      n = quad_layer_numbers[i];
    double   h = quad_layer_heights[i];
    double   g = quad_layer_growth[i];

    mesh.create_quad_layers(v1, v2, n, h, g);
  }

  /*------------------------------------------------------------------
  | Start meshing with respective algorithm
  ------------------------------------------------------------------*/
  if ( meshing_algorithm == "Paving" )
  {
    mesh.pave();
  }
  else if ( meshing_algorithm == "Tri-to-Quad")
  {
    mesh.triangulate();
    mesh.merge_triangles_to_quads();
  }
  else if ( meshing_algorithm == "Triangulation")
  {
    mesh.triangulate();
  }
  else
  {
    MSG("[ERROR] Invalid meshing algorithm: " << meshing_algorithm);
    return EXIT_FAILURE;
  }

  /*------------------------------------------------------------------
  | Mesh refinements
  ------------------------------------------------------------------*/
  for ( int i = 0; i < n_quad_refinements; ++i )
    mesh.refine_to_quads();

  /*------------------------------------------------------------------
  | Mesh smoothing
  ------------------------------------------------------------------*/
  Smoother smoother {};
  smoother.smooth( domain, mesh, 6, 0.5, 0.75, 0.95 );

  /*------------------------------------------------------------------
  | Export meshing
  ------------------------------------------------------------------*/
  mesh.assign_mesh_indices();
  std::cout << mesh;
  domain.export_size_function(xy_min, xy_max, 100, 100);

  return EXIT_SUCCESS;
}

