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
  ParaReader reader { argv[1] };

  auto para_size_fun = reader.new_string_parameter(
      "Element size:");
  auto para_element_type = reader.new_string_parameter(
      "Element type:");
  auto para_quad_layers = reader.new_double_list_parameter(
      "Add quad layers:");

  auto para_vertices = reader.new_double_list_parameter(
      "Define vertices:", "End vertices");
  auto para_fixed_vertices = reader.new_double_list_parameter(
      "Define fixed vertices:", "End fixed vertices");

  auto para_extr_bdry = reader.new_int_list_parameter(
      "Define exterior boundary:", "End exterior boundary");
  auto para_intr_bdry = reader.new_int_list_parameter(
      "Define interior boundary:", "End interior boundary");

  auto para_intr_rect_bdry = reader.new_double_list_parameter(
      "Define interior rectangular boundary:");
  auto para_intr_circ_bdry = reader.new_double_list_parameter(
      "Define interior circular boundary:");

  /*------------------------------------------------------------------
  | Query element size & initialize domain
  ------------------------------------------------------------------*/
  reader.query( para_size_fun );

  UserSizeFunction size_fun;

  if ( para_size_fun.found() )
  {
    DBG_MSG("USER DEFINED SIZE FUNCTION");
    size_fun = init_size_function( para_size_fun.value() );
  }
  else 
  {
    MSG("No mesh size definition provided.");
    return EXIT_SUCCESS;
  }

  /*------------------------------------------------------------------
  | Query boundary vertices and obtain domain extent
  ------------------------------------------------------------------*/
  reader.query( para_vertices );

  if ( !para_vertices.found() )
  {
    MSG("Invalid vertex definition.");
    return EXIT_SUCCESS;
  }

  // Estimate domain extent
  Vec2d xy_min {  DBL_MAX,  DBL_MAX };
  Vec2d xy_max { -DBL_MAX, -DBL_MAX };

  for ( size_t i = 0; i < para_vertices.rows(); ++i )
  {
    const double x = para_vertices.value(4*i);
    const double y = para_vertices.value(4*i + 1);

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
    const double x = para_vertices.value(4*i);
    const double y = para_vertices.value(4*i + 1);
    const double s = para_vertices.value(4*i + 2);
    const double r = para_vertices.value(4*i + 3);

    domain.add_vertex(x,y,s,r);
  }

  /*------------------------------------------------------------------
  | Query exterior boundary definitions
  ------------------------------------------------------------------*/
  int n_extr_bdry = 0;
  Vertices& vertices = domain.vertices();

  while( reader.query( para_extr_bdry ) )
  {
    ++n_extr_bdry;
    Boundary& b_ext = domain.add_boundary( BdryType::EXTERIOR );

    for ( size_t i = 0; i < para_extr_bdry.rows(); ++i )
    {
      size_t i1 = para_extr_bdry.value(3*i);
      size_t i2 = para_extr_bdry.value(3*i + 1);
      int    m  = para_extr_bdry.value(3*i + 2);

      DBG_MSG("ADD EDGE: ("<< i1 << "," << i2 << ") -> MARKER " << m );

      Vertex& v1 = vertices[i1];
      Vertex& v2 = vertices[i2];

      b_ext.add_edge( v1, v2, m );
    }
  }

  /*------------------------------------------------------------------
  | Stop proceeding in case of bad exterior boundaries
  | (currently only a single exterior boundary can be handled)
  ------------------------------------------------------------------*/
  if ( n_extr_bdry != 1 )
  {
    MSG("ERROR");
    MSG("--------------------------------------------------------");
    MSG("Invalid exterior boundary definition.");
    MSG("Either no or too many exterior boundaries are defined.");
    MSG("Currently, only one single exterior boundary definition ");
    MSG("is provided in TQMesh.");
    MSG("");
    MSG("Defined exterior boundaries: " << n_extr_bdry);
    return EXIT_SUCCESS;
  }

  /*------------------------------------------------------------------
  | Query interior boundary definitions
  ------------------------------------------------------------------*/
  while( reader.query( para_intr_bdry ) )
  {
    Boundary& b_int = domain.add_boundary( BdryType::INTERIOR );

    for ( size_t i = 0; i < para_intr_bdry.rows(); ++i )
    {
      size_t i1 = para_intr_bdry.value(3*i);
      size_t i2 = para_intr_bdry.value(3*i + 1);
      int m  = para_intr_bdry.value(3*i + 2);

      Vertex& v1 = vertices[i1];
      Vertex& v2 = vertices[i2];

      b_int.add_edge( v1, v2, m );
    }
  }

  while( reader.query( para_intr_rect_bdry ) )
  {
    Boundary& b_int = domain.add_boundary( BdryType::INTERIOR );

    int    m = static_cast<int>( para_intr_rect_bdry.value(0) );
    double x = para_intr_rect_bdry.value(1);
    double y = para_intr_rect_bdry.value(2);
    double w = para_intr_rect_bdry.value(3);
    double h = para_intr_rect_bdry.value(4);

    b_int.set_shape_rectangle( vertices, m, {x,y}, w, h );
  }

  while( reader.query( para_intr_circ_bdry ) )
  {
    Boundary& b_int = domain.add_boundary( BdryType::INTERIOR );

    int    m = static_cast<int>( para_intr_circ_bdry.value(0) );
    double x = para_intr_circ_bdry.value(1);
    double y = para_intr_circ_bdry.value(2);
    double r = para_intr_circ_bdry.value(3);
    double n = static_cast<int>( para_intr_circ_bdry.value(4) );

    b_int.set_shape_circle( vertices, m, {x,y}, r, n );
  }

  /*------------------------------------------------------------------
  | Query fixed vertex definitions
  ------------------------------------------------------------------*/
  reader.query( para_fixed_vertices );

  if ( para_fixed_vertices.found() )
  {
    for ( size_t i = 0; i < para_fixed_vertices.rows(); ++i )
    {
      const double x = para_fixed_vertices.value(4*i);
      const double y = para_fixed_vertices.value(4*i + 1);
      const double s = para_fixed_vertices.value(4*i + 2);
      const double r = para_fixed_vertices.value(4*i + 3);

      domain.add_fixed_vertex(x,y,s,r);
    }
  }

  /*------------------------------------------------------------------
  | Query quad layer definitions and store vertices for creation
  | of quad layers
  ------------------------------------------------------------------*/
  std::vector<std::pair<Vertex&,Vertex&>> quad_layer_vertices {};
  std::vector<int> quad_layer_numbers {};
  std::vector<double> quad_layer_heights {};
  std::vector<double> quad_layer_growth {};

  while( reader.query( para_quad_layers ) )
  {
    int i1 = static_cast<int>( para_quad_layers.value(0) );
    int i2 = static_cast<int>( para_quad_layers.value(1) );
    int n = static_cast<int>( para_quad_layers.value(2) );
    double h = para_quad_layers.value(3);
    double g = para_quad_layers.value(4);

    quad_layer_vertices.push_back( {vertices[i1], vertices[i2]} );
    quad_layer_numbers.push_back( n );
    quad_layer_heights.push_back( h );
    quad_layer_growth.push_back( g );

  }

  /*------------------------------------------------------------------
  | Run meshing
  ------------------------------------------------------------------*/
  Mesh mesh { domain, domain_size };

  for ( size_t i = 0; i < quad_layer_vertices.size(); ++i )
  {
    Vertex& v1 = quad_layer_vertices[i].first;
    Vertex& v2 = quad_layer_vertices[i].second;
    int n = quad_layer_numbers[i];
    double h = quad_layer_heights[i];
    double g = quad_layer_growth[i];

    mesh.create_quad_layers(v1, v2, n, h, g);
  }


  // Pick element type - use triangles by default
  std::string element_type = "Triangle";

  reader.query( para_element_type );

  if ( para_element_type.found() )
    element_type = para_element_type.value();

  if ( element_type == "Quadrilateral" )
  {
    mesh.pave();
  }
  else if ( element_type == "Tri->Quad")
  {
    mesh.triangulate();
    mesh.merge_triangles_to_quads();
  }
  else
  {
    mesh.triangulate();
  }

  // Smooth the grid
  mesh.smoothing(4, 0.9);

  /*------------------------------------------------------------------
  | Export meshing
  ------------------------------------------------------------------*/
  mesh.assign_mesh_indices();
  std::cout << mesh;
  domain.export_size_function(xy_min, xy_max, 100, 100);

  return EXIT_SUCCESS;
}

