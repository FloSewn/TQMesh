/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <vector>


#include "Vertex.h"
#include "Edge.h"
#include "Boundary.h"
#include "Domain.h"
#include "Mesh.h"
#include "Smoother.h"
#include "utils.h"

#include "ParaReader.h"
#include "Vec2.h"
#include "Helpers.h"
#include "Log.h"
#include "Container.h"

#include "size_function.h"


namespace TQMesh {

using namespace CppUtils;
using namespace TQAlgorithm;

/*--------------------------------------------------------------------
| Class for error handling
--------------------------------------------------------------------*/
class Invalid
{
public:
  Invalid(const std::string& msg){ error_message = msg; }
  const std::string& what() const { return error_message; }
private:
  std::string error_message;
};

/*********************************************************************
* A MeshGenerator entity  
*********************************************************************/
class MeshGenerator
{
public:
  /*------------------------------------------------------------------
  | Function for error handling
  ------------------------------------------------------------------*/
  void error(std::string msg) { throw Invalid{ msg }; }

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  MeshGenerator(int mesh_id ) 
  : mesh_id_ { mesh_id }
  {}

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  Mesh* mesh_ptr() 
  { 
    if ( mesh_ )
      return mesh_.get(); 
    return nullptr;
  }

  /*------------------------------------------------------------------
  | Read a new mesh from the input file and create it
  ------------------------------------------------------------------*/
  bool create_new_mesh(Mesh* base_mesh,
                       ParaReader& mesh_reader,
                       const std::string& prefix,
                       const std::string& format)
  {
    try
    {
      query_mandatory_parameters( mesh_reader );
      init_mesh_vertices( mesh_reader );
      init_mesh_domain( mesh_reader );
      init_domain_vertices( mesh_reader );
      init_exterior_boundary( mesh_reader );
      init_interior_boundaries( mesh_reader );
      init_fixed_vertices( mesh_reader );
      init_quad_layers( mesh_reader );
      init_element_color( mesh_reader );
      init_meshing_algorithm( mesh_reader );
      init_quad_refinements( mesh_reader );

      generate_mesh(base_mesh, prefix, format);
    }
    catch(const Invalid& inv)
    {
      LOG(ERROR) << inv.what();
      return false;
    }

    return true;

  } // MeshGenerator::create_new_mesh()

  /*------------------------------------------------------------------
  | Print out parameters 
  ------------------------------------------------------------------*/
  template<typename T>
  void print_parameter(ParaReader& reader, 
                       const std::string& name)
  {
    auto para = reader.get_parameter<T>(name);

    if (  para.type() == ParaType::scalar 
       || para.type() == ParaType::vector )
    {
      auto key = para.start_key();

      size_t n_values = para.columns();

      std::stringstream ss;

      for ( size_t i = 0; i < n_values; ++i )
      {
        ss << para.get_value(i);

        if ( i < n_values - 1 )
          ss << ", ";
      }

      LOG(INFO) << key << " " << ss.str();
      LOG(INFO) << "";

    }
    else
    {
      auto start_key = para.start_key();
      auto end_key = para.end_key();

      size_t n_rows = para.rows();
      size_t n_cols = para.columns();

      LOG(INFO) << start_key;

      for ( size_t j = 0; j < n_rows; ++j )
      {
        std::stringstream ss;

        for ( size_t i = 0; i < n_cols; ++i )
        {
          ss << para.get_value(i,j);

          if ( i < n_cols - 1 )
            ss << ", ";
        }
        
        LOG(INFO) << ss.str();
      }

      LOG(INFO) << end_key;
      LOG(INFO) << "";
    }

  } // MeshGenerator::print_parameter()


private:


  /*------------------------------------------------------------------
  | Query mandatory mesh parameters
  ------------------------------------------------------------------*/
  void query_mandatory_parameters(ParaReader& mesh_reader)
  {
    if ( !mesh_reader.query<std::string>("size_function") )
      error("Invalid size function definition for mesh " + mesh_id_);

    if ( !mesh_reader.query<double>("vertices") )
      error("Invalid definition of vertices for mesh " + mesh_id_);

  } // MeshGenerator::query_mandatory_parameters()

  /*------------------------------------------------------------------
  | Initialize the mesh vertices and estimate the domain extent
  ------------------------------------------------------------------*/
  void init_mesh_vertices(ParaReader& mesh_reader)
  {
    auto para_vertices = mesh_reader.get_parameter<double>("vertices");

    Vec2d xy_min {  DBL_MAX,  DBL_MAX };
    Vec2d xy_max { -DBL_MAX, -DBL_MAX };

    for ( size_t i = 0; i < para_vertices.rows(); ++i )
    {
      double x = para_vertices.get_value(0, i);
      double y = para_vertices.get_value(1, i);
      double s = para_vertices.get_value(2, i);
      double r = para_vertices.get_value(3, i);

      vertex_pos_.push_back( {x,y} );
      vertex_props_.push_back( {s, r} );

      // Estimate domain extents
      xy_min.x = MIN(xy_min.x, x);
      xy_min.y = MIN(xy_min.y, y);

      xy_max.x = MAX(xy_max.x, x);
      xy_max.y = MAX(xy_max.y, y);
    }

    const double dx = xy_max.x - xy_min.x;
    const double dy = xy_max.y - xy_min.y;
    domain_extent_  = ABS(10.0 * MAX(dx, dy));

    print_parameter<double>(mesh_reader, "vertices");

  } // MeshGenerator::init_mesh_vertices()

  /*------------------------------------------------------------------
  | Initialize the mesh domain
  ------------------------------------------------------------------*/
  void init_mesh_domain(ParaReader& mesh_reader)
  {
    // Initialize element size function
    UserSizeFunction size_fun = init_size_function( 
        mesh_reader.get_value<std::string>("size_function") 
    );

    // Initialize domain
    domain_ = std::make_unique<Domain>(size_fun, domain_extent_ );

    print_parameter<std::string>(mesh_reader, "size_function");

  } // MeshGenerator::init_mesh_domain()

  /*------------------------------------------------------------------
  | Initialize the domain's vertices 
  ------------------------------------------------------------------*/
  void init_domain_vertices(ParaReader& mesh_reader)
  {
    Domain& domain = *( domain_.get() );

    // Initialize domain vertices
    for ( size_t i = 0; i < vertex_pos_.size(); ++i )
    {
      const Vec2d& pos   = vertex_pos_[i];
      const Vec2d& props = vertex_props_[i];

      domain.add_vertex( pos.x, pos.y, props.x, props.y );
    }

  } // MeshGenerator::init_domain_vertices()

  /*------------------------------------------------------------------
  | Initialize the domain's exterior boundary 
  ------------------------------------------------------------------*/
  void init_exterior_boundary(ParaReader& mesh_reader)
  {
    ASSERT( domain_.get(), "DOMAIN NOT PROPERLY INITIALIZED" );

    Domain&   domain   = *( domain_.get() );
    Vertices& vertices = domain.vertices();

    // External boundary definition via direct edge placement
    mesh_reader.query<int>("extr_bdry");
    auto para_extr_bdry = mesh_reader.get_parameter<int>("extr_bdry");

    if ( para_extr_bdry.found() )
    {
      Boundary& b_ext    = domain.add_exterior_boundary();

      size_t n_vertices = vertices.size();

      for ( size_t i = 0; i < para_extr_bdry.rows(); ++i )
      {
        size_t i1 = para_extr_bdry.get_value(0, i);
        size_t i2 = para_extr_bdry.get_value(1, i);
        int    m  = para_extr_bdry.get_value(2, i);

        // Throw error if indices are larger than number of vertices
        if ( i1 > n_vertices || i2 > n_vertices )
          error("Invalid exterior boundary definition: " 
                "Some vertex index is larger that the number of "
                "provided input vertices.");

        Vertex& v1 = vertices[i1];
        Vertex& v2 = vertices[i2];

        b_ext.add_edge( v1, v2, m );
      }

      print_parameter<int>(mesh_reader, "extr_bdry");

      return;
    }

    // External boundary via rectangular boundary shape
    mesh_reader.query<double>("extr_bdry_rect");
    auto para_extr_bdry_rect
      = mesh_reader.get_parameter<double>("extr_bdry_rect");

    if ( para_extr_bdry_rect.found() )
    {
      Boundary& b_ext    = domain.add_exterior_boundary();
      
      print_parameter<double>(mesh_reader, "extr_bdry_rect");

      int    m = static_cast<int>( para_extr_bdry_rect.get_value(0) );
      double x = para_extr_bdry_rect.get_value(1);
      double y = para_extr_bdry_rect.get_value(2);
      double w = para_extr_bdry_rect.get_value(3);
      double h = para_extr_bdry_rect.get_value(4);

      b_ext.set_shape_rectangle( vertices, m, {x,y}, w, h );

      return;
    }

    // External boundary via circular boundary shape
    mesh_reader.query<double>("extr_bdry_circ");
    auto para_extr_bdry_circ
      = mesh_reader.get_parameter<double>("extr_bdry_circ");

    if ( para_extr_bdry_circ.found() )
    {
      Boundary& b_ext    = domain.add_exterior_boundary();
      
      print_parameter<double>(mesh_reader, "extr_bdry_circ");

      int    m = static_cast<int>( para_extr_bdry_circ.get_value(0) );
      double x = para_extr_bdry_circ.get_value(1);
      double y = para_extr_bdry_circ.get_value(2);
      double r = para_extr_bdry_circ.get_value(3);
      double n = static_cast<int>( para_extr_bdry_circ.get_value(4) );

      b_ext.set_shape_circle( vertices, m, {x,y}, r, n );

      return;
    }

    // External boundary via squared boundary shape
    mesh_reader.query<double>("extr_bdry_square");
    auto para_extr_bdry_square
      = mesh_reader.get_parameter<double>("extr_bdry_square");

    if ( para_extr_bdry_square.found() )
    {
      Boundary& b_ext    = domain.add_exterior_boundary();
      
      print_parameter<double>(mesh_reader, "extr_bdry_square");

      int    m = static_cast<int>( para_extr_bdry_square.get_value(0) );
      double x = para_extr_bdry_square.get_value(1);
      double y = para_extr_bdry_square.get_value(2);
      double h = para_extr_bdry_square.get_value(3);

      b_ext.set_shape_square( vertices, m, {x,y}, h );

      return;
    }

    // External boundary via triangular boundary shape
    mesh_reader.query<double>("extr_bdry_triangle");
    auto para_extr_bdry_triangle
      = mesh_reader.get_parameter<double>("extr_bdry_triangle");

    if ( para_extr_bdry_triangle.found() )
    {
      Boundary& b_ext    = domain.add_exterior_boundary();
      
      print_parameter<double>(mesh_reader, "extr_bdry_triangle");

      int    m = static_cast<int>( para_extr_bdry_triangle.get_value(0) );
      double x = para_extr_bdry_triangle.get_value(1);
      double y = para_extr_bdry_triangle.get_value(2);
      double h = para_extr_bdry_triangle.get_value(3);

      b_ext.set_shape_triangle( vertices, m, {x,y}, h );

      return;
    }


  } // MeshGenerator::init_exterior_boundary()

  /*------------------------------------------------------------------
  | Initialize the domain's interior boundaries
  ------------------------------------------------------------------*/
  void init_interior_boundaries(ParaReader& mesh_reader)
  {
    ASSERT( domain_.get(), "DOMAIN NOT PROPERLY INITIALIZED" );

    Domain&   domain   = *( domain_.get() );
    Vertices& vertices = domain.vertices();

    size_t n_vertices = vertices.size();

    // Query and initialize interior boundary definitions
    while( mesh_reader.query<int>("intr_bdry") )
    {
      Boundary& b_int = domain.add_interior_boundary();

      auto para_intr_bdry = mesh_reader.get_parameter<int>("intr_bdry");
      print_parameter<int>(mesh_reader, "intr_bdry");

      for ( size_t i = 0; i < para_intr_bdry.rows(); ++i )
      {
        size_t i1 = para_intr_bdry.get_value(0, i);
        size_t i2 = para_intr_bdry.get_value(1, i);
        int     m = para_intr_bdry.get_value(2, i);

        // Throw error if indices are larger than number of vertices
        if ( i1 > n_vertices || i2 > n_vertices )
          error("Invalid exterior boundary definition: " 
                "Some vertex index is larger that the number of "
                "provided input vertices.");

        Vertex& v1 = vertices[i1];
        Vertex& v2 = vertices[i2];

        b_int.add_edge( v1, v2, m );
      }
    }


    // Query and initialize interior boundary rectangles
    while( mesh_reader.query<double>("intr_bdry_rect") )
    {
      Boundary& b_int = domain.add_interior_boundary();

      auto para_intr_bdry_rect 
        = mesh_reader.get_parameter<double>("intr_bdry_rect");
      print_parameter<double>(mesh_reader, "intr_bdry_rect");

      int    m = static_cast<int>( para_intr_bdry_rect.get_value(0) );
      double x = para_intr_bdry_rect.get_value(1);
      double y = para_intr_bdry_rect.get_value(2);
      double w = para_intr_bdry_rect.get_value(3);
      double h = para_intr_bdry_rect.get_value(4);

      b_int.set_shape_rectangle( vertices, m, {x,y}, w, h );
    }


    // Query and initialize interior boundary circles
    while( mesh_reader.query<double>("intr_bdry_circ") )
    {
      Boundary& b_int = domain.add_interior_boundary();

      auto para_intr_bdry_circ 
        = mesh_reader.get_parameter<double>("intr_bdry_circ");
      print_parameter<double>(mesh_reader, "intr_bdry_circ");

      int    m = static_cast<int>( para_intr_bdry_circ.get_value(0) );
      double x = para_intr_bdry_circ.get_value(1);
      double y = para_intr_bdry_circ.get_value(2);
      double r = para_intr_bdry_circ.get_value(3);
      int    n = static_cast<int>( para_intr_bdry_circ.get_value(4) );

      b_int.set_shape_circle( vertices, m, {x,y}, r, n );
    }

    // Query and initialize interior boundary squares
    while( mesh_reader.query<double>("intr_bdry_square") )
    {
      Boundary& b_int = domain.add_interior_boundary();

      auto para_intr_bdry_square 
        = mesh_reader.get_parameter<double>("intr_bdry_square");
      print_parameter<double>(mesh_reader, "intr_bdry_square");

      int    m = static_cast<int>( para_intr_bdry_square.get_value(0) );
      double x = para_intr_bdry_square.get_value(1);
      double y = para_intr_bdry_square.get_value(2);
      double h = para_intr_bdry_square.get_value(3);

      b_int.set_shape_square( vertices, m, {x,y}, h );
    }

    // Query and initialize interior boundary triangle
    while( mesh_reader.query<double>("intr_bdry_triangle") )
    {
      Boundary& b_int = domain.add_interior_boundary();

      auto para_intr_bdry_triangle 
        = mesh_reader.get_parameter<double>("intr_bdry_triangle");
      print_parameter<double>(mesh_reader, "intr_bdry_triangle");

      int    m = static_cast<int>( para_intr_bdry_triangle.get_value(0) );
      double x = para_intr_bdry_triangle.get_value(1);
      double y = para_intr_bdry_triangle.get_value(2);
      double h = para_intr_bdry_triangle.get_value(3);

      b_int.set_shape_triangle( vertices, m, {x,y}, h );
    }

  } // MeshGenerator::init_interior_boundaries() 

  /*------------------------------------------------------------------
  | Initialize the domain's fixed vertices
  ------------------------------------------------------------------*/
  void init_fixed_vertices(ParaReader& mesh_reader)
  {
    ASSERT( domain_.get(), "DOMAIN NOT PROPERLY INITIALIZED" );

    Domain&   domain   = *( domain_.get() );

    if ( mesh_reader.query<double>("fixed_vertices") )
    {
      auto para_fixed_vertices 
        = mesh_reader.get_parameter<double>("fixed_vertices");
      print_parameter<double>(mesh_reader, "fixed_vertices");

      for ( size_t i = 0; i < para_fixed_vertices.rows(); ++i )
      {
        const double x = para_fixed_vertices.get_value(0, i);
        const double y = para_fixed_vertices.get_value(1, i);
        const double s = para_fixed_vertices.get_value(2, i);
        const double r = para_fixed_vertices.get_value(3, i);

        domain.add_fixed_vertex(x,y,s,r);
      }
    }

  } // MeshGenerator::init_fixed_vertices()

  /*------------------------------------------------------------------
  | Initialize quad layers 
  ------------------------------------------------------------------*/
  void init_quad_layers(ParaReader& mesh_reader)
  {
    ASSERT( domain_.get(), "DOMAIN NOT PROPERLY INITIALIZED" );

    Domain&   domain   = *( domain_.get() );
    Vertices& vertices = domain.vertices();

    while( mesh_reader.query<double>("quad_layers") )
    {
      auto para_quad_layers 
        = mesh_reader.get_parameter<double>("quad_layers");
      print_parameter<double>(mesh_reader, "quad_layers");

        int   i1 = static_cast<int>( para_quad_layers.get_value(0) );
        int   i2 = static_cast<int>( para_quad_layers.get_value(1) );
        int    n = static_cast<int>( para_quad_layers.get_value(2) );
        double h = para_quad_layers.get_value(3);
        double g = para_quad_layers.get_value(4);

        quad_layer_vertices_.push_back( {vertices[i1], vertices[i2]} );
        quad_layer_numbers_.push_back( n );
        quad_layer_heights_.push_back( h );
        quad_layer_growth_.push_back( g );
    }

  } // MeshGenerator::init_quad_layers()

  /*------------------------------------------------------------------
  | Initialize the meshing algorithm
  ------------------------------------------------------------------*/
  void init_meshing_algorithm(ParaReader& mesh_reader)
  {
    algorithm_ = "Triangulation";

    if ( mesh_reader.query<std::string>("algorithm") )
    {
      algorithm_ = mesh_reader.get_value<std::string>("algorithm");
      print_parameter<std::string>(mesh_reader, "algorithm");
    }

  } // MeshGenerator::init_meshing_algorithm()

  /*------------------------------------------------------------------
  | Initialize the number of quad refinements
  ------------------------------------------------------------------*/
  void init_quad_refinements(ParaReader& mesh_reader)
  {
    quad_refinements_ = 0;

    if ( mesh_reader.query<size_t>("quad_refinements") )
    {
      quad_refinements_ = mesh_reader.get_value<size_t>(
          "quad_refinements");
      print_parameter<size_t>(mesh_reader, "quad_refinements");
    }

  } // MeshGenerator::init_quad_refinements()

  /*------------------------------------------------------------------
  | Initialize the element color
  ------------------------------------------------------------------*/
  void init_element_color(ParaReader& mesh_reader)
  {
    element_color_ = 0;

    if ( mesh_reader.query<int>("elem_color") )
    {
      element_color_ = mesh_reader.get_value<int>( "elem_color" );
      print_parameter<int>(mesh_reader, "elem_color");
    }

  } // MeshGenerator::init_quad_refinements()


  /*------------------------------------------------------------------
  | Generate the mesh
  ------------------------------------------------------------------*/
  void generate_mesh(Mesh* base_mesh, 
                     const std::string& prefix, 
                     const std::string& format)
  {
    ASSERT( domain_.get(), "DOMAIN NOT PROPERLY INITIALIZED" );

    Domain&   domain   = *( domain_.get() );

    // Create the mesh
    mesh_ = std::make_unique<Mesh>(domain, mesh_id_, element_color_, 
                                   domain_extent_ );
    Mesh& mesh = *( mesh_.get() );

    // Connect mesh to base mesh
    if ( base_mesh )
      mesh.add_neighbor_mesh( *base_mesh );

    // Initialize the advancing front structure
    mesh.init_advancing_front();

    // Create quad layers
    for ( size_t i = 0; i < quad_layer_vertices_.size(); ++i )
    {
      Vertex& v1 = quad_layer_vertices_[i].first;
      Vertex& v2 = quad_layer_vertices_[i].second;
      int      n = quad_layer_numbers_[i];
      double   h = quad_layer_heights_[i];
      double   g = quad_layer_growth_[i];

      mesh.create_quad_layers(v1, v2, n, h, g);
    }

    // Start meshing 
    if ( algorithm_ == "Paving" )
    {
      mesh.pave();
    }
    else if ( algorithm_ == "Tri-to-Quad")
    {
      mesh.triangulate();
      mesh.merge_triangles_to_quads();
    }
    else if ( algorithm_ == "Triangulation")
    {
      mesh.triangulate();
    }
    else
    {
      error("Invalid meshing algorithm provided: " + algorithm_ );
    }

    if ( base_mesh )
      mesh.merge_neighbor_mesh( *base_mesh );

    // Apply mesh refinements
    for ( int i = 0; i < quad_refinements_; ++i )
      mesh.refine_to_quads();

    // Apply mesh smoothing
    Smoother smoother {};
    smoother.smooth( domain, mesh, 4, 0.5, 0.75, 0.95 );

    // Export the mesh
    std::string file_name = prefix + "_" + std::to_string(mesh_id_);

    if ( format == "VTU" || format == "vtu" )
    {
      LOG(INFO) << "Write mesh file to " << file_name << ".vtu";
      mesh.write_to_file( file_name + ".vtu", ExportType::vtu );
    }
    else if ( format == "TXT" || format == "txt" )
    {
      LOG(INFO) << "Write mesh file to " << file_name << ".txt";
      mesh.write_to_file( file_name + ".txt", ExportType::txt );
    }
    else
    {
      mesh.write_to_file( file_name, ExportType::cout );
    }

  } // MeshGenerator::generate_mesh()


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  int                                       mesh_id_;

  std::vector<Vec2d>                        vertex_pos_;
  std::vector<Vec2d>                        vertex_props_;

  double                                    domain_extent_;
  std::unique_ptr<Domain>                   domain_;

  std::vector<std::pair<Vertex&,Vertex&>>   quad_layer_vertices_ {};
  std::vector<int>                          quad_layer_numbers_  {};
  std::vector<double>                       quad_layer_heights_  {};
  std::vector<double>                       quad_layer_growth_   {};

  std::string                               algorithm_;
  size_t                                    quad_refinements_;
  int                                       element_color_;

  std::unique_ptr<Mesh>                     mesh_;

}; // MeshGenerator


/*********************************************************************
* The TQMesh main application
*********************************************************************/
class TQMeshApp
{
public:

  /*------------------------------------------------------------------
  | Function for error handling
  ------------------------------------------------------------------*/
  void error(std::string msg) { throw Invalid{ "[ERROR] " + msg}; }

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  TQMeshApp(const std::string& input_file)
  : reader_ { input_file }
  {
    init_para_reader();
  }

  /*------------------------------------------------------------------
  | Run the application
  ------------------------------------------------------------------*/
  bool run()
  {
    try
    {
      query_mandatory_parameters( );
    }
    catch(const Invalid& inv)
    {
      LOG(ERROR) << inv.what();
      return false;
    }

    ParaReader& mesh_reader = reader_.get_block("mesh_reader");
    int mesh_id = 0; 

    Mesh* base_mesh = nullptr;

    while( reader_.query( "mesh_reader" ) )
    {
      LOG(INFO) << "";
      LOG(INFO) << "============== "
                << "Create mesh " << mesh_id 
                << " ==============";

      meshes_.push_back( mesh_id );

      MeshGenerator& new_mesh = meshes_.back();

      bool success = new_mesh.create_new_mesh( base_mesh,
                                               mesh_reader, 
                                               output_prefix_, 
                                               output_format_ );

      if ( !success )
        return false;

      base_mesh = new_mesh.mesh_ptr();

      ++mesh_id;
    }


    return true;


  } // TQMeshApp::run()


private:
  /*------------------------------------------------------------------
  | Query mandatory mesh parameters
  ------------------------------------------------------------------*/
  void query_mandatory_parameters()
  {
    if ( !reader_.query<std::string>("output_prefix") )
      error("Invalid output file prefix" );

    if ( !reader_.query<std::string>("output_format") )
      error("Invalid mesh output format" );

    output_prefix_ = reader_.get_value<std::string>("output_prefix");
    output_format_ = reader_.get_value<std::string>("output_format");

  } // TQMeshApp::query_mandatory_parameters()


  /*------------------------------------------------------------------
  | Initialize file parameter 
  ------------------------------------------------------------------*/
  void init_para_reader()
  {
    reader_.new_scalar_parameter<std::string>(
        "output_prefix", "Prefix for mesh output file:");

    reader_.new_scalar_parameter<std::string>(
        "output_format", "Output format:");

    reader_.new_block_parameter(
        "mesh_reader", "Define mesh:", "End mesh");

    // Create sub-reader for mesh properties
    ParaReader& mesh_reader = reader_.get_block("mesh_reader");

    mesh_reader.new_matrix_parameter<double>(
        "vertices", "Define vertices:", "End vertices", 4);

    mesh_reader.new_scalar_parameter<std::string>(
        "size_function", "Element size:");

    mesh_reader.new_scalar_parameter<int>(
        "elem_color", "Element color:");

    mesh_reader.new_scalar_parameter<std::string>(
        "algorithm", "Meshing algorithm:");

    mesh_reader.new_scalar_parameter<size_t>(
        "quad_refinements", "Number of quad refinements:");

    mesh_reader.new_vector_parameter<double>(
        "quad_layers", "Add quad layers:", 5);

    mesh_reader.new_matrix_parameter<double>(
        "fixed_vertices", "Define fixed vertices:", "End fixed vertices", 4);

    mesh_reader.new_matrix_parameter<int>(
        "extr_bdry", "Define exterior boundary:", "End exterior boundary", 3);

    mesh_reader.new_vector_parameter<double>(
        "extr_bdry_rect", "Define exterior rectangular boundary:", 5);
    
    mesh_reader.new_vector_parameter<double>(
        "extr_bdry_circ", "Define exterior circular boundary:", 5);

    mesh_reader.new_vector_parameter<double>(
        "extr_bdry_square", "Define exterior squared boundary:", 4);

    mesh_reader.new_vector_parameter<double>(
        "extr_bdry_triangle", "Define exterior triangular boundary:", 4);

    mesh_reader.new_matrix_parameter<int>(
        "intr_bdry", "Define interior boundary:", "End interior boundary", 3);

    mesh_reader.new_vector_parameter<double>(
        "intr_bdry_rect", "Define interior rectangular boundary:", 5);
    
    mesh_reader.new_vector_parameter<double>(
        "intr_bdry_circ", "Define interior circular boundary:", 5);

    mesh_reader.new_vector_parameter<double>(
        "intr_bdry_square", "Define interior squared boundary:", 4);

    mesh_reader.new_vector_parameter<double>(
        "intr_bdry_triangle", "Define interior triangular boundary:", 4);

  } // TQMeshApp::init_para_reader()


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  ParaReader                 reader_;

  std::string                output_prefix_;
  std::string                output_format_;

  std::vector<MeshGenerator> meshes_;

}; // TQMeshApp

} // namespace TQMesh
