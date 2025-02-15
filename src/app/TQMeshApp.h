/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "TQMesh.h"
#include "size_function.h"

namespace TQMesh {

using namespace CppUtils;

/*********************************************************************
* The actual class to handle the generation of meshes 
*********************************************************************/
class MeshConstruction
{
public:
  using Vec2dPair = std::pair<Vec2d,Vec2d>;

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  MeshConstruction() {}

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

  /*------------------------------------------------------------------
  | Read a new mesh from the input file and create it
  ------------------------------------------------------------------*/
  bool construct_mesh(int mesh_id, ParaReader& mesh_reader)
  {
    mesh_id_ = mesh_id;

    query_mandatory_parameters( mesh_reader );

    init_mesh_vertices( mesh_reader );

    init_fixed_vertices( mesh_reader );

    init_mesh_domain( mesh_reader );

    init_domain_vertices( mesh_reader );

    init_exterior_boundary( mesh_reader );

    init_interior_boundaries( mesh_reader );

    init_fixed_edges( mesh_reader );

    init_quad_layers( mesh_reader );

    init_element_color( mesh_reader );

    init_meshing_algorithm( mesh_reader );

    init_refinement_parameters( mesh_reader );

    init_smoothing_parameters( mesh_reader );

    return generate_mesh();

  } // MeshConstruction::construct_mesh()

private:

  /*------------------------------------------------------------------
  | Generate the mesh
  ------------------------------------------------------------------*/
  bool generate_mesh()
  {
    ASSERT( domain_.get(), "MeshConstruction::generate_mesh: "
      "Domain has not been properly initialized." );

    Domain& domain = *( domain_.get() );

    // Create the mesh
    Mesh& mesh 
      = mesh_generator_.new_mesh( domain, mesh_id_, element_color_ );

    // Create quad layers
    for ( size_t i = 0; i < quad_layer_vertices_.size(); ++i )
    {
      const Vec2d& v1 = quad_layer_vertices_[i].first;
      const Vec2d& v2 = quad_layer_vertices_[i].second;
      int      n = quad_layer_numbers_[i];
      double   h = quad_layer_heights_[i];
      double   g = quad_layer_growth_[i];

      mesh_generator_.quad_layer_generation(mesh)
        .show_progress(true)
        .n_layers(n)
        .first_height(h)
        .growth_rate(g)
        .starting_position(v1)
        .ending_position(v2)
        .generate_elements();
    }

    // Start meshing 
    if ( algorithm_ == "Tri-to-Quad")
    {
      mesh_generator_.triangulation(mesh)
        .show_progress(true)
        .generate_elements();

      mesh_generator_.tri2quad_modification(mesh).modify();
    }
    else if ( algorithm_ == "Triangulation")
    {
      mesh_generator_.triangulation(mesh)
        .show_progress(true)
        .generate_elements();
    }
    else
    {
      throw_error("Invalid meshing algorithm provided: " + algorithm_ );
    }

    // Merge with other meshes
    if ( mesh_generator_.size() > 1 )
    {
      Mesh& other_mesh = mesh_generator_.mesh(0);
      ASSERT( &other_mesh != &mesh, "MeshConstruction::generate_mesh: "
        "Failed to access other mesh for merge operation.");
      mesh_generator_.merge_meshes(mesh, other_mesh);
      ASSERT( mesh_generator_.size() == 1, 
        "MeshConstruction::generate_mesh: "
        "Invalid total number of defined meshes.");
    }

    // Apply mesh refinements
    for ( std::size_t i = 0; i < quad_refinements_; ++i )
      mesh_generator_.quad_refinement(mesh).refine();

    // Apply mesh smoothing
    mesh_generator_.mixed_smoothing(mesh)
      .quad_layer_smoothing(smooth_quad_layers_)
      .smooth(smoothing_iterations_);

    // Finished progress bar requires newline
    LOG(INFO) << "\n";

    // Export the mesh
    if ( output_format_ == "VTU" || output_format_ == "vtu" )
    {
      std::string filename { output_prefix_ + ".vtu" };
      LOG(INFO) << "Write mesh file to " << filename;
      mesh_generator_.write_mesh(mesh, filename, MeshExportType::VTU);
    }
    else if ( output_format_ == "TXT" || output_format_ == "txt" )
    {
      std::string filename { output_prefix_ + ".txt" };
      LOG(INFO) << "Write mesh file to " << filename;
      mesh_generator_.write_mesh(mesh, filename, MeshExportType::TXT);
    }
    else
    {
      mesh_generator_.write_mesh(mesh, "DUMMY", MeshExportType::COUT);
    }

    // Check if meshing algorithm succeeded
    MeshChecker mesh_check { mesh, domain };
    if ( !mesh_check.check_completeness() )
    {
      LOG(ERROR) << "Failed to generate a valid mesh.";
      return false;
    }

    // Mesh generation succeeded 
    return true;

  } // MeshConstruction::generate_mesh()

  /*------------------------------------------------------------------
  | Initialize smoothing parameters
  ------------------------------------------------------------------*/
  void init_smoothing_parameters(ParaReader& mesh_reader)
  {
    smoothing_iterations_ = 2;
    smooth_quad_layers_   = false;

    if ( mesh_reader.query<size_t>("smoothing_iterations") )
    {
      smoothing_iterations_ = mesh_reader.get_value<size_t>("smoothing_iterations");
      print_parameter<size_t>(mesh_reader, "smoothing_iterations");
    }

    if ( mesh_reader.query<bool>("smooth_quad_layers") )
    {
      smooth_quad_layers_ = mesh_reader.get_value<bool>("smooth_quad_layers");
      print_parameter<bool>(mesh_reader, "smooth_quad_layers");
    }

  } // MeshConstruction::init_smoothing_parameters()

  /*------------------------------------------------------------------
  | Initialize the number of quad refinements
  ------------------------------------------------------------------*/
  void init_refinement_parameters(ParaReader& mesh_reader)
  {
    quad_refinements_ = 0;

    if ( mesh_reader.query<size_t>("quad_refinements") )
    {
      quad_refinements_ = mesh_reader.get_value<size_t>("quad_refinements");
      print_parameter<size_t>(mesh_reader, "quad_refinements");
    }

  } // MeshConstruction::init_refinement_parameters()

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

  } // MeshConstruction::init_meshing_algorithm()

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

  } // MeshConstruction::init_element_color()

  /*------------------------------------------------------------------
  | Initialize quad layers 
  ------------------------------------------------------------------*/
  void init_quad_layers(ParaReader& mesh_reader)
  {
    ASSERT( domain_.get(), "MeshConstruction::init_quad_layers: "
      "Domain has not been properly initialized." );

    quad_layer_vertices_.clear();
    quad_layer_numbers_.clear();
    quad_layer_heights_.clear();
    quad_layer_growth_.clear();

    while( mesh_reader.query<double>("quad_layers") )
    {
      auto para_quad_layers 
        = mesh_reader.get_parameter<double>("quad_layers");
      print_parameter<double>(mesh_reader, "quad_layers");

        double x1 = para_quad_layers.get_value(0);
        double y1 = para_quad_layers.get_value(1);

        double x2 = para_quad_layers.get_value(2);
        double y2 = para_quad_layers.get_value(3);

        int    n = static_cast<int>( para_quad_layers.get_value(4) );

        double h = para_quad_layers.get_value(5);
        double g = para_quad_layers.get_value(6);

        quad_layer_vertices_.push_back( { {x1,y1}, {x2,y2} } );
        quad_layer_numbers_.push_back( n );
        quad_layer_heights_.push_back( h );
        quad_layer_growth_.push_back( g );
    }

  } // MeshConstruction::init_quad_layers()

  /*------------------------------------------------------------------
  | Initialize the domain's fixed vertices
  ------------------------------------------------------------------*/
  void init_fixed_vertices(ParaReader& mesh_reader)
  {
    if ( mesh_reader.query<double>("fixed_vertices", true, -1.0) )
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

        vertex_pos_.push_back( {x,y} );
        vertex_props_.push_back( {s, r} );
        vertex_fixed_.push_back( true );
      }
    }

  } // MeshConstruction::init_fixed_vertices()

  /*------------------------------------------------------------------
  | Initialize the domain's fixed edges
  ------------------------------------------------------------------*/
  void init_fixed_edges(ParaReader& mesh_reader)
  {
    Domain&   domain   = *( domain_.get() );
    Vertices& vertices = domain.vertices();

    // Fixed edge definition via direct edge placement
    while( mesh_reader.query<int>("fixed_edges") )
    {
      size_t n_vertices = vertices.size();

      auto para_fixed_edge 
        = mesh_reader.get_parameter<int>("fixed_edges");

      print_parameter<int>(mesh_reader, "fixed_edges");

      for ( size_t i = 0; i < para_fixed_edge.rows(); ++i )
      {
        size_t i1 = para_fixed_edge.get_value(0, i);
        size_t i2 = para_fixed_edge.get_value(1, i);

        // Throw error if indices are larger than number of vertices
        if ( i1 >= n_vertices || i2 >= n_vertices )
          throw_error("Invalid fixed edge definition: " 
            "Some vertex index is larger that the number of "
            "provided input vertices.");

        Vertex& v1 = vertices[i1];
        Vertex& v2 = vertices[i2];

        domain.add_fixed_edge( v1, v2 );
      }
    }

  } // MeshConstruction::init_fixed_edges() 

  /*------------------------------------------------------------------
  | Initialize the domain's interior boundaries
  ------------------------------------------------------------------*/
  void init_interior_boundaries(ParaReader& mesh_reader)
  {
    ASSERT( domain_.get(), "MeshConstruction::init_interior_boundary: "
      "Domain has not been properly initialized." );

    Domain&   domain   = *( domain_.get() );
    Vertices& vertices = domain.vertices();

    // Query and initialize interior boundaries from edge definitions
    while( mesh_reader.query<int>("intr_bdry_edges") )
    {
      size_t n_vertices = vertices.size();

      Boundary& b_int = domain.add_interior_boundary();

      auto para_intr_bdry 
        = mesh_reader.get_parameter<int>("intr_bdry_edges");

      print_parameter<int>(mesh_reader, "intr_bdry_edges");

      for ( size_t i = 0; i < para_intr_bdry.rows(); ++i )
      {
        size_t i1 = para_intr_bdry.get_value(0, i);
        size_t i2 = para_intr_bdry.get_value(1, i);
        int     m = para_intr_bdry.get_value(2, i);

        // Throw error if indices are larger than number of vertices
        if ( i1 >= n_vertices || i2 >= n_vertices )
          throw_error("Invalid interior boundary definition: " 
            "Some vertex index is larger that the number of "
            "provided input vertices.");

        Vertex& v1 = vertices[i1];
        Vertex& v2 = vertices[i2];

        b_int.add_edge( v1, v2, m );
      }
    }


    // Query and initialize interior boundary coordinate definitions
    while( mesh_reader.query<double>("intr_bdry_coords", true, -1.0) )
    {
      Boundary& b_int = domain.add_interior_boundary();

      auto para_intr_bdry 
        = mesh_reader.get_parameter<double>("intr_bdry_coords");
      print_parameter<double>(mesh_reader, "intr_bdry_coords");

      std::vector<Vec2d> vertex_coords {};
      std::vector<Vec2d> vertex_props {};
      std::vector<int>   edge_colors {};

      for ( size_t i = 0; i < para_intr_bdry.rows(); ++i )
      {
        double x = para_intr_bdry.get_value(0, i);
        double y = para_intr_bdry.get_value(1, i);
        int m    = static_cast<int>( para_intr_bdry.get_value(2, i) );
        double s = para_intr_bdry.get_value(3, i);
        double r = para_intr_bdry.get_value(4, i);

        if ( m < 0 )
          throw_error("Invalid interior boundary definition");

        vertex_coords.push_back( {x,y} );
        vertex_props.push_back( {s,r} );
        edge_colors.push_back( m );
      }

      b_int.set_shape_from_coordinates(vertex_coords, edge_colors,
                                       vertex_props); 
    }


    // Query and initialize interior boundary rectangles
    while( mesh_reader.query<double>("intr_bdry_rect", true, -1.0) )
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
      double s = para_intr_bdry_rect.get_value(5);
      double r = para_intr_bdry_rect.get_value(6);

      b_int.set_shape_rectangle( m, {x,y}, w, h, s, r );
    }


    // Query and initialize interior boundary circles
    while( mesh_reader.query<double>("intr_bdry_circ", true, -1.0) )
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
      double ms = para_intr_bdry_circ.get_value(5);
      double mr = para_intr_bdry_circ.get_value(6);

      b_int.set_shape_circle( m, {x,y}, r, n, ms, mr );
    }

    // Query and initialize interior boundary squares
    while( mesh_reader.query<double>("intr_bdry_square", true, -1.0) )
    {
      Boundary& b_int = domain.add_interior_boundary();

      auto para_intr_bdry_square 
        = mesh_reader.get_parameter<double>("intr_bdry_square");
      print_parameter<double>(mesh_reader, "intr_bdry_square");

      int    m = static_cast<int>( para_intr_bdry_square.get_value(0) );
      double x = para_intr_bdry_square.get_value(1);
      double y = para_intr_bdry_square.get_value(2);
      double h = para_intr_bdry_square.get_value(3);
      double ms = para_intr_bdry_square.get_value(4);
      double mr = para_intr_bdry_square.get_value(5);

      b_int.set_shape_square( m, {x,y}, h, ms, mr );
    }

    // Query and initialize interior boundary triangle
    while( mesh_reader.query<double>("intr_bdry_triangle", true, -1.0) )
    {
      Boundary& b_int = domain.add_interior_boundary();

      auto para_intr_bdry_triangle 
        = mesh_reader.get_parameter<double>("intr_bdry_triangle");
      print_parameter<double>(mesh_reader, "intr_bdry_triangle");

      int    m = static_cast<int>( para_intr_bdry_triangle.get_value(0) );
      double x = para_intr_bdry_triangle.get_value(1);
      double y = para_intr_bdry_triangle.get_value(2);
      double h = para_intr_bdry_triangle.get_value(3);
      double ms = para_intr_bdry_triangle.get_value(4);
      double mr = para_intr_bdry_triangle.get_value(5);

      b_int.set_shape_triangle( m, {x,y}, h, ms, mr );
    }

    // Query and initialize interior boundaries from CSV file
    while( mesh_reader.query<std::string>("intr_bdry_csv") )
    {
      Boundary& b_int = domain.add_interior_boundary();

      b_int.set_shape_from_csv( 
          mesh_reader.get_value<std::string>("intr_bdry_csv") );

      print_parameter<std::string>(mesh_reader, "intr_bdry_csv");
    }

  } // MeshConstruction::init_interior_boundaries() 

  /*------------------------------------------------------------------
  | Initialize the domain's exterior boundary 
  ------------------------------------------------------------------*/
  void init_exterior_boundary(ParaReader& mesh_reader)
  {
    ASSERT( domain_.get(), "MeshConstruction::init_exterior_boundary: "
      "Domain has not been properly initialized." );

    Domain&   domain   = *( domain_.get() );
    Vertices& vertices = domain.vertices();

    // External boundary definition via direct edge placement
    if ( mesh_reader.query<int>("extr_bdry_edges") )
    {
      auto para_extr_bdry 
        = mesh_reader.get_parameter<int>("extr_bdry_edges");

      Boundary& b_ext = domain.add_exterior_boundary();

      size_t n_vertices = vertices.size();

      for ( size_t i = 0; i < para_extr_bdry.rows(); ++i )
      {
        size_t i1 = para_extr_bdry.get_value(0, i);
        size_t i2 = para_extr_bdry.get_value(1, i);
        int    m  = para_extr_bdry.get_value(2, i);

        // Throw error if indices are larger than number of vertices
        if ( i1 > n_vertices || i2 > n_vertices )
          throw_error("Invalid exterior boundary definition: " 
            "Some vertex index is larger that the number of "
            "provided input vertices.");

        Vertex& v1 = vertices[i1];
        Vertex& v2 = vertices[i2];

        b_ext.add_edge( v1, v2, m );
      }

      print_parameter<int>(mesh_reader, "extr_bdry_edges");

      return;
    }

    // External boundary  definition via direct edge coordinates 
    if ( mesh_reader.query<double>("extr_bdry_coords", true, -1.0) )
    {
      auto para_extr_bdry 
        = mesh_reader.get_parameter<double>("extr_bdry_coords");

      Boundary& b_ext = domain.add_exterior_boundary();

      std::vector<Vec2d> vertex_coords {};
      std::vector<Vec2d> vertex_props {};
      std::vector<int>   edge_colors {};

      for ( size_t i = 0; i < para_extr_bdry.rows(); ++i )
      {
        double x = para_extr_bdry.get_value(0, i);
        double y = para_extr_bdry.get_value(1, i);
        int m    = static_cast<int>( para_extr_bdry.get_value(2, i) );
        double s = para_extr_bdry.get_value(3, i);
        double r = para_extr_bdry.get_value(4, i);

        if ( m < 0 )
          throw_error("Invalid exterior boundary definition");

        vertex_coords.push_back( {x,y} );
        vertex_props.push_back( {s,r} );
        edge_colors.push_back( m );
      }

      b_ext.set_shape_from_coordinates(vertex_coords, edge_colors,
                                       vertex_props);

      print_parameter<double>(mesh_reader, "extr_bdry_coords");

      return;
    }


    // External boundary via rectangular boundary shape
    if ( mesh_reader.query<double>("extr_bdry_rect", true, -1.0) )
    {
      auto para_extr_bdry
        = mesh_reader.get_parameter<double>("extr_bdry_rect");

      Boundary& b_ext = domain.add_exterior_boundary();
      
      print_parameter<double>(mesh_reader, "extr_bdry_rect");

      int    m = static_cast<int>( para_extr_bdry.get_value(0) );
      double x = para_extr_bdry.get_value(1);
      double y = para_extr_bdry.get_value(2);
      double w = para_extr_bdry.get_value(3);
      double h = para_extr_bdry.get_value(4);
      double s = para_extr_bdry.get_value(5);
      double r = para_extr_bdry.get_value(6);

      b_ext.set_shape_rectangle( m, {x,y}, w, h, s, r );

      return;
    }

    // External boundary via circular boundary shape
    if ( mesh_reader.query<double>("extr_bdry_circ", true, -1.0) )
    {
      auto para_extr_bdry
        = mesh_reader.get_parameter<double>("extr_bdry_circ");

      Boundary& b_ext    = domain.add_exterior_boundary();
      
      print_parameter<double>(mesh_reader, "extr_bdry_circ");

      int    m = static_cast<int>( para_extr_bdry.get_value(0) );
      double x = para_extr_bdry.get_value(1);
      double y = para_extr_bdry.get_value(2);
      double r = para_extr_bdry.get_value(3);
      double n = static_cast<int>( para_extr_bdry.get_value(4) );
      double ms = para_extr_bdry.get_value(5);
      double mr = para_extr_bdry.get_value(6);

      b_ext.set_shape_circle( m, {x,y}, r, n, ms, mr );

      return;
    }

    // External boundary via squared boundary shape
    if ( mesh_reader.query<double>("extr_bdry_square", true, -1.0) )
    {
      auto para_extr_bdry
        = mesh_reader.get_parameter<double>("extr_bdry_square");

      Boundary& b_ext    = domain.add_exterior_boundary();
      
      print_parameter<double>(mesh_reader, "extr_bdry_square");

      int    m = static_cast<int>( para_extr_bdry.get_value(0) );
      double x = para_extr_bdry.get_value(1);
      double y = para_extr_bdry.get_value(2);
      double h = para_extr_bdry.get_value(3);
      double ms = para_extr_bdry.get_value(4);
      double mr = para_extr_bdry.get_value(5);

      b_ext.set_shape_square( m, {x,y}, h, ms, mr );

      return;
    }

    // External boundary via triangular boundary shape
    if ( mesh_reader.query<double>("extr_bdry_triangle", true, -1.0) )
    {
      auto para_extr_bdry
        = mesh_reader.get_parameter<double>("extr_bdry_triangle");

      Boundary& b_ext    = domain.add_exterior_boundary();
      
      print_parameter<double>(mesh_reader, "extr_bdry_triangle");

      int    m = static_cast<int>( para_extr_bdry.get_value(0) );
      double x = para_extr_bdry.get_value(1);
      double y = para_extr_bdry.get_value(2);
      double h = para_extr_bdry.get_value(3);
      double ms = para_extr_bdry.get_value(4);
      double mr = para_extr_bdry.get_value(5);

      b_ext.set_shape_triangle( m, {x,y}, h, ms, mr );

      return;
    }

    // External boundary definition from CSV file
    if ( mesh_reader.query<std::string>("extr_bdry_csv") )
    {
      Boundary& b_ext = domain.add_exterior_boundary();
      b_ext.set_shape_from_csv( 
          mesh_reader.get_value<std::string>("extr_bdry_csv") );

      print_parameter<std::string>(mesh_reader, "extr_bdry_csv");

      return;
    }

  } // MeshConstruction::init_exterior_boundary()

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

      if ( vertex_fixed_[i] )
        domain.add_fixed_vertex( pos.x, pos.y, props.x, props.y );
      else
        domain.add_vertex( pos.x, pos.y, props.x, props.y );
    }

  } // MeshConstruction::init_domain_vertices()

  /*------------------------------------------------------------------
  | Initialize the mesh domain
  ------------------------------------------------------------------*/
  void init_mesh_domain(ParaReader& mesh_reader)
  {
    TQMeshSetup::get_instance().set_quadtree_scale( domain_extent_ );
    // Initialize element size function
    UserSizeFunction size_fun = init_size_function( 
        mesh_reader.get_value<std::string>("size_function") 
    );

    // Initialize domain
    domain_ = std::make_unique<Domain>(size_fun );

    print_parameter<std::string>(mesh_reader, "size_function");

  } // MeshConstruction::init_mesh_domain()

  /*------------------------------------------------------------------
  | Initialize the mesh vertices and estimate the domain extent
  ------------------------------------------------------------------*/
  void init_mesh_vertices(ParaReader& mesh_reader)
  {
    // Initialize values
    vertex_pos_.clear();
    vertex_props_.clear();
    vertex_fixed_.clear();
    domain_extent_ = 0.0;

    // Vertices are given in the input file
    // -> Put vertex coordinates and properties into arrays 
    //    and set domain extents
    if ( mesh_reader.query<double>("vertices", true, -1.0) )
    {
      auto para_vertices 
        = mesh_reader.get_parameter<double>("vertices");

      print_parameter<double>(mesh_reader, "vertices");

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
        vertex_fixed_.push_back( false );

        // Estimate domain extents
        domain_extent_ = MAX(domain_extent_, ABS(x));
        domain_extent_ = MAX(domain_extent_, ABS(y));
      }

      // Double extent and enlarge slightly
      domain_extent_ *= 2.1;
    }
    // No input vertices given 
    // -> Set domain extent to large value
    else
    {
      domain_extent_  = TQ_MAX;
    }

  } // MeshConstruction::init_mesh_vertices()


  /*------------------------------------------------------------------
  | Query mandatory mesh parameters
  ------------------------------------------------------------------*/
  void query_mandatory_parameters(ParaReader& mesh_reader)
  {
    if ( !mesh_reader.query<std::string>("output_prefix") )
      throw_error("No output file prefix defined for mesh " + mesh_id_ );

    if ( !mesh_reader.query<std::string>("output_format") )
      throw_error("No output format defined for mesh " + mesh_id_ ); 

    output_prefix_ = mesh_reader.get_value<std::string>("output_prefix");
    output_format_ = mesh_reader.get_value<std::string>("output_format");

    print_parameter<std::string>(mesh_reader, "output_prefix");
    print_parameter<std::string>(mesh_reader, "output_format");

    if ( !mesh_reader.query<std::string>("size_function") )
      throw_error("Invalid size function definition for mesh " + mesh_id_);

  } // MeshConstruction::query_mandatory_parameters()


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  int                     mesh_id_;

  MeshGenerator           mesh_generator_ {};

  std::string             output_prefix_;
  std::string             output_format_;

  std::vector<Vec2d>      vertex_pos_;
  std::vector<Vec2d>      vertex_props_;
  std::vector<bool>       vertex_fixed_;

  double                  domain_extent_;
  std::unique_ptr<Domain> domain_;

  std::vector<Vec2dPair>  quad_layer_vertices_ {};
  std::vector<int>        quad_layer_numbers_  {};
  std::vector<double>     quad_layer_heights_  {};
  std::vector<double>     quad_layer_growth_   {};

  std::string             algorithm_;
  int                     element_color_;

  size_t                  quad_refinements_;

  size_t                  smoothing_iterations_;
  bool                    smooth_quad_layers_;


}; // MeshConstruction


/*********************************************************************
* The TQMesh main application
*********************************************************************/
class TQMeshApp
{
public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  TQMeshApp(const std::string& input_file)
  : reader_ { input_file }
  {
    init_parameter_file_reader();
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
    catch(const Error& inv)
    {
      LOG(ERROR) << inv.what();
      return false;
    }

    MeshConstruction mesh_construction {};
    int mesh_id = 0; 
    bool success = true;

    while( reader_.query( "mesh_reader" ) )
    {
      ParaReader& mesh_reader = reader_.get_block("mesh_reader");

      LOG(INFO) << "";
      LOG(INFO) << "============== " << "Create mesh " << mesh_id 
                << " ==============";

      success &= mesh_construction.construct_mesh(mesh_id, mesh_reader);

      ++mesh_id;
    }

    return success;

  } // TQMeshApp::run()


private:
  /*------------------------------------------------------------------
  | Query mandatory mesh parameters
  ------------------------------------------------------------------*/
  void query_mandatory_parameters()
  {

  } // TQMeshApp::query_mandatory_parameters()


  /*------------------------------------------------------------------
  | Initialize file parameter 
  ------------------------------------------------------------------*/
  void init_parameter_file_reader()
  {
    reader_.new_block_parameter(
        "mesh_reader", "Define mesh:", "End mesh");

    // Create sub-reader for mesh properties
    ParaReader& mesh_reader = reader_.get_block("mesh_reader");

    mesh_reader.new_scalar_parameter<std::string>(
        "output_prefix", "Output file prefix:");

    mesh_reader.new_scalar_parameter<std::string>(
        "output_format", "Output file format:");

    mesh_reader.new_matrix_parameter<double>(
        "vertices", "Define boundary vertices:", 
        "End boundary vertices", 4);

    mesh_reader.new_scalar_parameter<std::string>(
        "size_function", "Element size:");

    mesh_reader.new_scalar_parameter<int>(
        "elem_color", "Element color:");

    mesh_reader.new_scalar_parameter<std::string>(
        "algorithm", "Meshing algorithm:");

    mesh_reader.new_scalar_parameter<size_t>(
        "quad_refinements", "Number of quad refinements:");

    mesh_reader.new_scalar_parameter<size_t>(
        "smoothing_iterations", "Number of smoothing iterations:");

    mesh_reader.new_scalar_parameter<bool>(
        "smooth_quad_layers", "Smooth quad layers:");

    mesh_reader.new_vector_parameter<double>(
        "quad_layers", "Add quad layers:", 7);

    mesh_reader.new_matrix_parameter<double>(
        "fixed_vertices", "Define fixed vertices:", 
        "End fixed vertices", 4);

    mesh_reader.new_matrix_parameter<int>(
        "fixed_edges", "Define fixed edges:", 
        "End fixed edges", 2);


    mesh_reader.new_scalar_parameter<std::string>(
        "extr_bdry_csv", "Define exterior boundary from CSV file:");

    mesh_reader.new_scalar_parameter<std::string>(
        "intr_bdry_csv", "Define interior boundary from CSV file:");


    mesh_reader.new_matrix_parameter<int>(
        "extr_bdry_edges", "Define exterior boundary edges:", 
        "End exterior boundary edges", 3);

    mesh_reader.new_matrix_parameter<int>(
        "intr_bdry_edges", "Define interior boundary edges:", 
        "End interior boundary edges", 3);


    mesh_reader.new_matrix_parameter<double>(
        "extr_bdry_coords", "Define exterior boundary:", 
        "End exterior boundary", 5);

    mesh_reader.new_matrix_parameter<double>(
        "intr_bdry_coords", "Define interior boundary:", 
        "End interior boundary", 5);


    mesh_reader.new_vector_parameter<double>(
        "extr_bdry_rect", "Define exterior rectangular boundary:", 7);

    mesh_reader.new_vector_parameter<double>(
        "intr_bdry_rect", "Define interior rectangular boundary:", 7);

    
    mesh_reader.new_vector_parameter<double>(
        "extr_bdry_circ", "Define exterior circular boundary:", 7);

    mesh_reader.new_vector_parameter<double>(
        "intr_bdry_circ", "Define interior circular boundary:", 7);


    mesh_reader.new_vector_parameter<double>(
        "extr_bdry_square", "Define exterior squared boundary:", 6);

    mesh_reader.new_vector_parameter<double>(
        "intr_bdry_square", "Define interior squared boundary:", 6);


    mesh_reader.new_vector_parameter<double>(
        "extr_bdry_triangle", "Define exterior triangular boundary:", 6);

    mesh_reader.new_vector_parameter<double>(
        "intr_bdry_triangle", "Define interior triangular boundary:", 6);



  } // TQMeshApp::init_parameter_file_reader()


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  ParaReader                 reader_;

}; // TQMeshApp

} // namespace TQMesh
