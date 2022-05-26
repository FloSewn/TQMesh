

  auto para_size_fun = reader.new_string_parameter(
      "Element size:");
  auto para_meshing_algorithm = reader.new_string_parameter(
      "Algorithm:");
  auto para_quad_refinements = reader.new_int_parameter(
      "Quad-refinements:");
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
  | Initial query of some mandatory parameters
  ------------------------------------------------------------------*/
  reader.query( para_size_fun );
  reader.query( para_vertices );

  if ( !para_size_fun.found() )
  {
    MSG("[ERROR] No mesh size definition provided.");
    return EXIT_SUCCESS;
  }

  if ( !para_vertices.found() )
  {
    MSG("[ERROR] Invalid vertex definition.");
    return EXIT_SUCCESS;
  }


  /*------------------------------------------------------------------
  | Initialize element size function
  ------------------------------------------------------------------*/
  UserSizeFunction size_fun;
  size_fun = init_size_function( para_size_fun.value() );

  /*------------------------------------------------------------------
  | Query boundary vertices and obtain domain extent
  ------------------------------------------------------------------*/
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
    Boundary& b_ext = domain.add_exterior_boundary();

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
    Boundary& b_int = domain.add_interior_boundary();

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
    Boundary& b_int = domain.add_interior_boundary();

    int    m = static_cast<int>( para_intr_rect_bdry.value(0) );
    double x = para_intr_rect_bdry.value(1);
    double y = para_intr_rect_bdry.value(2);
    double w = para_intr_rect_bdry.value(3);
    double h = para_intr_rect_bdry.value(4);

    b_int.set_shape_rectangle( vertices, m, {x,y}, w, h );
  }

  while( reader.query( para_intr_circ_bdry ) )
  {
    Boundary& b_int = domain.add_interior_boundary();

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
      try 
      {
        const double x = para_fixed_vertices.value(4*i);
        const double y = para_fixed_vertices.value(4*i + 1);
        const double s = para_fixed_vertices.value(4*i + 2);
        const double r = para_fixed_vertices.value(4*i + 3);

        domain.add_fixed_vertex(x,y,s,r);
      }
      catch (const std::exception& e)
      {
        MSG("WARNING: Invalid fixed vertex input format.");
      }
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

  try
  {
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
  }
  catch( ... )
  {
    MSG("WARNING: Invalid quad layer input format.");
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
  std::string meshing_type = "Triangulation";

  reader.query( para_meshing_algorithm );

  if ( para_meshing_algorithm.found() )
    meshing_type = para_meshing_algorithm.value();

  if ( meshing_type == "Paving" )
  {
    mesh.pave();
  }
  else if ( meshing_type == "Tri-to-Quad")
  {
    mesh.triangulate();
    mesh.merge_triangles_to_quads();
  }
  else
  {
    mesh.triangulate();
  }

  /*------------------------------------------------------------------
  | Mesh refinements
  ------------------------------------------------------------------*/
  reader.query( para_quad_refinements );

  if ( para_quad_refinements.found() )
  {
    int n_refinements = para_quad_refinements.value();

    for ( int i = 0; i < n_refinements; ++i )
      mesh.refine_to_quads();
  }
  
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
