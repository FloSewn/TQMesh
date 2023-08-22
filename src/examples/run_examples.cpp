/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <iostream>
#include <vector>
#include <string>

#include "run_examples.h"
#include "Log.h"

/*********************************************************************
* Log utils
*********************************************************************/
using CppUtils::LOG;
using CppUtils::LogLevel::INFO;
using CppUtils::LogColor::GREEN;
using CppUtils::LogColor::RED;

/*********************************************************************
* The main test function
*********************************************************************/
int run_examples(const std::string& example)
{
  /*------------------------------------------------------------------
  | Print header
  ------------------------------------------------------------------*/
  LOG(INFO) << "";
  LOG(INFO) << "   -------------------------   ";
  LOG(INFO) << "   |   TQMesh - Examples   |   ";
  LOG(INFO) << "   -------------------------   ";
  LOG(INFO) << "";

  /*------------------------------------------------------------------
  | Run all examples
  ------------------------------------------------------------------*/
  if ( !example.compare("1") || 
       !example.compare("simple_triangular_mesh") )
  {
    LOG(INFO) << "Running example \"simple_triangular_mesh\"...";
    LOG(INFO) << "";
    simple_triangular_mesh();
  } 
  else if ( !example.compare("2") ||
            !example.compare("square_in_channel") )
  {
    LOG(INFO) << "Running example \"square_in_channel\"...";
    LOG(INFO) << "";
    square_in_channel();
  } 
  else if ( !example.compare("3") ||
            !example.compare("boundary_shapes") )
  {
    LOG(INFO) << "Running example \"boundary_shapes\"...";
    LOG(INFO) << "";
    boundary_shapes();
  }  
  else if ( !example.compare("4") ||
            !example.compare("fixed_vertices") )
  {
    LOG(INFO) << "Running example \"fixed_vertices\"...";
    LOG(INFO) << "";
    fixed_vertices();
  } 
  else if ( !example.compare("5") ||
            !example.compare("merge_meshes") )
  {
    LOG(INFO) << "Running example \"merge_meshes\"...";
    LOG(INFO) << "";
    merge_meshes();
  } 
  else if ( !example.compare("6") ||
            !example.compare("tqmesh_banner") )
  {
    LOG(INFO) << "Running example \"tqmesh_banner\"...";
    LOG(INFO) << "";
    tqmesh_banner();
  } /*
  else if ( !example.compare("7") )
  {
    LOG(INFO) << "Running example 7...";
    LOG(INFO) << "";
    run_example_7();
  }
  else if ( !example.compare("8") )
  {
    LOG(INFO) << "Running example 8...";
    LOG(INFO) << "";
    run_example_8();
  }
  else if ( !example.compare("9") )
  {
    LOG(INFO) << "Running example 9...";
    LOG(INFO) << "";
    run_example_9();
  }*/
  else
  {
    LOG(INFO) << "";
    LOG(INFO, RED) << "No example \"" << example  << "\" found";
    LOG(INFO) << "";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

} // run_examples()
