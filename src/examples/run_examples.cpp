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
  bool success = true;

  if ( !example.compare("1") || !example.compare("01") || 
       !example.compare("simple_triangular_mesh") )
  {
    LOG(INFO) << "Running example \"simple_triangular_mesh\"...";
    LOG(INFO) << "";
    success &= simple_triangular_mesh();
  } 
  else if ( !example.compare("2") || !example.compare("02") ||
            !example.compare("square_in_channel") )
  {
    LOG(INFO) << "Running example \"square_in_channel\"...";
    LOG(INFO) << "";
    success &= square_in_channel();
  } 
  else if ( !example.compare("3") || !example.compare("03") ||
            !example.compare("boundary_shapes") )
  {
    LOG(INFO) << "Running example \"boundary_shapes\"...";
    LOG(INFO) << "";
    success &= boundary_shapes();
  }  
  else if ( !example.compare("4") || !example.compare("04") ||
            !example.compare("fixed_vertices") )
  {
    LOG(INFO) << "Running example \"fixed_vertices\"...";
    LOG(INFO) << "";
    success &= fixed_vertices();
  } 
  else if ( !example.compare("5") || !example.compare("05") ||
            !example.compare("merge_meshes") )
  {
    LOG(INFO) << "Running example \"merge_meshes\"...";
    LOG(INFO) << "";
    success &= merge_meshes();
  } 
  else if ( !example.compare("6") || !example.compare("06") ||
            !example.compare("airfoil_from_csv") )
  {
    LOG(INFO) << "Running example \"airfoil_from_csv\"...";
    LOG(INFO) << "";
    success &= airfoil_from_csv();
  } 
  else if ( !example.compare("7") || !example.compare("07") ||
            !example.compare("multiple_meshes") )
  {
    LOG(INFO) << "Running example \"multiple_meshes\"...";
    LOG(INFO) << "";
    success &= multiple_meshes();
  } 
  else if ( !example.compare("8") || !example.compare("08") ||
            !example.compare("thin_fracture") )
  {
    LOG(INFO) << "Running example \"thin_fracture\"...";
    LOG(INFO) << "";
    success &= thin_fracture();
  } 
  else if ( !example.compare("9") || !example.compare("09") ||
            !example.compare("fixed_edges") )
  {
    LOG(INFO) << "Running example \"fixed_edges\"...";
    LOG(INFO) << "";
    success &= fixed_edges();
  } 
  else if ( !example.compare("10") ||
            !example.compare("tqmesh_banner") )
  {
    LOG(INFO) << "Running example \"tqmesh_banner\"...";
    LOG(INFO) << "";
    success &= tqmesh_banner();
  } 
  else
  {
    LOG(INFO) << "";
    LOG(INFO, RED) << "No example \"" << example  << "\" found";
    LOG(INFO) << "";
    return EXIT_FAILURE;
  }

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);

} // run_examples()
