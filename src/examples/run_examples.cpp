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
  if ( !example.compare("1") )
  {
    LOG(INFO) << "Running example 1...";
    LOG(INFO) << "";
    run_example_1();
  } 
  else if ( !example.compare("2") )
  {
    LOG(INFO) << "Running example 2...";
    LOG(INFO) << "";
    run_example_2();
  } 
  else if ( !example.compare("3") )
  {
    LOG(INFO) << "Running example 3...";
    LOG(INFO) << "";
    run_example_3();
  }  
  else if ( !example.compare("4") )
  {
    LOG(INFO) << "Running example 4...";
    LOG(INFO) << "";
    run_example_4();
  } /*
  else if ( !example.compare("5") )
  {
    LOG(INFO) << "Running example 5...";
    LOG(INFO) << "";
    run_example_5();
  }
  else if ( !example.compare("6") )
  {
    LOG(INFO) << "Running example 6...";
    LOG(INFO) << "";
    run_example_6();
  }
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
