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

#include "tests.h"
#include "Helpers.h"
#include "Testing.h"

/*********************************************************************
* Color text
*********************************************************************/
#define NC "\e[0m"
#define RED "\e[0;31m"
#define GRN "\e[0;32m"
#define CYN "\e[0;36m"
#define REDB "\e[41m"

/*********************************************************************
* The main test function
*********************************************************************/
int run_tests(const std::string& test_case)
{
  CppUtils::SimpleLogger TESTMSG(std::clog, "  ");

  /*------------------------------------------------------------------
  | Print header
  ------------------------------------------------------------------*/
  TESTMSG << "" << std::endl;
  TESTMSG << "   -------------------------   " << std::endl;
  TESTMSG << "   |  TQMesh - Test suite  |   " << std::endl;
  TESTMSG << "   -------------------------   " << std::endl;
  TESTMSG << "" << std::endl;

  /*------------------------------------------------------------------
  | Run all tests
  ------------------------------------------------------------------*/
  if ( !test_case.compare("Vertex") )
  {
    TESTMSG << "  Running tests for \"Vertex\" class..."
            << std::endl;
    run_tests_Vertex();
  }
  else if ( !test_case.compare("Triangle") )
  {
    TESTMSG << "  Running tests for \"Triangle\" class..."
            << std::endl;
    run_tests_Triangle();
  }
  else if ( !test_case.compare("Quad") )
  {
    TESTMSG << "  Running tests for \"Quad\" class..."
            << std::endl;
    run_tests_Quad();
  }
  else if ( !test_case.compare("Front") )
  {
    TESTMSG << "  Running tests for \"Front\" class..."
            << std::endl;
    run_tests_Front();
  }/*
  else if ( !test_case.compare("EdgeList") )
  {
    TESTMSG << "  Running tests for \"EdgeList\" class..."
            << std::endl;
    run_tests_EdgeList();
  }
  else if ( !test_case.compare("Boundary") )
  {
    TESTMSG << "  Running tests for \"Boundary\" class..."
            << std::endl;
    run_tests_Boundary();
  }
  else if ( !test_case.compare("SizeFunction") )
  {
    TESTMSG << "  Running tests for \"SizeFunction\" class..."
            << std::endl;
    run_tests_SizeFunction();
  }
  else if ( !test_case.compare("Smoother") )
  {
    TESTMSG << "  Running tests for \"Smoother\" class..."
            << std::endl;
    run_tests_Smoother();
  }
  else if ( !test_case.compare("Mesh") )
  {
    TESTMSG << "  Running tests for \"Mesh\" class..."
            << std::endl;
    run_tests_Mesh();
  }*/
  else
  {
    TESTMSG << std::endl;
    TESTMSG << RED "  No test case \"" << test_case 
            << "\" found" NC << std::endl;
    TESTMSG << std::endl;
    return EXIT_FAILURE;
  }

  /*------------------------------------------------------------------
  | Check for failed tests
  ------------------------------------------------------------------*/
  std::vector<CppUtils::TestData>& test_data 
    = CppUtils::TestDataSingleton::instance();

  bool   state = true;
  size_t error_count = 0;
  size_t total_tests = test_data.size();

  for (auto data : test_data )
  {
    if ( !data.state() )
    {
      ++error_count;
      TESTMSG << RED "[ERROR] Test (" << error_count 
              << "/" << total_tests << ") failed." NC << std::endl;
      TESTMSG << "        --> " << data << std::endl;
    }
    state &= data.state();
  }

  /*------------------------------------------------------------------
  | Succeess / fail
  ------------------------------------------------------------------*/
  TESTMSG << "" << std::endl;
  if (!state)
  {
    TESTMSG << RED "  --> (" << error_count << "/" 
            << total_tests << ") tests failed." NC  << std::endl;
  }
  else
  {
    TESTMSG << GRN "  --> (" << total_tests-error_count << "/" 
            << total_tests << ") tests succeeded." NC << std::endl;
  }
  TESTMSG << std::endl;
  TESTMSG << std::endl;

  if (!state)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

} // run_tests()

