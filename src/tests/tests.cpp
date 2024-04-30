/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include "tests.h"
#include "TQMesh.h"

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
int run_tests(const std::string& test_case)
{
  /*------------------------------------------------------------------
  | Print header
  ------------------------------------------------------------------*/
  LOG(INFO) << "";
  LOG(INFO) << "   -------------------------   ";
  LOG(INFO) << "   |  TQMesh - Test suite  |   ";
  LOG(INFO) << "   -------------------------   ";
  LOG(INFO) << "";

  /*------------------------------------------------------------------
  | Run all tests
  ------------------------------------------------------------------*/
  if ( !test_case.compare("Vertex") )
  {
    LOG(INFO) << "  Running tests for \"Vertex\" class...";
    run_tests_Vertex();
  }
  else if ( !test_case.compare("Triangle") )
  {
    LOG(INFO) << "  Running tests for \"Triangle\" class...";
    run_tests_Triangle();
  }
  else if ( !test_case.compare("Quad") )
  {
    LOG(INFO) << "  Running tests for \"Quad\" class...";
    run_tests_Quad();
  }
  else if ( !test_case.compare("Front") )
  {
    LOG(INFO) << "  Running tests for \"Front\" class...";
    run_tests_Front();
  }
  else if ( !test_case.compare("EdgeList") )
  {
    LOG(INFO) << "  Running tests for \"EdgeList\" class...";
    run_tests_EdgeList();
  }
  else if ( !test_case.compare("Boundary") )
  {
    LOG(INFO) << "  Running tests for \"Boundary\" class...";
    run_tests_Boundary();
  } 
  else if ( !test_case.compare("SizeFunction") )
  {
    LOG(INFO) << "  Running tests for \"SizeFunction\" class...";
    run_tests_SizeFunction();
  } 
  else if ( !test_case.compare("SmoothingStrategy") )
  {
    LOG(INFO) << "  Running tests for \"SmoothingStrategy\" class...";
    run_tests_SmoothingStrategy();
  }
  else if ( !test_case.compare("Mesh") )
  {
    LOG(INFO) << "  Running tests for \"Mesh\" class...";
    run_tests_Mesh();
  }
  else if ( !test_case.compare("MeshGenerator") )
  {
    LOG(INFO) << "  Running tests for \"MeshGenerator\" class...";
    run_tests_MeshGenerator();
  }
  else if ( !test_case.compare("MeshCleanup") )
  {
    LOG(INFO) << "  Running tests for \"MeshCleanup\" class...";
    run_tests_MeshCleanup();
  }
  else if ( !test_case.compare("ParaReader") )
  {
    LOG(INFO) << "  Running tests for \"ParaReader\" class...";
    run_tests_ParaReader();
  }
  else
  {
    LOG(INFO) << "";
    LOG(INFO, RED) << "  No test case \"" << test_case  << "\" found";
    LOG(INFO) << "";
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
      LOG(INFO,RED) << "[ERROR] Test (" << error_count 
              << "/" << total_tests << ") failed.";
      LOG(INFO) << "        --> " << data;
    }
    state &= data.state();
  }

  /*------------------------------------------------------------------
  | Succeess / fail
  ------------------------------------------------------------------*/
  LOG(INFO) << "";
  if (!state)
  {
    LOG(INFO, RED) << "  --> (" << error_count << "/" 
            << total_tests << ") tests failed." ;
  }
  else
  {
    LOG(INFO, GREEN) << "  --> (" << total_tests-error_count << "/" 
            << total_tests << ") tests succeeded.";
  }
  LOG(INFO) << "";
  LOG(INFO) << "";

  if (!state)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

} // run_tests()

/*********************************************************************
* Adjust logging output stream
* -> Log debug messages to specified output-file
*********************************************************************/
void adjust_logging_output_stream(const std::string& file)
{
  if (file == "COUT")
  {
    CppUtils::LOG_PROPERTIES.set_info_ostream( CppUtils::TO_COUT );
    CppUtils::LOG_PROPERTIES.set_debug_ostream( CppUtils::TO_COUT );
    return;
  }

  std::string source_directory { TQMESH_SOURCE_DIR };
  std::string filepath { source_directory + "/auxiliary/test_data/" + file };
  CppUtils::LOG_PROPERTIES.set_info_ostream( CppUtils::TO_FILE, filepath );
  CppUtils::LOG_PROPERTIES.set_debug_ostream( CppUtils::TO_FILE, filepath );

} // adjust_logging_output_stream()

