/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <TQMeshConfig.h>
#include "STLHeaders.h"
#include "CppUtils.h"

#include "tests.h"

using CppUtils::LOG_PROPERTIES;
using CppUtils::LOG;
using CppUtils::LogLevel::INFO;
using CppUtils::LogLevel::DEBUG;

/*********************************************************************
* The main function
*********************************************************************/
int main(int argc, char* argv[])
{
  LOG_PROPERTIES.set_level( DEBUG );
  LOG_PROPERTIES.show_header( true );
  LOG_PROPERTIES.use_color( true );
  LOG_PROPERTIES.set_info_header( "  " );
  LOG_PROPERTIES.set_debug_header( "# " );

  if ( argc < 2 )
  {
    LOG(INFO) << "";
    LOG(INFO) << "   -------------------------   ";
    LOG(INFO) << "   |  TQMesh - Test suite  |   ";
    LOG(INFO) << "   -------------------------   ";
    LOG(INFO) << "";
    LOG(INFO) << "Usage: " << argv[0] << " <Test-Case>";
    LOG(INFO) << "";
    LOG(INFO) << "";
    return EXIT_FAILURE;
  }

  std::string input { argv[1] };

  return run_tests( input );

  return EXIT_SUCCESS;
}
