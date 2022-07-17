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
#include <cstdlib>

#include "Log.h"

#include "run_examples.h"

using CppUtils::LOG_PROPERTIES;
using CppUtils::LOG;
using CppUtils::LogLevel::INFO;

/*********************************************************************
* The main function
*********************************************************************/
int main(int argc, char* argv[])
{
  LOG_PROPERTIES.set_level( INFO );
  LOG_PROPERTIES.set_info_header( "  " );
  LOG_PROPERTIES.set_debug_header( "# " );

  if ( argc < 2 )
  {
    LOG(INFO) << "";
    LOG(INFO) << "   -------------------------   ";
    LOG(INFO) << "   |   TQMesh - Examples   |   ";
    LOG(INFO) << "   -------------------------   ";
    LOG(INFO) << "";
    LOG(INFO) << "Usage: " << argv[0] << " <Example>";
    LOG(INFO) << "";
    LOG(INFO) << "";
    return EXIT_FAILURE;
  }

  std::string input { argv[1] };

  return run_examples( input );

  return EXIT_SUCCESS;
}
