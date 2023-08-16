/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <iostream>
#include <cstdlib>
#include <sstream>

#include <TQMeshConfig.h>

#include "Log.h"

#include "TQMeshApp.h"


using CppUtils::LOG_PROPERTIES;
using CppUtils::LOG;
using CppUtils::LogLevel::INFO;
using CppUtils::LogLevel::DEBUG;

/*********************************************************************
* Print TQMesh header
*********************************************************************/
static void print_header()
{
  LOG(INFO) << "";
  LOG(INFO) << "   - - - - - - - - - - - - - - - - - - -  ";
  LOG(INFO) << "  | TQMesh - A simple 2D mesh generator | ";
  LOG(INFO) << "   - - - - - - - - - - - - - - - - - - -  ";
  LOG(INFO) << "";
  LOG(INFO) << "  by Florian Setzwein";
  LOG(INFO) << "";
  LOG(INFO) << "  Version " 
            << TQMESH_VERSION_MAJOR << "." << TQMESH_VERSION_MINOR ;
  LOG(INFO) << "";

} // print_header()


using namespace CppUtils;
using namespace TQMesh;

int main(int argc, char* argv[])
{
  LOG_PROPERTIES.set_level( INFO );
  LOG_PROPERTIES.show_header( true );
  LOG_PROPERTIES.use_color( true );
  LOG_PROPERTIES.set_info_header( "  " );
  LOG_PROPERTIES.set_debug_header( "# " );
  LOG_PROPERTIES.set_error_header( "  [ERROR] " );
  LOG_PROPERTIES.set_warn_header( "  [WARNING] " );

  /*------------------------------------------------------------------
  | Handle command line arguments
  ------------------------------------------------------------------*/
  print_header();
  if ( argc < 2 )
  {
    LOG(INFO) << "Usage: ";
    LOG(INFO) <<  argv[0] << " \"Input-File\" > \"Output-Mesh\"" ;
    LOG(INFO) << "";
    return EXIT_SUCCESS;
  }

  try
  {
    TQMeshApp app { argv[1] };

    if ( !app.run() )
      return EXIT_FAILURE;
  }
  catch (const std::exception& e)
  {
    LOG(ERROR) << e.what();
    return EXIT_FAILURE;
  }


  return EXIT_SUCCESS;
}

