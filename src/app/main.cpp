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

#include "utils.h"

#include "TQMeshApp.h"

/*********************************************************************
* Print TQMesh header
*********************************************************************/
static void print_header()
{
  MSG("");
  MSG("   - - - - - - - - - - - - - - - - - - -  ");
  MSG("  | TQMesh - A simple 2D mesh generator | ");
  MSG("   - - - - - - - - - - - - - - - - - - -  ");
  MSG("");
  MSG("  by Florian Setzwein");
  MSG("");
  MSG("  Version " << TQMESH_VERSION_MAJOR << "." << TQMESH_VERSION_MINOR );
  MSG("");

} // print_header()


using namespace CppUtils;
using namespace TQMesh;

int main(int argc, char* argv[])
{
  /*------------------------------------------------------------------
  | Handle command line arguments
  ------------------------------------------------------------------*/
  print_header();
  if ( argc < 2 )
  {
    MSG("Usage: " );
    MSG( argv[0] << " \"Input-File\" > \"Output-Mesh\"" );
    MSG("");
    return EXIT_SUCCESS;
  }

  TQMeshApp app { argv[1] };

  if ( !app.run() )
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}

