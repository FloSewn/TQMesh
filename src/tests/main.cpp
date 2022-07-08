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

#include "Helpers.h"

#include "tests.h"


/*********************************************************************
* The main function
*********************************************************************/
int main(int argc, char* argv[])
{
  CppUtils::SimpleLogger TESTMSG(std::clog, "  ");

  if ( argc < 2 )
  {
    TESTMSG << "" << std::endl;
    TESTMSG << "   -------------------------   " << std::endl;
    TESTMSG << "   |  TQMesh - Test suite  |   " << std::endl;
    TESTMSG << "   -------------------------   " << std::endl;
    TESTMSG << "" << std::endl;
    TESTMSG << "Usage: " << argv[0] << " <Test-Case>" << std::endl;
    TESTMSG << "" << std::endl;
    TESTMSG << "" << std::endl;
    return EXIT_FAILURE;
  }

  std::string input { argv[1] };

  return run_tests( input );

  return EXIT_SUCCESS;
}
