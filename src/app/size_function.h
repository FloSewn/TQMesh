/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "utils.h"
#include "VecND.h"
#include "Domain.h"

using namespace CppUtils;
using namespace TQMesh::TQAlgorithm;

/********************************************************************
* Initialize the user defined size from a given input string
*******************************************************************/
UserSizeFunction init_size_function(const std::string& expr);
