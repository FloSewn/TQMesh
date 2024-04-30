/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "TQMesh.h"

namespace TQMesh {

using namespace CppUtils;

/********************************************************************* 
* Simple class for error handling
*********************************************************************/
class Error : public std::exception
{
public:
  Error(const std::string& msg) : message_ {msg} {}
  const char* what() const noexcept override { return message_.c_str(); }
private:
  std::string message_;
};

static inline void throw_error(std::string msg) { throw Error { msg }; }

} // namespace TQMesh
