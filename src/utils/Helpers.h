/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <stdlib.h>

namespace CppUtils {

/*********************************************************************
* Asserts messages
*********************************************************************/
#ifndef NDEBUG
static inline void assert_msg(bool cond, const std::string& msg, 
                              const std::string& filename, int line)
{
  if (!cond)
  {
    std::cerr << "[ERROR] " << msg << " (" << filename 
              << " - Line " << line << ")" << std::endl;
    assert(cond);
  }
}
#else
static inline void assert_msg(bool cond, const std::string& msg,
                              const std::string& filename, int line)
{}
#endif

#define ASSERT(cond, msg) assert_msg(cond, msg, __FILE__, __LINE__)

/*********************************************************************
* Terminate 
*********************************************************************/
static inline void TERMINATE(const std::string& msg)
{
  std::cerr << "[ERROR] " << msg << std::endl;
  exit(EXIT_FAILURE);
}

} // namespace CppUtils
