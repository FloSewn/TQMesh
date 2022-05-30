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

namespace CppUtils {

/*********************************************************************
* Asserts messages
*********************************************************************/
#ifndef NDEBUG
static inline void ASSERT(bool cond, const std::string& msg)
{
  if (!cond)
  {
    std::cerr << "[ERROR] " << msg << std::endl;
    assert(cond);
  }
}
#else
static inline void ASSERT(bool cond, const std::string& msg)
{}
#endif

/*********************************************************************
* A simple logging device. 
* It gives the possibility to choose the output, as well as a 
* a logging message at the start. 
* Must be used with std::endl - otherwise newlines are confusing
*
* Source: 
* -------
*   https://stackoverflow.com/questions/2212776/\
*   overload-handling-of-stdendl
*********************************************************************/
class SimpleLogger: public std::ostream
{

public:
  SimpleLogger(std::ostream& str, std::string key)
  : std::ostream(&buffer_)
  , buffer_(str, std::move(key))
  {}

private:

  /*------------------------------------------------------------------
  | A streambuffer that prefixes every line with a key
  ------------------------------------------------------------------*/
  class LogBuf: public std::stringbuf
  {
  public: 
    // Constructor
    LogBuf(std::ostream& str, std::string key) 
    : output_{str}
    , key_ {std::move(key)}
    {}

    // Destructor
    ~LogBuf() 
    {
      if (pbase() != pptr()) 
      {
        putOutput();
      }
    }

    // Syncronize stream with output
    virtual int sync() 
    {
      putOutput();
      return 0;
    }

    // This is were things happen
    void putOutput() 
    {
      output_ << key_ << str();
      str("");
      output_.flush();
    }
    
  private:
    std::ostream& output_;
    std::string key_;
  };

  LogBuf buffer_;

}; // SimpleLogger


} // namespace CppUtils
