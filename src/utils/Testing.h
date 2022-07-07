/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <vector>
#include <string>

namespace CppUtils {

/*********************************************************************
* This class provides as simple container to track the line,
* the function and the file, where a test is conducted.
*********************************************************************/
class TestData
{
public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  TestData(int line, const std::string& func, const std::string& file)
  {
    line_ = line;
    func_ = func;
    file_ = file;
  }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  void state(bool s) { state_ = s; }

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  bool state() const { return state_; }
  int  line() const { return line_; }
  const std::string& func() const { return func_; }
  const std::string& file() const { return file_; }

private:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  bool        state_ { true };
  int         line_  { 0 };
  std::string func_  {""};
  std::string file_  {""};

}; // TestData

/*********************************************************************
* TestData ostream overload 
*********************************************************************/
static std::ostream& operator<<(std::ostream& os, 
                                const TestData& t)
{
  return os << "Line " << t.line() << " in "
            << "function " << t.func() << "() in "
            << "file " << t.file();
}

/*********************************************************************
* This is where the test data will be stored
*********************************************************************/
class TestDataSingleton
{
private:
  TestDataSingleton();

public:
  static std::vector<TestData>& instance()
  {
    static std::vector<TestData> INSTANCE;
    return INSTANCE;
  }
};

/*********************************************************************
* The standard function to check for tests
*********************************************************************/
static inline bool eval_test(bool cond, 
                             int line,
                             const std::string& func,
                             const std::string& file)
{
  TestData sample { line, func, file };
  
  if (cond) 
    sample.state( true );
  else
    sample.state( false );

  TestDataSingleton::instance().push_back( sample );

  return cond;
}

#define CHECK(cond)                                        \
  do { eval_test( (cond), __LINE__, __func__, __FILE__);   \
  } while(false)


} // namespace CppUtils
