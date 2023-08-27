/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <string>

#include "utils.h"
#include "VecND.h"
#include "Domain.h"

#include "size_function.h"

#ifdef TQMESH_USE_EXPRTK
#include "exprtk.h"
#endif

/********************************************************************
* Static variables 
*******************************************************************/
#ifdef TQMESH_USE_EXPRTK
static double x_size_function, y_size_function;
static exprtk::symbol_table<double> symbol_table;
static exprtk::expression<double>   expression;
static exprtk::parser<double>       parser;
#endif

/********************************************************************
* Initialize the user defined size from a given input string
*******************************************************************/
UserSizeFunction init_size_function(const std::string& expr)
{
#ifdef TQMESH_USE_EXPRTK
  symbol_table.add_variable("x", x_size_function);
  symbol_table.add_variable("y", y_size_function);
  expression.register_symbol_table(symbol_table);
  parser.compile(expr, expression);

  UserSizeFunction f = [](const Vec2d& p) 
  { 
    x_size_function = p.x;
    y_size_function = p.y;
    return expression.value();  
  };

  return f;

#else

  double s = std::stod(expr);

  UserSizeFunction f = [s](const Vec2d& p) 
  { 
    return s;
  };

  return f;

#endif

} // init_size_function()



