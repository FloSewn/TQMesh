/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include "utils.h"
#include "Vec2.h"
#include "Domain.h"

#include "size_function.h"

#include "exprtk.h"

/********************************************************************
* Static variables 
*******************************************************************/
static double x_size_function, y_size_function;
static exprtk::symbol_table<double> symbol_table;
static exprtk::expression<double>   expression;
static exprtk::parser<double>       parser;

/********************************************************************
* Initialize the user defined size from a given input string
*******************************************************************/
UserSizeFunction init_size_function(const std::string& expr)
{
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

} // init_size_function()



