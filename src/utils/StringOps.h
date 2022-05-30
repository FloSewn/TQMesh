/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once
#include <iostream>
#include <sstream>
#include <vector>

namespace CppUtils {

/*********************************************************************
* A function to split a string by a delimiter.
* Put the results in a pre-constructed vector.
*
* Source: 
*   https://stackoverflow.com/questions/236129/how-do-i-iterate-\
*   over-the-words-of-a-string
*********************************************************************/
template<typename Out>
void split(const std::string& s, char delim, Out result, 
           bool remove_empty_elements=true)
{
  std::istringstream iss(s);
  std::string item;

  if (remove_empty_elements)
  {
    while (std::getline(iss, item, delim))
      if (!item.empty())
        *result++ = item;
  }
  else
  {
    while (std::getline(iss, item, delim))
      *result++ = item;
  }
}

/*********************************************************************
* A function to split a string by a delimiter.
* Put the results in a new vector.
*
* Source: 
*   https://stackoverflow.com/questions/236129/how-do-i-iterate-\
*   over-the-words-of-a-string
*********************************************************************/
std::vector<std::string> split(const std::string& s, char delim,
                               bool remove_empty_elements=true)
{
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems), remove_empty_elements);
  return elems;
}


/*********************************************************************
* Simple conversion of strings to integers or floats, depending 
* on the type
*
* Source:
*   https://stackoverflow.com/questions/60301788/template-friendly-\
*   string-to-numeric-in-c
*********************************************************************/
struct converter
{
  const std::string& x;
  template <typename T> operator T() { return 0; }
};

template <> converter::operator int() { return std::stoi(x); }
template <> converter::operator float() { return std::stof(x); }
template <> converter::operator double() { return std::stod(x); }
converter sto(const std::string& x) { return {x}; }


} // namespace CppUtils
