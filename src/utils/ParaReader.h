/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <memory>
#include <map>


namespace CppUtils {

/*--------------------------------------------------------------------
| Some local typedefs
--------------------------------------------------------------------*/
using string        = std::string;
using ofstream      = std::ofstream;
using ifstream      = std::ifstream;
using strVec        = std::vector<string>;
using istringstream = std::istringstream;

/*--------------------------------------------------------------------
| The default parameter parent class 
--------------------------------------------------------------------*/
class ParameterBase
{
public:
  ParameterBase(const string& key)
  : single_key_ { key }
  , start_key_ { "" }
  , end_key_ { "" }
  , line_start_ { 0 }
  , found_ { false } 
  , multi_lines_ { false }
  {}

  ParameterBase(const string& start, const string& end)
  : single_key_ { "" }
  , start_key_ { start }
  , end_key_ { end }
  , line_start_ { 0 }
  , found_ { false } 
  , multi_lines_ { true }
  {}

  virtual ~ParameterBase() = default;

  const string& single_key() const { return single_key_;   }
  string& single_key() { return single_key_;   }

  const string& start_key() const { return start_key_; }
  string& start_key() { return start_key_; }

  const string& end_key() const { return end_key_; }
  string& end_key() { return end_key_; }

  int line_to_start() const { return line_start_; }
  void line_to_start(int l) { line_start_ = l; }

  bool found() const { return found_; }
  void found(bool f) { found_ = f; }

  bool multi_line_definition() const { return multi_lines_; }

private:
  string single_key_;
  string start_key_;
  string end_key_;

  size_t line_start_;
  bool   found_;
  bool   multi_lines_ {false};
}; 

/*--------------------------------------------------------------------
| The derived template class for a specified parameter
--------------------------------------------------------------------*/
template <class T>
class Parameter : public ParameterBase
{
public:
  /*------------------------------------------------------------------
  | This constructor is used to create parameters that are defined
  | in a single line.
  | These parameters are queried with a given <key>, which is given
  | via the constructor. 
  | Additionally, the number <n> of parameters that must be given 
  | behind the query key, is set via the constructor.
  ------------------------------------------------------------------*/
  Parameter(const string& key, size_t n) 
    : ParameterBase(key)
    , ncol_{n}
    , nrow_{1}
    , values_(n,T{})
  {}

  /*------------------------------------------------------------------
  | This constructor is used to create parameter that are defined
  | over multiple lines
  | These parameters are queried with a given <start> and <end>
  | key, which is given via the constructor. 
  | Additionally, the number <n> of parameters that must be given 
  | in every line between the start and ending query keas 
  | are set via the constructor.
  ------------------------------------------------------------------*/
  Parameter(const string& start, const string& end, size_t n) 
    : ParameterBase(start, end)
    , ncol_{n}
    , nrow_{1}
    , values_(n,T{})
  {}

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  T get_value(size_t i) 
  { 
    if ( i >= ncol_*nrow_ )
      return values_[0];
    return values_[i]; 
  }

  T get_value(size_t i, size_t j) 
  { 
    size_t index = j * ncol_ + i;
    if ( index >= ncol_*nrow_ )
      return values_[0];
    return values_[index]; 
  }

  size_t columns() const { return ncol_; }
  size_t rows() const { return nrow_; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  void set_value(size_t i, T val) 
  { 
    if ( i >= ncol_*nrow_ ) 
      return;

    values_[i] = val; 
  }

  void set_value(size_t i, size_t j, T val) 
  { 
    if ( i >= ncol_ || j >= nrow_) 
      return;

    size_t index = i + j * ncol_;
    values_[index] = val; 
  }

  void add_row()
  {
    ++nrow_;
    for( size_t i = 0; i < ncol_; i++ )
      values_.push_back( T{} );
  }

private:
  // The dataset
  std::vector<T> values_;

  // Dimension of dataset
  size_t ncol_;
  size_t nrow_;

}; // Parameter



/*--------------------------------------------------------------------
| The parameter file reader 
--------------------------------------------------------------------*/
class ParaReader
{
  using ParameterList = std::vector<std::unique_ptr<ParameterBase>>; 
  using ParameterMap  = std::map<std::string, size_t>;

public:
  /*------------------------------------------------------------------
  | Class for error handling
  ------------------------------------------------------------------*/
  class Invalid
  {
  public:
    Invalid(const string& msg){ error_message = msg; }
    string& what() { return error_message; }
  private:
    string error_message;
  };

  /*------------------------------------------------------------------
  | Function for error handling
  ------------------------------------------------------------------*/
  void error(string msg) { throw Invalid{msg}; }

  /*------------------------------------------------------------------
  | Constructor
  | Read parameter file and store it in the file_content_ buffer, 
  | such that each line of the file is stored as a separate
  | string in it.
  ------------------------------------------------------------------*/
  ParaReader(const string& file_path) 
  {
    file_path_   = file_path;
    comment_     = "#";
    delimiter_   = ',';

    ifstream ifs {file_path_};
    if (!ifs) error("Can't open file: " + file_path);

    // General throw for bad reading of file
    ifs.exceptions(ifs.exceptions() | std::ios_base::badbit);

    // Read the input file
    string line;
    while (std::getline(ifs, line))
    {
      // Ignore commented lines
      string l_cpy = line;
      l_cpy.erase(std::remove_if(l_cpy.begin(), l_cpy.end(), isspace), 
                  l_cpy.end());
      if ( l_cpy[0] == '#' )
        continue;
      
      // Add line to content, if it not a comment line
      file_content_.push_back(line);
    }

    if (!ifs.eof()) 
      error("Failed to read entire file - "
            "Stopped reading in line " + 
            std::to_string(file_content_.size()));

  } // ParaReader::Constructor() 

  /*------------------------------------------------------------------
  | Create new scalar parameters to search for in a file
  | 
  | Arguments:
  | ----------
  |   name: The parameter name, under which it can be found in the
  |         ParaReader structure
  |   key: The parameter key, which is used to query the parameter
  |        in the input file
  ------------------------------------------------------------------*/
  template <typename T>
  void new_scalar_parameter(const string& name, const string& key)
  {
    // Check that parameter map does not yet contain given name
    if ( param_map_.count(name) )
      error("Multiple definitions for parmeter name \"" + name + "\"."); 

    // Empty parameter keys are not allowed
    if ( key.size() < 1 )
      error("Empty parameter queries are not allowed. Parameter: " + name);

    // Update parameter list
    param_list_.push_back( 
        std::make_unique<Parameter<T>>( key, 1 ) 
    );

    // Update parameter map
    param_map_[name] = param_list_.size() - 1;
  }

  /*------------------------------------------------------------------
  | Create new list parameters to search for in a file
  |
  | Arguments:
  | ----------
  |   name: The parameter name, under which it can be found in the
  |         ParaReader structure
  |   key: The parameter key, which is used to query the parameter
  |        in the input file
  |   n: The number of parameters that must be provided
  ------------------------------------------------------------------*/
  template <typename T>
  void new_list_parameter(const string& name, 
                          const string& key, size_t n)
  {
    // Check that parameter map does not yet contain given name
    if ( param_map_.count(name) )
      error("Multiple definitions for parmeter name \"" + name + "\"."); 

    // Empty parameter keys are not allowed
    if ( key.size() < 1 )
      error("Empty parameter queries are not allowed. Parameter: " + name);

    // Update parameter list
    param_list_.push_back( 
        std::make_unique<Parameter<T>>( key, n ) 
    );

    // Update parameter map
    param_map_[name] = param_list_.size() - 1;
  }

  /*------------------------------------------------------------------
  | Create new list parameters to search for in a file
  |
  | Arguments:
  | ----------
  |   name: The parameter name, under which it can be found in the
  |         ParaReader structure
  |   start: The parameter start key, which is used to query 
  |          the beginning of the parameter data in the input file
  |   end: The parameter end key, which is used to query 
  |          the ending of the parameter data in the input file
  |   n: The number of parameters that must be provided in each line
  ------------------------------------------------------------------*/
  template <typename T>
  void new_list_parameter(const string& name,
                          const string& start, const string& end, 
                          size_t n)
  {
    // Check that parameter map does not yet contain given name
    if ( param_map_.count(name) )
      error("Multiple definitions for parmeter name \"" + name + "\"."); 

    // Empty parameter keys are not allowed
    if ( start.size() < 1 || end.size() < 1 )
      error("Empty parameter queries are not allowed. Parameter: " + name);

    // Update parameter list
    param_list_.push_back( 
        std::make_unique<Parameter<T>>( start, end, n ) 
    );

    // Update parameter map
    param_map_[name] = param_list_.size() - 1;
  }

  /*------------------------------------------------------------------
  | Get a parameter reference from a given parameter name
  ------------------------------------------------------------------*/
  template <typename T>
  Parameter<T>& get_parameter(const string& name)
  {
    // Handle unknown parameter names
    if ( !param_map_.count(name) )
      error("No parameter with name \"" + name + "\" has been defined.");

    // Cast parameter 
    Parameter<T>* param_ptr 
      = dynamic_cast<Parameter<T>*>( 
          param_list_[ param_map_[name] ].get() 
    );

    Parameter<T>& param = *param_ptr;

    return param;

  } // ParaReader::get_parameter()


  /*------------------------------------------------------------------
  | Return the first value of a parameter with <name>.
  ------------------------------------------------------------------*/
  template <typename T>
  T get_value(const string& name)
  {
    Parameter<T>& param = get_parameter<T>( name );

    // Check if parameter is not found in the input file
    if ( !param.found() )
      error("Parameter with name \"" + name + "\" not found in input file.");
    
    return param.get_value(0);
  }

  /*------------------------------------------------------------------
  | Return the value of a parameter with <name>.
  ------------------------------------------------------------------*/
  template <typename T>
  T get_value(size_t i, const string& name)
  {
    Parameter<T>& param = get_parameter<T>( name );

    // Check if parameter is not found in the input file
    if ( !param.found() )
      error("Parameter with name \"" + name + "\" not found in input file.");
    
    return param.get_value(i);
  }

  /*------------------------------------------------------------------
  | Return the value of a parameter with <name> and return it.
  ------------------------------------------------------------------*/
  template <typename T>
  T get_value(size_t i, size_t j, const string& name)
  {
    Parameter<T>& param = get_parameter<T>( name );

    // Check if parameter is not found in the input file
    if ( !param.found() )
      error("Parameter with name \"" + name + "\" not found in input file.");
    
    return param.get_value(i,j);
  }

  /*------------------------------------------------------------------
  | Check if a parameter with a given name has been found in the file
  ------------------------------------------------------------------*/
  bool found(const string& name)
  {
    // Handle unknown parameter names
    if ( !param_map_.count(name) )
      error("No parameter with name \"" + name + "\" has been defined.");

    return param_list_[ param_map_[name] ].get()->found();
  }

  /*------------------------------------------------------------------
  | Query  parameters
  | -> Scan the file for the query key of a given parameter with 
  |    <name>. If query key is found, read the value and store it in
  |    the corresponding parameter.
  ------------------------------------------------------------------*/
  template <typename T>
  bool query(const string& name)
  {
    Parameter<T>& param = get_parameter<T>( name );

    // Scalar parameters
    if ( param.columns() == 1 && !param.multi_line_definition() )
      return query_scalar(param);
    // Single-line & multiple parameters
    else if ( param.columns() > 1 && !param.multi_line_definition() ) 
      return query_single_line(param);
    // Else: Multi-line parameters
    return query_multiple_lines(param);

  } // ParaReader::query()






private:

  /*------------------------------------------------------------------
  | Query scalar parameters
  | --> Parameter existence has been checked in calling function
  ------------------------------------------------------------------*/
  template <typename T>
  bool query_scalar(Parameter<T>& param)
  {
    string query = param.single_key();

    // Parameter query reached end of input file
    if ( param.line_to_start() >= file_content_.size() )
      return false;

    // Position of query in the line
    size_t line_pos = string::npos;

    // Index of line where query is found
    size_t line_index = string::npos;

    size_t index;

    // Search for query in file buffer line by line
    for (size_t cur_line = param.line_to_start(); 
         cur_line < file_content_.size(); 
         ++cur_line)
    {
      const string& line = file_content_[cur_line];

      // Search in all lines for query
      // Use last query definition that occures in a line
      size_t pos = 0;

      while ((index = line.find(query, pos)) != string::npos)
      {
        // Position is from next element of index
        line_pos   = index;
        line_index = cur_line;
        pos        = index + 1;

        // Set starting line for next query search
        param.line_to_start( cur_line + 1 );
      }

      // Stop if query has been found in current line
      if ( line_index != string::npos )
        break;
    }

    // Nothing found -> return empty string
    if (line_index == string::npos)
      return false;

    // Extract substring (all characters behind query)
    string sub_string = file_content_[line_index].substr(line_pos + query.size());
    
    // Remove everything behind the comment identifier
    sub_string = sub_string.substr(0, sub_string.find(comment_));

    // Remove all whitespace characters
    sub_string.erase(std::remove(sub_string.begin(), sub_string.end(), ' '), sub_string.end());

    // Write data to parameter
    try 
    {
      param.set_value( 0, string_to_single_value<T>(sub_string) );
      param.found( true );
    }
    catch (...)
    {
      return false;
    }

    return true;

  } // ParaReader::query_scalar()

  /*------------------------------------------------------------------
  | Query multiple parameters in a single line
  | --> Parameter existence has been checked in calling function
  ------------------------------------------------------------------*/
  template <typename T>
  bool query_single_line(Parameter<T>& param)
  {
    string query = param.single_key();

    // Parameter query reached end of input file
    if ( param.line_to_start() >= file_content_.size() )
      return false;

    // Position of query in the line
    size_t line_pos = string::npos;

    // Index of line where query is found
    size_t line_index = string::npos;

    size_t index;

    // Search for query in file buffer line by line
    for (size_t cur_line = param.line_to_start(); 
         cur_line < file_content_.size(); 
         ++cur_line)
    {
      const string& line = file_content_[cur_line];

      // Search in all lines for query
      // Use last query definition that occures in a line
      size_t pos = 0;

      while ((index = line.find(query, pos)) != string::npos)
      {
        // Position is from next element of index
        line_pos   = index;
        line_index = cur_line;
        pos        = index + 1;

        // Set starting line for next query search
        param.line_to_start( cur_line + 1 );
      }

      // Stop if query has been found in current line
      if ( line_index != string::npos )
        break;

    }

    // Nothing found -> return empty string
    if (line_index == string::npos)
      return false;

    // Extract substring (all characters behind query)
    string sub_string = file_content_[line_index].substr(line_pos + query.size());
    
    // Remove everything behind the comment identifier
    sub_string = sub_string.substr(0, sub_string.find(comment_));

    // Remove all whitespace characters
    sub_string.erase(std::remove(sub_string.begin(), sub_string.end(), ' '), sub_string.end());


    // Convert substring to a vector of type T
    std::vector<T> out;

    string s;
    std::stringstream ss(sub_string);

    // Split string at delimiter and put every sub-string
    // into "out" vector
    while(std::getline(ss, s, delimiter_))
    {
      // Remove parantheses
      s.erase(std::remove(s.begin(), s.end(), '('), s.end());
      s.erase(std::remove(s.begin(), s.end(), ')'), s.end());
      s.erase(std::remove(s.begin(), s.end(), '['), s.end());
      s.erase(std::remove(s.begin(), s.end(), ']'), s.end());

      try
      { 
        out.push_back( string_to_single_value<T>(s) );
      }
      catch (...)
      {
        continue;
      }
    }

    // Return false,  if not enough parameters
    if ( out.size() < param.columns() )
      return false;

    // write data into parameter object
    // --> consider only defined parameter size
    for ( size_t i = 0; i < param.columns(); ++i )
      param.set_value(i, out[i]);
    param.found( true );

    return true;

  } // ParaReader::query_single_line()

  /*------------------------------------------------------------------
  | Query list parameters over multiple lines
  ------------------------------------------------------------------*/
  template <typename T>
  bool query_multiple_lines(Parameter<T>& param)
  {
    string start_query = param.start_key();
    string end_query   = param.end_key();

    // Parameter query reached end of input file
    if ( param.line_to_start() >= file_content_.size() )
      return false;

    // Position of query in the line
    size_t line_start_pos = -1;
    size_t line_end_pos = -1;

    // Index of starting / ending line where query is found
    size_t line_start_index = -1;
    size_t line_end_index = -1;

    size_t index;

    // Search for starting query in file buffer line by line
    for (size_t cur_line = param.line_to_start(); 
         cur_line < file_content_.size(); 
         ++cur_line)
    {
      const string& line = file_content_[cur_line];

      // Search in all lines for starting query
      // use last query that occured in file
      size_t pos = 0;

      while ((index = line.find(start_query, pos)) != string::npos)
      {
        // Position is from next element of index
        line_start_pos   = index;
        line_start_index = cur_line;
        pos              = index + 1;

        // Set starting line for next query search
        param.line_to_start( cur_line + 1 );
      }

      // Stop if query has been found in current line
      if ( line_start_index != string::npos )
        break;
    }

    // Nothing found -> return empty string
    if (line_start_index == string::npos)
      return false;



    // Search for ending query in file buffer line by line
    for (size_t cur_line = param.line_to_start(); 
         cur_line < file_content_.size(); 
         ++cur_line)
    {
      const string& line = file_content_[cur_line];

      // Search in all lines for starting query
      // use last query that occured in file
      size_t pos = 0;

      while ((index = line.find(end_query, pos)) != string::npos)
      {
        // Position is from next element of index
        line_end_pos   = index;
        line_end_index = cur_line;
        pos            = index + 1;

        // Set starting line for next query search
        param.line_to_start( cur_line + 1 );
      }

      // Stop if query has been found in current line
      if ( line_end_index != string::npos )
        break;
    }

    // Nothing found -> return empty string
    if (line_end_index == string::npos)
      return false;

    std::vector<string> sub_strings {};

    for (int i = line_start_index+1; i < line_end_index; i++)
      sub_strings.push_back( file_content_[i] );


    // Convert substring to a vector of type T
    std::vector<T> out;

    for ( size_t i = 0; i < sub_strings.size(); ++i )
    {
      string s;
      std::stringstream ss( sub_strings[i] );

      // Split string at delimiter and put every sub-string
      // into "out" vector
      while(std::getline(ss, s, delimiter_))
      {
        // Remove parantheses
        s.erase(std::remove(s.begin(), s.end(), '('), s.end());
        s.erase(std::remove(s.begin(), s.end(), ')'), s.end());
        s.erase(std::remove(s.begin(), s.end(), '['), s.end());
        s.erase(std::remove(s.begin(), s.end(), ']'), s.end());

        try
        {
          out.push_back( string_to_single_value<T>(s) );
        }
        catch (...)
        {
          continue; 
        }
      }
    }

    // Check for correct shape of input data
    size_t n_rows = out.size() / param.columns();

    if ( out.size() != param.columns()*n_rows ) 
      return false;

    for ( size_t j = 0; j < n_rows; ++j )
    {
      for ( size_t i = 0; i < param.columns(); ++i )
      {
        size_t index = i + j * param.columns();
        param.set_value(index, out[index]);
      }

      if ( j < n_rows-1 )
        param.add_row();
    }

    param.found( true);


    return true;

  } // ParaReader::query_multiple_lines()


  /*------------------------------------------------------------------
  | Function to convert a string to a type T
  ------------------------------------------------------------------*/
  template <typename T> 
  T string_to_single_value (const string &str)
  {
    std::istringstream ss(str);
    T num;
    ss >> num;

    if ( ss.fail() )
      error("Failed to read input file");

    return num;
  }

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  string          file_path_;
  strVec          file_content_;
  string          comment_;
  char            delimiter_;

  ParameterList   param_list_; 
  ParameterMap    param_map_;



}; // ParaReader




} // namespace CppUtils
