/*
* This source file is part of the tqmesh library.  
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


/*--------------------------------------------------------------------
| A scalar parameter (defined in a single line)
--------------------------------------------------------------------*/
template <class T>
class ScalarParameter
{
  using string = std::string;

public:
  ScalarParameter(const string& key) 
    : key_{key}, line_start_{0}, found_{false} {}

  const T& value() const { return value_; }
  T& value() { return value_; }

  const string& key() const { return key_;   }
  string& key() { return key_;   }

  void value(T val) { value_ = val; };

  int line_to_start() const { return line_start_; }
  void line_to_start(int l) { line_start_ = l; }

  bool found() const { return found_; }
  void found(bool f) { found_ = f; }

private:
  string key_;
  size_t line_start_;
  bool found_;

  T value_;
};

/*--------------------------------------------------------------------
| A list parameter (may be defined over multiple lines)
--------------------------------------------------------------------*/
template <class T>
class ListParameter
{
  using string = std::string;

public:
  ListParameter(const string& start, const string& end) 
    : start_{start}, end_{end}, key_ {""}
    , line_start_ {0}, found_{false}
  { multi_lines_ = true; }

  ListParameter(const string& key) 
    : start_{""}, end_{""}, key_{key}
    , line_start_ {0}, found_{false}
  { multi_lines_ = false; }

  const T& value(size_t i) const { return values_[i]; }
  T& value(size_t i) { return values_[i]; }

  void values(std::vector<T> v) { values_ = std::move(v); }
  void add_value(T val) { values_.push_back( val ); }

  const string& start() const { return start_; }
  string& start() { return start_; }

  const string& end() const { return end_; }
  string& end() { return end_; }

  const string& key() const { return key_;   }
  string& key() { return key_;   }

  void columns(size_t i) { ncol_ = i; }
  size_t columns() const { return ncol_; }

  void rows(size_t i) { nrow_ = i; }
  size_t rows() const { return nrow_; }

  bool multiple_lines() const { return multi_lines_; }

  int line_to_start() const { return line_start_; }
  void line_to_start(int l) { line_start_ = l; }

  bool found() const { return found_; }
  void found(bool f) { found_ = f; }

private:
  string start_;
  string end_;
  string key_;
  size_t line_start_;
  bool   found_;

  std::vector<T> values_;

  // Dimension of dataset
  bool   multi_lines_ {false};
  size_t ncol_        { 0 };
  size_t nrow_        { 0 };

};


/*--------------------------------------------------------------------
| The parameter file reader 
--------------------------------------------------------------------*/
class ParaReader
{
  /*------------------------------------------------------------------
  | Some local typedefs
  ------------------------------------------------------------------*/
  using string        = std::string;
  using ofstream      = std::ofstream;
  using ifstream      = std::ifstream;
  using strVec        = std::vector<string>;
  using istringstream = std::istringstream;


public:

  /*------------------------------------------------------------------
  | Error handling
  ------------------------------------------------------------------*/
  class Invalid
  {
  public:
    Invalid(const string& msg){ error_message = msg; }
    string& what() { return error_message; }
  private:
    string error_message;
  };

  /*--------------------------------------------------------
  | Function for error handling
  --------------------------------------------------------*/
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
    if (!ifs) error("Can't open output file: " + file_path);

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
  }

  /*------------------------------------------------------------------
  | Create new scalar parameters to search for in a file
  ------------------------------------------------------------------*/
  ScalarParameter<bool>& 
  new_bool_parameter(const string& key)
  {
    ScalarParameter<bool> parameter {key};
    scalar_bool_params_.push_back(parameter);
    return scalar_bool_params_[scalar_bool_params_.size()-1];
  }

  ScalarParameter<int>& 
  new_int_parameter(const string& key)
  {
    ScalarParameter<int> parameter {key};
    scalar_int_params_.push_back(parameter);
    return scalar_int_params_[scalar_int_params_.size()-1];
  }

  ScalarParameter<double>& 
  new_double_parameter(const string& key)
  {
    ScalarParameter<double> parameter {key};
    scalar_dbl_params_.push_back(parameter);
    return scalar_dbl_params_[scalar_dbl_params_.size()-1];
  }

  ScalarParameter<string>& 
  new_string_parameter(const string& key)
  {
    ScalarParameter<string> parameter {key};
    scalar_str_params_.push_back(parameter);
    return scalar_str_params_[scalar_str_params_.size()-1];
  }

  /*------------------------------------------------------------------
  | Create new list parameters to search for in a file
  ------------------------------------------------------------------*/
  ListParameter<int>& 
  new_int_list_parameter(const string& key)
  {
    ListParameter<int> parameter {key};
    list_int_params_.push_back(parameter);
    return list_int_params_[list_int_params_.size()-1];
  }

  ListParameter<int>& 
  new_int_list_parameter(const string& start, const string& end)
  {
    ListParameter<int> parameter {start, end};
    list_int_params_.push_back(parameter);
    return list_int_params_[list_int_params_.size()-1];
  }

  ListParameter<double>& 
  new_double_list_parameter(const string& key)
  {
    ListParameter<double> parameter {key};
    list_dbl_params_.push_back(parameter);
    return list_dbl_params_[list_dbl_params_.size()-1];
  }

  ListParameter<double>& 
  new_double_list_parameter(const string& start, const string& end)
  {
    ListParameter<double> parameter {start, end};
    list_dbl_params_.push_back(parameter);
    return list_dbl_params_[list_dbl_params_.size()-1];
  }

  ListParameter<string>& 
  new_string_list_parameter(const string& key)
  {
    ListParameter<string> parameter {key};
    list_str_params_.push_back(parameter);
    return list_str_params_[list_str_params_.size()-1];
  }

  ListParameter<string>& 
  new_string_list_parameter(const string& start, const string& end)
  {
    ListParameter<string> parameter {start, end};
    list_str_params_.push_back(parameter);
    return list_str_params_[list_str_params_.size()-1];
  }

  /*------------------------------------------------------------------
  | Query scalar parameters
  ------------------------------------------------------------------*/
  template <typename T>
  bool query(ScalarParameter<T>& param)
  {
    string query = param.key();

    // Query check
    if (query.size() < 1) error("Invalid query");

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
      param.value( string_to_single_value<T>(sub_string) );
      param.found( true );
    }
    catch (...)
    {
      MSG("Failed to read parameter " << param.key() 
          << " in line " << param.line_to_start() 
          << " of input parameter file.");
      return false;
    }

    return true;
  }

  /*------------------------------------------------------------------
  | Query list parameters
  ------------------------------------------------------------------*/
  template <typename T>
  bool query(ListParameter<T>& param)
  {
    if ( param.multiple_lines() )
      return query_multiple_lines(param);
    else
      return query_single_line(param);
  }


private:

  /*------------------------------------------------------------------
  | Query list parameters over a single line
  ------------------------------------------------------------------*/
  template <typename T>
  bool query_single_line(ListParameter<T>& param)
  {
    string query = param.key();

    // Query check
    if (query.size() < 1) error("Invalid query");

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
        MSG("Failed to read parameter " << param.key() 
            << " in line " << param.line_to_start() 
            << " of input parameter file.");
        return false;
      }
    }

    // write data into parameter object
    param.values( out );
    param.rows( 1 );
    param.columns( out.size() );
    param.found( true );

    return true;

  } // ParaReader::query_single_line()

  /*------------------------------------------------------------------
  | Query list parameters over multiple lines
  ------------------------------------------------------------------*/
  template <typename T>
  bool query_multiple_lines(ListParameter<T>& param)
  {
    string start_query = param.start();
    string end_query   = param.end();

    // Query check
    if (start_query.size() < 1 || end_query.size() < 1) 
      error("Invalid query");

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

    size_t n_cols = 0;

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
          MSG("Failed to read parameter " << param.key() 
              << " in line " << param.line_to_start() 
              << " of input parameter file.");
          return false;
        }
      }

      // Check for correct shape of input data
      if ( n_cols == 0 )
      {
        n_cols = out.size();
        param.columns( n_cols );
      }
      else if ( n_cols != param.columns() )
        error("Invalid input data shape");

    }

    param.values( out );
    param.rows( sub_strings.size() );
    param.found( true );

    return true;

  } // ParaReader::query_multiple_lines()


  /*--------------------------------------------------------
  | Function to convert a string to a type T
  --------------------------------------------------------*/
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

  // Scalar parameters 
  std::vector<ScalarParameter<bool>>   scalar_bool_params_;
  std::vector<ScalarParameter<int>>    scalar_int_params_;
  std::vector<ScalarParameter<double>> scalar_dbl_params_;
  std::vector<ScalarParameter<string>> scalar_str_params_;

  // List parameters 
  std::vector<ListParameter<int>>      list_int_params_;
  std::vector<ListParameter<double>>   list_dbl_params_;
  std::vector<ListParameter<string>>   list_str_params_;



};
