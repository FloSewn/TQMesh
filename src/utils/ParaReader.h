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
using strVec_ptr    = std::shared_ptr<strVec>;


/*--------------------------------------------------------------------
| Basic parameter types
--------------------------------------------------------------------*/
enum class ParaType 
{ scalar, vector, matrix, block };

/*--------------------------------------------------------------------
| Delimiter char and comment marker
--------------------------------------------------------------------*/
static inline string ParaFileComment { "#" };
static inline char ParaFileDelimiter { ',' };

/*--------------------------------------------------------------------
| A simple container for query data
--------------------------------------------------------------------*/
struct QueryContainer
{
  QueryContainer() {};

  bool   found      { false };
  size_t line_index { string::npos };
  size_t line_pos   { string::npos };
  size_t query_size { 0 };
}; 


/*--------------------------------------------------------------------
| The parameter base class 
--------------------------------------------------------------------*/
class ParameterBase
{
public:
  ParameterBase(ParaType type)
  : type_        { type }
  , start_key_   { "" }
  , end_key_     { "" }
  , block_start_ { 0 }
  , block_end_   { 0 }
  {}

  ParameterBase(ParaType type, const string& key,
                size_t block_start, size_t block_end)
  : type_        { type }
  , start_key_   { key }
  , end_key_     { "" }
  , block_start_ { block_start }
  , block_end_   { block_end }
  {}

  ParameterBase(ParaType type,
                const string& start, const string& end,
                size_t block_start, size_t block_end)
  : type_        { type }
  , start_key_   { start }
  , end_key_     { end }
  , block_start_ { block_start }
  , block_end_   { block_end }
  {}

  virtual ~ParameterBase() = default;

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  ParaType type() const { return type_; }

  const string& key() const { return start_key_; }
  string& key() { return start_key_; }

  const string& start_key() const { return start_key_; }
  string& start_key() { return start_key_; }

  const string& end_key() const { return end_key_; }
  string& end_key() { return end_key_; }

  size_t block_start() const { return block_start_; }

  size_t block_end() const { return block_end_; }

  bool found() const { return found_; }

  int block_index() const { return block_index_; }

  bool reached_end() const { return block_index_ >= block_end_; }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  void block_start(size_t i) { block_start_ = i; }

  void block_end(size_t i) { block_end_ = i; }

  void found(bool f) { found_ = f; }

  void block_index(int l) { block_index_ = l; }

  /*------------------------------------------------------------------
  | Query the parameter -> Search in file content for its occurrence
  ------------------------------------------------------------------*/
  QueryContainer get_query_data(const string& query, 
                                const strVec& content) 
  {
    QueryContainer query_data {};

    // Parameter query reached end of input file
    if ( reached_end() )
      return query_data;

    // Position of query in the line
    size_t line_pos = string::npos;

    // Index of line where query is found
    size_t line_index = string::npos;

    size_t index;

    // Search for query in file buffer line by line
    for (size_t cur_line = block_start() + block_index(); 
         cur_line < content.size(); 
         ++cur_line)
    {
      const string& line = content[cur_line];

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
        block_index( cur_line + 1 - block_start() );
      }

      // Stop if query has been found in current line
      if ( line_index != string::npos )
        break;
    }

    // Nothing found -> return empty string
    if (line_index == string::npos)
      return query_data;

    // Copy data to query container
    query_data.found      = true;
    query_data.line_index = line_index;
    query_data.line_pos   = line_pos;
    query_data.query_size = query.size();

    return query_data;

  } // ParameterBase::get_query_data()





protected:
  ParaType type_;
  string   start_key_;
  string   end_key_;
  size_t   block_start_;
  size_t   block_end_;

  size_t   block_index_ { 0 };
  bool     found_      { false };


}; // ParameterBase


/*--------------------------------------------------------------------
| The derived template class for a specified parameter
--------------------------------------------------------------------*/
template <class T>
class Parameter : public ParameterBase
{
public:
  /*------------------------------------------------------------------
  | Constructor - Scalar / Vector
  ------------------------------------------------------------------*/
  Parameter(ParaType type, const string& key, size_t n,
            size_t block_start, size_t block_end) 
    : ParameterBase(type, key, block_start, block_end)
    , ncol_{n}
    , nrow_{1}
    , values_(n,T{})
  {}

  /*------------------------------------------------------------------
  | Constructor - Matrix 
  ------------------------------------------------------------------*/
  Parameter(ParaType type, 
            const string& start, const string& end, size_t n,
            size_t block_start, size_t block_end) 
    : ParameterBase(type, start, end, block_start, block_end)
    , ncol_{n}
    , nrow_{1}
    , values_(n,T{})
  {}

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  T get_value(size_t i) const
  { 
    if ( i >= ncol_*nrow_ )
      return values_[0];
    return values_[i]; 
  }

  T get_value(size_t i, size_t j) const
  { 
    size_t index = j * ncol_ + i;
    if ( index >= ncol_*nrow_ )
      return values_[0];
    return values_[index]; 
  }

  size_t columns() const { return ncol_; }
  size_t rows() const { return nrow_; }

  void resize() 
  { 
    nrow_ = 1;
    values_.resize( ncol_ ); 
  }

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
    if ( type_ != ParaType::matrix )
      return;

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
| A parameter block where data is stored
--------------------------------------------------------------------*/
class ParaBlock : public ParameterBase
{
  using ParaBlockList = std::vector<ParaBlock>;
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
  ------------------------------------------------------------------*/
  ParaBlock(const string& file_path)
  : ParameterBase( ParaType::block )
  { 
    // Read the file content
    content_ = read_content( file_path );

    // Set bounds of top ParaBlock
    block_start( 0 );
    block_end( content_->size() );
  }

  ParaBlock(const string& start, const string& end,
            size_t block_start, size_t block_end)
  : ParameterBase( ParaType::block, start, end, block_start, block_end )
  { }


  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  void set_content(strVec_ptr c) { content_ = c; }

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  const ParameterList& para_list() const { return para_list_; }
  ParameterList& para_list() { return para_list_; }

  /*------------------------------------------------------------------
  | Create new block parameter to search for in a file
  ------------------------------------------------------------------*/
  void new_block_parameter(const string& name,
                           const string& start, const string& end) 
  {
    // Check that parameter map does not yet contain given name
    if ( para_map_.count(name) )
      error("Multiple definitions for parmeter name \"" + name + "\"."); 

    // Empty parameter keys are not allowed
    if ( start.size() < 1 || end.size() < 1 )
      error("Empty parameter queries are not allowed. Parameter: " + name);

    ParaBlock* new_b = new ParaBlock { start, end, 
                                       block_start_, block_end_ };

    // Update parameter list
    para_list_.emplace_back( new_b );

    // Update parameter map
    para_map_[name] = para_list_.size() - 1;

    // Pass content to new block
    new_b->set_content( content_ );

  }

  /*------------------------------------------------------------------
  | Create new scalar parameters to search for in a file
  ------------------------------------------------------------------*/
  template <typename T>
  void new_scalar_parameter(const string& name, const string& key)
  {
    // Check that parameter map does not yet contain given name
    if ( para_map_.count(name) )
      error("Multiple definitions for parmeter name \"" + name + "\"."); 

    // Empty parameter keys are not allowed
    if ( key.size() < 1 )
      error("Empty parameter queries are not allowed. Parameter: " + name);

    // Update parameter list
    para_list_.emplace_back( 
      new Parameter<T> { ParaType::scalar, key, 1, 
                         block_start_, block_end_ } 
    );

    // Update parameter map
    para_map_[name] = para_list_.size() - 1;
  }

  /*------------------------------------------------------------------
  | Create new vector parameters to search for in a file
  ------------------------------------------------------------------*/
  template <typename T>
  void new_vector_parameter(const string& name, 
                            const string& key, size_t n)
  {
    // Check that parameter map does not yet contain given name
    if ( para_map_.count(name) )
      error("Multiple definitions for parmeter name \"" + name + "\"."); 

    // Empty parameter keys are not allowed
    if ( key.size() < 1 )
      error("Empty parameter queries are not allowed. Parameter: " + name);

    // Update parameter list
    para_list_.emplace_back( 
      new Parameter<T> { ParaType::vector, key, n,
                         block_start_, block_end_ } 
    );

    // Update parameter map
    para_map_[name] = para_list_.size() - 1;
  }

  /*------------------------------------------------------------------
  | Create new matrix parameters to search for in a file
  ------------------------------------------------------------------*/
  template <typename T>
  void new_matrix_parameter(const string& name,
                            const string& start, const string& end, 
                            size_t n)
  {
    // Check that parameter map does not yet contain given name
    if ( para_map_.count(name) )
      error("Multiple definitions for parmeter name \"" + name + "\"."); 

    // Empty parameter keys are not allowed
    if ( start.size() < 1 || end.size() < 1 )
      error("Empty parameter queries are not allowed. Parameter: " + name);

    // Update parameter list
    para_list_.emplace_back( 
      new Parameter<T> { ParaType::matrix, start, end, n,
                         block_start_, block_end_ } 
    );

    // Update parameter map
    para_map_[name] = para_list_.size() - 1;
  }

  /*------------------------------------------------------------------
  | Get a parameter reference from a given parameter name
  ------------------------------------------------------------------*/
  template <typename T>
  Parameter<T>& get_parameter(const string& name)
  {
    // Handle unknown parameter names
    if ( !para_map_.count(name) )
      error("No parameter with name \"" + name + "\" has been defined.");

    // Cast parameter 
    Parameter<T>* para_ptr = dynamic_cast<Parameter<T>*>(
        para_list_[ para_map_[name] ].get()
    );

    return *para_ptr;
  }

  /*------------------------------------------------------------------
  | Get a parameter reference from a given parameter name
  ------------------------------------------------------------------*/
  ParaBlock& get_block(const string& name)
  {
    // Handle unknown parameter names
    if ( !para_map_.count(name) )
      error("No block with name \"" + name + "\" has been defined.");

    // Cast parameter 
    ParaBlock* para_ptr = dynamic_cast<ParaBlock*>(
        para_list_[ para_map_[name] ].get()
    );

    return *para_ptr;
  }
   
  /*------------------------------------------------------------------
  | Return the first value of a parameter with <name>.
  ------------------------------------------------------------------*/
  template <typename T>
  T get_value(const string& name)
  {
    Parameter<T>& para = get_parameter<T>( name );

    // Check if parameter is not found in the input file
    if ( !para.found() )
      error("Parameter with name \"" + name + "\" not found in input file.");
    
    return para.get_value(0);
  }

  /*------------------------------------------------------------------
  | Return the value of a parameter with <name>.
  ------------------------------------------------------------------*/
  template <typename T>
  T get_value(size_t i, const string& name)
  {
    Parameter<T>& para = get_parameter<T>( name );

    // Check if parameter is not found in the input file
    if ( !para.found() )
      error("Parameter with name \"" + name + "\" not found in input file.");
    
    return para.get_value(i);
  }

  /*------------------------------------------------------------------
  | Return the value of a parameter with <name> and return it.
  ------------------------------------------------------------------*/
  template <typename T>
  T get_value(size_t i, size_t j, const string& name)
  {
    Parameter<T>& para = get_parameter<T>( name );

    // Check if parameter is not found in the input file
    if ( !para.found() )
      error("Parameter with name \"" + name + "\" not found in input file.");
    
    return para.get_value(i,j);
  }

  /*------------------------------------------------------------------
  | Check if a parameter with a given name has been found in the file
  ------------------------------------------------------------------*/
  bool found(const string& name)
  {
    // Handle unknown parameter names
    if ( !para_map_.count(name) )
      error("No parameter with name \"" + name + "\" has been defined.");

    return para_list_[ para_map_[name] ].get()->found();
  }

  /*------------------------------------------------------------------
  | Query  parameters
  ------------------------------------------------------------------*/
  template <typename T>
  bool query(const string& name)
  {
    // Account for empty content
    if (!content_)
      return false;

    // Query actual parameter
    Parameter<T>& para = get_parameter<T>( name );

    if ( para.type() == ParaType::scalar )
    {
      return query_scalar( para, *(content_) );
    }
    else if ( para.type() == ParaType::vector ) 
    {
      return query_vector( para, *(content_) );
    }
    else if ( para.type() == ParaType::matrix )
    {
      return query_matrix( para, *(content_) );
    }

    return false;

  } // ParaBlock::query()

  /*------------------------------------------------------------------
  | Query the block parameters
  ------------------------------------------------------------------*/
  bool query(const string& name)
  {
    // Account for empty content
    if (!content_)
      return false;

    // Query actual parameter
    ParaBlock& block = get_block( name );

    if ( block.block_end() < content_->size() )
    {
      block.block_start( block.block_end() );
      block.block_end( content_->size() );
      block.block_index( 0 );
    }
  
    auto start_data = block.get_query_data(block.start_key(), 
                                           *content_ );

    if ( !start_data.found )
      return false;

    auto end_data = block.get_query_data(block.end_key(), *content_ );

    if ( !end_data.found )
      return false;

    // Set start and end of this block
    block.block_start( start_data.line_index );
    block.block_end( end_data.line_index );

    // Set starts and ends of all parameters of this block
    for ( auto& para : block.para_list() )
    {
      para->block_start( block.block_start() );
      para->block_end( block.block_end() );
      para->block_index( 0 );
    }

    return true;

  } // ParaBlock::query() 

private:

  /*------------------------------------------------------------------
  | Query scalar parameters
  ------------------------------------------------------------------*/
  template <typename T>
  bool query_scalar(Parameter<T>& para, const strVec& content)
  {
    auto query_data = para.get_query_data(para.start_key(), content );

    if ( !query_data.found )
      return false;

    // Extract substring (all characters behind query)
    string sub_string = content[query_data.line_index].substr(
        query_data.line_pos + query_data.query_size
    );
    
    // Remove everything behind the comment identifier
    sub_string = sub_string.substr(0, sub_string.find(ParaFileComment));

    // Remove all whitespace characters
    sub_string.erase(
        std::remove(sub_string.begin(), sub_string.end(), ' '), 
        sub_string.end()
    );

    // Write data to parameter
    try 
    {
      para.set_value( 0, string_to_single_value<T>(sub_string) );
      para.found( true );
    }
    catch (...)
    {
      return false;
    }

    return true;

  } // ParaBlock::query_scalar()

  /*------------------------------------------------------------------
  | Query vector parameters 
  ------------------------------------------------------------------*/
  template <typename T>
  bool query_vector(Parameter<T>& para, const strVec& content)
  {
    auto query_data = para.get_query_data(para.start_key(), content );

    if ( !query_data.found )
      return false;

    // Extract substring (all characters behind query)
    string sub_string = content[query_data.line_index].substr(
        query_data.line_pos + query_data.query_size 
    );
    
    // Remove everything behind the comment identifier
    sub_string = sub_string.substr(0, sub_string.find(ParaFileComment));

    // Remove all whitespace characters
    sub_string.erase(
        std::remove(sub_string.begin(), sub_string.end(), ' '), 
        sub_string.end()
    );

    // Convert substring to a vector of type T
    std::vector<T> out;

    string s;
    std::stringstream ss(sub_string);

    // Split string at delimiter and put every sub-string
    // into "out" vector
    while(std::getline(ss, s, ParaFileDelimiter))
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
    if ( out.size() < para.columns() )
      return false;

    // write data into parameter object
    // --> consider only defined parameter size
    for ( size_t i = 0; i < para.columns(); ++i )
      para.set_value(i, out[i]);
    para.found( true );

    return true;

  } // ParaBlock::query_vector()

  /*------------------------------------------------------------------
  | Query matrix parameters 
  ------------------------------------------------------------------*/
  template <typename T>
  bool query_matrix(Parameter<T>& para, const strVec& content)
  {
    auto start_data = para.get_query_data(para.start_key(), content );

    if ( !start_data.found )
      return false;

    auto end_data = para.get_query_data(para.end_key(), content );

    if ( !end_data.found )
      return false;

    // Gather sub strings
    std::vector<string> sub_strings {};

    for (int i = start_data.line_index+1; i < end_data.line_index; i++)
      sub_strings.push_back( content[i] );

    // Convert substring to a vector of type T
    std::vector<T> out;

    for ( size_t i = 0; i < sub_strings.size(); ++i )
    {
      string s;
      std::stringstream ss( sub_strings[i] );

      // Split string at delimiter and put every sub-string
      // into "out" vector
      while(std::getline(ss, s, ParaFileDelimiter))
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
    size_t n_rows = out.size() / para.columns();

    if ( out.size() != para.columns()*n_rows ) 
      return false;

    // Resize the parameter values to initial size 
    para.resize();

    for ( size_t j = 0; j < n_rows; ++j )
    {
      for ( size_t i = 0; i < para.columns(); ++i )
      {
        size_t index = i + j * para.columns();
        para.set_value(index, out[index]);
      }

      if ( j < n_rows-1 )
        para.add_row();
    }

    para.found( true);


    return true;

  } // ParaBlock::query_matrix()


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
  | Get the content from a text file
  ------------------------------------------------------------------*/
  strVec_ptr read_content(const string& file_path)
  {
    strVec_ptr content = std::make_shared<strVec>();

    ifstream ifs {file_path};
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
      content->push_back(line);
    }

    if (!ifs.eof()) 
      error("Failed to read entire file - "
            "Stopped reading in line " + 
            std::to_string(content->size()));

    return content;
  }


  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  ParameterList   para_list_; 
  ParameterMap    para_map_;

  strVec_ptr      content_;

  size_t          line_end_ { 0 };


}; // ParaBlock

using ParaReader = ParaBlock;

} // namespace CppUtils
