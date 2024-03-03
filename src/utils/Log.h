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
#include <memory>

namespace CppUtils {

/*********************************************************************
* Log levels, stream types, colors
*********************************************************************/
enum LogLevel 
{ ERROR, WARNING, INFO, DEBUG };

enum OStreamType
{ TO_COUT, TO_CERR, TO_CLOG, TO_FILE };

enum LogColor {
  RED      = 31,
  GREEN    = 32,
  YELLOW   = 33,
  BLUE     = 34,
  PURPLE   = 35,
  CYAN     = 36,
  WHITE    = 37,
  DEFAULT  = 39,
};


/*********************************************************************
* Interface to create ostream unique_ptr, which gets properly 
* deleted.
*
* Reference:
* ----------
* -https://stackoverflow.com/questions/56521318/ostream-class-\
*  that-outputs-either-on-cout-or-on-a-file
*
*********************************************************************/
struct ConditionalDeleter
{
    bool must_delete;
    void operator()(std::ostream* os) const 
    { if (must_delete) delete os; }
};

using OStreamPtr = std::unique_ptr<std::ostream, ConditionalDeleter>;

static inline OStreamPtr create_stream(OStreamType type, const std::string& path="")
{
  switch( type ) {
    case TO_COUT:
      return OStreamPtr { &std::cout, ConditionalDeleter {false} };
    case TO_CERR: 
      return OStreamPtr { &std::cerr, ConditionalDeleter {false} };
    case TO_CLOG:
      return OStreamPtr { &std::clog, ConditionalDeleter {false} };
    case TO_FILE: 
      return OStreamPtr { new std::ofstream {path}, 
                          ConditionalDeleter {true} };
  }
  return OStreamPtr { &std::cout, ConditionalDeleter {false} };
}


/*********************************************************************
* The global logging properties
*********************************************************************/
class LogProperties
{
public:
  /*------------------------------------------------------------------
  | Default constructor
  ------------------------------------------------------------------*/
  LogProperties() 
  {
    error_os_ = create_stream( TO_COUT );
    warn_os_  = create_stream( TO_COUT );
    info_os_  = create_stream( TO_COUT );
    debug_os_ = create_stream( TO_COUT );
  }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  void set_level(LogLevel level) { level_ = level; }
  void show_header(bool show) { show_header_ = show; }
  void use_newline(bool nl) { use_newline_ = nl; }
  void use_color(bool c) { use_color_ = c; }

  void set_error_header(const std::string& msg) { error_header_ = msg; }
  void set_warn_header(const std::string& msg) { warn_header_ = msg; }
  void set_info_header(const std::string& msg) { info_header_ = msg; }
  void set_debug_header(const std::string& msg) { debug_header_ = msg; }

  void set_error_ostream(OStreamType type, const std::string& f="")
  { 
    error_os_      = create_stream( type, f ); 
    error_os_type_ = type;
  }
  void set_warn_ostream(OStreamType type, const std::string& f="")
  { 
    warn_os_      = create_stream( type, f ); 
    warn_os_type_ = type;
  }
  void set_info_ostream(OStreamType type, const std::string& f="")
  { 
    info_os_      = create_stream( type, f ); 
    info_os_type_ = type;
  }
  void set_debug_ostream(OStreamType type, const std::string& f="")
  { 
    debug_os_      = create_stream( type, f ); 
    debug_os_type_ = type;
  }

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  const LogLevel& level() const { return level_; }
  bool show_header() const { return show_header_; }
  bool use_newline() const { return use_newline_; }
  bool use_color() const { return use_color_; }

  const std::string& get_header(LogLevel level) const 
  {
    switch( level ) {
      case ERROR:   return error_header_; 
      case WARNING: return warn_header_; 
      case INFO:    return info_header_; 
      case DEBUG:   return debug_header_; 
    }
    return info_header_;
  }

  std::ostream& get_ostream(LogLevel level)
  {
    switch( level ) {
      case ERROR:   return *error_os_; 
      case WARNING: return *warn_os_; 
      case INFO:    return *info_os_; 
      case DEBUG:   return *debug_os_; 
    }
    return *info_os_;
  }

  OStreamType get_ostream_type(LogLevel level)
  {
    switch( level ) {
      case ERROR:   return error_os_type_; 
      case WARNING: return warn_os_type_; 
      case INFO:    return info_os_type_; 
      case DEBUG:   return debug_os_type_; 
    }
    return info_os_type_;
  }

  std::string get_color(LogLevel level)
  {
    std::string base = { "\033[" };
    switch( level ) {
      case ERROR:   
        base += std::to_string( error_col_ ); 
        break;
      case WARNING: 
        base += std::to_string( warn_col_ ); 
        break;
      case INFO:    
        base += std::to_string( info_col_ ); 
        break;
      case DEBUG:   
        base += std::to_string( debug_col_ ); 
        break;
      default:      
        base += std::to_string( LogColor::DEFAULT );
    }
    base += "m";
    return base;
  }


private:

  LogLevel    level_         = INFO;
  bool        show_header_   = true;
  bool        use_newline_   = true;
  bool        use_color_     = true;

  std::string error_header_  = "[ERROR] ";
  std::string warn_header_   = "[WARNING] ";
  std::string info_header_   = "[INFO] ";
  std::string debug_header_  = "[DEBUG] ";

  OStreamPtr error_os_ { nullptr };
  OStreamPtr warn_os_  { nullptr };
  OStreamPtr info_os_  { nullptr };
  OStreamPtr debug_os_ { nullptr };

  OStreamType error_os_type_ { TO_COUT };
  OStreamType warn_os_type_  { TO_COUT };
  OStreamType info_os_type_  { TO_COUT };
  OStreamType debug_os_type_ { TO_COUT };

  LogColor error_col_  { RED     };
  LogColor warn_col_   { YELLOW  };
  LogColor info_col_   { DEFAULT };
  LogColor debug_col_  { DEFAULT };

};

inline LogProperties LOG_PROPERTIES;

/*********************************************************************
* The interface for the actual SimpleLogger
*
* Reference:
* ----------
* -https://stackoverflow.com/questions/5028302/small-logger-class
*********************************************************************/
class LOG
{
public:
  /*------------------------------------------------------------------
  | Default constructor
  ------------------------------------------------------------------*/
  LOG() {}

  /*------------------------------------------------------------------
  | Constructror with log level specification
  ------------------------------------------------------------------*/
  LOG(LogLevel level)
  {
    level_ = level;

    if ( LOG_PROPERTIES.use_color() && 
         LOG_PROPERTIES.get_ostream_type(level_) != TO_FILE )
      operator<<( LOG_PROPERTIES.get_color( level_ ) );
    
    if ( LOG_PROPERTIES.show_header() )
      operator<<( LOG_PROPERTIES.get_header( level_ ) );
  }

  /*------------------------------------------------------------------
  | Constructror with log level and color specification
  ------------------------------------------------------------------*/
  LOG(LogLevel level, LogColor c)
  {
    level_ = level;

    if ( LOG_PROPERTIES.use_color() && 
         LOG_PROPERTIES.get_ostream_type(level_) != TO_FILE )
    {
      std::string color { "\033[" };
      color += std::to_string( c ); 
      color += "m";
      operator<<( color );
    }
    
    if ( LOG_PROPERTIES.show_header() )
      operator<<( LOG_PROPERTIES.get_header( level_ ) );
  }

  /*------------------------------------------------------------------
  | Destructor -> append new line, if property is set
  ------------------------------------------------------------------*/
  ~LOG() 
  {
    // Set default color
    if ( LOG_PROPERTIES.use_color() &&
         LOG_PROPERTIES.get_ostream_type(level_) != TO_FILE )
      operator<< ("\e[0m");

    // Append newline
    if ( opened_ && LOG_PROPERTIES.use_newline() )
      LOG_PROPERTIES.get_ostream( level_ ) << std::endl;

    // Reset 
    opened_ = false;
  }

  /*------------------------------------------------------------------
  | OStream operator
  ------------------------------------------------------------------*/
  template<class T>
  LOG& operator<<(const T& msg)
  {
    if ( level_ <= LOG_PROPERTIES.level() )
    {
      LOG_PROPERTIES.get_ostream( level_ ) << msg;
      opened_ = true;
    }
    return *this;
  }

private:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  bool     opened_ = false;
  LogLevel level_  = DEBUG;

}; // LOG

/*********************************************************************
* Additional debug macro
*********************************************************************/
#ifndef NDEBUG
#define DEBUG_LOG(str) \
  LOG(DEBUG) << str
#else
#define DEBUG_LOG(str) \
  do { } while (false)
#endif


} // namespace CppUtils
