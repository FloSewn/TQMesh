#***********************************************************
#
#***********************************************************
cmake_minimum_required(VERSION 3.10)

# --- C++ standard ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(TQMesh
        LANGUAGES CXX
        VERSION   1.4.0)

set( BIN ${CMAKE_SOURCE_DIR}/bin )
if(WIN32)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "$<0:>${CMAKE_SOURCE_DIR}/bin") # .exe and .dll
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "$<0:>${CMAKE_SOURCE_DIR}/lib") # .so and .dylib
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "$<0:>${CMAKE_SOURCE_DIR}/lib") # .lib and .a
endif()
# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# --- Compiler flags ---
if (CMAKE_BUILD_TYPE MATCHES "Release")
  if (WIN32)
    add_compile_definitions("_USE_MATH_DEFINES")
    if (MSVC)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bigobj")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL GNU)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreturn-type")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-compare")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTQMESH_USE_EXPRTK")
endif()

if (CMAKE_BUILD_TYPE MATCHES "Debug")
  if (WIN32)
    add_compile_definitions("_USE_MATH_DEFINES")
    if (MSVC)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bigobj")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Od")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreturn-type")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-compare")
  endif()
  if (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTQMESH_USE_EXPRTK")
endif()

# Add config file
configure_file(auxiliary/TQMeshConfig.h.in ${CMAKE_BINARY_DIR}/TQMeshConfig.h)
include_directories(${CMAKE_BINARY_DIR})

# Enable testing
# enable_testing()

# Directories
add_subdirectory( src/extern_libs )
add_subdirectory( src/utils )
add_subdirectory( src/algorithm )
add_subdirectory( src/tests )
add_subdirectory( src/examples )
add_subdirectory( src/app )

# Info
message(STATUS "CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_COMPILER_ID is ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION is ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")
