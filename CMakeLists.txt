#***********************************************************
#
#***********************************************************
cmake_minimum_required(VERSION 3.10)

# --- C++ standard ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(TQMesh
        LANGUAGES CXX
        VERSION   1.1)

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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bigobj")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
  else()
<<<<<<< HEAD
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
=======
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreturn-type")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-compare")
>>>>>>> fa0899f5faedbc3de2d30dba4c8c9fc7b7288940
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTQMESH_USE_EXPRTK")

  # --- For performance measuring only ---
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

if (CMAKE_BUILD_TYPE MATCHES "Debug")
  if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bigobj")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Od")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
<<<<<<< HEAD
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable")
=======
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreturn-type")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-compare")
>>>>>>> fa0899f5faedbc3de2d30dba4c8c9fc7b7288940
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTQMESH_USE_EXPRTK")
endif()

# Add config file
configure_file(auxiliary/TQMeshConfig.h.in ${CMAKE_BINARY_DIR}/TQMeshConfig.h)
include_directories(${CMAKE_BINARY_DIR})

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

#---------------------------------------------------------------------
# Testing
#---------------------------------------------------------------------
enable_testing()

add_test(NAME Vertex COMMAND ${CMAKE_SOURCE_DIR}/bin/run_tests "Vertex")
add_test(NAME Triangle COMMAND ${CMAKE_SOURCE_DIR}/bin/run_tests "Triangle")
add_test(NAME Quad COMMAND ${CMAKE_SOURCE_DIR}/bin/run_tests "Quad")
add_test(NAME Front COMMAND ${CMAKE_SOURCE_DIR}/bin/run_tests "Front")
add_test(NAME EdgeList COMMAND ${CMAKE_SOURCE_DIR}/bin/run_tests "EdgeList")
add_test(NAME Boundary COMMAND ${CMAKE_SOURCE_DIR}/bin/run_tests "Boundary")
add_test(NAME SizeFunction COMMAND ${CMAKE_SOURCE_DIR}/bin/run_tests "SizeFunction")
add_test(NAME Mesh COMMAND ${CMAKE_SOURCE_DIR}/bin/run_tests "Mesh")
add_test(NAME MeshGenerator COMMAND ${CMAKE_SOURCE_DIR}/bin/run_tests "MeshGenerator")
add_test(NAME MeshCleanup COMMAND ${CMAKE_SOURCE_DIR}/bin/run_tests "MeshCleanup")
add_test(NAME SmoothingStrategy COMMAND ${CMAKE_SOURCE_DIR}/bin/run_tests "SmoothingStrategy")


