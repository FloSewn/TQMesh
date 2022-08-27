#***********************************************************
#
#***********************************************************
cmake_minimum_required(VERSION 3.10)

# --- C++ standard ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Compiler flags ---
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG")

set( BIN ${CMAKE_SOURCE_DIR}/bin )

project(TQMesh
        LANGUAGES CXX
        VERSION   1.1)

# Add config file
configure_file(aux/TQMeshConfig.h.in ${CMAKE_BINARY_DIR}/TQMeshConfig.h)
include_directories(${CMAKE_BINARY_DIR})

add_subdirectory( src/extern_libs )
add_subdirectory( src/utils )
add_subdirectory( src/algorithm )
add_subdirectory( src/tests )
add_subdirectory( src/examples )
add_subdirectory( src/app )


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
add_test(NAME Smoother COMMAND ${CMAKE_SOURCE_DIR}/bin/run_tests "Smoother")
add_test(NAME Mesh COMMAND ${CMAKE_SOURCE_DIR}/bin/run_tests "Mesh")


