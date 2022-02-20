#***********************************************************
#
#***********************************************************
cmake_minimum_required(VERSION 3.10)
project(TQMesh
        LANGUAGES CXX
        VERSION   1.0)

# Set C++ flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wunused-variable -std=c++17 -O3 -DNDEBUG -pg")
set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wunused-variable -std=c++17 -O3 -DNDEBUG")
#set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++17 -O0 -g -pg")
#set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wunused-variable -std=c++17 -O0 -g")


set( BIN ${CMAKE_SOURCE_DIR}/bin )


add_subdirectory( src/utils )
add_subdirectory( src/algorithm )
add_subdirectory( src/tests )
add_subdirectory( src/app )
