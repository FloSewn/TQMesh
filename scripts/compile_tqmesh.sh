#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "$0 <compiler> <build type> <cores>"
    exit 1
fi


COMPILER=$1
BUILDTYPE=$2
CORES=$3

mkdir -p build

cd build

cmake .. -DCMAKE_CXX_COMPILER=${COMPILER} -DCMAKE_BUILD_TYPE=${BUILDTYPE}
#make install -j${CORES}

mingw32-make install -j${CORES}

cd ..

