#!/bin/bash

set -e

#cd build
#ctest

./bin/run_tests Vertex
./bin/run_tests Triangle
./bin/run_tests Quad
./bin/run_tests Front
./bin/run_tests EdgeList
./bin/run_tests Boundary
./bin/run_tests SizeFunction
./bin/run_tests SmoothingStrategy
./bin/run_tests Mesh
./bin/run_tests MeshGenerator
./bin/run_tests MeshCleanup
./bin/run_tests ParaReader
./bin/run_tests MeshChecker

#cd ..

if [ $? -eq 0 ]; then
  echo ALL TESTS PASSED
else
  echo SOME TESTS FAILED
fi
