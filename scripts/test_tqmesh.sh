#!/bin/bash

set -e

<<<<<<< HEAD
cd build
ctest
=======
#cd build
#ctest

./bin/run_tests Vertex
./bin/run_tests Triangle
./bin/run_tests Quad
./bin/run_tests Front
./bin/run_tests EdgeList
./bin/run_tests Boundary
./bin/run_tests SizeFunction
./bin/run_tests Mesh
./bin/run_tests MeshGenerator
./bin/run_tests MeshCleanup
./bin/run_tests SmoothingStrategy

#./../bin/TQMesh ../input/Example_1.para
#./../bin/TQMesh ../input/Example_2.para
#./../bin/TQMesh ../input/Example_3.para
#./../bin/TQMesh ../input/Example_4.para
#./../bin/TQMesh ../input/Example_5.para

#./../bin/run_examples 1
#./../bin/run_examples 2
#./../bin/run_examples 3
#./../bin/run_examples 4
#./../bin/run_examples 5

>>>>>>> fa0899f5faedbc3de2d30dba4c8c9fc7b7288940
cd ..

if [ $? -eq 0 ]; then
  echo ALL TESTS PASSED
else
  echo SOME TESTS FAILED
fi
