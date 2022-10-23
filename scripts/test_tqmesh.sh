#!/bin/bash

set -e

cd build
ctest

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

cd ..

if [ $? -eq 0 ]; then
  echo ALL TESTS PASSED
else
  echo SOME TESTS FAILED
fi
