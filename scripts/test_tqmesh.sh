#!/bin/bash

cd build
ctest
cd ..

if [ $? -eq 0 ]; then
  echo ALL TESTS PASSED
else
  echo SOME TESTS FAILED
fi
