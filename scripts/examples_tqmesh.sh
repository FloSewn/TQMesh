#!/bin/bash

set -e

./bin/run_examples 01
./bin/run_examples 02
./bin/run_examples 03
./bin/run_examples 04
./bin/run_examples 05
./bin/run_examples 06
./bin/run_examples 07
./bin/run_examples 08
./bin/run_examples 09

if [ $? -eq 0 ]; then
  echo ALL TESTS PASSED
else
  echo SOME TESTS FAILED
fi
