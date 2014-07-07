#!/bin/bash

set -e

cd $(dirname $0);

TEST_EXE=$(ls ../ext/_result/*/*/bin/signature_contest | head -n 1)
ECC_LIBRARY="../build/libgost_ecc.so"

exec $TEST_EXE "${ECC_LIBRARY}" "test2.dat"
