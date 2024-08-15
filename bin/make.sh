#!/bin/bash

CURRENT_DIR="`pwd`"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
BUILD_DIR=${PROJECT_DIR}"/build_cmake"

mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# Make (and build, if necessary)
if [[ "$#" -eq 0 ]] ; then
    if [[ ! -f "Makefile" ]] ; then
        cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release ..
    fi
    make -j
fi

# Return to original directory
cd ${CURRENT_DIR}
