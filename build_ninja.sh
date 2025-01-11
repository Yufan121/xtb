#!/bin/bash

# Check if the first argument is 'debug'
if [ "$1" == "debug" ]; then
    BUILD_TYPE="Debug"
    C_FLAGS="-g"
    FORTRAN_FLAGS="-g"
else
    BUILD_TYPE="Release"
    C_FLAGS="-O3"
    FORTRAN_FLAGS="-O3"
fi

# Load the necessary module and run the build command
module load cmake/3.27.7
cmake -B_build -S. -GNinja -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_C_FLAGS_$BUILD_TYPE="$C_FLAGS" -DCMAKE_Fortran_FLAGS_$BUILD_TYPE="$FORTRAN_FLAGS" -DCMAKE_INSTALL_PREFIX=/scratch/pawsey0799/yx7184/xtb_install
ninja -C _build
ninja -C _build install

