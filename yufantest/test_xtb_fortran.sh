#!/bin/bash

# Check number of arguments
if [ $# -lt 5 ]; then
    echo "Usage: $0 <input_file> <xtb_path> <--output> <--fe|--vib> <executable>"
    echo "  <input_file>    : Input coordinate file"
    echo "  <xtb_path>      : Path to xtb installation"
    echo "  <--output>      : Optional flag to show output"
    echo "  <--fe|--vib>    : Calculation type (free energy or vibrational)"
    echo "  <executable>    : Path to xtb executable"
    exit 1
fi

# Set the XTBPATH environment variable
export XTBPATH=$2

# Obtain name without extension for $1
filename=$(basename -- "$1")
# Remove filename.engrad if it exists
rm -f $filename.engrad
rm -f vibspectrum
exec=$5

# Set default parallel number to 4 if not provided
parallel=${6:-4}

# Check if the flag output is enabled
if [ "$3" == "--output" ]; then
    if [ "$4" == "--fe" ]; then
        $5 $1 --sp --grad --gfn 2 --parallel $parallel --norestart --acc 0.0001 --iterations 1000
    fi
    if [ "$4" == "--vib" ]; then
        $5 $1 --hess --gfn 2 --parallel $parallel --norestart --acc 0.0001 --iterations 1000
    fi
    # if [ "$4" == "--nci" ]; then
    # fi
else
    if [ "$4" == "--fe" ]; then
        $5 $1 --sp --grad --gfn 2 --parallel $parallel --norestart --acc 0.0001 --iterations 1000 > /dev/null 2>&1
    fi
    if [ "$4" == "--vib" ]; then
        $5 $1 --hess --gfn 2 --parallel $parallel --norestart --acc 0.0001 --iterations 1000 > /dev/null 2>&1
    fi
    # if [ "$4" == "--nci" ]; then
    # fi
fi

