#!/bin/bash
find . -type f -name "*.txt" ! -name "config.txt" -exec rm -r {} +
rm -rf *.csv
mpirun -np 1 /media/pc-14-2/Data/LinearElasticity/github_chenghau/LinearElasticity/cmake-build-debug_2d/sle-kt

