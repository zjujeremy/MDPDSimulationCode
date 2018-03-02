#! /bin/bash

rm ../../preprocess
ifort -O3 -heap-arrays *.f90 -o ../../preprocess
echo '***compile successfully -> preprocess*^_^'