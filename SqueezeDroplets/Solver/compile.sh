#! /bin/bash

rm ../squeeze_solver
ifort -O3 -heap-arrays *.f90 -o ../squeeze_solver
echo '***compile successfully -> squeeze_solver*^_^'
