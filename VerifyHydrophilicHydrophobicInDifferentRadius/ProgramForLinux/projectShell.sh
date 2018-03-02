#! /bin/bash
rm ./preprocess
rm ./solver

ifort -O3 -heap-arrays ../PreProcess/*.f90 -o ./preprocess
echo '***compile successfully -> preprocess*^_^'
ifort -O3 -heap-arrays ../Solver/*.f90 -o ./solver
echo '***compile successfully -> squeeze_solver*^_^'

chmod a+x ./preprocess ./solver
echo '--------------start preprocess--------------'
./preprocess
echo '--------------start solver--------------'
./solver