#! /bin/bash
rm ../ProgramForLinux/preprocess
ifort -O3 -heap-arrays *.f90 -o ../ProgramForLinux/preprocess
echo '***compile successfully -> preprocess*^_^'