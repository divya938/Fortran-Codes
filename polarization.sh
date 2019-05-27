#!/bin/bash

a="`grep -w "@" -c -v DP_along_stack1.xvg `"
b=$((a))

# In awk n is line number
awk 'NR==n{$1=a}1' n=7 a=$b polarization.inp > tmp && mv tmp polarization.inp

#lineNumber=3
#newValue=$b 
#sed -i ${lineNumber}'s/\S\+/'"${newValue}"'/1' movie.inp

gfortran polarization.f90 -o polarization.x
./polarization.x polarization.inp

exit
