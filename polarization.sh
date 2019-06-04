#!/bin/bash


#..............USAGE...............#
#./polarization.sh 0.3 "'../dump.efield_Z_0.3'"
#..............USAGE...............#

a="`grep -w "@" -c -v DP_$1.xvg `"
b=$((a))
dumpfile=$2

# In awk n is line number
awk 'NR==n{$1=a}1' n=7 a=$b polarization.inp > tmp && mv tmp polarization.inp
awk 'NR==n{$1=a}1' n=6 a=$1 polarization.inp > tmp && mv tmp polarization.inp
awk 'NR==n{$1=a}1' n=5 a=$2 polarization.inp > tmp && mv tmp polarization.inp

#lineNumber=3
#newValue=$b 
#sed -i ${lineNumber}'s/\S\+/'"${newValue}"'/1' movie.inp

gfortran polarization.f90 -o polarization.x
./polarization.x polarization.inp

exit
