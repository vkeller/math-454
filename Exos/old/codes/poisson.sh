#!/bin/bash  
module purge
module load intel/15.0.0

icc poisson.c -o poisson

./poisson > poisson.out


for i in *.pgm
do
    convert $i `basename $i .pgm`.png
done

rm *.pgm

grep l2 poisson.out | sed 's/l2=//' > l2.out
