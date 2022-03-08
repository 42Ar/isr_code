#!/bin/bash

mkdir -p multi_builds

for n in {2..10}; do
    for ((m=1; m < n; m++)) do
        echo building $n $m
	    g++ code.c -Wall -march=native -fopenmp -O3 -o multi_builds/brute_${n}_${m} -DN_VAL=$n -DM_VAL=$m
    done
done
