#!/bin/bash

mkdir -p multi_run_logs

for n in {2..10}; do
    for ((m=1; m < n; m++)) do
        echo running $n $m
	    ./multi_builds/brute_${n}_${m} > multi_run_logs/log_${n}_${m}
    done
done
