
brute: code.c
	g++ code.c -Wall -march=native -O3 -fopenmp -o brute -DN_VAL=8 -DM_VAL=6 -DMAX_L_VAL=34 -DMIN_L_VAL=34 #-DIMPL_ENUM #-DIMPL_RUN_STACK

code.s: code.c
	g++ code.c -Wall -S -masm=intel -march=native -fopenmp -O3 -DN_VAL=8 -DM_VAL=6 -DMAX_L_VAL=36 -DMIN_L_VAL=36 -DCOUNT -DIMPL_ENUM

