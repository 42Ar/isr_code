
brute: code.c
	g++ code.c -Wall -march=native -O3 -fopenmp -o brute -DN_VAL=9 -DM_VAL=7 -DMAX_L_VAL=56 -DMIN_L_VAL=50 -DBIT_128 #-DIMPL_ENUM #-DIMPL_RUN_STACK

code.s: code.c
	g++ code.c -Wall -fdump-tree-optimized -S -masm=intel -march=native -fopenmp -g3 -DN_VAL=8 -DM_VAL=6 -DMAX_L_VAL=34 -DMIN_L_VAL=34 #-DIMPL_ENUM #-DIMPL_RUN_STACK
