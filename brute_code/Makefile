
brute: code.c
	g++ code.c -Wall -march=native -O3 -fopenmp -o brute -DN_VAL=8 -DM_VAL=6 -DMAX_L_VAL=36 -DMIN_L_VAL=36 -DIMPL_ENUM

code.s: code.c
	g++ code.c -Wall -S -masm=intel -march=native -fopenmp -O3

