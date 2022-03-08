
brute: code.c
	g++ code.c -Wall -march=native -fopenmp -O3 -o brute -DN_VAL=5 -DM_VAL=4

code.s: code.c
	g++ code.c -Wall -S -masm=intel -march=native -fopenmp -O3 -DN_VAL=5 -DM_VAL=4

