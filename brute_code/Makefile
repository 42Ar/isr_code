
brute: code.c
	g++ code.c -Wall -march=native -fopenmp -O3 -o brute

code.s: code.c
	g++ code.c -Wall -S -masm=intel -march=native -fopenmp -O3

