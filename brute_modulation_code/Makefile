
brute: brute.cpp
	g++ brute.cpp -Wall -march=native -O3 -o brute

brute.s: brute.cpp
	g++ brute.cpp -Wall -march=native -O3 -fdump-tree-optimized -S -masm=intel
