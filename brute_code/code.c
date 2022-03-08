#include <iostream>
#include <chrono>

using namespace std;

static const unsigned int N = N_VAL, M = M_VAL, maxL = 60; // we assume M < N
static const unsigned int minL = (max(N, M + 1) + 1)/2*2;

typedef uint64_t code_t;

/*
 * print code in +- format
 */
void print(unsigned int L, code_t c){
    for(unsigned int i = 0; i < L; i++){
        if((c >> (L - i - 1)) & 1){
            cout << "-";
        }else{
            cout << "+";
        }
    }
    cout << endl;
}

/*
 * print code in python array format
 */
void pyprint(unsigned int L, code_t c){
    cout << "[";
    for(unsigned int i = 0; i < L; i++){
        if((c >> (L - i - 1)) & 1){
            cout << "-1";
        }else{
            cout << "+1";
        }
        if(i != L - 1){
            cout << ", ";
        }
    }
    cout << "]" << endl;
}

/*
 * Checks if the given code is valid. If yes, it is printed and the program is
 * exited.
 */
__attribute__((optimize("unroll-loops"), always_inline)) inline void check(unsigned int L, code_t c){
    code_t mask = (code_t(1) << L) - 1;
    c = (c << L) | c;
    // when checking the machine code make sure that these loops are unrolled
    // and that this function is inlined
    for(unsigned int delta = 1; delta <= M; delta++){ // 181 machine instructions for code check at N = 8 and M = 4
        code_t c0 = c ^ (c >> delta);
        for(unsigned int k = delta; k < N; k++){
            code_t c1 = c >> k;
            code_t c2 = c >> (delta + k);
            // the following machine instruction calculates the sum over the bitvector
            // if it is not supported, it must be implemented in software which
            // will come with a large performance impact
            if(__builtin_popcountll((c1 ^ c2 ^ c0) & mask) != L/2){
                return; // code not valid
            }
        }
    }
    // if we get here, the code is valid
    cout << "solution (L=" << L << "): ";
    print(L, c);
    pyprint(L, c);
    exit(0);
}

int main(){
    // loop over bitvector length
    for(unsigned int L = minL; L <= maxL; L += 2){
        cout << "L: " << L << endl;
        auto start = chrono::system_clock::now();
        // note that this parallelized loop acually checks only a quarter of the codes
        #pragma omp parallel for
        for(code_t c = code_t(1) << (L - 1); code_t(c) < (code_t(1) << L); c += 2){
            check(L, c);
        }
        auto end = chrono::system_clock::now();
        int elapsed = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        if(elapsed > 0){
            cout << "codes per sec " << 500*((code_t(1) << L) - (code_t(1) << (L - 1)))/elapsed << "e6" << endl;
        }
    }
}
