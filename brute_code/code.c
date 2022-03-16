#include <iostream>
#include <chrono>

using namespace std;

#ifndef N_VAL
#error "define N_VAL"
#endif
#ifndef M_VAL
#error "define M_VAL"
#endif
#ifndef MAX_L_VAL
#define "define MAX_L_VAL"
#endif
static const unsigned int N = N_VAL, M = M_VAL, maxL = MAX_L_VAL; // we assume M < N
#ifdef MIN_L_VAL
static const unsigned int minL = MIN_L_VAL;
#else
static const unsigned int minL = (max(N, M + 1) + 1)/2*2;
#endif

#ifndef BIT_128
typedef unsigned __int128 code_t;
#else
typedef uint64_t code_t;
#endif

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
    #pragma omp critical 
    {
        cout << "solution_text," << L << "," << N << "," << M << ",";
        print(L, c);
        cout << "solution_py," << L << "," << N << "," << M << ",";
        pyprint(L, c);
        exit(0);
    }
}

int main(){
#ifdef IMPL_ENUM
    cout << "running enumerate impl" << endl;
    // loop over bitvector length
    for(unsigned int L = minL; L <= maxL; L += 2){
        cout << "L: " << L << endl;
        // note that this parallelized loop acually checks only a quarter of the codes
        #pragma omp parallel for
        for(code_t c = code_t(1) << (L - 1); code_t(c) < (code_t(1) << L); c += 2){
            check(L, c);
        }
    }
#else
    typedef unsigned int uint;
    typedef unsigned int run_len_t;

    cout << "running anti-cyclic impl" << endl;
    for(uint L = minL; L <= maxL; L += 2){
        cout << "L: " << L << endl;
        code_t alternating = 0;
        for(uint i = 0; i < sizeof(code_t)*8; i += 2){
            alternating += (code_t(1) << i);
        }
        #pragma omp parallel for num_threads(16)
        for(uint leading = L - 1; leading >= 2; leading--){ // leading = 1 not supported
            run_len_t sol[L - 2];
            uint c = L - leading;
            code_t b = (code_t(1) << (leading - 1)) - 1; // start with leading sequence one bit shorter
            for(uint i = 0; i < L - leading; i++){
                sol[i] = 1;
            }
            uint runs = 0;
            b = b << (c + 1);
            c += 1;
            #pragma omp critical 
            {
                cout << "leading " << leading << endl;
            }
            while(1){
                b = ((b >> c) << c) | (alternating >> (sizeof(code_t)*8 - 1 - c + (runs&1)));
                runs += c - 2; // all added runs should have length 1
                if((runs&1) == 0){
                    check(L, b);
                }
                if(runs == 0){
                    break;
                }
                c = sol[runs];
                sol[runs] = 1;
                code_t prev = sol[runs - 1];
                while(prev == leading){
                    runs--;
                    if(runs == 0){
                        goto out;
                    }
                    c += prev;
                    sol[runs] = 1; // prev = 1
                    prev = sol[runs - 1];
                }
                sol[runs - 1] = prev + 1; // prev += 1
            }
            out:;
        }
    }
#endif
}
