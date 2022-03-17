#include <iostream>
#include <chrono>
#include <omp.h>

using namespace std;

#define NUM_THREADS 6

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
        if((c >> (L - i - 1))&1){
            cout << "1";
        }else{
            cout << "0";
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
    omp_set_num_threads(NUM_THREADS);
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
#elif IMPL_RUN_STACK
    typedef unsigned int uint;
    typedef unsigned int run_len_t;

    // this algorithm checks all cyclic-unique combinations which start with 1 except 101010101...
    // note that by convention the constant sequence is not anti-cyclic
    cout << "running anti-cyclic impl" << endl;
#ifdef COUNT
    uint64_t counter = 0;
#endif
    code_t end, b, code, alternating = 0;
    for(uint i = 0; i < sizeof(code_t)*8; i += 2){
        alternating += (code_t(1) << i);
    }
    for(uint L = minL; L <= maxL; L += 2){
        cout << "L: " << L << endl;
        for(uint leading = L/2 - 1; leading >= 2; leading--){
            uint L_algo = L - leading;
            run_len_t sol[L_algo - 2];
            uint c = L_algo - leading;
            b = (code_t(1) << (leading - 1)) - 1; // start with leading sequence one bit shorter
            for(uint i = 0; i < L_algo - leading; i++){
                sol[i] = 1;
            }
            uint runs = 0;
            b = b << (c + 1);
            c += 1;
            cout << "leading " << leading << endl;
            while(1){
                b = ((b >> c) << c) | (alternating >> (sizeof(code_t)*8 - 1 - c + (runs&1)));
                runs += c - 2; // all added runs should have length 1
                c = sol[runs];
                if((runs&1) == 0){
                    code = code_t(1) << (c - 1);
                    if(code == 1){
                        code = 2;
                    }
                    end = (code_t(1) << leading) - 2;
                }else{
                    code = 0;
                    end = (code_t(1) << leading) - 1 - (code_t(1) << (c - 1));
                }
#ifdef COUNT
                counter += ((end - code) >> 1) + 1;
#endif
                code |= (b << leading);
                end |= (b << leading);
                for(code_t c = code; c <= end; c += 2){
                    check(L, c);
                }
                if(__builtin_expect(runs == 0, 0)){
                    break;
                }
                sol[runs] = 1;
                code_t prev = sol[runs - 1];
                while(__builtin_expect(prev == leading, 0)){
                    runs--;
                    if(__builtin_expect(runs == 0, 0)){
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
        cout << "remaining" << endl;
        code = ((code_t(1) << (L/2)) - 1) << L/2;
        end = (code_t(1) << L) - 2;
#ifdef COUNT
        counter += ((end - code) >> 1) + 1;
#endif
        //#pragma omp parallel for
        for(code_t c = code; c <= end; c += 2){
            check(L, c);
        }
#ifdef COUNT
        cout << "counter: " << counter << endl;
#endif
    }
#else
    typedef unsigned int uint;
    typedef unsigned int sub_code_t;

    // this algorithm checks all cyclic-unique combinations which start with 1 except 101010101...
    // note that by convention the constant sequence is not anti-cyclic
    cout << "running stacked loop impl" << endl;
#ifdef COUNT
    uint64_t counter = 0;
#endif
    for(uint L = minL; L <= maxL; L += 2){
        cout << "L: " << L << endl;
        #pragma omp parallel
        {
            code_t code, end_code, prev_code;
            int split_counter = 0;
            for(uint leading = L/2 - 1; leading >= 2; leading--){
                //cout << "leading " << leading << endl;
                uint remaining = L - leading - 1;
                uint seqs = remaining/leading;
                if(seqs == 1){
                    // can not split this one, so split here
                    if(split_counter % NUM_THREADS != omp_get_thread_num()){
                        split_counter += 1;
                        continue;
                    }
                    split_counter += 1;
                }
                // make work pieces as large as possible but not to large
                uint split_pos = (seqs - 1)*leading; // zero is not supported, must be a multiple of leading
                if(leading < 5 && split_pos > leading){
                    split_pos -= leading;
                }
                uint buffer_seq_len = remaining%leading;
                sub_code_t mask = (sub_code_t(1) << leading) - 1;
                code = mask;
                code = code << (1 + buffer_seq_len); // we start with leading sequence followed by the zero and the buffer seq
                sub_code_t buffer_end_exclusive = 1 << buffer_seq_len;
                for(sub_code_t buffer_seq = 0; buffer_seq < buffer_end_exclusive; buffer_seq++){
                    uint start = 0;
                    if(buffer_seq == 0){
                        start = sub_code_t(1) << buffer_seq_len; 
                    }else if((buffer_seq&1) == 0){
                        start = sub_code_t(1) << (__builtin_ctz(buffer_seq) - 1);
                    }
                    code = ((((code >> buffer_seq_len) << buffer_seq_len) | buffer_seq) << leading) | start;
                    uint cur_pos = leading*(seqs - 1);
                    uint goup = 0;
                    do{
                        prev_code = code >> leading;
                        sub_code_t end = mask;
                        if((prev_code&1) == 1){
                            end -= sub_code_t(1) << (__builtin_ctz(~sub_code_t(prev_code)) - 1);
                        }
                        if(__builtin_expect(cur_pos == 0, 0)){
                            if((code&1) == 1){
                                code += 1; // round up if uneven, as all final codes must end with 0
                            }
                            end_code = (prev_code << leading) | end; // might end with 1, but we are using <= in the loop and not !=
#ifdef COUNT
                            #pragma omp atomic
                            counter += ((end_code - code) >> 1) + 1;
#endif
                            for(code_t c = code; c <= end_code; c += 2){
                                check(L, c);
                            }
                            cur_pos += leading;
                            code = prev_code;
                            goup = 1;
                        }else{
increase:
                            if(__builtin_expect((sub_code_t(code)&mask) == end, 0)){
                                // goup must be 1 here
                                code = prev_code;
                                cur_pos += leading;
                            }else{
                                code += goup;
                                goup = 0;
                                if(__builtin_expect(cur_pos == split_pos, 0)){
                                    if(split_counter % NUM_THREADS != omp_get_thread_num()){
                                        goup = 1;
                                        split_counter += 1;
                                        goto increase;
                                    }
                                    split_counter += 1;
                                }
                                sub_code_t start = 0;
                                if((code&1) == 0){
                                    start = start | (sub_code_t(1) << (__builtin_ctz(sub_code_t(code)) - 1));
                                }
                                code = (code << leading) | start;
                                cur_pos -= leading;
                            }
                        }
                    }while(cur_pos < leading*seqs);
                }
            }
            code = ((code_t(1) << (L/2)) - 1) << L/2;
            end_code = (code_t(1) << L) - 2;
#ifdef COUNT
            if(omp_get_thread_num() == 0){
                #pragma omp atomic
                counter += ((end_code - code) >> 1) + 1;
            }
#endif
            #pragma omp for
            for(code_t c = code; c <= end_code; c += 2){
                check(L, c);
            }
        } // parallel section
    } // loop over L
#ifdef COUNT
    cout << "counter: " << counter << endl;
#endif
#endif
}
