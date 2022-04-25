#include <iostream>
#include <limits>

using namespace std;

typedef uint64_t code_t;
typedef int32_t code_len_t;

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

__attribute__((optimize("unroll-loops"))) int main(){
    for(code_len_t L = 1; L <= 48; L += 2){
        // we only need to check uneven length codes
        int best = numeric_limits<int>::max();
        code_t cc;
        for(code_t c = 0; c < (code_t(1) << L); c++){
restart_outer:
            for(code_len_t i = 1; i < L - 1; i += 2){
                if(2*__builtin_popcountll(((c >> i) ^ c) & ((code_t(1) << (L - i)) - 1)) != L - i){
                    c++;
                    goto restart_outer;
                }
            }
            int a = 0;
            for(code_len_t i = 2; i < L; i += 2){
                a += abs(2*int(__builtin_popcountll(((c >> i) ^ c) & ((code_t(1) << (L - i)) - 1))) - (L - i));
            }
            if(a < best){
                cc = c;
                best = a;
            }
        }
        if(best != numeric_limits<int>::max()){
            cout << "L=" << L << ", a=" << best - (L - 1)/2 << ": ";
            pyprint(L, cc);
        }else{
            cout << "no sol for " << L << endl;
        }
    }
    return 0;
}


