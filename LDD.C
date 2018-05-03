#include <stdlib.h>
#include <math.h>
#include <limits>
#include "ligra.h"
#include "sequence.h"
#include "utils.h"
#include <iostream>
#include <omp.h>

using namespace std;
double exp(double lambda){
    double u = (double) rand() / (double) RAND_MAX;
    return -log(u) / lambda;
}

struct LDD_UPDATE_F{
    double *S;
    double *C;
    LDD_UPDATE_F(double *_S, double *_C) : S(_S), C(_C){}
    inline bool update(uintE u, uintE v) {
        if(C[2 * v + 1] ==  numeric_limits<double>::max()){
            C[2 * v] = minF<double>()(S[u], C[2 * v]);
        }
        return false;
    }
    inline bool updateAtomic(uintE u, uintE v) {
        if(C[2 * v + 1] ==  numeric_limits<double>::max()){
            writeMin(C + 2 * v, S[u]);
        }
        return false;
    }
    inline bool cond(uintE u){ return true; }
};

struct LDD_CHECK_F{
    double *S;
    double *C;
    LDD_CHECK_F(double *_S, double *_C) : S(_S), C(_C){}
    inline bool update(uintE u, uintE v) {
        if(C[2 * v] ==  S[u]){
            C[2 * v + 1] = u;
            return true;
        }
        return false;
    }
    inline bool updateAtomic(uintE u, uintE v) {
        if(C[2 * v] ==  S[u]){
            C[2 * v + 1] = u;
            return true;
        }
        return false;
    }
    inline bool cond(uintE u){ true; }
};


template <class vertex>
void Compute(graph<vertex>& GA, commandLine P){
    long n = GA.n;
    long m = GA.m;
    double beta = P.getOptionDoubleValue("-b", .005);
    setWorkers(1);
    double *E = newA(double, n);
    double *S = newA(double, n);
    double *C = newA(double, 2 * n);
    double MAX_VALUE = numeric_limits<double>::max                            (); 
    // cout << getWorkers() << endl;

    {parallel_for(long i = 0; i < n; i++) E[i] = exp(beta);}
    double maxE = sequence::reduce(E, n, maxF<double>());
    {parallel_for(long i = 0; i < n; i++) {
        S[i] = maxE - E[i]; 
    }}
    {parallel_for(long i = 0; i < 2 * n; i++) C[i] = MAX_VALUE;}
    
    long num_processed = 0;
    long round = 1;
    vertexSubset F(n);
    while(num_processed < n){
        // cout << "Round "  << round << endl << "Processed: " << num_processed << endl;
        F.toDense();
        
        bool *active = newA(bool, n);
        {parallel_for(long i = 0; i < n; i++) active[i] = F.isIn(i) || (S[i] < round && C[2 * i] == MAX_VALUE);}
        vertexSubset F_new(n, active);
        num_processed += F_new.size();
        F.del();
        edgeMap(GA, F_new, LDD_UPDATE_F(S, C));
        vertexSubset temp = edgeMap(GA, F_new, LDD_CHECK_F(S, C));
        // cout << "TEMP : " << &temp << endl;
        // cout << "F    : " << &F << endl;
        // cout << "F_new: " << &F_new << endl;
        F_new.del();
        F = temp;
        // cout << "TEMP : " << &temp << endl;
        // cout << "F    : " << &F << endl;
        // cout << "F_new: " << &F_new << endl;
        round ++;
    }
}

