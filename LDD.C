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
    double MAX_VALUE; 

    LDD_UPDATE_F(double *_S, double *_C, double _MAX_VALUE) : S(_S), C(_C), MAX_VALUE(_MAX_VALUE){}
    inline bool update(uintE u, uintE v) {
        if(C[2 * v + 1] == MAX_VALUE){
            C[2 * v] = minF<double>()(S[u], C[2 * v]);
        }
        return false;
    }
    inline bool updateAtomic(uintE u, uintE v) {
        if(C[2 * v + 1] == MAX_VALUE){
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
    inline bool cond(uintE u){ return true; }
};

struct CUT_EDGES_F {
    uintE* IDs;
    long &cut_edges;
    CUT_EDGES_F(uintE* _IDs, long& _cut_edges) : 
        IDs(_IDs), cut_edges(_cut_edges) {}

    inline bool update(uintE s, uintE d){ //Update function writes min ID
        if(IDs[s] != IDs[d] && s < d){
            cut_edges ++;
        }
        return false;
    }
    inline bool updateAtomic (uintE s, uintE d) { //atomic Update
        if(IDs[s] != IDs[d] && s < d){
            writeAdd(&cut_edges, 1l);
        }
        return false;
    }
    inline bool cond (uintE d) { return cond_true(d); } //does nothing
};

struct CC_Vertex_F {
    uintE* IDs;
    CC_Vertex_F(uintE* _IDs) :
        IDs(_IDs) {}
    inline bool operator () (uintE i) {
        if(IDs[i] == IDs[IDs[i]]){
            return false;
        } else{
            if(i == IDs[IDs[i]]){
                IDs[i] = min(IDs[i], i);
            } else{
                IDs[i] = IDs[IDs[i]];
            }
            return true;
        }
    }
};

void contract(long n, uintE* IDs){

    bool* frontier = newA(bool,n);
    {parallel_for(long i=0;i<n;i++) frontier[i] = 1;} 
    vertexSubset Frontier(n,n,frontier); //initial frontier contains all vertices
    while(!Frontier.isEmpty()){ //iterate until IDS converge
        // cout << Frontier.size() << endl;
        vertexSubset output = vertexFilter(Frontier,CC_Vertex_F(IDs));
        Frontier.del();
        Frontier = output;
    }
    Frontier.del();
}

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P){
    long n = GA.n;
    long m = GA.m;
    double beta = P.getOptionDoubleValue("-e", .005);
    int procs = P.getOptionIntValue("-p", 4);

    setWorkers(procs);
    double *E = newA(double, n);
    double *S = newA(double, n);
    double *C = newA(double, 2 * n);
    double MAX_VALUE = numeric_limits<double>::max(); 
    cout << "N: " << n << " M: " << m << " Num Workers: " << getWorkers() << endl;

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
        edgeMap(GA, F_new, LDD_UPDATE_F(S, C, MAX_VALUE));
        vertexSubset temp = edgeMap(GA, F_new, LDD_CHECK_F(S, C));
        F_new.del();
        F = temp;
        round ++;
    }
    // cout << "ROUND: " << round << endl;
    uintE *IDs = newA(uintE, n);
    {parallel_for(long i = 0; i < n; i++){
        if(C[2 * i + 1] == MAX_VALUE){
            IDs[i] = i;
        } else{
            IDs[i] = (uintE) C[2 * i + 1];
        }
    }}
    contract(n, IDs);
    // for(int i =0; i < n; i++){
    //     cout << "vertex: " << i << " cluster: " << IDs[i] << " START: " << S[i]<< endl;
    // }
    bool *uniq = newA(bool, n);
    bool *verts = newA(bool, n);
    {parallel_for(long i = 0; i < n; i++){
        uniq[i] = 0;
    }}
    {parallel_for(long i = 0; i < n; i++){
        uniq[IDs[i]] = 1;
        verts[i] = 1;
    }}
    vertexSubset clusters(n, uniq);
    long num_clusters = clusters.size();
    clusters.del();

    long cut_edges = 0;
    vertexSubset allVerts(n, verts);
    edgeMap(GA, allVerts, CUT_EDGES_F(IDs, cut_edges));

    cout << "Beta: " << beta << " # Clusters: " << num_clusters << " # cut-edges: " << cut_edges << endl; 
}

