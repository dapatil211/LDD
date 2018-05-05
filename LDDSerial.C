#include <stdlib.h>
#include <math.h>
#include <limits>
#include "ligra.h"
#include "sequence.h"
#include "utils.h"
#include <iostream>
#include <omp.h>
#include <vector>
#include<algorithm>
#include<unordered_set>
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

    inline bool update(uintE s, uintE d){ 
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
        vertexSubset output = vertexFilter(Frontier,CC_Vertex_F(IDs));
        Frontier.del();
        Frontier = output;
    }
    Frontier.del();
}

template <class vertex>
struct VOLUME_F{
    long &degree_sum;
    graph<vertex>& GA;
    VOLUME_F(long &_sum, graph<vertex>& _GA):
        degree_sum(_sum), GA(_GA) {}
    
    void operator () (uintE i){
        degree_sum += GA.V[i].getOutDegree();
    }
};

struct SET_F{
    unordered_set<uintE> &remove;
    SET_F(unordered_set<uintE> &_remove):
        remove(_remove) {}
    
    void operator () (uintE i){
        !remove.count(i);
    }
};
struct BOUNDARY_F{
    uintE* balls;
    long &outGoing;
    BOUNDARY_F(long& _outGoing, uintE* _balls):
        balls(_balls), outGoing(_outGoing){}
    
    inline bool update(uintE s, uintE d){ 
        outGoing++;
        return false;
    }

    inline bool updateAtomic (uintE s, uintE d) { //atomic Update
        return update(s, d);
    }

    inline bool cond (uintE d) { return balls[d] == -1; } //does nothing
};

struct BFS_F{
    uintE *balls;
    uintE ballCenter;
    // long &numInBall;
    unordered_set<uintE> &verticesInBall;
    BFS_F(uintE* _balls, uintE _ballCenter, unordered_set<uintE>& _verticesInBall):
        balls(_balls), ballCenter(_ballCenter), verticesInBall(_verticesInBall){}
    
    inline bool update(uintE s, uintE d){ 
        balls[d] = ballCenter;
        // numInBall++;
        verticesInBall.erase(d);
        return true;
    }

    inline bool updateAtomic (uintE s, uintE d) { //atomic Update
        return update(s, d);
    }

    inline bool cond (uintE d) { return balls[d] == -1; } //does nothing
};
timer t1, t2, t3, t4;
template <class vertex>
void growBall(graph<vertex>& GA, uintE* balls, uintE ballCenter, double beta, unordered_set<uintE>& all){
    vertexSubset F(GA.n, ballCenter);
    balls[ballCenter] = ballCenter;
    // set<uintE>verticesInBall(1, ballCenter);
    all.erase(ballCenter);
    long degrees = 0;
    long outGoing = 0;
    // long numInBall = 1;
    vertexMap(F, VOLUME_F<vertex>(degrees, GA));
    edgeMap(GA, F, BOUNDARY_F(outGoing, balls));
    // cout << "CENTER: " << ballCenter << " OUTGOING: " << outGoing << " DEGREES: "<< degrees << endl;
    while(outGoing > beta * degrees){
        // t1.start();
        vertexSubset output = edgeMap(GA, F, BFS_F(balls, ballCenter, all));
        // t1.stop();
        F.del();
        F = output;
        long frontierDegrees = 0;
        outGoing = 0;
        // t2.start();
        vertexMap(F, VOLUME_F<vertex>(frontierDegrees, GA));
        // t2.stop();
        // t3.start();
        edgeMap(GA, F, BOUNDARY_F(outGoing, balls));
        // t3.stop();
        degrees += frontierDegrees;
        // cout << "CENTER: " << ballCenter << " OUTGOING: " << outGoing << " DEGREES: "<< degrees << endl;
    }
    F.del();
    // return verticesInBall;
}

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P){
    long n = GA.n;
    long m = GA.m;
    double beta = P.getOptionDoubleValue("-e", .005);

    setWorkers(1);
    cout << endl << "SEQUENTIAL LDD" << endl;
    cout << "N: " << n << " M: " << m  << endl;
    timer t;
    t.start();
    uintE *balls = newA(uintE, n);
    for(int i = 0; i < n; i ++){
        balls[i] = -1;
    }

    // long numActive = n;
    // uintE *active = newA(uintE, n);

    unordered_set<uintE> all;
    for(long i = 0; i < n; i++){
        // active[i] = 1;
        all.insert(i);
    }

    // vertexSubset F(n, active);
    while(!all.empty()){
        // cout << all.size() << endl;
        // int ballCenterIdx = rand() % numActive;
        // uintE ballCenter = 0;
        // while(ballCenterIdx != -1){
        //     if(balls[ballCenter] == -1){
        //         ballCenterIdx -= 1;
        //     }
        //     ballCenter ++;
        // }
        // F.toSparse();
        // uintE ballCenter = F.vtx[ballCenterIdx];
        // cout << ballCenterIdx << " " << numActive << " " << ballCenter << endl;
        // t4.start();
        auto it = all.begin();
        // t4.stop();
        // cout << *it << endl;
        growBall(GA, balls, *it, beta, all);
        // cout << verticesInBall.size()<<endl;
        // for(int i = 0; i < n; i ++){
        //     cout << active[i] << " ";
        // }
        // sort(verticesInBall.rbegin(), verticesInBall.rend());
        // cout << verticesInBall.size() << endl << endl;
        // long ind = numActive - 1;
        // for(long i = 0; i < verticesInBall.size(); i++){
        //     while(active[ind] != verticesInBall[i]){
        //         ind --;
        //     }
        //     active[ind] = active[numActive - 1];
        //     numActive --;
        //     // cout << verticesInBall[i] << endl;
        // }
        // sort(active, active + numActive);
        // for(int i = 0; i < n; i ++){
        //     cout << active[i] << " ";
        // }
        // cout << endl;
    }
    t.stop();
    // cout << t1.total() << endl;
    // cout << t2.total() << endl;
    // cout << t3.total() << endl;
    // cout << t4.total() << endl;

    bool *uniq = newA(bool, n);
    bool *verts = newA(bool, n);
    for(long i = 0; i < n; i++){
        uniq[i] = 0;
        // cout << balls[i] << endl;
    }

    for(long i = 0; i < n; i++){
        uniq[balls[i]] = 1;
        verts[i] = 1;
    }
    vertexSubset clusters(n, uniq);
    long num_clusters = clusters.size();
    clusters.del();
    
    long cut_edges = 0;
    vertexSubset allVerts(n, verts);
    edgeMap(GA, allVerts, CUT_EDGES_F(balls, cut_edges));
    allVerts.del();

    cout << "Beta: " << beta << " # Clusters: " << num_clusters << " # cut-edges: " << cut_edges << endl; 
    t.reportTotal("Sequential Computation Time: "); 
    free(balls);
}

