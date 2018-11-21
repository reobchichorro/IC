#include <iostream>
#include <vector>
#include "chain.h"
#include <omp.h>
using namespace std;

typedef double T;

int main() {
    int resolution = 50;

    int n; cin>>n;
    vector<Edge<T> > map0(n);
    for(int i=0; i<n; i++)
        cin >> map0[i];

    cin>>n;
    vector<Edge<T> > map1(n);
    for(int i=0; i<n; i++)
        cin >> map1[i];

    vector<vector<pair<vector<Edge<T> >,vector<Edge<T> > > > > grid(resolution);
    for(int i=0; i<resolution; i++)
        grid[i].resize(resolution);

    #pragma omp parallel for
    for(int i=0; i<map0.size(); i++) {
        
    }

    return 0;
}