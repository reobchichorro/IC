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

    vector<pair<int,int> > cell0(map0.size());
    #pragma omp parallel for
    for(int i=0; i<map0.size(); i++) {
        //determine cell for each edge of map0
    }

    for(int i=0; i<map0.size(); i++) {
        //grid[cell0[i].first][cell0[i].second].first.push_back()
        //Insert edges of map0 in grid
    }

    vector<pair<int,int> > cell1(map1.size());
    #pragma omp parallel for
    for(int i=0; i<map1.size(); i++) {
        //determine cell for each edge of map1
    }

    for(int i=0; i<map1.size(); i++) {
        //grid[cell1[i].first][cell1[i].second].second.push_back()
        //Insert edges of map1 in grid
    }

    return 0;
}