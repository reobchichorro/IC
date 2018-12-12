#include <iostream>
#include <vector>
#include <fstream>
#include "chain.h"
#include <omp.h>
using namespace std;

typedef double T;
typedef Edge<T> E;
typedef vector<E> ve;
typedef pair<T,T> ptt;
typedef pair<int,int> pii;

void getX(const ve& map0, const ve& map1, ptt& x){
    x.first = map0[0].getInit().getX();
    x.second = map0[0].getInit().getX();
    //#pragma omp parallel for (posso paralelizar isso? se sim, vale a pena?)
    for(int i=1; i<map0.size(); i++) {
        if(map0[i].getInit().getX() < x.first)
            x.first = map0[i].getInit().getX();
        else if(map0[i].getInit().getX() > x.second)
            x.second = map0[i].getInit().getX();
    }
    for(int i=0; i<map1.size(); i++) {
        if(map1[i].getInit().getX() < x.first)
            x.first = map1[i].getInit().getX();
        else if(map1[i].getInit().getX() > x.second)
            x.second = map1[i].getInit().getX();
    }
}

void getY(const ve& map0, const ve& map1, ptt& y){
    y.first = map0[0].getInit().getY();
    y.second = map0[0].getInit().getY();
    //#pragma omp parallel for (posso paralelizar isso? se sim, vale a pena?)
    for(int i=1; i<map0.size(); i++) {
        if(map0[i].getInit().getY() < y.first)
            y.first = map0[i].getInit().getY();
        else if(map0[i].getInit().getY() > y.second)
            y.second = map0[i].getInit().getY();
    }
    for(int i=0; i<map1.size(); i++) {
        if(map1[i].getInit().getY() < y.first)
            y.first = map1[i].getInit().getY();
        else if(map1[i].getInit().getY() > y.second)
            y.second = map1[i].getInit().getY();
    }
}

void findCell(const ptt& x, const ptt& y, const ptt& cellSize, pii& cell, const int resolution, const E& edgeFromMap) {
    for(int i=0; i<resolution; i++)
        for(int j=0; j<resolution; j++) {
            cell = make_pair(i,j);
            if(edgeFromMap.getInit().getX() - (x.first + i*cellSize.first) > cellSize.first)
                break;
            if(edgeFromMap.getInit().getY() - (y.first + j*cellSize.second) <= cellSize.second)
                return;
        }
}   // Uma versao usando busca binaria seria melhor

int main() {
    int resolution = 10;
    ifstream m0("map0.txt"); ifstream m1("map1.txt"); 
    
    int n; m0>>n;
    ve map0(n);
    for(int i=0; i<n; i++)
        m0 >> map0[i];

    m1>>n;
    ve map1(n);
    for(int i=0; i<n; i++)
        m1 >> map1[i];

    vector<vector<pair<ve,ve > > > grid(resolution);
    for(int i=0; i<resolution; i++)
        grid[i].resize(resolution);

    ptt x; getX(map0,map1,x); x.second += 0.5; // x.first e' o limite inferior da coordenada x da grade e x.second e' o limite superior.
    ptt y; getY(map0,map1,y); y.second += 0.5; // Analogo para y.

    cerr << x.first << " " << y.first << "\n";
    cerr << x.second << " " << y.second << "\n";
    vector<pii> cell0(map0.size());
    ptt cellSize; cellSize.first = (x.second - x.first)/resolution; cellSize.second = (y.second - y.first)/resolution;
    cerr << cellSize.first << " x " << cellSize.second << "\n\n";

    // #pragma omp parallel for
    for(int i=0; i<map0.size(); i++) {
        //determine cell for each edge of map0
        findCell(x,y,cellSize,cell0[i],resolution,map0[i]);
        cerr << map0[i].getInit().getX() << " " << map0[i].getInit().getY() << " - " << cell0[i].first << " " << cell0[i].second << "\n";
    }

    // for(int i=0; i<map0.size(); i++) {
    //     //grid[cell0[i].first][cell0[i].second].first.push_back()
    //     //Insert edges of map0 in grid
    // }
    cerr << "\n";
    vector<pii> cell1(map1.size());
    // #pragma omp parallel for
    for(int i=0; i<map1.size(); i++) {
        //determine cell for each edge of map1
        findCell(x,y,cellSize,cell1[i],resolution,map1[i]);
        cerr << map1[i].getInit().getX() << " " << map1[i].getInit().getY() << " - " << cell1[i].first << " " << cell1[i].second << "\n";
    }

    // for(int i=0; i<map1.size(); i++) {
    //     //grid[cell1[i].first][cell1[i].second].second.push_back()
    //     //Insert edges of map1 in grid
    // }

    return 0;
}