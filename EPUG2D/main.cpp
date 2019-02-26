#include <iostream>
#include <vector>
#include <fstream>
#include "chain.h"
#include <omp.h>
using namespace std;

typedef double T;
typedef Vertex<T> vtx;
typedef Edge<T> E;
typedef vector<vtx> vv;
typedef vector<E> ve;
typedef pair<T,T> ptt;
typedef pair<int,int> pii;

void read(vv& points, ve& map0, ve& map1) {
    ifstream m0("map0.txt"); ifstream m1("map1.txt"); 
    
    int p0; m0 >> p0;
    int p1; m1 >> p1;
    points.resize(p0+p1);
    for(int i=0; i<p0; i++)
        m0 >> points[i];
    for(int i=0; i<p1; i++)
        m1 >> points[p0+i];

    int n;
    m0>>n;
    map0.resize(n);
    for(int i=0; i<n; i++)
        m0 >> map0[i];

    m1>>n;
    map1.resize(n);
    for(int i=0; i<n; i++)
        m1 >> map1[i];
}

void getX(const vv& points, ptt& x) {
    x.first = points[0].getX();
    x.second = points[0].getX();
    //#pragma omp parallel for (posso paralelizar isso? se sim, vale a pena?)
    for(int i=1; i<points.size(); i++) {
        if(points[i].getX() < x.first)
            x.first = points[i].getX();
        else if(points[i].getX() > x.second)
            x.second = points[i].getX();
    }
}

void getY(const vv& points, ptt& y) {
    y.first = points[0].getY();
    y.second = points[0].getY();
    //#pragma omp parallel for (posso paralelizar isso? se sim, vale a pena?)
    for(int i=1; i<points.size(); i++) {
        if(points[i].getY() < y.first)
            y.first = points[i].getY();
        else if(points[i].getY() > y.second)
            y.second = points[i].getY();
    }
}

void findCell(const ptt& x, const ptt& y, const ptt& cellSize, vector<vector<vv> >& pointGrid, const vtx& point) {
    int row = (int)( (point.getX() - x.first)/cellSize.first );
    int column = (int)( (point.getY() - y.first)/cellSize.second );
    pointGrid[row][column].push_back(point);
    cerr << point.getX() << " " << point.getY() << " - " << row << " " << column << "\n";
}   

int main() {
    int resolution = 10;
    vv points; ve map0, map1;
    read(points,map0,map1);
    //

    // vector<vector<pair<ve,ve > > > grid(resolution);
    vector<vector<vv> > pointGrid(resolution);
    for(int i=0; i<resolution; i++)
        pointGrid[i].resize(resolution);
    // for(int i=0; i<resolution; i++)
    //     grid[i].resize(resolution);

    ptt x; getX(points,x); x.second += 0.5; // x.first e' o limite inferior da coordenada x da grade e x.second e' o limite superior.
    ptt y; getY(points,y); y.second += 0.5; // Analogo para y.

    cerr << x.first << " " << y.first << "\n";
    cerr << x.second << " " << y.second << "\n";
    ptt cellSize; cellSize.first = (x.second - x.first)/resolution; cellSize.second = (y.second - y.first)/resolution;
    cerr << cellSize.first << " x " << cellSize.second << "\n\n";

    // #pragma omp parallel for
    for(int i=0; i<points.size(); i++) {
        //determine cell for each point
        findCell(x,y,cellSize,pointGrid,points[i]);
    }

    // for(int i=0; i<map0.size(); i++) {
    //     //grid[cell0[i].first][cell0[i].second].first.push_back()
    //     //Insert edges of map0 in grid
    // }
    cerr << "\n";

    // for(int i=0; i<map1.size(); i++) {
    //     //grid[cell1[i].first][cell1[i].second].second.push_back()
    //     //Insert edges of map1 in grid
    // }

    return 0;
}