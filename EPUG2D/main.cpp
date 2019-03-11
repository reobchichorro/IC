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

void read(vv& points, ve& map0, ve& map1, int& p0, int& p1) {
    ifstream m0("map0.txt"); ifstream m1("map1.txt"); 
    
    m0 >> p0; m1 >> p1;
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

void findCell(const ptt& x, const ptt& y, const ptt& cellSize, vector<pii>& cellOfPoint, const vtx& point, const int i) {
    int row = (int)( (point.getX() - x.first)/cellSize.first );
    int column = (int)( (point.getY() - y.first)/cellSize.second );
    cellOfPoint[i] = make_pair(row,column);

    cerr << "Point " << i << " and its cell: (" << point.getX() << "," << point.getY() << ") -> (" << row << "," << column << ")\n";
}   

void insertEdge(vector<vector<pair<ve,ve > > >& grid, const vector<pii>& cellOfPoint, const ve& map, const bool whichMap) {
    int initialPoint, finalPoint;
    pii x, y;
    
    for(int i=0; i<map.size(); i++) {
        initialPoint = map[i].getInit();
        finalPoint = map[i].getFin();

        x = make_pair(cellOfPoint[initialPoint].first,cellOfPoint[finalPoint].first);
        if(x.first > x.second) swap(x.first,x.second);

        y = make_pair(cellOfPoint[initialPoint].second,cellOfPoint[finalPoint].second);
        if(y.first > y.second) swap(y.first,y.second);

        cerr << "Edge " << i << ": Points " << initialPoint << " - " << finalPoint << " / (" << x.first << "," << y.first << ") - (" << x.second << "," << y.second << ")\n";

        for(int horizontal=x.first; horizontal<x.second; horizontal++) {
            for(int vertical=y.first; vertical<y.second; vertical++) {
                if(whichMap)
                    grid[horizontal][vertical].second.push_back(map[i]);
                else
                    grid[horizontal][vertical].first.push_back(map[i]);
            }    
        }
    }
}

int main() {
    int resolution = 10;
    int p0; int p1;
    vv points; ve map0, map1;

    read(points,map0,map1,p0,p1); //Reading points and edges of both maps

    vector<vector<pair<ve,ve > > > grid(resolution); //Creating uniform grid and setting its size
    for(int i=0; i<resolution; i++)
        grid[i].resize(resolution);

    vector<pii> cellOfPoint0(p0); // Creating vectors to store which cell a point is on the grid
    vector<pii> cellOfPoint1(p1);


    ptt x; getX(points,x); x.second += 0.5; // x.first is the lower limit of the grid's x-coordinate x.second is the upper limit.
    ptt y; getY(points,y); y.second += 0.5; // Same for y.

    cerr << "Limits of the grid:\n";
    cerr << x.first << " " << y.first << "\n";
    cerr << x.second << " " << y.second << "\n\n";
    
    // Calculating the cell size of the grid
    ptt cellSize; cellSize.first = (x.second - x.first)/resolution; cellSize.second = (y.second - y.first)/resolution;
    cerr << "Cell Size: " << cellSize.first << " x " << cellSize.second << "\n\n";

    // #pragma omp parallel for
    for(int i=0; i<p0; i++) {
        //determining cell for each point of map0
        findCell(x,y,cellSize,cellOfPoint0,points[i],i);
    }
    for(int i=0; i<p1; i++) {
        //determining cell for each point of map1
        findCell(x,y,cellSize,cellOfPoint1,points[i+p0],i);
    }

    cerr << "\nInserting edges \n";

    //Inserting edges into grid
    insertEdge(grid, cellOfPoint0, map0, 0);
    insertEdge(grid, cellOfPoint1, map1, 1);

    cerr << "\n";

    return 0;
}