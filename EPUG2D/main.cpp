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
    for(int i=0; i<p0; i++) {
        m0 >> points[i];
        points[i].setOriginLabels(-1,-1);
    }
    for(int i=0; i<p1; i++) {
        m1 >> points[p0+i];
        points[p0+i].setLabel(p0+i);
        points[i].setOriginLabels(-2,-2);
    }

    int n;
    m0>>n;
    map0.resize(n);
    for(int i=0; i<n; i++)
        m0 >> map0[i];

    m1>>n;
    map1.resize(n);
    for(int i=0; i<n; i++) {
        m1 >> map1[i];
        map1[i].setVs(map1[i].getInit()+p0, map1[i].getFin()+p0);
    }
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

        for(int horizontal=x.first; horizontal<=x.second; horizontal++) {
            for(int vertical=y.first; vertical<=y.second; vertical++) {
                //cerr << horizontal << " " << vertical << "\n";
                if(whichMap)
                    grid[horizontal][vertical].second.push_back(map[i]);
                else
                    grid[horizontal][vertical].first.push_back(map[i]);
            }    
        }
    }
}

int orientation(const vtx& p, const vtx& q, const vtx& r) {
    double determinant = (p.getX() - q.getX()) * (q.getY() - r.getY()) - (q.getX() - r.getX()) * (p.getY() - q.getY());
    if(determinant == 0.0) return 0;
    return (determinant > 0)? 1: -1;
}

//Checks if two edges intersect. Does not include degenerate cases (handled by SoS).
int checkIntersection(const vv& points, const E& a, const E& b) {
    cerr << "in function\n";
    if(
        (orientation(points[a.getInit()], points[a.getFin()], points[b.getInit()]) != orientation(points[a.getInit()], points[a.getFin()], points[b.getFin()]) )
        &&
        (orientation(points[b.getInit()], points[b.getFin()], points[a.getInit()]) != orientation(points[b.getInit()], points[b.getFin()], points[a.getFin()]) )
    ) {
        cerr << "Edges " << a.getLabel() << " and " << b.getLabel() << " intersect.\n";
        return 1
    }
    else {
        cerr << "Edges " << a.getLabel() << " and " << b.getLabel() << " do not intersect.\n";
        return -1;
    }
}

//Returns the intersection point of the LINES, but doesn't verify if this point in actually IN THE EDGES (could be outside).
void edgeIntersection(vtx& intersection, const E& e0, const E& e1, const vv& points) {
    T a=points[e0.getInit()].getX(), b=points[e0.getInit()].getY(), c=points[e0.getFin()].getX(), d=points[e0.getFin()].getY();
    T e=points[e1.getInit()].getX(), f=points[e1.getInit()].getY(), g=points[e1.getFin()].getX(), h=points[e1.getFin()].getY();
    T bigFractionNumerator   = ( (g-e)*(f-b) - (h-f)*(e-a) );
    T bigFractionDenominator = ( (g-e)*(d-b) - (c-a)*(h-f) );
    cerr << a << "," << b << " - " << c << "," << d << "\n" << e << "," << f << " " << g << "," << h << "\n";
    if(bigFractionNumerator == 0) cerr << "Numerator 0\n";
    if(bigFractionDenominator == 0) cerr << "Denominator 0 - parallel lines(I believe), SoS handles it so it never happens\n";
    else bigFractionNumerator /= bigFractionDenominator;
    //cerr << "bigFraction: " << bigFractionNumerator << "/" << bigFractionDenominator <<  "\n";

    intersection = vtx(points.size(),(c-a)*bigFractionNumerator+a,(d-b)*bigFractionNumerator+b,e0.getLabel(),e1.getLabel());
    return;
}

/*
"Pseudocode" for the Point Location function:
Given a query vertex P of map A, find which region/face of map B contains P.

for j = yOfCellOfP to  j <= maxY
    for each edge e of map B in cell [xOfCellOfP][j]:
        if a semi-infinite vertical ray from P hits edge E, then P is inside the region/face that the ray eas when it hit E
        (how to make this check? i'm guessing something using the x-coordinates of the bounding points of E)
        (also, in the first cell (the one that contains P), we have to check the y-coordinates of E *, as E could be BELOW P, 
        but as it is in the same cell, it is being checked, but it shouldn't)
        * doing this check: if the y-coordinate of both points of the edge are higher than P's y-coordinate, 
        then we should check if the ray hits the edge;
        if the y-coordinate of both points is below P, than we shouldn't check if the ray hits the edge.
        if one point of the edge is above P and the other is below, then check P's orientation in relation to the edge.

*/

int main() {
    int resolution = 10;
    int p0; int p1;
    vv points; ve map0, map1;

    read(points,map0,map1,p0,p1); //Reading points and edges of both maps

    vector<vector<pair<ve,ve > > > grid(resolution); //Creating uniform grid and setting its size
    for(int i=0; i<resolution; i++)
        grid[i].resize(resolution);

    vector<pii> cellOfPoint(p0+p1); // Creating a vector to store which cell a point is on the grid
    //vector<pii> cellOfPoint1(p1);


    ptt x; getX(points,x); x.second += 0.5; // x.first is the lower limit of the grid's x-coordinate x.second is the upper limit.
    ptt y; getY(points,y); y.second += 0.5; // Same for y.

    cerr << "Limits of the grid:\n";
    cerr << x.first << " " << y.first << "\n";
    cerr << x.second << " " << y.second << "\n\n";
    
    // Calculating the cell size of the grid
    ptt cellSize; cellSize.first = (x.second - x.first)/resolution; cellSize.second = (y.second - y.first)/resolution;
    cerr << "Cell Size: " << cellSize.first << " x " << cellSize.second << "\n\n";

    // #pragma omp parallel for
    for(int i=0; i<p0+p1; i++) {
        //determining cell for each point
        findCell(x,y,cellSize,cellOfPoint,points[i],i);
    }

    cerr << "\nInserting edges \n";

    //Inserting edges into grid
    insertEdge(grid, cellOfPoint, map0, 0);
    insertEdge(grid, cellOfPoint, map1, 1);

    cerr << "\n";

    vtx intersectionPoint;
    for(int i=0; i<grid.size(); i++) {
        for(int j=0; j<grid[i].size(); j++) {
            for(int k=0; k<grid[i][j].first.size(); k++) {
                for(int l=0; l<grid[i][j].second.size(); l++) {
                    if(checkIntersection(points, grid[i][j].first[k], grid[i][j].second[l]) == 1) {
                        edgeIntersection(intersectionPoint,grid[i][j].first[k],grid[i][j].second[l],points);
                        points.push_back(intersectionPoint);
                    }
                }
            }
        }
    }

    //Testing case for intersection function
    checkIntersection(points, map0[0], map0[1]);
    edgeIntersection(intersectionPoint,map0[0],map0[1],points); //Edge 0's dinal point is the same as edge 1's initial point.
    cout << "Intersection: (" << intersectionPoint.getX() << "," << intersectionPoint.getY() << ")\n"; //Thus, cout should print this point (point 1 in the points vector).

    return 0;
}