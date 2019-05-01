#include <bits/stdc++.h>
using namespace std;

struct Vertex {
    double x;
    double y;
};

bool operator<(const Vertex& a, const Vertex& b) {
    if(a.x < b.x)
        return true;
    return a.y < b.y;
}

struct Edge1 {
    int label;
    Vertex p0;
    Vertex p1;
    int fPos,fNeg;
};

bool operator<(const Edge1& a, const Edge1& b) {
    if(a.p0 < b.p0)
        return true;
    else if (a.p1 < b.p1)
        return true;
    else if(a.fPos < b.fPos)
        return true;
    return a.fNeg < b.fNeg;
}

struct Edge2 {
    int label;
    int p0;
    int p1;
    int fPos,fNeg;
};

int trueSearch(const vector<Vertex>& points, const Vertex& p0, int l, int r) {
    if (r >= l) { 
        int mid = l + (r - l) / 2; 

        if(!(points[mid] < p0 || p0 < points[mid]))
            return mid; 
  
        if(p0 < points[mid]) 
            return trueSearch(points, p0, l, mid - 1); 
  
        return trueSearch(points, p0, mid + 1, r); 
    }
    cerr << p0.x << " " << p0.y  << "\n";
    return -1; 
} 

int binarySearch(const vector<Vertex>& points, const Vertex& p0) {
    return trueSearch(points,p0,0,points.size()-1);
}

int main() {
    int l,n,v0,v1,fPos,fNeg;
    int nPoints=0,nEdges=0;
    vector<Vertex> points; set<Edge1> edges;
    set<Vertex> pointsSet;
    while(cin >> l >> n >> v0 >> v1 >> fPos >> fNeg) {
        Vertex point0,point1;
        Edge1 edge = {-1,{0,0},{0,0},fPos,fNeg};

        cin >> point1.x >> point1.y;
        pointsSet.insert(point1);
        edge.p0 = point1;
        for(int i=1; i<n-1; i++) {
            point0 = point1;
            cin >> point1.x >> point1.y;
            pointsSet.insert(point1);
            edge.p1 = point1;
            edges.insert(edge);
            edge.p0 = point1;
        }
        point0 = point1;
        cin >> point1.x >> point1.y;
        edge.p0 = point0; edge.p1 = point1;
        pointsSet.insert(point1);
        edges.insert(edge);
    }

    cout << pointsSet.size() << "\n";
    points.resize(pointsSet.size());
    int i=0;
    for(auto &c:pointsSet) {
        points[i].x=c.x; points[i].y=c.y;
        cout << i << " " << c.x << " " << c.y << "\n";
        i++;
    }
    i=0; Edge2 out;
    cout << edges.size() << "\n";
    for(auto it = edges.begin(); it != edges.end(); it++, i++) {
        out.label = i;
        out.fPos = (*it).fPos;
        out.fNeg = (*it).fNeg;
        out.p0 = binarySearch(points,(*it).p0);
        out.p1 = binarySearch(points,(*it).p1);
        cout << out.label << " " << out.p0 << " " << out.p1 << " " << out.fPos << " " << out.fNeg << "\n";
    }
    return 0;
}