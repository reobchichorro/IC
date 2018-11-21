/*

*/

#ifndef CHAIN_H
#define CHAIN_H
#include <iostream>
#include <vector>
//#include <algorithm>
using namespace std;

template <class T>
istream &operator>>(istream &, Vertex<T> &c);

template <class T>
class Vertex {
    public:
    Vertex(): label(-1), x(0), y(0) {}
    Vertex(const int label_, const T& x_, const T& y_) : label(label_), x(x_), y(y_) {}
    Vertex(const Vertex& other) : label(-1), x(0), y(0) { *this = other; }
	Vertex<T> &operator=(const Vertex<T>& other);
    friend istream &operator>> <T>(istream &, Vertex<T> &c);


    private:
    T x, y;
    int label;
};

template <class T>
Vertex<T>& Vertex<T>::operator=(const Vertex<T>& other) {
    if(&other==this) 
        return *this;
    label = other.label;
    x = other.x; y = other.y;
    return *this;
}

template <class T>
istream &operator>>(istream &is, Vertex<T> &c) {
    is >> c.label >> c.x >> c.y;
	return is;
}


template <class T>
istream &operator>>(istream &, Edge<T> &c);

template <class T>
class Edge {
    public:
    Edge(): label(-1), vInitial(), vFinal(), fPos(-1), fNeg(-1) {}
    Edge(const int label_, const Vertex<T>& vInitial_, const Vertex<T> vFinal_, const int fPos_, const int fNeg_) : label(label_), vInitial(vInitial_), vFinal(vFinal_), fPos(fPos_), fNeg(fNeg_) {}
    Edge(const Edge& other) : label(-1), vInitial(), vFinal(), fPos(-1), fNeg(-1) { *this = other; }
	Edge<T> &operator=(const Edge<T>& other);
    friend istream &operator>> <T>(istream &, Edge<T> &c);


    private:
    int label;
    Vertex<T> vInitial, vFinal;
    int fPos, fNeg;
};

template <class T>
Edge<T>& Edge<T>::operator=(const Edge<T>& other) {
    if(&other==this) 
        return *this;
    label = other.label;
    vInitial = other.vInitial; vFinal = other.vFinal; 
    fPos = other.fPos; fNeg = other.fNeg;
    return *this;
}

template <class T>
istream &operator>>(istream &is, Edge<T> &c) {
    is >> c.label >> c.vInitial >> c.vFinal >> c.fPos >> c.fNeg;
	return is;
}

/*
template <class T>
class Chain {
    public:
    Chain(): label(-1), nEdges(0), vInitial(-1), vFinal(-1), fPos(-1), fNeg(-1) { vertices(0); }
    Chain(const int label, const int nEdges, const int vInitial, const int vFinal, const int fPos, const int fNeg, const vector<pair<T,T> >& vertices);
    Chain(const Chain& other) : label(-1), nEdges(0), vInitial(-1), vFinal(-1), fPos(-1), fNeg(-1) { vertices(0); *this = other; }
	Chain<T> &operator=(const Chain<T>& other);



    private:
    int label, nEdges;
    int vInitial, vFinal;
    int fPos, fNeg;
    vector<Vertex<T> > vertices;

};

template <class T>
Chain<T>& Chain<T>::operator=(const Chain<T>& other) {
    if(&other==this) 
        return *this;
    label = other.label; nEdges = other.nEdges; 
    vInitial = other.vInitial; vFinal = other.vFinal; 
    fPos = other.fPos; fNeg = other.fNeg;
    vertices = other.vertices;
    return *this;
}
*/

#endif