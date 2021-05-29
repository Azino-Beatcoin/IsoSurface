#ifndef TRIANGULATION
#define TRIANGULATION

#include <list>
#include <vector>
#include <cstdio>
// #include "R2Graph.h"
#include "r2geom.h"
#include "R3Graph.h"

// using namespace R3Graph;

class Triangulation {
public:
    class Triangle {
    public:
        int indices[3];
        //??? int adjacentTriangles[3];

        Triangle(int i0, int i1, int i2) {
            indices[0] = i0;
            indices[1] = i1;
            indices[2] = i2;
        }

        int& operator[](int i) {
            return indices[i];
        }

        int operator[](int i) const {
            return indices[i];
        }

        void invert() {
            int tmp = indices[1];
            indices[1] = indices[2];
            indices[2] = tmp;
        }
    };

    class Vertex {
    public:
        R3Graph::R3Point point;
        R3Graph::R3Vector normal;
        std::list<Triangle> adjacentTriangles;

        Vertex(
            const R3Graph::R3Point& p = R3Graph::R3Point(),
            const R3Graph::R3Vector& n = R3Graph::R3Vector(0., 1., 0.)
        ):
            point(p),
            normal(n)
        {}
    };

    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
    R3Graph::R3Box box;

public:
    Triangulation():
        vertices(),
        triangles(),
        box()
    {}

    void clear();
    void computeFramingBox();
    void orientate();
    R3Graph::R3Point center() const;

    bool save(const char *path) const;
    bool load(const char *path);
};

#endif
