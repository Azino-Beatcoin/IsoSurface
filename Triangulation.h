#ifndef TRIANGULATION
#define TRIANGULATION

#include <list>
#include <vector>
#include <set>
#include <map>
#include <cstdio>
#include <cassert>
// #include "R2Graph.h"
#include "r2geom.h"
#include "R3Graph.h"
#include <QString>

using namespace R3Graph;

class Triangulation {
public:
    class Triangle {
    public:
        int indices[3];

        Triangle(int i0 = (-1), int i1 = (-1), int i2 = (-1)) {
            int minIndex = 0;
            int minValue = i0;
            if (i1 < minValue) {
                minIndex = 1;
                minValue = i1;
            }
            if (i2 < minValue) {
                minIndex = 2;
                minValue = i2;
            }

            int ind[3];
            ind[0] = i0; ind[1] = i1; ind[2] = i2;
            int i = minIndex;
            indices[0] = minValue;
            ++i;
            if (i >= 3)
                i = 0;
            indices[1] = ind[i];
            ++i;
            if (i >= 3)
                i = 0;
            indices[2] = ind[i];
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

        bool isAdjacent(const Triangle& t) const {
            // For each edge of first triangle
            for (int i = 0; i < 3; ++i) {
                int v0 = indices[i];
                int v1;
                if (i < 2)
                    v1 = indices[i+1];
                else
                    v1 = indices[0];

                // For each edge of second triangle
                for (int j = 0; j < 3; ++j) {
                    int w0 = t.indices[j];
                    if (w0 != v0 && w0 != v1)
                        continue;
                    int w1;
                    if (j < 2)
                        w1 = t.indices[j+1];
                    else
                        w1 = t.indices[0];
                    if (
                        (v0 == w0 && v1 == w1) ||
                        (v0 == w1 && v1 == w0)
                    )
                        return true;
                } // end for (j...
            } // end for (i...
            return false;
        }

        bool operator==(const Triangle& t) const {
            return (
                indices[0] == t.indices[0] &&
                indices[1] == t.indices[1] &&
                indices[2] == t.indices[2]
            );
        }

        bool operator!=(const Triangle& t) const {
            return !operator==(t);
        }

        bool operator<(const Triangle& t) const {
            return (
                indices[0] < t.indices[0] || (
                    indices[0] == t.indices[0] && (
                        indices[1] < t.indices[1] || (
                            indices[1] == t.indices[1] &&
                            indices[2] < t.indices[2]
                        )
                    )
                )
            );
        }

        bool operator<=(const Triangle& t) const {
            return (
                indices[0] < t.indices[0] || (
                    indices[0] == t.indices[0] && (
                        indices[1] < t.indices[1] || (
                            indices[1] == t.indices[1] &&
                            indices[2] <= t.indices[2]
                        )
                    )
                )
            );
        }

        bool operator>(const Triangle& t) const {
            return !operator<=(t);
        }

        bool operator>=(const Triangle& t) const {
            return !operator<(t);
        }
    };

    class Vertex {
    public:
        R3Point point;
        R3Vector normal;

        Vertex(
            const R3Point& p = R3Point(),
            const R3Vector& n = R3Vector(0., 1., 0.)
        ):
            point(p),
            normal(n)
        {}
    };

    class Edge {
    public:
        int vertIdx[2]; // In ascending order

        Edge(int v0 = (-1), int v1 = (-1)) {
            if (v0 <= v1) {
                vertIdx[0] = v0; vertIdx[1] = v1;
            } else {
                vertIdx[0] = v1; vertIdx[1] = v0;
            }
        }

        bool operator==(const Edge& e) const {
            return (
                vertIdx[0] == e.vertIdx[0] &&
                vertIdx[1] == e.vertIdx[1]
            );
        }

        bool operator!=(const Edge& e) const {
            return !operator==(e);
        }

        bool operator<(const Edge& e) const {
            return (
                vertIdx[0] < e.vertIdx[0] ||
                (vertIdx[0] == e.vertIdx[0] &&
                vertIdx[1] < e.vertIdx[1])
            );
        }

        bool operator<=(const Edge& e) const {
            return (
                vertIdx[0] < e.vertIdx[0] ||
                (vertIdx[0] == e.vertIdx[0] &&
                vertIdx[1] <= e.vertIdx[1])
            );
        }

        bool operator>(const Edge& e) const {
            return !operator<=(e);
        }

        bool operator>=(const Edge& e) const {
            return !operator<(e);
        }
    };

    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
    R3Box box;

    class AdjacentTriangles {
    public:
        int adjacentTriangles[3];

        void clear() {
            adjacentTriangles[0] = (-1);
            adjacentTriangles[1] = (-1);
            adjacentTriangles[2] = (-1);
        }

        AdjacentTriangles() { clear(); }

        int numberOfTriangles() const {
            if (adjacentTriangles[0] < 0)
                return 0;
            else if (adjacentTriangles[1] < 0)
                return 1;
            else if (adjacentTriangles[2] < 0)
                return 2;
            else
                return 3;
        }

        int size() const { return numberOfTriangles(); }

        void add(int t) {
            for (int i = 0; i < 3; ++i) {
                if (adjacentTriangles[i] < 0) {
                    adjacentTriangles[i] = t;
                    return;
                }
                if (adjacentTriangles[i] == t) {
                    // Already added
                    return;
                }
            }
            assert(false);  // More than 3 adjacent triangles - cannot be so
        }

        void push_back(int t) { add(t); }
    };

    typedef std::vector<int> TrianglesOfVertex;
    typedef std::vector<int> TrianglesOfEdge;
    typedef std::set<int> VertexStar; // Vertices incident to this vertex

    // typedef std::vector<int> VertexRing; // Vertices incident to this vertex
    //                                     // in the ring order
    // For border edges, a vertex ring may be non-closed
    class VertexRing: public std::vector<int> {
    public:
        bool borderVertex;
    public:
        bool isBorderVertex() const { return borderVertex; }
        VertexRing(size_t n = 0):
            std::vector<int>(n),
            borderVertex(false)
        {}
    };

    mutable bool adjacentTrianglesCalculated;
    mutable std::vector<AdjacentTriangles> *adjacentTriangles;
    mutable bool trianglesOfVerticesCalculated;
    mutable std::vector<TrianglesOfVertex> trianglesOfVertices;

    mutable std::vector<VertexStar> starsOfVertices;
    mutable std::vector<VertexRing> ringsOfVertices;
    mutable std::map<Edge, TrianglesOfEdge> trianglesOfEdges;

    typedef std::list<int> LinkedComponent;
    mutable bool linkedComponentsCalculated;
    mutable std::vector<LinkedComponent>* linkedComponents;

public:
    Triangulation():
        vertices(),
        triangles(),
        box(),
        adjacentTrianglesCalculated(false),
        adjacentTriangles(0),
        trianglesOfVerticesCalculated(false),
        trianglesOfVertices(),
        linkedComponentsCalculated(false),
        linkedComponents(0)
    {}

    Triangulation(const Triangulation& t):
        vertices(t.vertices),
        triangles(t.triangles),
        box(t.box),
        adjacentTrianglesCalculated(t.adjacentTrianglesCalculated),
        adjacentTriangles(0),
        trianglesOfVerticesCalculated(t.trianglesOfVerticesCalculated),
        trianglesOfVertices(t.trianglesOfVertices),
        starsOfVertices(t.starsOfVertices),
        ringsOfVertices(t.ringsOfVertices),
        linkedComponentsCalculated(false),
        linkedComponents(0)
    {
        if (t.adjacentTrianglesCalculated) {
            assert(t.adjacentTriangles != 0);
            adjacentTriangles = new std::vector<AdjacentTriangles>;
            *adjacentTriangles = *(t.adjacentTriangles);
        }
    }

    ~Triangulation() {
        delete adjacentTriangles;
        delete linkedComponents;
    }

    Triangulation& operator=(const Triangulation& t) {
        vertices = t.vertices;
        triangles = t.triangles;
        box = t.box;
        delete adjacentTriangles; adjacentTriangles = 0;
        adjacentTrianglesCalculated = t.adjacentTrianglesCalculated;
        if (t.adjacentTrianglesCalculated) {
            assert(t.adjacentTriangles != 0);
            adjacentTriangles = new std::vector<AdjacentTriangles>;
            *adjacentTriangles = *(t.adjacentTriangles);
        }
        trianglesOfVerticesCalculated = t.trianglesOfVerticesCalculated;
        trianglesOfVertices = t.trianglesOfVertices;
        starsOfVertices = t.starsOfVertices;
        ringsOfVertices = t.ringsOfVertices;
        linkedComponentsCalculated = false;
        delete linkedComponents; linkedComponents = 0;
        return *this;
    }

    void clear();
    void computeFramingBox();
    void orientate();
    void computeNormals();
    R3Point center() const;

    void refine();      // Remove double vertices

    void defineAdjacentTriangles() const;
    void defineTrianglesOfVertices() const;
    int defineLinkedComponents() const;
    void clearLinkedComponents() const {
        delete linkedComponents; linkedComponents = 0;
        linkedComponentsCalculated = false;
    }

    void invalidateAdjacentTriangles() const {
        adjacentTrianglesCalculated = false;
        delete adjacentTriangles; adjacentTriangles = 0;
    }

    void invalidateTrianglesOfVertices() const {
        trianglesOfVerticesCalculated = false;
        trianglesOfVertices.clear();
        starsOfVertices.clear();
        ringsOfVertices.clear();
        trianglesOfEdges.clear();
    }

    void copyMaximalComponent(Triangulation& t) const;

    // Point is inside the 3D model
    bool contains(const R3Point& t) const;

    bool save(const char *path) const;
    bool load(const char *path);
    bool saveVRML(
        const QString& path,  double pixel_mm = 0.5
    ) const;
    bool saveSTLbinary(
        const QString& path,  double pixel_mm = 0.5
    ) const;
    bool saveSTLascii(
        const QString& path,  double pixel_mm = 0.5
    ) const;
    bool saveOBJ(       // Wavefront OBJ format
        const QString& path,  double pixel_mm = 0.5
    ) const;

    void computeVertexRing(
        int vertexIdx, VertexRing& vertexRing
    ) const;

    void checkEdges(
        int& numManifoldEdges,
        int& numBorderEdges,
        int& numNonManifoldEdges
    ) const;

    void annihilateSmallTriangles(
        int& removedEdges,
        int& removedCompleteTriangles,
        int& numRemovedTriangles,
        double eps = 0.1
    );

    // Return a maximal edge length after simplification
    double simplify(
        int& removedVerticesNum,
        int& removedEdgesNum,
        int& removedTrianglesNum,
        double eps = 0.01
    );

    /*double simplify2(
        int& removedEdges,
        int& numRemovedTriangles,
        double eps = 0.1
    );*/

    double mySimplify(
        int& removedVerticesNum,
        int& removedEdgesNum,
        int& removedTrianglesNum,
        double eps = 0.01
    );

    std::vector<int> findPath(int startVert, int endVert, std::set<int> stopList = std::set<int>()) const;
    std::set<int> findInterVertices(int startVert, std::set<int> stopList = std::set<int>()) const;

    void cotangentLaplaceSmoothing(
        double lambda = 0.330   // May be negative for inflation/Taubin smooth
    );

    void uniformLaplaceSmoothing(
        double lambda = 0.330   // May be negative for inflation/Taubin smooth
    );

    // Smoothing with cotangent Laplacian and Taubin shrink-inflate sequence
    void taubinSmoothing(
        int iterations = 1,
        double lambda = 0.330,
        double mu = 0.331,
        bool useCotangentLaplace = false
    );

    R3Vector cotangentLaplace(int vertexIdx) const;
    R3Vector uniformLaplace(int vertexIdx) const;
};

#endif
