#include "Triangulation.h"
#include "triFunc.h"
#include <string>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

using namespace R3Graph;

void plot2DFunction(
    double (*f)(const R2Point&),
    const R2Rectangle& rect,
    int numX, int numY,
    Triangulation& triangulation
) {
    triangulation.clear();
    double x0 = rect.left();
    double y0 = rect.bottom();
    double dx = rect.width() / numX;
    double dy = rect.height() / numY;

    // Array of vertices
    for (int iy = 0; iy <= numY; ++iy) {
        for (int ix = 0; ix <= numX; ++ix) {
            double x = x0 + dx*ix;
            double y = y0 + dy*iy;
            R2Point p(x, y);
            triangulation.vertices.push_back(
                Triangulation::Vertex(
                    R3Point(x, y, f(p)),
                    gradient(f, p)
                )
            );
        }
    }
    // Array of triangles
    for (int iy = 0; iy < numY; ++iy) {
        for (int ix = 0; ix < numX; ++ix) {
            //  p2=(ix, iy+1)  p3=(ix+1, iy+1)
            //  p0=(ix, iy)    p1=(ix+1, iy)
            // Two triangles:
            //     (p0, p1, p3),
            //     (p0, p3, p2)
            int p0 = iy*(numX + 1) + ix;
            int p1 = p0 + 1;
            int p2 = p0 + numX + 1;
            int p3 = p2 + 1;
            triangulation.triangles.push_back(
                Triangulation::Triangle(p0, p1, p3)
            );
            triangulation.triangles.push_back(
                Triangulation::Triangle(p0, p3, p2)
            );
        }
    }
    triangulation.computeFramingBox();
}

static const double GRADIENT_DX = 1e-3;
static const double GRADIENT_DY = 1e-3;

R3Vector gradient(
    double (*f)(const R2Point&),
    const R2Point& p
) {
    R2Vector dx(GRADIENT_DX, 0.);
    R2Vector dy(0., GRADIENT_DY);
    double df_dx = (f(p + dx) - f(p - dx))/(2.* GRADIENT_DX);
    double df_dy = (f(p + dy) - f(p - dy))/(2.* GRADIENT_DY);
    // F(x, y, z) = f(x, y) - z
    R3Vector grad(-df_dx, -df_dy, 1.);
    grad.normalize();
    return grad;
}
