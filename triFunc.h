#ifndef TRI_FUNC_H
#define TRI_FUNC_H

// #include "R2Graph.h"
#include "r2geom.h"
#include "R3Graph.h"
#include "Triangulation.h"

void plot2DFunction(
    double (*f)(const R2Point&),
    const R2Rectangle& rect,
    int numX, int numY,
    Triangulation& triangulation
);

R3Graph::R3Vector gradient(
    double (*f)(const R2Point&),
    const R2Point& p
);

#endif
