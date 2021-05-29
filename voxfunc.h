#ifndef VOXEL_FUNC_H
#define VOXEL_FUNC_H

#include "R3Graph.h"
// #include "voxelset.h"

class Voxel;

// const double XMIN = 0., XMAX = 512.;
// const double YMIN = 0., YMAX = 512.;
// const double ZMIN = 0., ZMAX = 512.;
const double XMIN = (-8.), XMAX = 8.;
const double YMIN = (-8.), YMAX = 8.;
const double ZMIN = (-8.), ZMAX = 8.;

const int NUM_VOXELS_X = 64;
const int NUM_VOXELS_Y = 64;
const int NUM_VOXELS_Z = 64;

double voxelFunc(const Voxel& voxel);
double spaceFunc(const R3Graph::R3Point& p);

double distanceToCurve(
    const R3Graph::R3Point& p,
    R3Graph::R3Point (*curve)(double t),
    double t0 = 0., double t1 = 1.,
    double* tmin = 0
);

#endif
