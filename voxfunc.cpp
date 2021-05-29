#include "R3Graph.h"
#include "voxelset.h"
#include "voxfunc.h"

using namespace R3Graph;

static R3Point curve(double t);

double torusFunc(const R3Point& p);
static R3Point curve(double t);

double torusFunc(const R3Point& p) {
    // Slanted torus
    R3Vector n(-2., -1., 1);
    n.normalize();
    R3Vector ez = n;

    R3Vector ey(1., 2., 4.);
    ey.normalize();

    R3Vector ex = ey.vectorProduct(ez);
    ex.normalize();

    R3Point o(0., 0., 0.);  // Center
    double R = 4.;
    double r = 1.;

    R3Vector v = p - o;
    R3Vector proj_v = ex*(v*ex) + ey*(v*ey);
    double l = v.length();
    if (l < (R - r)/2.)
        return (-(R - r)/2.);
    proj_v.normalize();
    R3Point t = o + proj_v*R;
    double d = t.distance(p);
    return (r - d);
    // double d = (p - t).length2();
    // return (r*r - d);
}

double voxelFunc(const Voxel& voxel) {
    double dx = (XMAX - XMIN)/double(NUM_VOXELS_X);
    double dy = (YMAX - YMIN)/double(NUM_VOXELS_Y);
    double dz = (ZMAX - ZMIN)/double(NUM_VOXELS_Z);

    double dx2 = dx/2.;
    double dy2 = dy/2.;
    double dz2 = dz/2.;

    R3Point p(
        XMIN + dx2 + voxel.point.x*dx,
        YMIN + dy2 + voxel.point.y*dy,
        ZMIN + dz2 + voxel.slice*dz
    );

    return spaceFunc(p);
}

double distanceToCurve(
    const R3Graph::R3Point& p,
    R3Graph::R3Point (*curve)(double),
    double t0 /* = 0. */, double t1 /* = 1. */,
    double* tmin /* = 0 */
) {
    double eps = 1e-6;      // Presizion
    double dt = 0.05;       // Initial step
    double dmin = 1e+30;    // Infinity
    double res = 1e+30;
    double t_min = (t0 + t1)/2.;
    double tLower = t0;
    double tUpper = t1;
    double tMinus = tLower + dt, tPlus = tUpper - dt;
    while (tUpper - tLower > eps) {
        double t = tLower;
        dmin = 1e+30;       // Infinity
        while (t <= tUpper) {
            R3Point q = (*curve)(t);
            double d = p.distance(q);
            if (d < dmin) {
                dmin = d;
                tMinus = t - dt;
                tPlus = t + dt;
                if (d < res) {
                    res = d;
                    t_min = t;
                }
            }
            t += dt;
        }
        tLower = tMinus;
        if (tLower < t0)
            tLower = t0;
        tUpper = tPlus;
        if (tUpper > t1)
            tUpper = t1;
        dt /= 8.;
    }
    if (tmin != 0) {
        *tmin = t_min;
    }
    return dmin;
}

static R3Point curve(double t) {
    // z(t) = (-(t-0.25)^3)*12./0.4375 + 5.57142857142857
    // y(t) = (-(t-0.25)^3)*16./1.75 + 1.85714285714286
    // x(t) = 8 - t*16
    double t3 = (t + 0.1*sin(10.*t)- 0.25);
    t3 = (-t3*t3*t3);
    return R3Point(
        7. - t*14.,
        -t*t*t3*12./1.75 + 1.85714285714286,
        t3*12./0.4375 + 5.57142857142857
    );
}

double spaceFunc(const R3Point& p) {
    /*...
    double tmin;
    double d = distanceToCurve(
        p, &curve, -10., 11., &tmin
    );
    double r = 0.75;
    return (r - d);
    ...*/
    return torusFunc(p);
}
