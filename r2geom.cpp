#include <cassert>
#include "r2geom.h"
#include "gauss.h"

// Not used any more...
static const int NUM_SMOOTHING_WEIGHTS = 19;
static const double SMOOTHING_WEIGHTS[NUM_SMOOTHING_WEIGHTS] = {
    0.0303539138078867,
    0.166155727317393,
    0.744658307092434,
    2.73237224472926,
    8.20849986238988,
    20.1896517994655,
    40.6569659740599,
    67.0320046035639,
    90.4837418035959,
    100.000000000000,
    90.4837418035959,
    67.0320046035639,
    40.6569659740599,
    20.1896517994655,
    8.20849986238988,
    2.73237224472926,
    0.744658307092434,
    0.166155727317393,
    0.0303539138078867
};
static const double SMOOTHING_WEIGHTS_SUM = 560.488808472044;

R2Vector R2Contour::ray_directions[MAX_INTERPOLATION_RAYS];
bool R2Contour::ray_directions_calculated = false;

void R2Contour::calculate_ray_directions() {
    if (!ray_directions_calculated) {
        double dAlpha = 2.*PI/double(MAX_INTERPOLATION_RAYS);
        double alpha = 0.;
        for (int i = 0; i < MAX_INTERPOLATION_RAYS; ++i) {
            ray_directions[i] = R2Vector(
                cos(alpha), sin(alpha)
            );
            alpha += dAlpha;
        }
        ray_directions_calculated = true;
    }
}

static double smoothing_rate(int contour_perimeter) {
    const double a = 3.696;
    const double b = (-12.8);
    const double c = 0.03;
    // f(x, a, b, c) = a/(x - b) + c
    // f(4, a, b, c) = 0.25
    // f(40, a, b, c) = 0.1
    // f(infinity, a, b, c) = 0.03
    double x = double(contour_perimeter);
    return (a/(x - b) + c);
}

void R2Contour::computeSmoothedContour(R2Contour& smoothedContour) const {
    calculate_ray_directions();
    if (!rayDistancesCalculated)
        calculateRayDistances();
    std::vector<double> distances = rayDistances;
    R2Point center = centroid();
    R2Contour res(MAX_INTERPOLATION_RAYS);
    int s = int(size());
    double smoothingRate = smoothing_rate(s);
    double sigma = double(MAX_INTERPOLATION_RAYS)*smoothingRate/6.;
    if (sigma < 1.)
        sigma = 1.;
    if (sigma > double(MAX_INTERPOLATION_RAYS)/(4.*6.))
        sigma = double(MAX_INTERPOLATION_RAYS)/(4.*6.);
    std::vector<double> smoothingWeights;
    int numWeights;
    defineSmoothingWeights(
        sigma,
        numWeights, smoothingWeights
    );
    assert((numWeights & 1) != 0); // numWeights must be odd
    assert(numWeights < MAX_INTERPOLATION_RAYS/2);

    int i0 = MAX_INTERPOLATION_RAYS -
        (numWeights - 1)/2;  // -5
    for (int ray = 0; ray < MAX_INTERPOLATION_RAYS; ++ray) {
        R2Vector s(0., 0.);
        int i = i0;
        for (int k = 0; k < numWeights; ++k) {
            s += ray_directions[i]*
                (rayDistances[i]*smoothingWeights.at(k));
            ++i;
            if (i >= MAX_INTERPOLATION_RAYS)
                i -= MAX_INTERPOLATION_RAYS;
        }
        ++i0;
        if (i0 >= MAX_INTERPOLATION_RAYS)
            i0 -= MAX_INTERPOLATION_RAYS;
        res[ray] = center + s;
    }
    smoothedContour = res;
}

R2Contour::R2Contour(const I2Contour& c):
    contourType(c.contourType),
    rayDistancesCalculated(false)
{
    resize(c.size());
    for (size_t i = 0; i < c.size(); ++i) {
        at(i) = R2Point(double(c.at(i).x), double(c.at(i).y));
    }
}

R2Contour& R2Contour::operator=(const I2Contour& c) {
    contourType = c.contourType;
    rayDistancesCalculated = false;
    resize(c.size());
    for (size_t i = 0; i < c.size(); ++i) {
        at(i) = R2Point(double(c.at(i).x), double(c.at(i).y));
    }
    return *this;
}

I2Contour::I2Contour(const R2Contour& c)
{
    *this = c;
}

I2Contour& I2Contour::operator=(const R2Contour& c) {
    contourType = c.contourType;
    resize(c.size());
    if (c.size() == 0)
        return *this;

    resize(c.size());
    I2Point p0(
        int(round(c.at(0).x)),
        int(round(c.at(0).y))
    );
    at(0) = p0;
    int n = 1;
    for (size_t i = 1; i < c.size(); ++i) {
        I2Point p(
            int(round(c.at(i).x)),
            int(round(c.at(i).y))
        );
        if (p != p0) {
            at(n) = p;
            ++n;
            p0 = p;
        }
    }
    if (n < int(c.size()))
        resize(n);
    return *this;
}

double I2Contour::internalArea() const {
    return (area() - perimeter()/2. + 1.);
}

bool intersectStraightLines(
    const R2Point& p, const R2Vector& v,    // First line
    const R2Point& q, const R2Vector& w,    // Second line
    R2Point& intersection                   // Result
) {
    assert(
        v != R2Vector(0., 0.) &&
        w != R2Vector(0., 0.)
    );
    if (
        v == R2Vector(0., 0.) ||
        w == R2Vector(0., 0.)
    )
        return false;
    R2Vector n = v.normal();
    // Point s on the second line:
    //     s = q + w*t
    // Intersection: scalar product equals zero
    // (s - p, n) == 0
    // (q + w*t - p, n) = 0
    // t = (p - q, n) / (w, n)
    double wn = w*n;
    double qpn = (p - q)*n;

    if (fabs(wn) <= R2_EPSILON) {
        // Parallel lines
        if (fabs(qpn) <= R2_EPSILON) {
            // Equal lines
            intersection = p;
            return true;
        }
        return false;
    }
    double t = qpn / wn;
    intersection = q + w*t;
    return true;
}

bool intersectLineSegments(
    const R2Point& p0, const R2Point& p1,   // First line segment
    const R2Point& q0, const R2Point& q1,   // Second line segmant
    R2Point& intersection                   // Result
) {
    R2Point inter;
    if (!intersectStraightLines(
        p0, p1 - p0,
        q0, q1 - q0,
        inter
    ))
        return false;

    R2Vector n = (p1 - p0).normal();
    R2Vector m = (q1 - q0).normal();
    if (
        ((p0 - q0)*m) * ((p1 - q0)*m) <= 0. &&
        ((q0 - p0)*n) * ((q1 - p0)*n) <= 0.
    ) {
        intersection = inter;
        return true;
    } else {
        return false;
    }
}

bool intersectLineSegmentAndLine(
    const R2Point& p0, const R2Point& p1,   // Line segment
    const R2Point& q, const R2Vector& v,    // Straight line
    R2Point& intersection                   // Result
) {
    R2Point inter;
    if (!intersectStraightLines(
        p0, (p1 - p0),
        q, v,
        inter
    ))
        return false;

    R2Vector n = v.normal();
    if (
        ((p0 - q)*n) * ((p1 - q)*n) <= 0.
    ) {
        intersection = inter;
        return true;
    } else {
        return false;
    }
}

// Interpolate contours that have star-like shapes.
// 0 <= t <= 1. For t=0 the result equals *this,
//              for t=1 the result equals c.
R2Contour R2Contour::starInterpolation(const R2Contour& c, double t) const {
    calculate_ray_directions();
    if (!rayDistancesCalculated)
        calculateRayDistances();
    if (!c.rayDistancesCalculated)
        c.calculateRayDistances();
    double t1 = 1. - t;
    R2Point center0 = centroid();
    R2Point center1 = c.centroid();
    R2Point center = center0 + (center1 - center0)*t;
    R2Contour res(MAX_INTERPOLATION_RAYS);
    for (int ray = 0; ray < MAX_INTERPOLATION_RAYS; ++ray) {
        double dist = rayDistances[ray]*t1 + c.rayDistances[ray]*t;
        res[ray] = center + ray_directions[ray]*dist;
    }
    return res;
}

static const double DALPHA = 2.*PI/double(MAX_INTERPOLATION_RAYS);

void R2Contour::calculateRayDistances() const {
    assert(orientation() >= 0);
    R2Vector xAxis(1., 0);
    // rayDistances.resize(MAX_INTERPOLATION_RAYS);
    rayDistances.assign(MAX_INTERPOLATION_RAYS, 0.);

    R2Point c = centroid();
    for (int i = 0; i < int(size()); ++i) {
        int j = i + 1;
        if (j >= int(size()))
            j = 0;
        const R2Point& p0 = at(i);
        const R2Point& p1 = at(j);
        R2Vector v0 = p0 - c;
        R2Vector v1 = p1 - c;
        R2Vector v = p1 - p0;
        R2Vector n = v0.normal();
        if (v*n < 0.)
            continue;   // Skip an edge with backward direction

        double alpha0 = xAxis.angle(v0);
        double alpha1 = xAxis.angle(v1);
        if (alpha0 < 0.)
            alpha0 += 2.*PI;
        if (alpha1 < 0.)
            alpha1 += 2.*PI;
        int ray0 = int(alpha0/DALPHA);
        int ray1 = int(alpha1/DALPHA);
        if (ray1 < MAX_INTERPOLATION_RAYS-1)
            ++ray1;
        int ray = ray0;
        while (true) {
            R2Point q;
            if (intersectLineSegmentAndLine(
                p0, p1,
                c, ray_directions[ray],
                q
            )) {
                double d = (q - c)*ray_directions[ray]; // distance to center
                assert (d >= 0.);
                if (d > rayDistances.at(ray)) {
                    rayDistances.at(ray) = d;
                }
            }
            if (ray == ray1)
                break;
            ++ray;
            if (ray == MAX_INTERPOLATION_RAYS)
                ray = 0;
        }
    }
    rayDistancesCalculated = true;
}

bool I2Contour::canAdd(const I2Point& p) const {
    if (size() == 0)
        return true;
    if (p == back())
        return false;
    R2Point q0(back().x, back().y);
    R2Point q1(p.x, p.y);
    for (int i = 0; i < int(size() - 2); ++i) {
        R2Point p0(at(i).x, at(i).y);
        R2Point p1(at(i+1).x, at(i+1).y);
        R2Point t;
        if (intersectLineSegments(p0, p1, q0, q1, t))
            return false;
    }
    return true;
}

bool I2Contour::canClose() const {
    if (size() < 3)
        return false;
    R2Point q0(back().x, back().y);
    R2Point q1(front().x, front().y);
    for (int i = 1; i < int(size() - 2); ++i) {
        R2Point p0(at(i).x, at(i).y);
        R2Point p1(at(i+1).x, at(i+1).y);
        R2Point t;
        if (intersectLineSegments(p0, p1, q0, q1, t))
            return false;
    }
    return true;
}
