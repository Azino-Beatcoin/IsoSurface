#include "calc.h"
#include <cassert>
#include <string>
#include <map>
#include <cmath>
#include <cfenv>
#include <cerrno>

#ifndef M_LN2
# define M_LN2          0.69314718055994530942  /* log_e 2 */
#endif

#ifndef M_LN10
# define M_LN10         2.30258509299404568402  /* log_e 10 */
#endif

using namespace std;

static int powmod(int b, int e, int m);

typedef double (*StdFunc)(bool& defined, double x, double y, double z);

class StdFuncDsc {      // Standard function descriptor
public:
    StdFunc func;       // Interface
    int numArgs;        // Can be 1, 2, 3

    StdFuncDsc(StdFunc f = 0, int nArgs = 1):
        func(f),
        numArgs(nArgs)
    {}
};

// Interfaces for mathematical functions
static double StdSin(bool& defined, double x, double, double) {
    defined = true;
    return sin(x);
}

static double StdCos(bool& defined, double x, double, double) {
    defined = true;
    return cos(x);
}

static double StdSqrt(bool& defined, double x, double, double) {
    if (x < 0.) {
        defined = false;
        return 0.;
    }
    defined = true;
    return sqrt(x);
}

static double StdExp(bool& defined, double x, double, double) {
    defined = true;
    errno = 0;
    feclearexcept(FE_ALL_EXCEPT);
    double res = exp(x);
    if (
        errno != 0 ||
        fetestexcept(
            FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW
        ) != 0
    ) {
        defined = false; res = 0.;
    }
    return res;
}

static double StdLog(bool& defined, double x, double, double) {
    if (x <= 0.) {
        defined = false;
        return 0.;
    }
    defined = true;
    errno = 0;
    feclearexcept(FE_ALL_EXCEPT);
    double res = log(x);
    if (
        errno != 0 ||
        fetestexcept(
            FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW
        ) != 0
    ) {
        defined = false; res = 0.;
    }
    return res;
}

static double StdLog2(bool& defined, double x, double, double) {
    if (x <= 0.) {
        defined = false;
        return 0.;
    }
    defined = true;
    errno = 0;
    feclearexcept(FE_ALL_EXCEPT);
    double res = log(x);
    if (
        errno != 0 ||
        fetestexcept(
            FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW
        ) != 0
    ) {
        defined = false; res = 0.;
    }
    return res/M_LN2;
}

static double StdLog10(bool& defined, double x, double, double) {
    if (x <= 0.) {
        defined = false;
        return 0.;
    }
    defined = true;
    errno = 0;
    feclearexcept(FE_ALL_EXCEPT);
    double res = log(x);
    if (
        errno != 0 ||
        fetestexcept(
            FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW
        ) != 0
    ) {
        defined = false; res = 0.;
    }
    return res/M_LN10;
}

static double StdTan(bool& defined, double x, double, double) {
    defined = true;
    errno = 0;
    feclearexcept(FE_ALL_EXCEPT);
    double res = tan(x);
    if (
        errno != 0 ||
        fetestexcept(
            FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW
        ) != 0
    ) {
        defined = false; res = 0.;
    }
    return res;
}

static double StdAtan(bool& defined, double x, double, double) {
    defined = true;
    return atan(x);
}

static double StdAtan2(bool& defined, double y, double x, double) {
    if (fabs(x) <= 0. && fabs(y) <= 0.) {
        defined = false;
        return 0.;
    }
    defined = true;
    return atan2(y, x);
}

static double StdPowmod(bool& defined, double z, double y, double x) {
    int m = abs(int(z));
    int e = int(y);
    int b = int(x);
    if (m == 0 || e < 0) {
        defined = false;
        return 0.;
    }
    defined = true;
    return double(powmod(b, e, m));
}

static double StdAbs(bool& defined, double x, double, double) {
    defined = true;
    return fabs(x);
}

static double StdMax(bool& defined, double y, double x, double) {
    defined = true;
    return (x >= y? x : y);
}

static double StdMin(bool& defined, double y, double x, double) {
    defined = true;
    return (x <= y? x : y);
}

static map<string, StdFuncDsc> standardFunctionMap;
static bool standardFunctionMapInitialized = false;

static const pair<string, StdFuncDsc> standardFunctions[] = {
    { string("sin"), StdFuncDsc(&StdSin, 1) },
    { string("cos"), StdFuncDsc(&StdCos, 1) },
    { string("exp"), StdFuncDsc(&StdExp, 1) },
    { string("sqrt"), StdFuncDsc(&StdSqrt, 1) },
    { string("log"), StdFuncDsc(&StdLog, 1) },
    { string("log2"), StdFuncDsc(&StdLog2, 1) },
    { string("log10"), StdFuncDsc(&StdLog10, 1) },
    { string("tan"), StdFuncDsc(&StdTan, 1) },
    { string("atan"), StdFuncDsc(&StdAtan, 1) },
    { string("abs"), StdFuncDsc(&StdAbs, 1) },
    { string("fabs"), StdFuncDsc(&StdAbs, 1) },
    { string("atan2"), StdFuncDsc(&StdAtan2, 2) },
    { string("max"), StdFuncDsc(&StdMax, 2) },
    { string("min"), StdFuncDsc(&StdMin, 2) },
    { string("powmod"), StdFuncDsc(&StdPowmod, 3) },
    { string(), StdFuncDsc(0, 0) }
};

static void initializeStandardFunctionMap() {
    if (standardFunctionMapInitialized)
        return;
    const pair<string, StdFuncDsc>* f = standardFunctions;
    while (!(f->first.empty())) {
        standardFunctionMap[f->first] = f->second;
        ++f;
    }
}

// Calc theCalc;

Calc::Calc():
    actarray(),
    stack(256)
{
}

Calc::~Calc() {
}

void Calc::addAction(const Calc::action& act) {
    actarray.push_back(act);
}

const Calc::action& Calc::getAction(int i) const {
    return actarray[i];
}

void Calc::initialize() {
    actarray.clear();
    stack.clear();
}

double Calc::calculate(
    bool& defined,
    double x /* = 0. */, 
    double y /* = 0. */, 
    double z /* = 0. */
) {
    bool def = true;
    stack.clear();
    // printf("Calculate\n");
    if (!standardFunctionMapInitialized)
        initializeStandardFunctionMap();

    for (int i = 0; def && i < (int) actarray.size(); i++) {
        double v1, v2, v3;
        const action& act = getAction(i);
        switch (act.op) {
            case Double_const:
                // printf("const %lf, stack = %p\n", act.value, stack);
                stack.push_back(act.value);
                break;
            case Plus:
                // printf("Plus\n");
                v2 = stack.back();
                stack.pop_back();
                v1 = stack.back();
                stack.pop_back();
                stack.push_back(v1 + v2);
                break;
            case Minus:
                v2 = stack.back();
                stack.pop_back();
                v1 = stack.back();
                stack.pop_back();
                stack.push_back(v1 - v2);
                break;
            case Mult:
                v2 = stack.back();
                stack.pop_back();
                v1 = stack.back();
                stack.pop_back();
                stack.push_back(v1 * v2);
                break;
            case Div:
                v2 = stack.back();
                stack.pop_back();
                v1 = stack.back();
                stack.pop_back();
                if (fabs(v2) > CALC_EPS) {
                    stack.push_back(v1 / v2);
                } else {
                    stack.push_back(0.);
                    def = false;
                }
                break;
            case Pow:
                v2 = stack.back();
                stack.pop_back();
                v1 = stack.back();
                stack.pop_back();
                if (v1 >= CALC_EPS) {
                    stack.push_back(pow(v1, v2));
                } else if (fabs(v1) <= CALC_EPS) {
                    if (v2 >= 0.) {
                        stack.push_back(pow(v1, v2));
                    } else {
                        // Zero division
                        stack.push_back(0.);
                        def = false;
                    }
                } else {
                    assert(v1 < 0.);
                    // Only integer power is defined
                    if (fabs(round(v2) - v2) <= 0.) {
                        stack.push_back(pow(v1, v2));
                    } else {
                        stack.push_back(0.);
                        def = false;
                    }
                }
                break;
            case Uminus:
                v1 = stack.back();
                stack.pop_back();
                stack.push_back(-v1);
                break;
            case X_value:
                stack.push_back(x);
                break;
            case Y_value:
                stack.push_back(y);
                break;
            case Z_value:
                stack.push_back(z);
                break;
            case StdFunc:
                {
                    if (
                        standardFunctionMap.find(act.name) ==
                        standardFunctionMap.end()
                    ) {
                        /*...
                        printf(
                            "Unknown function: %s\n", act.name.c_str()
                        );
                        ...*/
                        stack.push_back(0.);
                        break;
                    }
                    const StdFuncDsc& fDsc = standardFunctionMap[act.name];
                    double res = 0.;
                    v3 = stack.back();
                    stack.pop_back();
                    if (fDsc.numArgs > 1) {
                        v2 = stack.back();
                        stack.pop_back();
                    }
                    if (fDsc.numArgs > 2) {
                        v1 = stack.back();
                        stack.pop_back();
                    }
                    if (fDsc.numArgs == 1) {
                        res = fDsc.func(def, v3, 0., 0.);
                    } else if (fDsc.numArgs == 2) {
                        res = fDsc.func(def, v2, v3, 0.);
                    } else {
                        assert(fDsc.numArgs == 3);
                        res = fDsc.func(def, v1, v2, v3);
                    }
                    stack.push_back(res);
                    break;
                }
            case No_operation:
                break;
        }
    } // end for

    defined = def;
    if (def) {
        return stack.back();
    } else {
        return 0.;
    }
}

static int powmod(int b, int e, int m) { // Fast power algorithm
    int p = 1;
    while (e > 0) {
        // Invariant: p*b^e == const
        if ((e & 1) == 0) {     // e is even
            e /= 2;
            b = (b*b)%m;
        } else {                // e is odd
            --e;
            p = (p*b)%m;
        }
    }
    return p;
}
