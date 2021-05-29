//
// The class Calc
// stores the reverse polish notation of an expression
// and calculates its value using a stack of real numbers
//
#ifndef CALC_H
#define CALC_H

#include <vector>
#include <string>
#include <cmath>

const double CALC_EPS = 1e-15;

class Calc {
public:
    enum Operation {   // Operations
        Plus,
        Minus,
        Mult,
        Div,
        Pow,
        Uminus,
        Double_const,
        X_value,
        Y_value,
        Z_value,
        StdFunc,
        No_operation
    };

    class action {
    public:
        enum Operation op;
        double value;
        std::string name; // name of standard function

        action():
            op(No_operation),
            value(0.),
            name()
        {}

        action(Operation oper, double v = 0., std::string n = std::string()):
            op(oper),
            value(v),
            name(n)
        {}
    };

private:
    // Array of actions (Reverse polish notation of expression)
    std::vector<action> actarray;

    std::vector<double> stack;

public:
    Calc();
    ~Calc();

    void addAction(const action& act);
    const action& getAction(int i) const;
    const action& operator[](int i) const { return getAction(i); }
    int getActarraySize() const { return (int) actarray.size(); }
    void initialize();
    void clear() { initialize(); }
    int size() const { return int(actarray.size()); }

    // Calcualte an expression in point (x, y, z)
    // Returns: value of the function in a point x
    // Out: defined == the function is defined in this point
    double calculate(
        bool& defined, 
        double x = 0., double y = 0., double z = 0.
    );
};

// extern Calc theCalc;

#endif /* CALC_H */
