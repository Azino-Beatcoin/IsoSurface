// Recursive parser for the grammar of
// arithmetic formulas
//
// S -> F $
// F -> TE
// E -> empty | +TE | -TE
// E ->
// T -> MG
// G -> empty | *MG | /MG
// M -> Base EndBase | -M
// EndBase -> empty | POW M
// Base -> Const | X | Y | Z | (F) | Func
// Func -> NAME(Args)
// Args -> empty | ArgList
// ArgList -> F ArgsTail
// ArgsTail -> empty | , F ArgsTail
//
//
//
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cctype>
#include <cmath>

#include <string>

#include "lparser.h"
#include "calc.h"

using namespace lparser;

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#    define M_E 2.7182818284590452354
#endif

// void parseS(StreamPosition& pos, Calc& calc);        // Parse a line
static void parseF(StreamPosition& pos, Calc& calc);    // Parse a formula
static void parseT(StreamPosition& pos, Calc& calc);    // Parse a term,
static void parseE(StreamPosition& pos, Calc& calc);    //   end of term
static void parseM(StreamPosition& pos, Calc& calc);    // Parse a factor,
static void parseG(StreamPosition& pos, Calc& calc);    //   end of factor
static void parseBase(StreamPosition& pos, Calc& calc); // Parse base of power
static void parseEndBase(StreamPosition& pos, Calc& calc); // end of base
static void parseFunc(StreamPosition& pos, Calc& calc); // Parse function
static void parseArgs(StreamPosition& pos, Calc& calc); //   function arguments
static void parseArgsTail(StreamPosition& pos, Calc& calc);

static void syntaxError(const char* txt = 0);

// Token types
static const int CONSTANT = 0;
static const int NAME = 1;
static const int LPAR = 2;
static const int RPAR = 3;
static const int PLUS = 4;
static const int MINUS = 5;
static const int MULT = 6;
static const int DIV = 7;
static const int CM = 8;        // Comma ","
static const int X = 9;         // x
static const int Y = 10;        // y
static const int Z = 11;        // z
static const int POW = 12;
static const int ENDLINE = 13;
static const int UNDEFINED = 14;

static int nextToken(LexValue* v, StreamPosition& pos);
static void skipToken(StreamPosition& pos);
static void readToken(StreamPosition& pos);

int nextToken(LexValue* v, StreamPosition& pos) {
    if (!pos.tokenRead)
        readToken(pos);
    assert(pos.tokenRead);
    if (v != 0)
        *v = pos.value;
    return pos.token;
}

void skipToken(StreamPosition& pos) {
    if (!pos.tokenRead)
        readToken(pos);
    pos.position = pos.tokenBeg + pos.tokenLength;
    pos.tokenRead = false;
}

static void readToken(StreamPosition& pos) {
    pos.tokenRead = true;

    // Skip space
    while (
        pos.offset < pos.streamLength &&
        pos.currentChar() != 0 && isspace(pos.currentChar())
    ) {
        ++(pos.position);
        ++(pos.offset);
    }

    if (pos.offset >= pos.streamLength || pos.currentChar() == 0) {
        pos.token = ENDLINE;
        pos.tokenLength = 1;
        pos.value.value = 0.;
        pos.value.name = std::string("\n");
        return;
    }

    pos.tokenBeg = pos.position;
    pos.tokenLength = 1;

    const char* curpos = pos.position;
    int nextChar = pos.currentChar();
    int offset = pos.offset;
    ++curpos; ++offset;

    switch (nextChar) {
    case '(':
        pos.token = LPAR; break;
    case ')':
        pos.token = RPAR; break;
    case '+':
        pos.token = PLUS; break;
    case '-':
        pos.token = MINUS; break;
    case '*':
        if (offset < pos.streamLength && pos.charAt(curpos) == '*') {
            // "**" is power
            pos.token = POW;
            pos.tokenLength = 2;
        } else {
            pos.token = MULT;
        }
        break;
    case '^':
        pos.token = POW; break;
    case '/':
        pos.token = DIV; break;
    case ',':
        pos.token = CM; break;
    case '$':
        pos.token = ENDLINE; break;
    default:
        if (isalpha(nextChar)) {

            // Name
            pos.token = NAME;
            while (
                offset < pos.streamLength && pos.charAt(curpos) != 0 &&
                (
                    isalpha(pos.charAt(curpos)) ||
                    isdigit(pos.charAt(curpos))
                )
            ) {
                ++curpos; ++offset; ++(pos.tokenLength);
            }
            assert(curpos - pos.tokenBeg == pos.tokenLength);
            if (pos.tokenLength > NAME_MAXLEN - 1)
                pos.tokenLength = NAME_MAXLEN - 1;
            pos.value.value = 0.;
            pos.value.name = std::string(pos.tokenBeg, pos.tokenLength);
            if (pos.value.name == "x") {
                pos.token = X;
            } else if (pos.value.name == "y") {
                pos.token = Y;
            } else if (pos.value.name == "z") {
                pos.token = Z;
            } else if (pos.value.name == "pi") {
                pos.token = CONSTANT;
                pos.value.value = M_PI;
            } else if (pos.value.name == "e") {
                pos.token = CONSTANT;
                pos.value.value = M_E;
            }

        } else if (isdigit(nextChar)) {

            // Numeric constant
            pos.token = CONSTANT;
            while (
                offset < pos.streamLength && pos.charAt(curpos) != 0 &&
                isdigit(pos.charAt(curpos))
            ) {
                ++curpos; ++offset; ++(pos.tokenLength);
            }
            if (offset < pos.streamLength && pos.charAt(curpos) == '.') {
                // Skip decimal point
                ++curpos; ++offset; ++(pos.tokenLength);
                while (
                    offset < pos.streamLength && pos.charAt(curpos) != 0 &&
                    isdigit(pos.charAt(curpos))
                ) {
                    ++curpos; ++offset; ++(pos.tokenLength);
                }
            }
            assert(curpos - pos.tokenBeg == pos.tokenLength);
            if (pos.tokenLength > NAME_MAXLEN - 1) {
                pos.tokenLength = NAME_MAXLEN - 1;
            }

            pos.value.name = std::string(pos.tokenBeg, pos.tokenLength);
            pos.value.value = atof(pos.value.name.c_str());
        } else {
            pos.token = UNDEFINED;
            pos.value.name = std::string(pos.tokenBeg, pos.tokenLength);
            pos.value.value = 0.;
        }
        break;
    } // end switch
}

void lparser::parseS(StreamPosition& pos, Calc& calc) {
    calc.initialize();
    parseF(pos, calc);
    if (nextToken(0, pos) != ENDLINE)
        syntaxError("end of formula");
}

void parseF(StreamPosition& pos, Calc& calc) {
    // printf("parseF\n");

    parseT(pos, calc);
    parseE(pos, calc);
}

void parseE(StreamPosition& pos, Calc& calc) {
    // printf("parseE\n");

    int next = nextToken(0, pos);
    if (next == PLUS) {
        skipToken(pos);
        parseT(pos, calc);
        calc.addAction(Calc::action(Calc::Plus));
        parseE(pos, calc);
    } else if (next == MINUS) {
        skipToken(pos);
        parseT(pos, calc);
        calc.addAction(Calc::action(Calc::Minus));
        parseE(pos, calc);
    } else if (next == RPAR || next == CM || next == ENDLINE) {
        // Nothing to do
    } else {
        syntaxError("Illegal end of sum");
    }
}

void parseT(StreamPosition& pos, Calc& calc) {
    // printf("parseT\n");

    int next = nextToken(0, pos);
    if (next == MINUS) {
        skipToken(pos);
        parseT(pos, calc);
        calc.addAction(Calc::action(Calc::Uminus));
    } else {
        parseM(pos, calc);
        parseG(pos, calc);
    }
}

void parseG(StreamPosition& pos, Calc& calc) {
    // printf("parseG\n");

    int next = nextToken(0, pos);
    if (next == MULT) {
        skipToken(pos);
        parseM(pos, calc);
        calc.addAction(Calc::action(Calc::Mult));
        parseG(pos, calc);
    } else if (next == DIV) {
        skipToken(pos);
        parseM(pos, calc);
        calc.addAction(Calc::action(Calc::Div));
        parseG(pos, calc);
    } else if (
        next == PLUS || next == MINUS ||
        next == RPAR || next == CM ||
        next == ENDLINE
    ) {
        // Nothing to do
    } else {
        /*
        printf(
            "Next token: %d \"%s\"\n",
            next, pos.position
        );
        */
        syntaxError("Illegal end of term");
    }
}

void parseM(StreamPosition& pos, Calc& calc) {
    int next = nextToken(0, pos);
    if (next == MINUS) {
        // Unary minus
        skipToken(pos);
        parseM(pos, calc);
        calc.addAction(Calc::action(Calc::Uminus));
    } else {
        parseBase(pos, calc);
        parseEndBase(pos, calc);
    }
}

void parseBase(StreamPosition& pos, Calc& calc) { // Base of power
    LexValue lval;
    int next = nextToken(&lval, pos);
    if (next == CONSTANT) {
        skipToken(pos);
        calc.addAction(
            Calc::action(Calc::Double_const, lval.value)
        );
    } else if (next == X) {
        skipToken(pos);
        calc.addAction(Calc::action(Calc::X_value));
    } else if (next == Y) {
        skipToken(pos);
        calc.addAction(Calc::action(Calc::Y_value));
    } else if (next == Z) {
        skipToken(pos);
        calc.addAction(Calc::action(Calc::Z_value));
    } else if (next == NAME) {
        parseFunc(pos, calc);
    } else if (next == LPAR) {
        skipToken(pos);
        parseF(pos, calc);
        if (nextToken(0, pos) != RPAR)
            syntaxError();
        skipToken(pos);
    } else {
        syntaxError("Illegal factor");
    }
}

void parseEndBase(StreamPosition& pos, Calc& calc) {
    int next = nextToken(0, pos);
    if (next == POW) {
        skipToken(pos);
        parseM(pos, calc);
        calc.addAction(Calc::action(Calc::Pow));
        parseG(pos, calc);
    } else if (
        next == MULT || next == DIV ||
        next == PLUS || next == MINUS ||
        next == RPAR || next == CM ||
        next == ENDLINE
    ) {
        // Nothing to do
    } else {
        /*
        printf(
            "Next token: %d \"%s\"\n",
            next, pos.position
        );
        */
        syntaxError("Illegal end of power");
    }
}

static void parseFunc(StreamPosition& pos, Calc& calc) {
    // printf("parseFunc\n");

    LexValue v;
    int t = nextToken(&v, pos);
    assert(t == NAME);
    if (t != NAME)
        syntaxError("Incorrect function name");
    std::string funcName = v.name;
    skipToken(pos);
    t = nextToken(0, pos);
    if (t != LPAR)
        syntaxError("Function: no \"(\"");
    skipToken(pos);
    parseArgs(pos, calc);
    t = nextToken(0, pos);
    if (t != RPAR)
        syntaxError("Function: no \")\"");
    skipToken(pos);
    calc.addAction(
        Calc::action(
            Calc::StdFunc, 0., funcName
        )
    );
}

static void parseArgs(StreamPosition& pos, Calc& calc) {
    // printf("parseArgs\n");

    int t = nextToken(0, pos);
    if (t == RPAR)
        return;
    parseF(pos, calc);
    parseArgsTail(pos, calc);
}

static void parseArgsTail(StreamPosition& pos, Calc& calc) {
    // printf("parseArgsTail\n");

    int t = nextToken(0, pos);
    if (t == RPAR)
        return;
    if (t != CM)
        syntaxError("Function args: no \",\"");
    skipToken(pos);
    parseF(pos, calc);
    parseArgsTail(pos, calc);
}

static void syntaxError(const char* txt) {
#   ifdef PARSER_TEST
    if (txt == 0)
        printf("Syntax error.\n");
    else
        printf("Syntax error: %s.\n", txt);
#   endif

    throw ParseException(txt);
}

void lparser::initializeParser(const char* line, StreamPosition& pos) {
    pos.initialize(line);
}

bool lparser::parseFormula(
    const char* line,
    Calc& calc,
    std::string& errorMessage
) {
    StreamPosition pos;
    initializeParser(line, pos);

    try {
        parseS(pos, calc);
        errorMessage = "";
    } catch (const ParseException& e) {
        // Syntax error
        errorMessage = e.what();
        return false;
    }
    return true;
}

#ifdef PARSER_TEST

int main() {
    char line[1024];
    int line_len = 0;
    double x = 10., y = 20., z = 30.;

    /*
    printf("Input the values of x, y, z:\n");
    // scanf("%lf%lf%lf\n", &x);
    fgets(line, 1022, stdin);
    // x = atof(line);
    sscanf(line, "%lf%lf%lf", &x, &y, &z);
    */

    Calc calc;

    printf("Input a formula to evaluate:\n");

    while (true) {

        fgets(line, 1022, stdin);
        line_len = int(strlen(line));
        if (line_len > 0 && line[line_len-1] == '\n') {
            line[line_len-1] = 0; --line_len;
        }
        if (line_len > 0 && line[line_len-1] == '\r') {
            line[line_len-1] = 0; --line_len;
        }
        if (line_len == 0)
            break;

        std::string errorMessage;
        bool parsingDone = parseFormula(
            line, calc, errorMessage
        );
        if (parsingDone) {
            // Calculate the value of the expression, show the result
            bool defined;
            double res = calc.calculate(defined, x, y, z);
            if (defined) {
                printf("= %g\n", res);
            } else {
                printf("= undefined\n");
            }
        }
    }
    return 0;
}

#endif /* PARSER_TEST */
