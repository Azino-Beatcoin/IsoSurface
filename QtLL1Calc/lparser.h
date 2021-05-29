#ifndef LPARSER_H
#define LPARSER_H

#include <string>

class Calc;

namespace lparser {

typedef const char* char_ptr;

// Exception in parsing
class ParseException {
public:
    const char* reason;

    ParseException(const char* cause = 0):
        reason(cause)
    {}

    const char* what() const {
        if (reason == 0)
            return "syntax error";
        else
            return reason;
    }
};

// Maximal length of tokens in input stream
const int NAME_MAXLEN = 128;

class LexValue {
public:
    double value;
    std::string name;

    LexValue():
        value(0.),
        name()
    {
    }

    void clear() {
        value = 0.;
        name.clear();
    }
};

class StreamPosition {
public:
    char_ptr initialPosition;
    char_ptr position;
    bool tokenRead;
    int token;
    LexValue value;
    char_ptr tokenBeg;  // Beginning of a token
    int tokenLength;    // Length of a token
    int offset;         // offset from beginning of line
    int streamLength;   // Total length of a stream

    StreamPosition(const char* line = 0):
        initialPosition(line),
        position(line),
        tokenRead(false),
        token(0),
        value(),
        tokenBeg(0),
        tokenLength(0),
        offset(0),
        streamLength(0)
    {}

    void initialize(const char* line = 0) {
        initialPosition = line;
        position = line;
        tokenRead = false;
        token = 0;
        value.clear();
        tokenBeg = line;
        tokenLength = 0;
        offset = 0;
        streamLength = 0;
        if (line != 0)
            streamLength = int(strlen(line));
    }

    int currentChar() const {
        if (offset < 0 || offset >= streamLength)
            return 0;
        int c = (int(*position) & 255);
        return c;
    }

    int charAt(const char* ptr) const {
        if (ptr < initialPosition || ptr >= (initialPosition + streamLength))
            return 0;
        int c = (int(*ptr) & 255);
        return c;
    }
};

bool parseFormula(
    const char* line,
    Calc& calc,
    std::string& errorMessage
);

void initializeParser(
    const char* line,
    StreamPosition& pos
);

void parseS(StreamPosition& pos, Calc& calc); // Parse a line

} // end of namespace lparser

#endif
