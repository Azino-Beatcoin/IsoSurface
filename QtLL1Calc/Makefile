CC = g++ $(CFLAGS)
CFLAGS =

lparser: lparser.cpp lparser.h calc.h calc.cpp
	$(CC) -o lparser -D PARSER_TEST lparser.cpp calc.cpp

clean:
	rm -f lparser *.o
