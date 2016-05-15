CC      =g++
CFLAGS  =-c -Wall -pedantic -std=c++11

all: func.o TMatrix.o temp.o
	$(CC) $^ -o solve

temp.o: temp.cpp
	$(CC) $(CFLAGS) $<

TMatrix.o: TMatrix.cpp
	$(CC) $(CFLAGS) $<

func.o: func.cpp
	$(CC) $(CFLAGS) $<

clean:
	rm -rf *.o *~ temp *.log result
