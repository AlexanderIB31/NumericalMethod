CC      =g++
CFLAGS  =-c -Wall -pedantic -std=c++11
#LFLAGS  =

all: func.o TMatrix.o TVector.o TSolve.o solve1.o
	$(CC) $^ -o solve

solve1.o: solve1.cpp
	$(CC) $(CFLAGS) $<

TSolve.o: TSolve.cpp
	$(CC) $(CFLAGS) $<

TVector.o: TVector.cpp
	$(CC) $(CFLAGS) $<

TMatrix.o: TMatrix.cpp
	$(CC) $(CFLAGS) $<

#TCommonClass.o: TCommonClass.cpp
#	$(CC) $(CFLAGS) $<

func.o: func.cpp
	$(CC) $(CFLAGS) $<

clean:
	rm -rf *.o *~ solve
