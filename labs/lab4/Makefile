CC = g++
CFLAGS = -std=c++11 -Wall -pedantic
DEP = dependences
OBJ = obj
DELFILES = *.o *~ solve *.log res*
DELETE = rm -f

all: solve4_1.cpp $(OBJ)/euler.o $(OBJ)/runge-kutta.o
	$(CC) $(CFLAGS) $^ -o solve

$(OBJ)/euler.o: $(DEP)/TMethodEuler.cpp
	$(CC) -c $(CFLAGS) $^ -o $(OBJ)/euler.o

$(OBJ)/runge-kutta.o: $(DEP)/TMethodRungeKutta.cpp
	$(CC) -c $(CFLAGS) $^ -o $(OBJ)/runge-kutta.o

clean:
	$(DELETE) $(DELFILES) $(OBJ)/*
