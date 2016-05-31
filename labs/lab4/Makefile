CC = g++
CFLAGS = -std=c++11 -Wall -pedantic
DEP = dependences
OBJ = obj
DELFILES = *.o *~ solve *.log res*
DELETE = rm -f

ifeq ($(arg), 1)
all: solve4_1.cpp $(OBJ)/euler.o $(OBJ)/runge-kutta.o $(OBJ)/adams.o
	$(CC) $(CFLAGS) $^ -o solve

$(OBJ)/euler.o: $(DEP)/TMethodEuler.cpp
	$(CC) -c $(CFLAGS) $^ -o $(OBJ)/euler.o

$(OBJ)/runge-kutta.o: $(DEP)/TMethodRungeKutta.cpp
	$(CC) -c $(CFLAGS) $^ -o $(OBJ)/runge-kutta.o

$(OBJ)/adams.o: $(DEP)/TMethodAdams.cpp
	$(CC) -c $(CFLAGS) $^ -o $(OBJ)/adams.o	
endif

ifeq ($(arg), 2)
all: $(OBJ)/_tmp.o $(OBJ)/_TSolve.o $(OBJ)/_TVector.o $(OBJ)/_TMatrix.o $(OBJ)/_func.o
	$(CC) $^ -o solve

$(OBJ)/_tmp.o: solve4_2.cpp
	$(CC) $(CFLAGS) -c $^ -o $(OBJ)/_tmp.o

$(OBJ)/_TSolve.o: $(DEP)/TSolve.cpp
	$(CC) $(CFLAGS) -c $^ -o $(OBJ)/_TSolve.o

$(OBJ)/_TVector.o: $(DEP)/TVector.cpp
	$(CC) $(CFLAGS) -c $^ -o $(OBJ)/_TVector.o

$(OBJ)/_TMatrix.o: $(DEP)/TMatrix.cpp
	$(CC) $(CFLAGS) -c $^ -o $(OBJ)/_TMatrix.o

$(OBJ)/_func.o: $(DEP)/func.cpp
	$(CC) $(CFLAGS) -c $^ -o $(OBJ)/_func.o
endif
clean:
	$(DELETE) $(DELFILES) $(OBJ)/*
