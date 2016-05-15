CC = g++
CFLAGS = -std=c++11 -Wall -pedantic
DEP = dependens
OBJ = obj
DELFILES = *.o *~ solve *.log res* solution* plotData*
DELETE = rm -f

ifeq ($(arg), 1)
all: solve3_1.cpp
	$(CC) $(CFLAGS) $^ -o solve
else ifeq ($(arg), 2)
all: $(OBJ)/_tmp.o $(OBJ)/_TSolve.o $(OBJ)/_TVector.o $(OBJ)/_TMatrix.o $(OBJ)/_func.o
	$(CC) $^ -o solve

$(OBJ)/_tmp.o: solve3_2.cpp
	$(CC) $(CFLAGS) -c $^ -o $(OBJ)/_tmp.o

$(OBJ)/_TSolve.o: $(DEP)/TSolve.cpp
	$(CC) $(CFLAGS) -c $^ -o $(OBJ)/_TSolve.o

$(OBJ)/_TVector.o: $(DEP)/TVector.cpp
	$(CC) $(CFLAGS) -c $^ -o $(OBJ)/_TVector.o

$(OBJ)/_TMatrix.o: $(DEP)/TMatrix.cpp
	$(CC) $(CFLAGS) -c $^ -o $(OBJ)/_TMatrix.o

$(OBJ)/_func.o: $(DEP)/func.cpp
	$(CC) $(CFLAGS) -c $^ -o $(OBJ)/_func.o
else ifeq($(arg), 3)
all: $(OBJ)/_tmp.o $(OBJ)/_TSolve.o $(OBJ)/_TVector.o $(OBJ)/_TMatrix.o $(OBJ)/_func.o
	$(CC) $^ -o solve

$(OBJ)/_tmp.o: solve3_3.cpp
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
