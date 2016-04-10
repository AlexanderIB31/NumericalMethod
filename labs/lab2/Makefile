CC = g++
CFLAGS = -c -std=c++11 -Wall -pedantic
DEP = dependens
OBJ = obj
DELFILES = *.o *~ solve *.log res*
DELETE = rm -rf

ifeq ($(arg), 2)
all: $(OBJ)/_tmp.o $(OBJ)/_TSolve.o $(OBJ)/_TVector.o $(OBJ)/_TMatrix.o $(OBJ)/_func.o $(OBJ)/func.o
	$(CC) $^ -o solve

$(OBJ)/_tmp.o: solve2_2.cpp
	$(CC) $(CFLAGS) $^ -o $(OBJ)/_tmp.o

$(OBJ)/_TSolve.o: $(DEP)/TSolve.cpp
	$(CC) $(CFLAGS) $^ -o $(OBJ)/_TSolve.o

$(OBJ)/_TVector.o: $(DEP)/TVector.cpp
	$(CC) $(CFLAGS) $^ -o $(OBJ)/_TVector.o

$(OBJ)/_TMatrix.o: $(DEP)/TMatrix.cpp
	$(CC) $(CFLAGS) $^ -o $(OBJ)/_TMatrix.o

$(OBJ)/_func.o: $(DEP)/func.cpp
	$(CC) $(CFLAGS) $^ -o $(OBJ)/_func.o

$(OBJ)/func.o: func.cpp
	$(CC) $(CFLAGS) $^ -o $(OBJ)/func.o
else
all: $(OBJ)/_tmp.o $(OBJ)/func.o
	$(CC) $^ -o solve

$(OBJ)/_tmp.o: solve2_1.cpp
	$(CC) $(CFLAGS) $^ -o $(OBJ)/_tmp.o

$(OBJ)/func.o: func.cpp
	$(CC) $(CFLAGS) $^ -o $(OBJ)/func.o
endif

clean:
	$(DELETE) $(DELFILES)

prepare:
	$(DELETE) $(OBJ) $(DEP)
	mkdir $(OBJ)
	mkdir $(DEP)
	cp ~/src/GitHub/NumericalMethod/labs/lab1/TSolve.* ./$(DEP)/
	cp ~/src/GitHub/NumericalMethod/labs/lab1/TMatrix.* ./$(DEP)/
	cp ~/src/GitHub/NumericalMethod/labs/lab1/TVector.* ./$(DEP)/
	cp ~/src/GitHub/NumericalMethod/labs/lab1/func.* ./$(DEP)/

