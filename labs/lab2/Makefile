CC = g++
CFLAGS = -std=c++11 -Wall -pedantic
DEP = dependens/
OBJ = obj/
DELFILES = *.o *~ solve *.log res* $(OBJ)/*
DELETE = rm -rf

ifeq ($(arg), 2)
all: $(OBJ)_tmp.o $(OBJ)_TSolve.o $(OBJ)_TVector.o $(OBJ)_TMatrix.o $(OBJ)_func.o $(OBJ)func.o
	$(CC) $^ -o solve

$(OBJ)_tmp.o: solve2_2.cpp
	$(CC) -c $(CFLAGS) $^ -o $(OBJ)_tmp.o

$(OBJ)_TSolve.o: $(DEP)TSolve.cpp
	$(CC) -c $(CFLAGS) $^ -o $(OBJ)_TSolve.o

$(OBJ)_TVector.o: $(DEP)TVector.cpp
	$(CC) -c $(CFLAGS) $^ -o $(OBJ)_TVector.o

$(OBJ)_TMatrix.o: $(DEP)TMatrix.cpp
	$(CC) -c $(CFLAGS) $^ -o $(OBJ)_TMatrix.o

$(OBJ)_func.o: $(DEP)func.cpp
	$(CC) -c $(CFLAGS) $^ -o $(OBJ)_func.o

$(OBJ)func.o: func.cpp
	$(CC) -c $(CFLAGS) $^ -o $(OBJ)func.o
else
all: $(OBJ)_tmp.o $(OBJ)func.o
	$(CC) $^ -o solve
$(OBJ)_tmp.o: solve2_2.cpp
	$(CC) -c $(CFLAGS) $^ -o $(OBJ)_tmp.o

$(OBJ)func.o: func.cpp
	$(CC) -c $(CFLAGS) $^ -o $(OBJ)func.o
endif

clean:
	$(DELETE) $(DELFILES)

prepare:
	$(DELETE) $(OBJ)
	mkdir $(OBJ)

