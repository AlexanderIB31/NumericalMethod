CC 		= g++
CFLAGS 	= -std=c++11 -Wall -pedantic -I/opt/AMDAPPSDK-2.9-1/include 
LFLAGS 	= -L/opt/AMDAPPSDK-2.9-1/lib/x86_64 -lOpenCL

all: testProg.cpp
	$(CC) $(CFLAGS) $< -o solve $(LFLAGS)

clean:
	rm -rf *.o *~ solve *.log res*
