CC 		= g++
CFLAGS 	= -std=c++11 -Wall -pedantic -I/opt/AMDAPPSDK-3.0/include -L/opt/AMDAPPSDK-3.0/lib/x86_64/

solve: testProg.cpp
	$(CC) $(CFLAGS) $^ -o solve -lOpenCL

clean:
	rm -rf *.o *~ solve *.log res*
