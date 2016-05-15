#include "TMatrix.h"
#include <ctime>

using namespace std;

int main() {
	static const int sz = 1000;
	clock_t start = clock();
	double** g1 = (double**)malloc(sizeof(double*) * sz);
	double** g2 = (double**)malloc(sizeof(double*) * sz);
	for (int i = 0; i < sz; ++i) {
		g2[i] = g1[i] = (double*)malloc(sizeof(double) * sz);
		for (int j = 0; j < sz; ++j)
			g2[i][j] = g1[i][j] = i * sz + j;
	}
	TMatrix A(g1, sz, sz);
	TMatrix B(g2, sz, sz);
	A = A * B;
	A = B * A;
	clock_t end = clock();
	cout << (end - start) << endl;
	return 0;
}
