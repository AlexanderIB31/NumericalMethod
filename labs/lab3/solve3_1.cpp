#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <cmath>

using namespace std;

void PolynomialLagrange() {
	auto Y = [](double x) -> double { return tan(x); };
	auto Product = [](const vector<double>& X, double curX, size_t pos) -> double {
		double temp = 1.0;
		for (size_t i = 0; i < 4; ++i) {
			if (i == pos) continue;
			temp *= (curX - X[i]) / (X[pos] - X[i]);
		}
		return temp;
	};

	vector<double> X;
	// task a)
	
	X.push_back(0);
	X.push_back(1.0 * M_PI / 8);
	X.push_back(2.0 * M_PI / 8);
	X.push_back(3.0 * M_PI / 8);
	
	//task b)
	/*
	X.push_back(0);
	X.push_back(1.0 * M_PI / 8);
	X.push_back(1.0 * M_PI / 3);
	X.push_back(3.0 * M_PI / 8);
	*/
	double perfectX = 3.0 * M_PI / 16;
	double sum = 0.0;
	for (size_t i = 0; i < 4; ++i) {
		sum += Y(X[i]) * Product(X, perfectX, i);	
	}
	cout << "L(" << perfectX << ") = " << sum << endl;
	cout << "y(" << perfectX << ") = " << Y(perfectX) << endl;
	cout << "|L(" << perfectX << ") - y(" << perfectX << ")| = " << fabs(Y(perfectX) - sum) << endl;
}

void PolynomialNewtoon() {
	auto Y = [](double x) -> double { return sin(M_PI * x / 6); };
	function<double(const vector<double>&, double, size_t, size_t)> FuncY; 
	FuncY = [Y, &FuncY](const vector<double>& X, double curX, size_t start, size_t end) -> double {
		if (start == end) {
			return Y(X[start]);
		}
		else {
			return (FuncY(X, curX, start, end - 1) - FuncY(X, curX, start + 1, end)) / (X[start] - X[end]);
		}
	};
	auto Product = [](const vector<double>& X, double curX, size_t sz) -> double {
		double prod = 1.0;
		for (size_t i = 0; i < sz; ++i)
			prod *= (curX - X[i]);
		return prod;
	};

	vector<double> X;
	// task a)
	
	X.push_back(0);
	X.push_back(1.0 * M_PI / 8);
	X.push_back(2.0 * M_PI / 8);
	X.push_back(3.0 * M_PI / 8);
	
	//task b)
	/*
	X.push_back(0);
	X.push_back(1.0 * M_PI / 8);
	X.push_back(1.0 * M_PI / 3);
	X.push_back(3.0 * M_PI / 8);
	*/

	double perfectX = 1.5;
	double sum = 0.0;
	for (size_t i = 1; i < 4; ++i) {
		sum += FuncY(X, perfectX, 0, i) * Product(X, perfectX, i);
	}
	cout << "P(" << perfectX << ") = " << sum << endl;
	cout << "y(" << perfectX << ") = " << Y(perfectX) << endl;
	cout << "|P(" << perfectX << ") - y(" << perfectX << ")| = " << fabs(Y(perfectX) - sum) << endl;

}

int main() {
	cout << "Polynomial Lagrange" << endl;
	PolynomialLagrange();
	cout << "****************************************" << endl;
	cout << "Polynomial Newtoon" << endl;
	PolynomialNewtoon();
	cout << "****************************************" << endl;
	return 0;
}
