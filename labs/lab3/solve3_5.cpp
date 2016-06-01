#include <iostream>
#include <fstream>
#include <functional>
#include <cmath>

using namespace std;

int main(int argc, char* argv[]) {
	if (argc != 2) {
		cerr << "Error: arg is incorrect" << endl;
		exit(-1);
	}

	string dataFile = argv[1];
	double h1, h2, X1, Xn;
	
	ifstream in(dataFile, ios::in);
	in >> X1 >> Xn >> h1 >> h2;
	in.close();

	auto Y = [](double x) -> double {
		return x / pow(3 * x + 4, 3.0);
	};

	auto secondDerivateY = [](double x) -> double {
		return 18 * (3 * x - 4) / pow(3 * x + 4, 5.0);
	};

	auto fourthDerivateY = [](double x) -> double {
		return 3240 * (3 * x - 8) / pow(3  * x + 4, 7.0);
	};

	double M2_h1 = 0.0, M2_h2 = 0.0, M4_h1 = 0.0, M4_h2 = 0.0;

	for (double cur = X1; cur <= Xn; cur += h1) {
		M2_h1 = max(M2_h1, abs(secondDerivateY(cur)));
		M4_h1 = max(M4_h1, abs(fourthDerivateY(cur)));
	}
	for (double cur = X1; cur <= Xn; cur += h2) {
		M2_h2 = max(M2_h2, abs(secondDerivateY(cur)));
		M4_h2 = max(M4_h2, abs(fourthDerivateY(cur)));		
	}

	auto rectangleMethod = [&Y](double startX, double endX, double step) -> double {		
		double res = 0.0;
		for (double cur = startX + step; cur <= endX; cur += step) {
			res += Y(0.5 * (2 * cur - step));
		}
		return step * res;
	};

//	auto estimateRectangleMethod = [](double startX, double endX, double step, double M) -> double {
//		return step * step * M * (endX - startX) / 24;
//	};

	auto trapezoidalMethod = [&Y](double startX, double endX, double step) -> double {
		double res = 0.0;
		for (double cur = startX + step; cur <= endX; cur += step) {
			res += Y(cur) + Y(cur - step);
		}
		return 0.5 * step * res;
	};

//	auto estimateTrapezoidalMethod = [](double startX, double endX, double step, double M) -> double {
//		return step * step * M * (endX - startX) / 12;
//	};

	auto SimpsonMethod = [&Y](double startX, double endX, double step) -> double {
		double res = 0.0;
		res += Y(startX) + Y(endX);
		for (double cur = startX + step; cur < endX; cur += 2 * step) {
			res += 4 * Y(cur);
		}
		for (double cur = startX + 2 * step; cur < endX; cur += 2 * step) {
			res += 2 * Y(cur);
		}
		return step * res / 3;
	};

//	auto estimateSimpsonMethod = [](double startX, double endX, double step, double M) -> double {
//		return (endX - startX) * pow(step, 4.0) * M / 180;
//	};

	double k = max(h1, h2) / min(h1, h2);

	cout << "Rectangle method (h1): " << rectangleMethod(X1, Xn, h1) << endl;
	cout << "Rectangle method (h2): " << rectangleMethod(X1, Xn, h2) << endl;
	cout << "Runge-Romberg estimate: ";
	if (h1 >= h2) cout << fabs(rectangleMethod(X1, Xn, h2) + (rectangleMethod(X1, Xn, h2) - rectangleMethod(X1, Xn, h1)) / (pow(k, 2.0) - 1.0)) << endl;
	else cout << fabs(rectangleMethod(X1, Xn, h1) + (rectangleMethod(X1, Xn, h1) - rectangleMethod(X1, Xn, h2)) / (pow(k, 2.0) - 1.0)) << endl;
	cout << "Trapezoidal method (h1): " << trapezoidalMethod(X1, Xn, h1) << endl;
	cout << "Trapezoidal method (h2): " << trapezoidalMethod(X1, Xn, h2) << endl;
	cout << "Runge-Romberg estimate: ";	
	if (h1 >= h2) cout << fabs(trapezoidalMethod(X1, Xn, h2) + (trapezoidalMethod(X1, Xn, h2) - trapezoidalMethod(X1, Xn, h1)) / (pow(k, 2.0) - 1.0)) << endl;
	else cout << fabs(trapezoidalMethod(X1, Xn, h1) + (trapezoidalMethod(X1, Xn, h1) - trapezoidalMethod(X1, Xn, h2)) / (pow(k, 2.0) - 1.0)) << endl;
	cout << "Simpson method (h1): " << SimpsonMethod(X1, Xn, h1) <<  endl;	
	cout << "Simpson method (h2): " << SimpsonMethod(X1, Xn, h2) <<  endl;	
	cout << "Runge-Romberg estimate: ";	
	if (h1 >= h2) cout << fabs(SimpsonMethod(X1, Xn, h2) + (SimpsonMethod(X1, Xn, h2) - SimpsonMethod(X1, Xn, h1)) / (pow(k, 4.0) - 1.0)) << endl;
	else cout << fabs(SimpsonMethod(X1, Xn, h1) + (SimpsonMethod(X1, Xn, h1) - SimpsonMethod(X1, Xn, h2)) / (pow(k, 4.0) - 1.0)) << endl;	
	return 0;
}