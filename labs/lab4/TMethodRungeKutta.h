#pragma once
#include <cmath>
#include <fstream>
#include <iostream>
#include <functional>
#include <iomanip>

using namespace std;

class TMethodRungeKutta {
public:
	TMethodRungeKutta(const string& outFile);
	~TMethodRungeKutta();
	void ToSolve();
private:	
	double _FuncExpression(double x, double y) {
		return 2.0 * cos(x) - y;
		//return (y + x) * (y + x);
	};
	double _FuncY(double x) {
		return exp(x * x) + exp(x * sqrt(2.0)) + exp(-x * sqrt(2.0));
		//return tan(x) - x;
	};
	static constexpr double _firstDerivativeFuncY0 = 0.0;
	static constexpr double _funcY0 = 3.0;//0.0;
	static constexpr double _a = 0.0;
	static constexpr double _b = 1.0;//0.5;
	static constexpr double _h = 0.1;
	ofstream out;
	ofstream log;
};