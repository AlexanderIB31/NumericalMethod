#pragma once
#include <cmath>
#include <fstream>
#include <iostream>
#include <functional>
#include <iomanip>

using namespace std;

class TMethodEuler {
public:
	TMethodEuler(const string& outFile);
	~TMethodEuler();
	void RungeRomberg();
	void ToSolve();
private:	
	double _FuncExpression(double x, double y) {
		return 2.0 * y + 4.0 * x * x * exp(x * x);
		//return (y + x) * (y + x);
	};
	double _FuncY(double x) {
		return exp(x * x) + exp(x * sqrt(2.0)) + exp(-x * sqrt(2.0));
		//return tan(x) - x;
	};
	static constexpr double _funcZ0 = 0.0;
	static constexpr double _funcY0 = 3.0;//0.0;
	static constexpr double _a = 0.0;
	static constexpr double _b = 1.0;//0.5;
	static constexpr double _h = 0.1;
	ofstream out;
	ofstream log;
};