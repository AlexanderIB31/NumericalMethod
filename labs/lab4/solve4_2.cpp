#include <iostream>
#include <functional>
#include <vector>
#include <utility>
#include <cmath>
#include <fstream>
#include "dependences/TSolve.h"

using namespace std;

// TODO: write method Runge-Romberg

pair<double, double> ToSolveRungeKutta(double a,
									   double b,
									   double h,
									   double funcY0,
									   double funcZ0,
									   function<double(double, double)> FuncExpression,
									   const char* fName) {
	/* Runge-Kutta method */	
	bool flag = false;
	ofstream out;
	if (!!fName) {
		flag = true;
		out.open(fName, ios::out | ios::app);
	}
	vector<double> X, Y, Z;
	double x0 = a, y0 = funcY0, z0 = funcZ0;
	size_t N = (b - a) / h;
	double K1, K2, K3, K4, L1, L2, L3, L4;
	double y_k = y0;
	double x_k = x0;
	double z_k = z0;
	double deltaZ, deltaY;
	for (size_t k = 0; k < N; ++k) {
		L1 = h * FuncExpression(x_k, y_k);
		K1 = h * z_k;		

		K2 = h * (z_k + 0.5 * L1);
		L2 = h * FuncExpression(x_k + 0.5 * h, y_k + 0.5 * K1);

		K3 = h * (z_k + 0.5 * L2);
		L3 = h * FuncExpression(x_k + 0.5 * h, y_k + 0.5 * K2);

		K4 = h * (z_k + L3);
		L4 = h * FuncExpression(x_k + h, y_k + K3);

		deltaY = 1.0 / 6.0 * (K1 + 2.0 * K2 + 2.0 * K3 + K4);
		deltaZ = 1.0 / 6.0 * (L1 + 2.0 * L2 + 2.0 * L3 + L4);

		X.push_back(x_k);
		Y.push_back(y_k);
		Z.push_back(z_k);

		x_k += h;
		z_k += deltaZ;
		y_k += deltaY;		
	}	
	X.push_back(x_k);
	Y.push_back(y_k);
	Z.push_back(z_k);
	if (flag) {
		auto FuncY = [](double x) -> double {
			return 1.0 / x + 1.0;
		};
		vector <double> tmp;
		out.precision(5);
		out << "x : ";
		for (size_t i = 0; i < X.size(); ++i) {
			out << fixed << X[i] << "\t";
			tmp.push_back(FuncY(X[i]));
		}
		out << endl;
		out << "y : ";
		for (size_t i = 0; i < Y.size(); ++i) {
			out << fixed << Y[i] << "\t";
		}
		out << endl;	
		out << "y*: ";
		for (size_t i = 0; i < Y.size(); ++i) {
			out << fixed << tmp[i] << "\t";
		}
		out << endl;			
		out.close();
	}
	return make_pair(Y.back(), Z.back());
}

void ShooterMethod(const char* f) {
	double a = 1.0, b = 2.0, h = 0.1;
	double n1 = b, n2 = b - h, n;
	double eps = 0.0001;
	double funcZ0 = -1.0;
	size_t j = 0;
	bool flag = false;
	vector<double> X, Y;
	auto FuncExpression = [](double x, double y) -> double {
		return 2.0 * y / (x * x * (x + 1.0));
	};
	auto calc = [](double y, double dy) -> double {
		return 2.0 * y - 4.0 * dy;
	};
	auto calcFi = [&calc](const pair<double, double>& p) -> double {
		return calc(p.first, p.second) - 4.0;
	};
	auto checkRightBorder = [&calcFi](const pair<double, double>& p, double eps) -> bool {
		return fabs(calcFi(p)) < eps;
	};

	string fn = f;
	fn += "ShooterMethod";
	ofstream out(fn, ios::out);
	out.precision(5);
	
	out << "\tn\ty\tFi" << endl;
	for (; !flag; j += 2) {
		auto tmp1 = ToSolveRungeKutta(a, b, h, n1, funcZ0, FuncExpression, NULL);
		auto tmp2 = ToSolveRungeKutta(a, b, h, n2, funcZ0, FuncExpression, NULL);			
		out << fixed << j << "\t" << n1 << "\t" << tmp1.first << "\t" << calcFi(tmp1) << endl;
		out << fixed << j + 1 << "\t" << n2 << "\t" << tmp2.first << "\t" << calcFi(tmp2) << endl;
		n = n2 - (n2 - n1) 
			/ (calcFi(tmp2) 
				- calcFi(tmp1)) 
			* calcFi(tmp2);
		n1 = n2;
		n2 = n;			
		auto tmp3 = ToSolveRungeKutta(a, b, h, n, funcZ0, FuncExpression, NULL);
		flag = checkRightBorder(tmp3, eps);
	}
	auto tmp = ToSolveRungeKutta(a, b, h, n2, funcZ0, FuncExpression, NULL);
	out << fixed << j << "\t" << n2 << "\t" << tmp.first << "\t" << calcFi(tmp) << endl;
	out.close();
	ToSolveRungeKutta(a, b, h, n2, funcZ0, FuncExpression, fn.c_str());		
}

void FiniteDifferenceMethod (const char* f) {
	double a = 1.0, b = 2.0, h = 0.1, ya, yb;
	size_t N = fabs(b - a) / h;
	vector<double> X, Y;	
	for (size_t i = 0; i <= N; ++i) {
		X.push_back(a + h * i);
	}
	auto qx = [](double x) -> double {
		return -2.0 / (x * x * (x + 1.0));
	};
	auto funcY = [](double x) -> double {
		return 1.0 / x + 1;
	};

	string fn = f;	
	string inFile = "dependences/inputData"; 
	string outFile = "dependences/outputData";

	ofstream out(fn + "FiniteDifferenceMethod", ios::out);
	ofstream output(outFile, ios::out);
	
	output << N - 1 << endl;
	output << -1.0 + h * h * qx(X[0]) << " " << 1.0 << " " << 0.0 << " " << -h << endl;
	for (size_t i = 1; i < N - 2; ++i) {
		output << 1.0 << " " << -2.0 + h * h * qx(X[i]) << " " << 1.0 << " " << 0.0 << endl;
	}
	output << 1.0 << " " << -2.0 + h * h * qx(X[N - 1]) - 4.0 / (2.0 * h - 4.0) << " " << 0.0 << " " << -4.0 * h / (2.0 * h - 4.0) << endl; 
	output.close();
	
	TSolve solution(outFile, inFile);
	if (!!solution.ToSolveByTripleDiagMatrix()) {
		cerr << "Error: Troubles with method TripleDiagMatrix!" << endl;
		exit(-1);
	}	
	
	ifstream in(inFile, ios::in);
	double temp;
	while (in >> temp) {
		Y.push_back(temp);
	}
	in.close();
	ya = Y[0] + h;
	yb = 4.0 * (h - Y.back()) / (2.0 * h - 4.0);
	out << "\tx\ty\ty*" << endl;
	out.precision(5);
	out << fixed << 0 << "\t" << X[0] << "\t" << ya << "\t" << funcY(X[0]) << endl;
	for (size_t i = 0; i < Y.size(); ++i) {
		out << fixed << i + 1 << "\t" << X[i + 1] << "\t" << Y[i] << "\t" << funcY(X[i + 1]) << endl;
	}
	out << fixed << N << "\t" << X[N] << "\t" << yb << "\t" << funcY(X[N]) << endl;
	out.close();
	output.close();
}

pair<double, double> tempSolutionRK(double a,
								    double b,
								    double h,
								    double funcY0,
								    double funcZ0,
								    function<double(double, double)> FuncExpression,
								    vector<double>* v) {
	/* Runge-Kutta method */	
	bool flag = false;
	if (!!v) {
		flag = true;
	}
	vector<double> X, Y, Z;
	double x0 = a, y0 = funcY0, z0 = funcZ0;
	size_t N = (b - a) / h;
	double K1, K2, K3, K4, L1, L2, L3, L4;
	double y_k = y0;
	double x_k = x0;
	double z_k = z0;
	double deltaZ, deltaY;
	for (size_t k = 0; k < N; ++k) {
		L1 = h * FuncExpression(x_k, y_k);
		K1 = h * z_k;		

		K2 = h * (z_k + 0.5 * L1);
		L2 = h * FuncExpression(x_k + 0.5 * h, y_k + 0.5 * K1);

		K3 = h * (z_k + 0.5 * L2);
		L3 = h * FuncExpression(x_k + 0.5 * h, y_k + 0.5 * K2);

		K4 = h * (z_k + L3);
		L4 = h * FuncExpression(x_k + h, y_k + K3);

		deltaY = 1.0 / 6.0 * (K1 + 2.0 * K2 + 2.0 * K3 + K4);
		deltaZ = 1.0 / 6.0 * (L1 + 2.0 * L2 + 2.0 * L3 + L4);

		X.push_back(x_k);
		Y.push_back(y_k);
		Z.push_back(z_k);

		x_k += h;
		z_k += deltaZ;
		y_k += deltaY;		
	}	
	X.push_back(x_k);
	Y.push_back(y_k);
	Z.push_back(z_k);
	if (flag) {
		for (size_t i = 0; i < Y.size(); ++i) {
			v->push_back(Y[i]);
		}
	}
	return make_pair(Y.back(), Z.back());
}


vector<double> RungeRombergForShooterMethod(double step) {
	double a = 1.0, b = 2.0, h = step;
	double n1 = b, n2 = b - h, n;
	double eps = 0.0001;
	double funcZ0 = -1.0;
	size_t j = 0;
	bool flag = false;
	vector<double> X, Y;
	auto FuncExpression = [](double x, double y) -> double {
		return 2.0 * y / (x * x * (x + 1.0));
	};
	auto calc = [](double y, double dy) -> double {
		return 2.0 * y - 4.0 * dy;
	};
	auto calcFi = [&calc](const pair<double, double>& p) -> double {
		return calc(p.first, p.second) - 4.0;
	};
	auto checkRightBorder = [&calcFi](const pair<double, double>& p, double eps) -> bool {
		return fabs(calcFi(p)) < eps;
	};

	for (; !flag; j += 2) {
		auto tmp1 = tempSolutionRK(a, b, h, n1, funcZ0, FuncExpression, NULL);
		auto tmp2 = tempSolutionRK(a, b, h, n2, funcZ0, FuncExpression, NULL);			
		n = n2 - (n2 - n1) 
			/ (calcFi(tmp2) 
				- calcFi(tmp1)) 
			* calcFi(tmp2);
		n1 = n2;
		n2 = n;			
		auto tmp3 = tempSolutionRK(a, b, h, n, funcZ0, FuncExpression, NULL);
		flag = checkRightBorder(tmp3, eps);
	}
	auto tmp = tempSolutionRK(a, b, h, n2, funcZ0, FuncExpression, NULL);
	tempSolutionRK(a, b, h, n2, funcZ0, FuncExpression, &Y);	
	return Y;
}

vector<double> RungeRombergForFiniteDifferenceMethod(double step) {
	double a = 1.0, b = 2.0, h = step, ya, yb;
	size_t N = fabs(b - a) / h;
	vector<double> X, Y;	
	for (size_t i = 0; i <= N; ++i) {
		X.push_back(a + h * i);
	}
	auto qx = [](double x) -> double {
		return -2.0 / (x * x * (x + 1.0));
	};
	auto funcY = [](double x) -> double {
		return 1.0 / x + 1;
	};

	string inFile = "dependences/inputData"; 
	string outFile = "dependences/outputData";

	ofstream output(outFile, ios::out);
	
	output << N - 1 << endl;
	output << -1.0 + h * h * qx(X[0]) << " " << 1.0 << " " << 0.0 << " " << -h << endl;
	for (size_t i = 1; i < N - 2; ++i) {
		output << 1.0 << " " << -2.0 + h * h * qx(X[i]) << " " << 1.0 << " " << 0.0 << endl;
	}
	output << 1.0 << " " << -2.0 + h * h * qx(X[N - 1]) - 4.0 / (2.0 * h - 4.0) << " " << 0.0 << " " << -4.0 * h / (2.0 * h - 4.0) << endl; 
	output.close();
	
	TSolve solution(outFile, inFile);
	if (!!solution.ToSolveByTripleDiagMatrix()) {
		cerr << "Error: Troubles with method TripleDiagMatrix!" << endl;
		exit(-1);
	}	
	
	ifstream in(inFile, ios::in);
	double temp;
	while (in >> temp) {
		Y.push_back(temp);
	}
	in.close();
	ya = Y[0] + h;
	yb = 4.0 * (h - Y.back()) / (2.0 * h - 4.0);
	vector<double> result;
	result.push_back(ya);
	for (size_t i = 0; i < Y.size(); ++i) {
		result.push_back(Y[i]);
	}
	result.push_back(yb);
	return result;
}

int main(int argc, char* argv[]) {
	if (argc != 2) {
		cerr << "Error! Incorrect number of arguments" << endl;
		exit(-1);
	}
	ShooterMethod(argv[1]);
	{
		string fn = argv[1];
		fn += "Runge-RombergForShooter";
		ofstream out(fn, ios::out);

		double a = 1.0, b = 2.0, h = 0.1;
		double N = fabs(b - a) / h;
		vector<double> Y1, Y2, X;
		for (size_t i = 0; i <= N; ++i) {
			X.push_back(a + h * i);
		}
		Y1 = RungeRombergForShooterMethod(h);		
		Y2 = RungeRombergForShooterMethod(h / 2);
		out.precision(5);
		out << "\tx\ty\ty*\terr" << endl;
		for (size_t i = 0, j = 0; i < X.size() && j < Y2.size(); ++i, j += 2) {
			out << fixed << i << "\t" << X[i] << "\t" << Y1[i] << "\t" << Y2[j] + (Y2[j] - Y1[i]) / 15.0 << "\t" << fabs(Y1[i] - Y2[j] - (Y2[j] - Y1[i]) / 15.0)  << endl;
		}
		out.close();
	}
	FiniteDifferenceMethod(argv[1]);
	{
		string fn = argv[1];
		fn += "Runge-RombergForFiniteDifferenceMethod";
		ofstream out(fn, ios::out);

		double a = 1.0, b = 2.0, h = 0.1;
		double N = fabs(b - a) / h;
		vector<double> Y1, Y2, X;
		for (size_t i = 0; i <= N; ++i) {
			X.push_back(a + h * i);
		}
		Y1 = RungeRombergForFiniteDifferenceMethod(h);		
		Y2 = RungeRombergForFiniteDifferenceMethod(h / 2);
		out.precision(5);
		out << "\tx\ty\ty*\terr" << endl;
		for (size_t i = 0, j = 0; i < X.size() && j < Y2.size(); ++i, j += 2) {
			out << fixed << i << "\t" << X[i] << "\t" << Y1[i] << "\t" << Y2[j] + (Y2[j] - Y1[i]) << "\t" << fabs(Y1[i] - Y2[j] - (Y2[j] - Y1[i]))  << endl;
		}
		out.close();		
	}
	return 0;
}