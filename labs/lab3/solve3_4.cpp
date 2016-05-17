#include <iostream>
#include <functional>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char* argv[]) {
	auto firstDerivative2 = [](const vector<double>& Y, const vector<double>& X, int pos, double x) -> double {
		if (pos >= Y.size() - 2) return 0;
		return (Y[pos + 1] - Y[pos]) / (X[pos + 1] - X[pos]) + 
			((Y[pos + 2] - Y[pos + 1]) / (X[pos + 2] - X[pos + 1]) - (Y[pos + 1] - Y[pos]) / (X[pos + 1] - X[pos])) 
				/ (X[pos + 2] - X[pos]) * (2 * x - X[pos] - X[pos + 1]);
	};
	auto secondDerivative2 = [](const vector<double>& Y, const vector<double>& X, int pos, double x) -> double {
		if (pos >= Y.size() - 2) return 0;
		return 2 * ((Y[pos + 2] - Y[pos + 1]) / (X[pos + 2] - X[pos + 1]) 
				- (Y[pos + 1] - Y[pos]) / (X[pos + 1] - X[pos])) 
			/ (X[pos + 2] - X[pos]);
	};

	if (argc != 2) {
		cerr << "Error: arg is incorrect" << endl;
		exit(-1);
	}

	string dataFile = argv[1];
	vector<double> X, Y;
	size_t N;
	double temp, perfectX;
	ifstream in(dataFile, ios::in);

	in >> N;
	for (size_t i = 0; i < N; ++i) {
		in >> temp;
		X.push_back(temp);
	}
	for (size_t i = 0; i < N; ++i) {
		in >> temp;
		Y.push_back(temp);
	}
	in >> perfectX;
	in.close();

	for (size_t i = 0; i < N - 1; ++i) {
		if (X[i] <= perfectX && perfectX <= X[i + 1]) {
			cout << "F'(x) = " << firstDerivative2(Y, X, i, perfectX) << endl;
			cout << "F\"(x) = " << secondDerivative2(Y, X, i, perfectX) << endl;
			break;
		}
	}
	return 0;
}