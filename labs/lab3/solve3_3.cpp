#include "./dependens/TSolve.h"
#include <fstream>
#include <functional>

using namespace std;

static size_t powMNK = 2;

int ToInt(const string& str) {
	int res = 0;
	for (int i = str.length() - 1; i >= 0 ; --i) {
		res = res * 10 + (str[i] - '0');
	}
	return res;
}

int main(int argc, char* argv[]) {
	if (argc != 3) {
		cerr << "Error: arg is incorrect" << endl;
		exit(-1);
	}
	string dataFile = argv[1];
	powMNK = max(0, ToInt(argv[2]));
	ifstream in(dataFile, ios::in);

	size_t N;
	vector<double> X, Y;
	double temp;
	
	in >> N;
	for (size_t i = 0; i < N; ++i) {
		in >> temp;
		X.push_back(temp);
	}
	for (size_t i = 0; i < N; ++i) {
		in >> temp;
		Y.push_back(temp);
	}
	in.close();

	TMatrix matrixFI(N, powMNK + 1, Zero);
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < powMNK + 1; ++j) {
			matrixFI[i][j] = pow(X[i], 1.0 * j);
		}
	}
	TVector vecY(Y);
	string inputData 	= "./dependens/inputData";
	string outputData 	= "./dependens/outputData";
	ofstream input(inputData, ios::out);
	ifstream output(outputData, ios::in);
	TSolve solution(inputData, outputData);
 
	TMatrix writeMatrix = matrixFI.Rotate() * matrixFI;
	TVector writeVector = (matrixFI.Rotate() * vecY.Rotate()).Rotate();
	input << powMNK + 1 << endl;
	for (size_t j = 0; j < writeMatrix.GetSizeRow(); ++j) {
		for (size_t i = 0; i < writeMatrix.GetSizeCol(); ++i) {
			input << writeMatrix[j][i] << " ";
		}
		input << writeVector[j] << endl;
	}
	if (!!solution.ToSolveByGauss()) {
		cerr << "Error: troubles with LUP method, check input data" << endl;
		exit(-1);
	}
	input.close();
	output.close();

	vector<double> aproxVec(powMNK + 1);
	in.open(outputData);

	for (size_t i = 0; i < aproxVec.size(); ++i) {
		in >> aproxVec[i];
	}
	in.close();

	auto CalcApproxFunction = [](double x, const vector<double>& arr) -> double {
		double res = 0.0;
		for (size_t i = 0; i < arr.size(); ++i)
			res += arr[i] * pow(x, 1.0 * i);
		return res;
	};

	string resultFile = "./solution3_3";
	ofstream o(resultFile, ios::out);

	o << "X: ";
	for (size_t i = 0; i < N; ++i) {
		o << X[i] << " ";
	}
	o << endl;
	o << "Y: ";
	for (size_t i = 0; i < N; ++i) {
		o << Y[i] << " ";
	}
	o << endl;
	o << "Approx koef: ";
	double sumOfError = .0;
	for (size_t i = 0; i < N; ++i) {
		o << CalcApproxFunction(X[i], aproxVec) << " ";
		sumOfError += (CalcApproxFunction(X[i], aproxVec) - Y[i]) * (CalcApproxFunction(X[i], aproxVec) - Y[i]);
	}	
	o << endl;
	o << "Summ of errors: " << sumOfError << endl;
	o.close();

	string pointFiles = "./plotDataPoint";
	string approxFiles = "./plotDataApprox" + to_string(powMNK);
	ofstream outPoint(pointFiles, ios::out);
	ofstream outApprox(approxFiles, ios::out);

	for (size_t i = 0; i < N; ++i) {
		outPoint << X[i] << " " << Y[i] << endl;
	}
	outPoint.close();
	double deltaX = 0.01;
	for (size_t i = 0; i < N - 1; ++i) {
		for (double cur = X[i]; cur <= X[i + 1]; cur += deltaX) {
			outApprox << cur << " " << CalcApproxFunction(cur, aproxVec) << endl;
		}
	}
	outApprox.close();
	return 0;
}
