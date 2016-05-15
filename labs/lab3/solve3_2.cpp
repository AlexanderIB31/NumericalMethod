#include "./dependens/TSolve.h"
#include <map>
#include <utility>

using namespace std;

static const size_t powerSplane = 3;

struct TData {
	pair<double, double> pSection;
	vector<double> vKoef;
};

int main(int argc, char* argv[]) {
	if (argc != 2) {
		cerr << "Error: count of args is incorrect" << endl;
		exit(-1);
	}	
	string dataFile = argv[1];
	ifstream in(dataFile, ios::in);

	size_t N;
	vector<double> X, Y;
	vector<TData> table;
	double perfectX;
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
	in >> perfectX;
	in.close();	

	for (size_t i = 1; i < N; ++i) {
		TData temp;
		temp.pSection = make_pair(X[i - 1], X[i]);
		table.push_back(temp);
	}

	string inputFile 	= "./dependens/inputData";
	string outputFile 	= "./dependens/outputData";

	ofstream out(inputFile, ios::out);
	N--;	
	out << N - 1 << endl;
	out << 2 * (X[2] - X[0]) << " " << (X[2] - X[1]) << " " << 0.0 << " " << 3 * ((Y[2] - Y[1]) / (X[2] - X[1]) - (Y[1] - Y[0]) / (X[1] - X[0])) << endl;
	for (size_t i = 3; i <= N - 1; ++i) {
		out << (X[i - 1] - X[i - 2]) << " " << 2 * (X[i] - X[i - 2]) << " " << (X[i] - X[i - 1]) << " " << 3 * ((Y[i] - Y[i - 1]) / (X[i] - X[i - 1]) - (Y[i - 1] - Y[i - 2]) / (X[i - 1] - X[i - 2])) << endl;
	}
	out << 0.0 << " " << (X[N - 1] - X[N - 2]) << " " << 2 * (X[N] - X[N - 2]) << " " << 3 * ((Y[N] - Y[N - 1]) / (X[N] - X[N - 1]) - (Y[N - 1] - Y[N - 2]) / (X[N - 1] - X[N - 2])) << endl;
	out.close();
	TSolve sol(inputFile, outputFile);
	if (!!sol.ToSolveByTripleDiagMatrix()) {
		cerr << "Error: Troubles with method TripleDiagMatrix!" << endl;
		exit(-1);
	}
	in.close();
	in.open(outputFile);
	vector<double> vecC;
	vecC.push_back(0.0);
	while (in >> temp) {
		vecC.push_back(temp);
	}
	for (size_t i = 1; i < N; ++i) {
		vector<double> t(powerSplane + 1);
		t[0] = Y[i - 1];
		t[1] = (Y[i] - Y[i - 1]) / (X[i] - X[i - 1]) - 1.0 / 3 * (X[i] - X[i - 1]) * (vecC[i + 1] + 2 * vecC[i]);
		t[2] = vecC[i - 1];
		t[3] = (vecC[i + 1] - vecC[i]) / (3 * (X[i] - X[i - 1]));
		table[i - 1].vKoef = t;	
 	}
 	vector<double> vecTemp;
 	vecTemp.push_back(Y[N - 1]);
 	vecTemp.push_back((Y[N] - Y[N - 1]) / (X[N] - X[N - 1]) - 2.0 / 3 * (X[N] - X[N - 1]) * vecC[N]);
 	vecTemp.push_back(vecC[N]);
 	vecTemp.push_back(-vecC[N] / (3 * (X[N] - X[N - 1])));
 	table[N - 1].vKoef = vecTemp;
/*
	for (size_t i = 1; i < N - 1; i += 2) {	
		TSolve solution(inputFile, outputFile);
		ofstream input(inputFile, ios::out);		
		
		double x0 = X[i - 1];
		double x1 = X[i];
		double x2 = X[i + 1];
		double y0 = Y[i - 1];
		double y1 = Y[i];
		double y2 = Y[i + 1];
		
		input << 2 * (powerSplane + 1) << endl;
		for (size_t k = 0; k < powerSplane + 1; ++k) {
			input << pow(x0, k * 1.0) << " ";
		}
		for (size_t k = 0; k < powerSplane + 1; ++k) {
			input << 0.0 << " ";
		}
		input << y0 << endl;

		for (size_t k = 0; k < powerSplane + 1; ++k) {
			input << pow(x1, k * 1.0) << " ";
		}
		for (size_t k = 0; k < powerSplane + 1; ++k) {
			input << 0.0 << " ";
		}
		input << y1 << endl;

		for (size_t k = 0; k < powerSplane + 1; ++k) {
			input << 0.0 << " ";
		}
		for (size_t k = 0; k < powerSplane + 1; ++k) {
			input << pow(x1, k * 1.0) << " ";
		}
		input << y1 << endl;

		for (size_t k = 0; k < powerSplane + 1; ++k) {
			input << 0.0 << " ";
		}
		for (size_t k = 0; k < powerSplane + 1; ++k) {
			input << pow(x2, k * 1.0) << " ";
		}
		input << y2 << endl;

		input << 0.0 << " " << 1.0 << " " << 2.0 * x1 << " " << 3.0 * x1 * x1 << " ";
		input << 0.0 << " " << -1.0 << " " << -2.0 * x1 << " " << -3.0 * x1 * x1 << " ";
		input << 0.0 << endl;

		input << 0.0 << " " << 0.0 << " " << 2.0 << " " << 6.0 * x1 << " ";
		input << 0.0 << " " << 0.0 << " " << -2.0 << " " << -6.0 * x1 << " ";
		input << 0.0 << endl;

		input << 0.0 << " " << 0.0 << " " << 2.0 << " " << 6.0 * X[0] << " ";
		for (size_t k = 0; k < powerSplane + 1; ++k) {
			input << 0.0 << " ";
		}
		input << 0.0 << endl;

		for (size_t k = 0; k < powerSplane + 1; ++k) {
			input << 0.0 << " ";
		}
		input << 0.0 << " " << 0.0 << " " << -2.0 << " " << -6.0 * X[4] << " ";
		input << 0.0 << endl;
		
		if (!!solution.ToSolveByGauss()) {
			cerr << "Error: LUP method dropped exception"<< endl;
			exit(-1);
		}		
		input.close();		
		ifstream output(outputFile, ios::in);
		double temp;
		TData dataTemp;
		for (size_t k = 0; k < powerSplane + 1; ++k) {
			output >> temp;
			table[i - 1].vKoef.push_back(temp);
		}
		dataTemp.vKoef.clear();
		for (size_t k = 0; k < powerSplane + 1; ++k) {
			output >> temp;
			table[i].vKoef.push_back(temp);
		}
	}
	for (size_t i = 0; i < table.size(); ++i) {
		if (perfectX >= table[i].pSection.first && perfectX <= table[i].pSection.second) {
			double valueInPoint = .0;
			double a = table[i].pSection.first;
			for (size_t j = 0; j < table[i].vKoef.size(); ++j) {
				valueInPoint += table[i].vKoef[j] * pow(perfectX - a, j * 1.0);
			}
			cout << "f(" << perfectX << ") = " << valueInPoint << endl;
			break;
		}
		else continue;
	}
	*/
	// draw graphics
	double deltaX = 0.01;	
	for (size_t i = 0; i < table.size(); ++i) {
		double start = table[i].pSection.first;
		double end = table[i].pSection.second;
		ofstream out("plotData" + to_string(i), ios::out);
		for (double cur = start; cur <= end; cur += deltaX) {
			double f = .0;
			for (size_t k = 0; k < table[i].vKoef.size(); ++k) {
				f += table[i].vKoef[k] * pow(cur - start, k * 1.0);
			}
			out << cur << " " << f << endl;
		}
		cout << "[" << table[i].pSection.first << "; " << table[i].pSection.second << "]\t";
		for (size_t j = 0; j < table[i].vKoef.size(); ++j)
			cout << table[i].vKoef[j] << " ";
		cout << endl;
	}
	
	return 0;
}
