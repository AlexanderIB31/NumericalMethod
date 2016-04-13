#include "func.h"
#include "./../lab1/TSolve.h"

// System of functions:
//
// ((x1) ^ 2 + a ^ 2) * x2 - a ^ 3 = 0
// (x1 - a/2) ^ 2 + (x2 - a/2) ^ 2 - a ^ 2 = 0
//

#define DEBUG 1


class BadFileExcept : public exception {
private:
    const char* _fName;
public:
    BadFileExcept(const char* file) {
        _fName = file;
    }

    virtual const char* what() const throw() {
        return strcat("Can`t open file: ", _fName);
    }
};

vector<double> GetAverageVector(const vector<pair<double, double> >& v) {
    vector<double> res;
    for (auto i : v)
        res.push_back((i.first + i.second) / 2);
    return res;
}

class TSolveSystem {
public:
    TSolveSystem(const string& from, const string& to) {
        in.open(from.c_str(), ios::in);
        out.open(to.c_str(), ios::out);
        if (!in.is_open()) {
            throw new BadFileExcept(from.c_str());
        }
        if (!out.is_open()) {
            throw new BadFileExcept(to.c_str());
        }
        in >> EPS;
        while (1) {
            double a, b;
            in >> a >> b;
            if (in.eof())
                break;
            _borders.push_back(make_pair(a, b));
        }
        
        auto f1 = [](double x1, double x2, double a) -> double { 
            return (x1 * x1 + a * a) * x2 - a * a * a; };
        auto f2 = [](double x1, double x2, double a) -> double {
            return pow(x1 - a/2, 2.0) + pow(x2 - a/2, 2.0) - a * a; };
            
        _systemFunctions.push_back(f1); 
        _systemFunctions.push_back(f2); 
    }

    ~TSolveSystem() { in.close(); out.close(); _systemFunctions.clear(); }

    void ToSolveBySimpleIter() {
        ofstream log("solveSimpleIter.log", ios::out);

#ifdef DEBUG
        cout << __LINE__ << endl;
#endif        
        TVector x_0(GetAverageVector(_borders)); 
#ifdef DEBUG
        cout << __LINE__ << endl;
#endif        

        auto CalcJakobian_F = [](double x1, double x2, double a) -> double { 
            return max(abs(2 * x2 * x1), 
                        max(abs(2 * x2 - a), 
                            max(abs(2 * x1 - a), abs(x1 * x1 + a * a)))); 
        };

        double J_0 = CalcJakobian_F(x_0[0], x_0[1], A);
#ifdef DEBUG
        cout << "J_0 = " << J_0 << endl;
#endif        
       
        auto CalcJakobian_fi = [=](double x1, double x2, double a) -> double {
            return max(abs(1.0 - (2 * x2 * x1) / J_0), 
                    max(abs((x1 * x1 + a * a) / J_0), 
                        max(abs((2 * x1 - a) / J_0), abs(1.0 - (2 * x2 - a) / J_0))));
        };

        vector<double (*)(double, double, double)> _systemFI;

        auto fi_1 = [&](double x1, double x2, double a) -> double {
            return x1 - _systemFunctions[0](x1, x2, A) / J_0; 
        };
        
        auto fi_2 = [&](double x1, double x2, double a) -> double {
            return x2 - _systemFunctions[1](x1, x2, A) / J_0; 
        };

        double q;
        try {
            q = CalcJakobian_fi(x_0[0], x_0[1], A);
        }
        catch (const out_of_range& e) {
            cerr << e.what() << endl;
        }
        cout << "q = " << q << endl;        

        auto norm = [](const TVector& x) -> double { 
            double res = 0;
            for (int i = 0; i < x.GetSize(); ++i)
                res += x[i] * x[i];
            return sqrt(res); 
        };

        double eps_k = EPS;    
        TVector solution(x_0.GetSize());

#ifdef DEBUG
        cout << __LINE__ << endl;
#endif        
        for (int iter = 0; iter < 1000; ++iter) {
            solution[0] = fi_1(x_0[0], x_0[1], A);
            solution[1] = fi_2(x_0[0], x_0[1], A);
            log << "x[" << iter << "] = (" << x_0[0] << ", " << x_0[1] << "); ";
            log << "x[" << iter + 1 << "] = (" << solution[0] << ", " << solution[1] << "); ";
            eps_k = q / (1.0 - q) * norm(solution - x_0);
            log << "eps_" << iter << " = " << eps_k << endl;
            if (eps_k < EPS)
                break;
            x_0 = solution;
        }
        log.close();
        out << "solution = (" << solution[0] << ", " << solution[1] << ")" << endl;
    }

    void ToSolveByNewtoon() {
        ofstream log("solveNewtoon.log", ios::out);
#ifdef DEBUG
        cout << __LINE__ << endl;
#endif        
        TVector x_0(GetAverageVector(_borders)); 
		TVector solution(x_0);
	
		string inputStr = "./dependens/inputData";
		string outputStr = "./dependens/outputData";
		size_t szOfMatrix = 2;
		auto print = [&] (double x1, double x2, double a, ofstream& o) {
			o << (2 * x2 * x1) << " " << (x1 * x1 + a * a) << " " << (x1 * x1 + a * a) * x2 - a * a * a << endl;
			o << (2 * x1 - a) << " " << (2 * x2 - a) << " " << (pow(x1 - a/2, 2.0) + pow(x2 - a/2, 2.0) - a * a) << endl;
		};

		auto findMax = [] (const vector<double>& v) -> double {
			double res = abs(v.back());
			for (size_t i = 0; i < v.size(); ++i)
				res = max(res, abs(v[i]));
			return res;
		};

		double eps_k = EPS + 1;
#ifdef DEBUG
        cout << __LINE__ << endl;
#endif        
	    for (size_t iter = 0; iter < 1000; iter++) {
		    log << "x_[" << iter << "] = (" << x_0[0] << ", " << x_0[1] << "); ";
			ofstream o(inputStr.c_str(), ios::out);
		    o << szOfMatrix << endl;
			print(x_0[0], x_0[1], A, o);
		    o.close();
#ifdef DEBUG
        cout << __LINE__ << endl;
#endif        
			TSolve solve(inputStr, outputStr);
			if (!solve.ToSolveByGauss()) {
				ifstream in(outputStr.c_str(), ios::in);
	            vector<double> tempVec;
				double tmp;
				while (1) {
				    in >> tmp;
				    if (in.eof())
				        break;
	                tempVec.push_back(tmp);
				}
				eps_k = findMax(tempVec);
				solution = TVector(tempVec);
				in.close();
		    }
	 	    else {
	           cerr << "Some troubles..." << endl;
		       exit(-1);
	 	    }
#ifdef DEBUG
        cout << __LINE__ << endl;
#endif        
		    solution = x_0 + solution;
		    log << "x_[" << iter + 1 << "] = (" << solution[0] << ", " << solution[1] << "); ";
		    log << "eps[k] = " << eps_k << endl;
	        if (eps_k >= EPS) {
				x_0 = solution;
	        }
	        else break;
	    }
        log.close();
        out << "solve: ";
		solution.Print(out);
    }
private:
    ifstream in;
    ofstream out, log;
    
    vector<double (*)(double, double, double)> _systemFunctions; // F(x) = 0
    vector<pair<double, double> > _borders;
};


int main(int argc, char* argv[]) {
    string pathFrom, pathTo;
    try {
        if (argc != 3) throw out_of_range("You try access to index out of range.");
        pathFrom = argv[1];
        pathTo = argv[2];
    }
    catch (const out_of_range& e) {
        cerr << e.what() << endl;
        return -1;
    }
    try {
        TSolveSystem solve(pathFrom, pathTo);
//        solve.ToSolveBySimpleIter();
        solve.ToSolveByNewtoon();
    }
    catch (const exception& e) {
        cerr << e.what() << endl;
        return -1;
    }
    return 0;
}
