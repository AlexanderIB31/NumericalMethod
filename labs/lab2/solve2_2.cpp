#include "func.h"
#include "./../lab1/TSolve.h"

// System of functions:
//
// ((x1) ^ 2 + a ^ 2) * x2 - a ^ 3 = 0
// (x1 - a/2) ^ 2 + (x2 + a/2) ^ 2 - a ^ 2 = 0
//


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
        _a = 0.0;
        _b = 1.0;

        auto f1 = [](double x1, double x2, double a) -> double { 
            return (x1 * x1 + a * a) * x2 - a * a * a; };
        auto f2 = [](double x1, double x2, double a) -> double {
            return pow(x1 - a/2, 2.0) + pow(x2 + a/2, 2.0) - a * a; };
            
        _systemFunctions.push_back(f1); 
        _systemFunctions.push_back(f2); 
    }

    ~TSolveSystem() { in.close(); out.close(); _systemFunctions.clear(); }

    void ToSolveBySimpleIter() {
        ofstream log("solveSimpleIter.log", ios::out);
        _a = 0.0;
        _b = 0.5;
        double x_0 = (_a + _b) / 2; 
        double q = max(abs(Dfi(_b)), abs(Dfi(_a))); // as function monotonously
        double eps_k = EPS;
        double solve = _b + 1.0;
        for (int iter = 0; iter < 1000; ++iter) {
            solve = fi(x_0);
            log << "x_" << iter << " = " << x_0 << "; ";
            log << "x_" << iter + 1 << " = " << solve << "; ";
            eps_k = q / (1 - q) * abs(solve - x_0);
            log << "eps_" << iter << " = " << eps_k << endl;
            if (eps_k < EPS)
                break;
            x_0 = solve;
        }
        log.close();
        out << "solve = " << solve << endl;
    }

    void ToSolveByNewtoon() {
        ofstream log("solveNewtoon.log", ios::out);
        _a = 0.0;
        _b = 1.0;
        double x_0 = (_a + _b) / 2;
        if (F(x_0) * DDF(x_0) <= 0) { // must check this 
            x_0 = _a;
            if (F(x_0) * DDF(x_0) <= 0)
                x_0 = _b;
        }
        log << "x_0 = " << x_0 << endl;
        double solve = _b + 1.0;
        for (size_t iter = 0; iter < 1000; iter++) {
            double x_new = x_0 - F(x_0) / DF(x_0);
            log << "x_" << iter << " = " << x_new << endl;
            if (abs(F(x_new)) < EPS && abs(x_new - x_0) < EPS) {
                solve = x_new;
                break;
            }
            else {
                x_0 = x_new;
            }
        }
        log.close();
        out << "solve: " << solve << endl;
    }
private:
    ifstream in;
    ofstream out, log;
    
    double _a, _b;
    vector<double (*)(double, double, double)> _systemFunctions; // F(x) = 0
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
        solve.ToSolveBySimpleIter();
//        solve.ToSolveByNewtoon();
    }
//    catch (const BadFileExcept& e) {
    catch (const exception& e) {
        cerr << e.what() << endl;
        return -1;
    }
    return 0;
}
