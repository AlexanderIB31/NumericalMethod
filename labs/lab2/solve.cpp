#include "func.h"

//extern double EPS;

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

class TSolve {
public:
    TSolve(const string& from, const string& to) {
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
    }

    ~TSolve() { in.close(); out.close(); }

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
        TSolve solve(pathFrom, pathTo);
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
