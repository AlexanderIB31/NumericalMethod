#pragma once
#include "TVector.h"

struct TFRFF { // TypeForReadFromFile
    TMatrix matr;
    TVector vec;    

    TFRFF(const TMatrix& m) {
        matr.SetLink(m);
    }

    TFRFF(const TMatrix& m, const TVector& v) {
        matr.SetLink(m);
        vec.SetLink(v);
    }
};

class TSolve {
private:
    TMatrix _matrA;
    TVector _vecB;
    TVector _vecX;
   
    ofstream log;
    ofstream output;
    ifstream input;    

    string pathFrom, pathTo;

    double _t(const TMatrix& m) const;

    TVector _solveAx_is_b(const TMatrix& mL, const TMatrix& mU, const TVector& vB);
    void _findSolve(double P, double Q, int n, TVector& x, ofstream& log);
    TFRFF* _readFromFile(const string& path, TypeProblem tP);
    void _writeToFile(const string& path);
    void _clear();
public:
    TSolve(const string& pathFrom, const string& pathTo);
    
    int ToSolveByGauss();    
    int ToSolveBySimpleIterations();
    int ToSolveByZeydel();
    int ToSolveByTripleDiagMatrix();
    int ToSolveByRotateMethod();
    int ToSolveByQR();
};

