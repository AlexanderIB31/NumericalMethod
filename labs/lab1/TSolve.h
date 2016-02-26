#pragma once
#include "TMatrix.h"

class TSolve {
private:
    TMatrix _matrA;
    TVector _vecB;
    TVector _vecX;
    
    ofstream out;
    ifstream in;    

    TVector _Solve(const TMatrix& mL, const TMatrix& mU, const TVector& vB);
    void _findSolve(double P, double Q, int n, TVector& x, ofstream& log);
    
public:
    TSolve(const string& pathFrom, const string& pathTo);
    
    void ToSolveBySimpleIterations();
    void ToSolveByZeydel();
    void ToSolveByGauss();    
    void ToSolveByTripleDiagMatrix();
};

