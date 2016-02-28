#include "TSolve.h"

TFRFF* TSolve::_readFromFile(const string& path, TypeProblem tP) {
    input.open(path.c_str(), ifstream::in);
    if (!input.is_open()) {
        cerr << "Incorrect file name" << endl;
        return 0; 
    }
    int tmpSz = -1;
    input >> tmpSz;
    TMatrix matr(tmpSz, (tP == TripleDiagMatrix ? 3 : tmpSz), Zero);
    TVector vec(tmpSz);
    for (int i = 0; i < tmpSz; ++i) {
        for (int j = 0; j < matr.GetSizeCol(); ++j) {
            input >> matr[i][j];
        }
        input >> vec[i];
    }
    if (tP == Zeydel || tP == SimpleIter) {
        for (int i = 0; i < tmpSz; ++i) {
            for (int j = 0; j < tmpSz; ++j)
                if (i != j)
                    matr[i][j] /= -1.0 * matr[i][i];
            vec[i] /= matr[i][i];
            matr[i][i] = 0.0;
        }
    }
    else if (tP != Gauss && tP != TripleDiagMatrix) {
        cerr << "This type of problem was not found!" << endl;
    }
    return new TFRFF(matr, vec);
}

void TSolve::_writeToFile(const string& path) {
    output.open(path.c_str(), ofstream::out);
    _matrA.Print(output, "A");
    _vecB.Print(output, "B");
    _vecX.Print(output, "X");
}

void TSolve::_clear() {
    _matrA.Clear();
    _vecB.Clear();
    _vecX.Clear();
    input.close();
    output.close();
    log.close();
}

TVector TSolve::_solveAx_is_b(const TMatrix& mL, 
                            const TMatrix& mU, 
                            const TVector& vB) {
    int tmpSz = vB.GetSize();    
	TVector z(tmpSz);
    TVector x(tmpSz);
	z[0] = vB[0];	
    for (int i = 1; i < tmpSz; ++i) {
        double tmpVal = 0.0;
        for (int j = 0; j < i; ++j)
            tmpVal += mL[i][j] * z[j];
        z[i] = vB[i] - tmpVal;
    }    
    x[tmpSz - 1] = z[tmpSz - 1] / mU[tmpSz - 1][tmpSz - 1];    
    for (int i = tmpSz - 1; i >= 0; --i) {
        double tmpVal = 0.0;
        for (int j = i + 1; j < tmpSz; ++j) {
            tmpVal += mU[i][j] * x[j];
        }
        x[i] = (z[i] - tmpVal) / mU[i][i];
    }
    z.Clear();
    return x;
}


void TSolve::_findSolve(double P, double Q, int n, TVector& x, ofstream& log) {
    if (n == x.GetSize()) {
        x[n - 1] = Q;
        log << "P_" << n << " = " << 0 << endl;
        log << "Q_" << n << " = " << Q << endl;
    } else {
        try {
//            _findSolve(-_matrA[n][n + 1] / (_matrA[n][n - 1] * P + _matrA[n][n]),
//                        (_vecB[n] - _matrA[n][n - 1] * Q) / (_matrA[n][n - 1] * P 
//                            + _matrA[n][n]),
//                        n + 1, x, log);
            _findSolve(-_matrA[n][2] / (_matrA[n][0] * P + _matrA[n][1]),
                        (_vecB[n] - _matrA[n][0] * Q) / (_matrA[n][0] * P + _matrA[n][1]),
                        n + 1, x, log);
            log << "P_" << n << " = " << P << endl;
            log << "Q_" << n << " = " << Q << endl;
            x[n - 1] = P * x[n] + Q;           
        }
        catch (const out_of_range& e) {
            cerr << e.what() << endl;
        }
    }        
}

TSolve::TSolve(const string& pFrom, const string& pTo) {
    pathFrom = pFrom;
    pathTo = pTo;
    output.precision(3);
    log.precision(3);
}

int TSolve::ToSolveBySimpleIterations() {
    TFRFF* tmpRead = _readFromFile(pathFrom, SimpleIter);
    _matrA.SetLink(tmpRead->matr);
    _vecB.SetLink(tmpRead->vec);
    if (abs(_matrA.GetNorm()) >= 1.0) {
        cerr << "Error!!!" << endl;
        return 1;
    }
    log.open("solve1SimpleIter.log", ofstream::out);
    log << "|Method Simple Iterations| by Alexander Bales 80-308" << endl << endl;
    log << "|A| = " << _matrA.GetNorm() << endl << endl;
    _vecX = _vecB;    
    TVector vecRes = _vecX + _matrA * _vecX;    
    _vecB.Print(log, "B");    
    _matrA.Print(log, "A");    
    double tmp = abs(_matrA.GetNorm());
    double eps_k = 1.0 * tmp / (1.0 - tmp) * (vecRes - _vecX).GetNorm();	
    for (int i = 1; eps_k > eps; ++i) {
        _vecX = vecRes;
        log << "x_" << i - 1 << " = (";
        _vecX.Print(log);
        log << "); ";	    
        vecRes = _vecB + _matrA * _vecX;
        log << "x_" << i << " = (";
        vecRes.Print(log);
        log << "); ";
        eps_k = tmp / (1.0 - tmp) * (vecRes - _vecX).GetNorm();
        log << "eps(" << i << ") = " << eps_k << endl;
    }    
    for (int i = 0; i < vecRes.GetSize(); ++i) {
        output << vecRes[i] << endl;
    }
    _writeToFile(pathTo);
    _clear();   
    return 0;
} 

int TSolve::ToSolveByZeydel() {    
    TFRFF* tmpRead = _readFromFile(pathFrom, Zeydel);
    _matrA.SetLink(tmpRead->matr);
    _vecB.SetLink(tmpRead->vec);
    //if (abs(_matrA.GetNorma()) >= 1.0) cout << "Error!!![3]" << endl;
    log.open("solve1Zeydel.log", ofstream::out);
    log << "|Method Zeydel| by Alexander Bales 80-308" << endl << endl;
    //log << "|L| = " << _matrA.GetNorma() << endl;
    _vecX = _vecB;    
    TVector curVec;
    curVec = _vecB;    
    for (int i = 0; i < _vecB.GetSize(); ++i) {
        double tmp = _vecB[i];
        for (int j = 0; j < _vecB.GetSize(); ++j) {
            if (j == i)
                continue;
            else 
                tmp += _matrA[i][j] * curVec[j];
        }
        curVec[i] = tmp;
    }    
    _vecB.Print(log, "B");          
    _matrA.Print(log, "A");
    double tmp1 = abs(_vecB.GetNorm()); 
    double tmp2 = abs(_matrA.GetNorm());
    double eps_k = 1.0 * tmp1 / (1.0 - tmp2) * (curVec - _vecX).GetNorm();	
    for (int i = 1; eps_k > eps; ++i) {
        _vecX = curVec;
        log << "x_" << i - 1 << " = (";
        _vecX.Print(log);
        log << "); ";	    
        for (int k = 0; k < _vecB.GetSize(); ++k) {
            double tmp = _vecB[k];
            for (int j = 0; j < _vecB.GetSize(); ++j) {
                if (j == k)
                    continue;
                else
                    tmp += _matrA[k][j] * curVec[j];
            }
            curVec[k] = tmp;
        }
        log << "x_" << i << " = (";
        curVec.Print(log);
        log << "); ";
        eps_k = tmp1 / (1.0 - tmp2) * (curVec - _vecX).GetNorm();
        log << "eps(" << i << ") = " << eps_k << endl;
    }    
    _vecX = curVec;
    curVec.Clear();    
    _writeToFile(pathTo);
    _clear();
    return 0;
} 

int TSolve::ToSolveByGauss() {                 
    TFRFF* tmpRead = _readFromFile(pathFrom, Gauss);
    _matrA.SetLink(tmpRead->matr);
    _vecB.SetLink(tmpRead->vec);
    ofstream log("solve1Gauss.log", ios::out);
    log << "|Method Gauss (LUP)| by Alexander Bales 80-308" << endl << endl;    
    TMatrix L(_matrA.GetSizeRow(), _matrA.GetSizeCol(), Identity); 
    TMatrix U(_matrA);           
    int cntSwitchRowsAndColumns = 0;
    int posPrecc = -1;
    int tmpSz = min(U.GetSizeRow(), U.GetSizeCol());    
    for (int i = 0; i < tmpSz; ++i) {
        int posMax = U.FindPosMaxElemInColumn(i);           
        if (tmpSz - 1 != i) {
            U.SwapRows(posMax, i);            
            _matrA.SwapRows(posMax, i);
            _vecB.Swap(posMax, i); 
            if (posMax != i)
                cntSwitchRowsAndColumns++;                
            if (posPrecc != -1) {
                L.SwapRows(posMax, i);
                L.SwapColumns(posMax, i);                            
            }               
        }
                                
        for (int j = i + 1; j < tmpSz; ++j) {
            double koef = - U[j][i] / U[i][i];
            U[j][i] = 0.0;
            for (int k = i + 1; k < tmpSz; ++k) {
                U[j][k] = U[j][k] + koef * U[i][k];
            }
            L[j][i] = -koef;
        }
        posPrecc = posMax;
    }

    _vecX.SetLink(_solveAx_is_b(L, U, _vecB));	    
    L.Print(log, "L");
    U.Print(log, "U");    
    output << "det(A) = ";
    double detA = 1.0;
    for (int i = 0; i < tmpSz; ++i)
        detA *= U[i][i];
    output << pow(-1.0, 1.0 * cntSwitchRowsAndColumns) * detA << endl << endl;              
    TMatrix reverseA(_matrA.GetSizeRow(), _matrA.GetSizeCol(), Zero);
    for (int i = 0; i < tmpSz; ++i) {
        TVector vec(_matrA.GetSizeRow());
        TVector calcVec;        
        vec[i] = 1.0;        
        calcVec.SetLink(_solveAx_is_b(L, U, vec));		        
        reverseA.AssignColumn(calcVec, i);
        calcVec.Clear();
        vec.Clear();
    }
    _writeToFile(pathTo); 
    reverseA.Print(output, "A^(-1)");          
    TMatrix check(_matrA * reverseA);
    check.Print(log, "A * A^(-1)");    
    L.Clear();
    U.Clear();
    reverseA.Clear();
    check.Clear();
    _clear();
    return 0;
}

int TSolve::ToSolveByTripleDiagMatrix() {
    TFRFF* tmpRead = _readFromFile(pathFrom, TripleDiagMatrix);
    _matrA.SetLink(tmpRead->matr);
    _vecB.SetLink(tmpRead->vec);
    ofstream log("solve1TripleDiagMatrix.log", ios::out);
    log << "|Method TripleDiagMatrix| by Alexander Bales 80-308" << endl << endl;
    double P, Q;
    try {
        P = -_matrA[0][1] / _matrA[0][0];
        Q = _vecB[0] / _matrA[0][0];            
    }
    catch (const out_of_range& e) {
        cerr << "Out of range: " << e.what() << endl;
    }
    _vecX = _vecB;    
    _findSolve(P, Q, 1, _vecX, log);        
    _writeToFile(pathTo);
    _clear();
    return 0;
}
