#include "TSolve.h"

TFRFF* TSolve::_readFromFile(const string& path, TypeProblem tP) {
    input.open(path.c_str(), ifstream::in);
    if (!input.is_open()) {
        cerr << "Incorrect file name" << endl;
        return NULL; 
    }
    int tmpSz = -1;
    input >> tmpSz;
    TMatrix matr(tmpSz, (tP == TripleDiagMatrix ? 3 : tmpSz), Zero);
    if (tP == Rotate) {
         for (int i = 0; i < tmpSz; ++i) {
            for (int j = 0; j < matr.GetSizeCol(); ++j) 
                input >> matr[i][j];
        }
        return new TFRFF(matr);
    }
    else {
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
}

double TSolve::_t(const TMatrix& m) const {
    double res = 0.0;
    for (int i = 0; i < m.GetSizeRow(); i++)
        for (int j = i + 1; j < m.GetSizeCol(); j++)
            res += m[i][j] * m[i][j];
    return sqrt(res);
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
    if (tmpRead == NULL)
        return -1;
    _matrA.SetLink(tmpRead->matr);
    _vecB.SetLink(tmpRead->vec);
    if (abs(_matrA.GetNorm()) >= 1.0) {
        cerr << "Error!!!" << endl;
        return -1;
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
    delete tmpRead; 
    return 0;
} 

int TSolve::ToSolveByZeydel() {    
    TFRFF* tmpRead = _readFromFile(pathFrom, Zeydel);
    if (tmpRead == NULL)
        return -1;
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
    delete tmpRead; 
    return 0;
} 

int TSolve::ToSolveByGauss() {                 
    TFRFF* tmpRead = _readFromFile(pathFrom, Gauss);
    if (tmpRead == NULL)
        return -1;
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
    delete tmpRead; 
    return 0;
}

int TSolve::ToSolveByTripleDiagMatrix() {
    TFRFF* tmpRead = _readFromFile(pathFrom, TripleDiagMatrix);
    if (tmpRead == NULL)
        return -1;
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
    delete tmpRead; 
    return 0;
}

int TSolve::ToSolveByRotateMethod() {
    TFRFF* tmpRead = _readFromFile(pathFrom, Rotate);
    if (tmpRead == NULL)
        return -1;
    _matrA.SetLink(tmpRead->matr);
    ofstream log("solve1RotateMethod.log", ios::out);
    log << "|Method Rotate| by Alexander Bales 80-308" << endl << endl;
    TMatrix rotateMatr( _matrA.GetSizeRow(), 
                        _matrA.GetSizeCol(),
                        Identity );
    TMatrix A(_matrA);
    TMatrix OwnVectors(rotateMatr);
    //ofstream urs("tmp.log", ios::out);
    int cnt = 0;
    while (_t(A) > eps) {
        cout << _t(A) << endl;
        log << "|" << cnt++ + 1 << " iteration|" << endl;
        log << "********************************************" << endl;
        A.Print(log, "A");
        pair<int, int> pos = A.FindPosMaxNotDiagElem();
        int i = pos.first, j = pos.second;
        log << "maxPos = (" << i + 1 << "; " << j + 1 << ");" << endl << endl;
        double angel = 0.0;
        try {
            angel = A[i][i] == A[j][j] ? M_PI / 4 :
                            .5 * atan(2 * A[i][j] /
                                    (A[i][i] - A[j][j]));
        }
        catch (const out_of_range& e) {
            cerr << "Out of range: " << e.what() << endl;
        }
        rotateMatr[i][i] = rotateMatr[j][j] = cos(angel);
        rotateMatr[i][j] = -sin(angel);
        rotateMatr[j][i] = -rotateMatr[i][j];
        rotateMatr.Print(log, "rotateMatr");
        A = rotateMatr.Rotate() * A * rotateMatr;
        OwnVectors = OwnVectors * rotateMatr;
        rotateMatr[i][i] = rotateMatr[j][j] = 1.0;
        rotateMatr[i][j] = rotateMatr[j][i] = 0.0;
        OwnVectors.Print(log, "curMultiplyRotateMatrix");
        log << "############################################" << endl << endl;
    }
    _writeToFile(pathTo);
    for (int i = 0; i < min(A.GetSizeRow(), A.GetSizeCol()); i++)
        output << "l_" << i + 1 << " = " << A[i][i] << endl;
    output << endl;
    for (int i = 0; i < OwnVectors.GetSizeCol(); i++) {
        output << "x_" << i << " = (";
        for (int j = 0; j < OwnVectors.GetSizeRow(); j++) {
            output << OwnVectors[j][i];
            if (j != OwnVectors.GetSizeRow() - 1)
                output << ", ";
        }
        output << ");" << endl;
    }
    _clear();
    A.Clear();
    OwnVectors.Clear();
    rotateMatr.Clear(); 
    return 0;
}

int TSolve::ToSolveByQR() {
     TFRFF* tmpRead = _readFromFile(pathFrom, Rotate);
     if (tmpRead == NULL)
         return -1;
    _matrA.SetLink(tmpRead->matr);
    TMatrix A(_matrA);
    ofstream log("solve1QR.log", ios::out);
    log << "|Method QR| by Alexander Bales 80-308" << endl << endl;
    int curSz = min(A.GetSizeRow(), A.GetSizeCol());
    bool check = false;
    TMatrix E(curSz, curSz, Identity);
    vector<pair<pair<double, double>, char> > tmpVecOwnValues(A.GetSizeRow(), make_pair(make_pair(0.0, 0.0), 255));

    while (!check) {
        TMatrix Q(E);
        for (int i = 0; i < curSz; i++) {
            TVector v(A.GetSizeRow());
            vector<double> tmpVec;
            for (int j = i; j < A.GetSizeRow(); j++)
                tmpVec.push_back(A[j][i]);
            v[i] = A[i][i] + sign(A[i][i]) * EvklidNorm(tmpVec);
            for (int j = i + 1; j < A.GetSizeRow(); j++) {
                v[j] = A[j][i];        
            }
            TMatrix H(E - (2.0 / ((v * v.Rotate())[0][0]) * (v.Rotate() * v)) );
            Q = Q * H;
            A = H * A;
            H.Clear();
        }
        A = A * Q;
        Q.Clear();
        check = true;
        for (int i = 0; check && i < A.GetSizeCol(); i++) {
            for (int j = i + 2; check && j < A.GetSizeRow(); j++)
                if (((int)((abs(A[j][i]) + eps / 2) * 100)) / 100.0 > eps) {
                    check = false;
                    break;
                }
            if (check && i + 1 < A.GetSizeRow() 
                    && ((int)((abs(A[i + 1][i]) + eps / 2) * 100)) / 100.0 > eps) { // than find complex
                vector<double> vec;
                vec.push_back(1.0);
                vec.push_back(-(A[i][i] + A[i + 1][i + 1]));
                vec.push_back(A[i][i] * A[i + 1][i + 1] - 
                            A[i + 1][i] * A[i][i + 1]);
                pair<pair<double, double>, char> res = 
                    solveQuadEquation(vec);
                if (res.second != 1) {
                    cerr << "Trouble!!! " << __LINE__ << endl;    
                }
                if (tmpVecOwnValues[i].second != 255) {
                    res.second = -1;
                    tmpVecOwnValues[i] = res;
                    res.second = -2;
                    tmpVecOwnValues[i + 1] = res;
                }
                else if (tmpVecOwnValues[i].first != res.first){
                    cerr << "Trouble!!! " << __LINE__ << endl;          
                }
                i++;
            }
        }   
        A.Print(log);
        log << endl;
    }
    E.Clear();
    A.Print(log);
    _writeToFile(pathTo);
    for (int i = 0; i < A.GetSizeRow(); i++) {        
        if (i + 1 < A.GetSizeRow() && ((int)((abs(A[i + 1][i]) + eps / 2) * 100)) / 100.0 > eps) {
            vector<double> vec;
            vec.push_back(1.0);
            vec.push_back(-(A[i][i] + A[i + 1][i + 1]));
            vec.push_back(A[i][i] * A[i + 1][i + 1] - 
                        A[i + 1][i] * A[i][i + 1]);
            pair<pair<double, double>, char> res = 
                solveQuadEquation(vec);

            if (res.second != 1) {
                cerr << "Troubles..." << endl;
            }
            else if (tmpVecOwnValues[i].first == res.first) {
                output << "own value[" << i + 1 << "] = " 
                    << res.first.first << " + " << res.first.second 
                    << "i" << endl;
                output << "own value[" << i + 2 << "] = "
                    << res.first.first << " - " << res.first.second
                    << "i" << endl;
                i++;
            }
        }
        else {
            tmpVecOwnValues[i] = make_pair(make_pair(A[i][i], 0.0), 2);
            output << "own value[" << i + 1 << "] = " << A[i][i] << endl;
        }
    }
    _clear();
    return 0;
}
