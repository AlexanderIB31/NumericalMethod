#pragma once
#include "func.h"

class TMatrix {
protected:
    double** _vec;
    int _szRow, _szCol;    
public:
    TMatrix();
	TMatrix(double** v, int szR, int szC);
//    TMatrix(int sz);    
//    TMatrix(int sz, double** v);    
    TMatrix(int szR, int szC, TypeMatrix typeMatrix);   
    TMatrix(const TMatrix& M);
//    TMatrix(ifstream& fstr, TVector& vec, TypeProblem typeProblem);
                
    TMatrix operator* (const TMatrix& m) const;
    TMatrix operator+ (const TMatrix& m) const;
    TMatrix operator- (const TMatrix& m) const;
    double* operator[] (int pos) const;    
    TMatrix& operator= (const TMatrix& matr);

	int FindPosMaxElemInColumn(int col) const;
    pair<int, int> FindPosMaxNotDiagElem() const;
    int FindDiagElemIsNotNull(int pos) const;
	void Print() const;
    void Print(ofstream& o) const;    
    virtual void Print(ofstream& o, const string& str) const;
    void SetLink(const TMatrix& matr);     
    void AssignColumn(const TMatrix& vec, int pos);
    TMatrix Rotate() const;
    
    int GetSizeRow() const;
    int GetSizeCol() const;
    double** GetVec() const;
    double GetNorm() const;
    
    void Clear();    
    void SwapRows(int pos1, int pos2) const;    
    void SwapColumns(int pos1, int pos2) const;    
};

