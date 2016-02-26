#pragma once
#include "TVector.h"

class TMatrix : public TCommonClass {
public:
	TMatrix();
    TMatrix(int sz);    
    TMatrix(int sz, double** v);    
    TMatrix(int sz, TypeMatrix typeMatrix);   
    TMatrix(const TMatrix& M);
    TMatrix(ifstream& fstr, TVector& vec, TypeProblem typeProblem);
                
    double* operator[] (int pos) const;    
    TVector operator* (const TVector& vec);
    TMatrix operator* (const TMatrix& matr);    
    TMatrix& operator= (const TMatrix& matr);

	int FindPosMaxElemInColumn(int col);
    int FindDiagElemIsNotNull(int pos);
	void Print();
    void Print(ofstream& o);    
    void Print(ofstream& o, const string& str);
    void SetLink(const TMatrix& matr);     
    void AssignColumn(const TVector& vec, int pos);
    
    virtual void SetSize(int sz);    
    virtual double GetNorma() const;
    virtual void Clear();    
    virtual void SwapRows(int pos1, int pos2);    
    virtual void SwapColumns(int pos1, int pos2);    
};

