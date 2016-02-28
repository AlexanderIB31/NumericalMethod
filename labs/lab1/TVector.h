#pragma once
#include "TMatrix.h"

class TVector : public TMatrix {   
public:
    TVector();
	TVector(double** v, int sz);    
	TVector(int sz);    
    TVector(const TMatrix& m);
//    TVector(int sz, double val);
      
    //using TMatrix::operator+;
    //using TMatrix::operator-;
    //using TMatrix::operator=;
    //TVector operator + (const TVector& vec);
    //TVector operator - (const TVector& vec);    
    //TVector& operator= (const TVector& vec);
    double operator[] (int pos) const;
    double& operator[] (int pos);

    int GetSize() const;
    void Print() const;
    void Print(ofstream& o) const;
    virtual void Print(ofstream& o, const string& str) const; 
    void Swap(int pos1, int pos2) const;
};

