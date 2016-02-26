#pragma once
#include "TCommonClass.h"

class TVector : public TCommonClass {   
public:
	TVector();    
	TVector(int sz);    
    TVector(int sz, double val);
        
    TVector operator + (const TVector& vec);
    TVector operator - (const TVector& vec);    
    TVector& operator= (const TVector& vec);
    double operator[] (int pos) const;
    double& operator[] (int pos);

    void SetLink(const TVector& vec);

	void Print();
    void Print(ofstream& o);
    void Print(ofstream& o, const string& str);        
    
    virtual void SetSize(int sz);    
    virtual double GetNorma() const;
    virtual void Clear();
    virtual void SwapRows(int pos1, int pos2);    
    virtual void SwapColumns(int pos1, int pos2);
};

