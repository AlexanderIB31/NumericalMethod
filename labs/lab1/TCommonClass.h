#pragma once
#include "func.h"

class TCommonClass {
protected:
    double** _vec;
    int _sz;
    
public: 
    TCommonClass();
    
    int GetSize() const;
    double** GetVec() const;
    
    virtual void SetSize(int sz) = 0;    
    virtual double GetNorma() const = 0;
    virtual void Clear() = 0;
    virtual void SwapRows(int pos1, int pos2) = 0;
    virtual void SwapColumns(int pos1, int pos2) = 0;    
};
