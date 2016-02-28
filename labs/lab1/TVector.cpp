#include "TVector.h"

TVector::TVector() : TMatrix() { _szRow = 1; }
    
TVector::TVector(double** v, int sz) : TMatrix(v, sz, 1) { }
//    TMatrix::TMatrix(v, sz, 1);
//}

TVector::TVector(int sz) {
    _szCol = sz;
    _szRow = 1;
    _vec = new double*[_szRow];
    _vec[0] = new double[_szCol];
	for (int i = 0; i < _szCol; ++i)
		_vec[0][i] = 0.0;
}

TVector::TVector(const TMatrix& m) {
    _szRow = 1;
    _szCol = m.GetSizeCol();
    _vec = m.GetVec();   
}
/*
TVector::TVector(int sz, double val) {
    this->_sz = sz;
    this->_vec = new double*[1];
    this->_vec[0] = new double[sz];
    for (int i = 0; i < sz; ++i) 
        this->_vec[0][i] = val;
}    
*/
//TVector TVector::operator+ (const TVector& vec) {
//    TVector res(this->_sz);
//    for (int i = 0; i < this->_sz; ++i) {
//        res[i] = this->_vec[0][i] + vec[i];
//    }
//    return res;     
//}

//TVector TVector::operator- (const TVector& vec) {
//    TVector res(this->_sz);
//    for (int i = 0; i < this->_sz; ++i) {
//        res[i] = this->_vec[0][i] - vec[i];
//    }
//    return res;     
//}

//TVector& TVector::operator= (const TVector& vec) {
//    if (!this->_sz) {
//        this->_sz = vec.GetSize();
//        this->_vec = new double*[1];
//        this->_vec[0] = new double[this->_sz];
//    }
//    else if (this->_sz != vec.GetSize()) {
//        cout << "Error!!![1]" << endl;
//    }
//    for (int i = 0; i < this->_sz; ++i) {
//        this->_vec[0][i] = vec[i];
//    }
//    return *this;
//}   

double TVector::operator[] (int pos) const {
    try {
        return _vec[0][pos];
    }
    catch (const out_of_range& e) {
        cerr << "Out of range: " << e.what() << endl;
    }
}

double& TVector::operator[] (int pos) {
    try {
        return this->_vec[0][pos];
    }
    catch (const out_of_range& e) {
        cerr << "Out of ragne: " << e.what() << endl;
    }
}

int TVector::GetSize() const {
    return _szCol;
}
/*
void TVector::SetLink(const TVector& vec) { 
    // be carefull
    if (!!this->_sz) {
        this->Clear();
    }
    this->_sz = vec.GetSize();     
    this->_vec = vec.GetVec(); 
}
*/
void TVector::Print() const {
    TMatrix::Print();
//    for (int i = 0; i < _szCol; ++i) {
//        cout << _vec[0][i] << (i == _szCol - 1 ? "" : ", ");
//    }
}    
/*
void TVector::Print(ofstream& o) {
    for (int i = 0; i < this->_sz; ++i) {
        o << this->_vec[0][i] << (i == this->_sz - 1 ? "" : ", ");
    }
}    
*/
void TVector::Print(ofstream& o) const {
    TMatrix::Print(o);
}

void TVector::Print(ofstream& o, const string& str) const {
    o << "Vector " + str + " = (";
    for (int i = 0; i < _szCol; ++i) {
        o << _vec[0][i] << (i == _szCol - 1 ? "" : ", ");
    }
    o << ")" << endl << endl;
}
/*
void TVector::Clear() {
    if (this->_sz > 0) {
        delete[] this->_vec[0];
        delete[] this->_vec;
        this->_sz = 0;
        this->_vec = NULL;
    }
}

double TVector::GetNorma() const {
    double tmp = -1.0;
    for (int i = 0; i < this->_sz; ++i)
        tmp = max(tmp, abs(this->_vec[0][i]));
    return tmp;
}

void TVector::SetSize(int sz) { 
    if (!!this->_sz) {
        delete[] this->_vec[0];
        delete[] this->_vec;
    }
    this->_sz = sz; 
    this->_vec = new double*[1];
    this->_vec[0] = new double[sz];
    for (int i = 0; i < sz; ++i)
        this->_vec[0][i] = 0.0;
}

void TVector::SwapRows(int pos1, int pos2) {
    if (pos1 < 0 || pos2 < 0)
        cout << "Error!!![7]" << endl;
    
    if (pos1 == pos2) return;
    
    double tmp = this->_vec[0][pos1];
    this->_vec[0][pos1] = this->_vec[0][pos2];
    this->_vec[0][pos2] = tmp;
}
*/
void TVector::Swap(int pos1, int pos2) const {
    try {
        if (pos1 < 0 || pos2 < 0 || pos1 >= _szCol || pos1 >= _szCol)
            throw out_of_range("Trying swap element out of range[TVector]");
    }
    catch (const out_of_range& e) {
        cerr << "Out of range: " << e.what() << endl;
    }
    
    if (pos1 == pos2) return;
    
    double tmp = _vec[0][pos1];
    _vec[0][pos1] = _vec[0][pos2];
    _vec[0][pos2] = tmp;    
}
