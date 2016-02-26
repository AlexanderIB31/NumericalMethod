#include "TVector.h"
    
TVector::TVector() {
	this->_sz = 0;
	this->_vec = NULL;
}

TVector::TVector(int sz) {
    this->_sz = sz;
    this->_vec = new double*[1];
    this->_vec[0] = new double[sz];
	for (int i = 0; i < sz; ++i)
		this->_vec[0][i] = 0.0;
}   

TVector::TVector(int sz, double val) {
    this->_sz = sz;
    this->_vec = new double*[1];
    this->_vec[0] = new double[sz];
    for (int i = 0; i < sz; ++i) 
        this->_vec[0][i] = val;
}    

TVector TVector::operator+ (const TVector& vec) {
    TVector res(this->_sz);
    for (int i = 0; i < this->_sz; ++i) {
        res[i] = this->_vec[0][i] + vec[i];
    }
    return res;     
}

TVector TVector::operator- (const TVector& vec) {
    TVector res(this->_sz);
    for (int i = 0; i < this->_sz; ++i) {
        res[i] = this->_vec[0][i] - vec[i];
    }
    return res;     
}

TVector& TVector::operator= (const TVector& vec) {
    if (!this->_sz) {
        this->_sz = vec.GetSize();
        this->_vec = new double*[1];
        this->_vec[0] = new double[this->_sz];
    }
    else if (this->_sz != vec.GetSize()) {
        cout << "Error!!![1]" << endl;
    }
    for (int i = 0; i < this->_sz; ++i) {
        this->_vec[0][i] = vec[i];
    }
    return *this;
}   

double TVector::operator[] (int pos) const {
    return this->_vec[0][pos];
}

double& TVector::operator[] (int pos) {
    return this->_vec[0][pos];
}

void TVector::SetLink(const TVector& vec) { 
    // be carefull
    if (!!this->_sz) {
        this->Clear();
    }
    this->_sz = vec.GetSize();     
    this->_vec = vec.GetVec(); 
}

void TVector::Print() {
    for (int i = 0; i < this->_sz; ++i) {
        cout << this->_vec[0][i] << (i == this->_sz - 1 ? "" : ", ");
    }
}    

void TVector::Print(ofstream& o) {
    for (int i = 0; i < this->_sz; ++i) {
        o << this->_vec[0][i] << (i == this->_sz - 1 ? "" : ", ");
    }
}    

void TVector::Print(ofstream& o, const string& str) {
    o << str + " = (";
    for (int i = 0; i < this->_sz; ++i) {
        o << this->_vec[0][i] << (i == this->_sz - 1 ? "" : ", ");
    }
    o << ")" << endl << endl;
}

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

void TVector::SwapColumns(int pos1, int pos2) {
    if (pos1 < 0 || pos2 < 0)
        cout << "Error!!![7]" << endl;
    
    if (pos1 == pos2) return;
    
    double tmp = this->_vec[0][pos1];
    this->_vec[0][pos1] = this->_vec[0][pos2];
    this->_vec[0][pos2] = tmp;    
}
