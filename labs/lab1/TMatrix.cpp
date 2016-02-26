#include "TMatrix.h"

TMatrix::TMatrix() {
	this->_sz = 0;
	this->_vec = NULL;
}

TMatrix::TMatrix(int sz) {
    this->_sz = sz;
    this->_vec = new double*[sz];
    for (int i = 0; i < sz; ++i) {
        this->_vec[i] = new double[sz];
        for (int j = 0; j < sz; ++j) {
            this->_vec[i][j] = 0.0;
        }
    }
}

TMatrix::TMatrix(int sz, double** v) {
    this->_sz = sz;
    this->_vec = new double*[sz];
    for (int i = 0; i < sz; ++i) {
        this->_vec[i] = new double[sz];
        for (int j = 0; j < sz; ++j) {
            this->_vec[i][j] = v[i][j];
        }
    }
}

TMatrix::TMatrix(int sz, TypeMatrix typeMatrix) {
    this->_sz = sz;
    this->_vec = new double*[sz];
    for (int i = 0; i < sz; ++i) {
        this->_vec[i] = new double[sz];
        for (int j = 0; j < sz; ++j)
            this->_vec[i][j] = 0.0;
    }
    
    if (typeMatrix == Identity) {
        for (int i = 0; i < this->_sz; ++i)
            this->_vec[i][i] = 1.0;
    }
}    

TMatrix::TMatrix(const TMatrix& M) {
    this->_sz = M.GetSize();
    this->_vec = new double*[this->_sz];
    for (int i = 0; i < this->_sz; ++i) {
        this->_vec[i] = new double[this->_sz];
        for (int j = 0; j < this->_sz; ++j)
            this->_vec[i][j] = M[i][j];
    }
}

TMatrix::TMatrix(ifstream& fstr, TVector& vec, TypeProblem typeProblem) {
    fstr >> this->_sz;
    this->_vec = new double*[this->_sz];
    vec.SetLink(TVector(this->_sz));
	
	for (int i = 0; i < this->_sz; ++i) {
        this->_vec[i] = new double[this->_sz];
    }
    
    for (int i = 0; i < this->_sz; ++i) {
        for (int j = 0; j < this->_sz; ++j) {
            fstr >> this->_vec[i][j];
        }
        fstr >> vec[i];
    }
	
    if (typeProblem == Zeydel || typeProblem == SimpleIter) {
        for (int i = 0; i < this->_sz; ++i) {
            double div = this->_vec[i][i]; 
            for (int j = 0; j < this->_sz; ++j) {
                this->_vec[i][j] = (i == j ? 0 : -1.0 * this->_vec[i][j] / div);
            }
            vec[i] /= div;
        }
    }
}
    
void TMatrix::SetLink(const TMatrix& matr) { 
    if (!!this->_sz) {
        this->Clear();
    }
    this->_sz = matr.GetSize(); 
    this->_vec = matr.GetVec(); 
}
    
double* TMatrix::operator[] (int pos) const {
    return this->_vec[pos];
}

TVector TMatrix::operator* (const TVector& vec) {
    TVector res(this->_sz);
    
    for (int i = 0; i < this->_sz; ++i) {
        double tmp = 0.0;
        for (int j = 0; j < this->_sz; ++j) {
            tmp += this->_vec[i][j] * vec[j];
        }
        res[i] = tmp;
    }
    return res;
}

TMatrix TMatrix::operator* (const TMatrix& matr) {
    if (matr.GetSize() != this->_sz) {
        cout << "Error!!![5]" << endl;
    }
    
    TMatrix res(matr.GetSize(), Normal);
                
    for (int i = 0; i < this->_sz; ++i) {
        for (int j = 0; j < this->_sz; ++j) {
            double tmp = 0.0;
            for (int k = 0; k < this->_sz; ++k) {
                tmp += this->_vec[j][k] * matr[k][i];
            }
            res[j][i] = tmp;
        }
    }
    
    return res;
}

TMatrix& TMatrix::operator= (const TMatrix& matr) {
	if (this->_sz != matr.GetSize()) {
	this->_sz = matr.GetSize();
	}
	this->_vec = new double*[this->_sz];
	for (int i = 0; i < this->_sz; ++i) {
		this->_vec[i] = new double[this->_sz];
		for (int j = 0; j < this->_sz; ++j)
			this->_vec[i][j] = matr[i][j];
	}
	return *this;	
}

void TMatrix::Print() {
    for (int i = 0; i < this->_sz; ++i) {
        for (int j = 0; j < this->_sz; ++j)
            cout << this->_vec[i][j] << " ";
        cout << endl;
    }
}    

void TMatrix::Print(ofstream& o) {
    for (int i = 0; i < this->_sz; ++i) {
        for (int j = 0; j < this->_sz; ++j)
            o << this->_vec[i][j] << " ";
        o << endl;
    }
}    

void TMatrix::Print(ofstream& o, const string& str) {
    o << "Matrix " + str + ":" << endl;
    for (int i = 0; i < this->_sz; ++i) {
        for (int j = 0; j < this->_sz; ++j)
            o << this->_vec[i][j] << " ";
        o << endl;
    }
    o << endl << endl;
}

int TMatrix::FindPosMaxElemInColumn(int col) {
	double res = -1.0;
	int pos = -1;
	for (int i = col; i < this->_sz; ++i) {
		if (res < abs(this->_vec[i][col])) {
			res = abs(this->_vec[i][col]);
			pos = i;
		}
	}
	return pos;
}

int TMatrix::FindDiagElemIsNotNull(int pos) {
    for (int i = pos + 1; i < this->_sz; ++i) {
        if (!!this->_vec[i][i])
            return i;
    }
	return -1;
}

void TMatrix::AssignColumn(const TVector& vec, int pos) {
	for (int i = 0; i < this->_sz; ++i) {
		this->_vec[i][pos] = vec[i];			
	}
}    

void TMatrix::SwapRows(int pos1, int pos2) {
    if (pos1 < 0 || pos2 < 0)
        cout << "Error!!![4]" << endl;
    
    if (pos1 == pos2) return;
    
    double* tmp = this->_vec[pos1];
    this->_vec[pos1] = this->_vec[pos2];
    this->_vec[pos2] = tmp;
}

void TMatrix::SwapColumns(int pos1, int pos2) {
    if (pos1 < 0 || pos2 < 0)
        cout << "Error!!![4]" << endl;
    
    if (pos1 == pos2) return;
    
    for (int i = 0; i < this->_sz; ++i) {
        double tmp = this->_vec[i][pos1];
        this->_vec[i][pos1] = this->_vec[i][pos2];
        this->_vec[i][pos2] = tmp;
    }
}

void TMatrix::SetSize(int sz) { 
    if (!!this->_sz) {
        this->Clear();
    }
    this->_sz = sz; 
    this->_vec = new double*[sz];
	for (int i = 0; i < sz; ++i)
		this->_vec[i] = new double[sz];
}

void TMatrix::Clear() {
    if (this->_sz > 0) {
        for (int i = 0; i < this->_sz; ++i)
            delete[] this->_vec[i];
        delete[] this->_vec;
        this->_vec = NULL;
        this->_sz = 0;
    }
}   

double TMatrix::GetNorma() const {
   double tmp = 0.0;
   for (int i = 0; i < this->_sz; ++i)
       tmp += abs(this->_vec[0][i]);
   for (int i = 1; i < this->_sz; ++i) {
       double t = 0.0;
       for (int j = 0; j < this->_sz; ++j) {
           t += abs(this->_vec[i][j]);
       }
       tmp = max(tmp, t);
   }
   return tmp;
}
