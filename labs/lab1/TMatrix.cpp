#include "TMatrix.h"

TMatrix::TMatrix() {
    _szRow = _szCol = 0;
    _vec = NULL;
}

TMatrix::TMatrix(double** v = NULL,
                int szR = 0, int szC = 0) {
	_szRow = szR;
    _szCol = szC;
	_vec = v;
}

TMatrix::TMatrix(int szR, int szC, TypeMatrix tM = Zero) {
    _szRow = szR;
    _szCol = szC;
    _vec = new double*[szR];
    for (int i = 0; i < szR; ++i) {
        _vec[i] = new double[szC];
        for (int j = 0; j < szC; ++j) {
            _vec[i][j] = 0.0;
        }
    }
    if (tM == Identity) {
        for (int i = 0; i < min(szR, szC); ++i)
            _vec[i][i] = 1.0;
    }
}

TMatrix::TMatrix(const TMatrix& M) {
    _szRow = M.GetSizeRow();
    _szCol = M.GetSizeCol();
    _vec = new double*[_szRow];
    for (int i = 0; i < _szRow; ++i) {
        _vec[i] = new double[_szCol];
        for (int j = 0; j < _szCol; ++j)
            _vec[i][j] = M[i][j];
    }
}

void TMatrix::SetLink(const TMatrix& matr) { 
    if (!!_szRow && !!_szCol) {
        Clear();
    }
    _szRow = matr.GetSizeRow(); 
    _szCol = matr.GetSizeCol();
    _vec = matr.GetVec(); 
}
    
double* TMatrix::operator[] (int pos) const {
    try {
        return _vec[pos];
    }
    catch (const out_of_range& e) {
        cerr << "Out of range: " << endl; 
    }
}

TMatrix TMatrix::operator* (const TMatrix& matr) const {
    try {
        if (matr.GetSizeRow() != _szCol)
            throw incomp_matrix();
    }
    catch (const incomp_matrix& e) { // must add class incomp_matrix : public exception
        cerr << e.what() << endl;
        return TMatrix(_szRow, _szCol, Zero);
    }
    
    TMatrix res(matr.GetSizeRow(), matr.GetSizeCol(), Zero);                for (int i = 0; i < _szRow; ++i) {
        for (int j = 0; j < _szCol; ++j) {
            double tmp = 0.0;
            for (int k = 0; k < _szRow; ++k) {
                tmp += _vec[j][k] * matr[k][i];
            }
            res[j][i] = tmp;
        }
    }
    return res;
}

TMatrix& TMatrix::operator= (const TMatrix& matr) {
	if (!!_szRow && !!_szCol) {
	    for (int i = 0; i < _szRow; ++i)
            delete[] _vec[i];
        delete[] _vec;
    }
	_szRow = matr.GetSizeRow();
    _szCol = matr.GetSizeCol();
    _vec = new double*[_szRow];
	for (int i = 0; i < _szRow; ++i) {
		_vec[i] = new double[_szCol];
		for (int j = 0; j < _szCol; ++j)
			_vec[i][j] = matr[i][j];
	}
	return *this;	
}

void TMatrix::Print() const {
    for (int i = 0; i < _szRow; ++i) {
        for (int j = 0; j < _szCol; ++j)
            cout << _vec[i][j] << " ";
        cout << endl;
    }
}    

void TMatrix::Print(ofstream& o) const {
    for (int i = 0; i < _szRow; ++i) {
        for (int j = 0; j < _szCol; ++j)
            o << _vec[i][j] << " ";
        o << endl;
    }
}    

void TMatrix::Print(ofstream& o, const string& str) const {
    o << "Matrix " + str + ":" << endl;
    for (int i = 0; i < _szRow; ++i) {
        for (int j = 0; j < _szCol; ++j)
            o << _vec[i][j] << " ";
        o << endl;
    }
    o << endl;
}

int TMatrix::FindPosMaxElemInColumn(int col) const {
    try {
        if (col < 0)
            throw out_of_range("Trying access to element out of range");        
    }
    catch (const out_of_range& e) {
        cerr << "Out of range: " << e.what() << endl;
        return 0;
    }
	double res = -1.0;
	int pos = 0;
	for (int i = col; i < _szRow; ++i) {
		if (res < abs(_vec[i][col])) {
			res = abs(_vec[i][col]);
			pos = i;
		}
	}
	return pos;
}

int TMatrix::FindDiagElemIsNotNull(int pos) const {
    for (int i = pos + 1; i < min(_szRow, _szCol); ++i) {
        if (!!_vec[i][i])
            return i;
    }
	return -1;
}

void TMatrix::AssignColumn(const TMatrix& vec, int pos) {
    try {
        if (pos < 0 || pos >= _szCol)
            throw out_of_range("Trying access to element out of range"); 
    }
    catch (const out_of_range e) {
        cerr << "Out of range: " << e.what() << endl;
        return;
    }

	for (int i = 0; i < _szRow; ++i) {
		_vec[i][pos] = vec[0][i];			
	}
}    

void TMatrix::SwapRows(int pos1, int pos2) const {
    try {
        if (pos1 < 0 || pos2 < 0 || 
            pos1 >= _szRow || pos2 >= _szRow)
            throw out_of_range("Trying access to element out of range"); 
    }
    catch (const out_of_range& e) {
        cerr << "Out of range: " << e.what() << endl;
        return;
    }

    if (pos1 == pos2) return;
    
    double* tmp = _vec[pos1];
    _vec[pos1] = _vec[pos2];
    _vec[pos2] = tmp;
}

void TMatrix::SwapColumns(int pos1, int pos2) const {
    try {
        if (pos1 < 0 || pos2 < 0 || 
            pos1 >= _szCol || pos2 >= _szCol)
            throw out_of_range("Trying access to element out of range"); 
    }
    catch (const out_of_range& e) {
        cerr << "Out of range: " << e.what() << endl;
        return;
    }
   
    if (pos1 == pos2) return;
    
    for (int i = 0; i < _szRow; ++i) {
        double tmp = _vec[i][pos1];
        _vec[i][pos1] = _vec[i][pos2];
        _vec[i][pos2] = tmp;
    }
}
/*
void TMatrix::SetSize(int sz) { 
    if (!!this->_sz) {
        this->Clear();
    }
    this->_sz = sz; 
    this->_vec = new double*[sz];
	for (int i = 0; i < sz; ++i)
		this->_vec[i] = new double[sz];
}
*/
void TMatrix::Clear() {
    if (_szRow > 0) {
        for (int i = 0; i < _szRow; ++i)
            delete[] _vec[i];
        delete[] _vec;
        _vec = NULL;
        _szRow = _szCol = 0;
    }
}   

double TMatrix::GetNorm() const {
   double tmp = 0.0;
   for (int i = 0; i < _szCol; ++i)
       tmp += abs(_vec[0][i]);
   for (int i = 1; i < _szRow; ++i) {
       double t = 0.0;
       for (int j = 0; j < _szCol; ++j) {
           t += abs(_vec[i][j]);
       }
       tmp = max(tmp, t);
   }
   return tmp;
}

int TMatrix::GetSizeRow() const {
    return _szRow;
}

int TMatrix::GetSizeCol() const {
    return _szCol;
}

double** TMatrix::GetVec() const {
    return _vec;
}
