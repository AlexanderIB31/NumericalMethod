#include "TCommonClass.h"

TCommonClass::TCommonClass() {
    this->_sz = 0;
    this->_vec = NULL;
}

int TCommonClass::GetSize() const {
    return this->_sz;
}

double** TCommonClass::GetVec() const {
    return this->_vec;
}
