#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <exception>
#include <stdexcept>
#include <utility>

using namespace std;

const static double eps = 1e-3;

enum TypeProblem {  Gauss, 
                    Zeydel, 
                    SimpleIter, 
                    TripleDiagMatrix,
                    Rotate  };
enum TypeMatrix { Zero, Identity };

int strToInt(char* s);


class incomp_matrix : public exception {
public:
    virtual const char* what() const throw() {
        return "Matrix are incompable";
    }
};

class bad_file_name : public exception {
public:
    virtual const char* what() const throw() {
        return "File was not found";
    }
};
