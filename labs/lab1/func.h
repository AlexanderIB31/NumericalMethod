#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <exception>
#include <stdexcept>

using namespace std;

const static double eps = 1e-5;

enum TypeProblem { Gauss, Zeydel, SimpleIter, TripleDiagMatrix };
enum TypeMatrix { Zero, Identity };

int strToInt(char* s);


class incomp_matrix : public exception {
public:
    virtual const char* what() const throw() {
        return "Matrix are incompable";
    }
};
