#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <exception>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace std;

static double eps = 0.01;

enum TypeProblem {  Gauss, 
                    Zeydel, 
                    SimpleIter, 
                    TripleDiagMatrix,
                    Rotate,
                    QR  };

enum TypeMatrix { Zero, Identity };

int strToInt(const char* s);
double EvklidNorm(const vector<double>& v);
int sign(double a);
pair<pair<double, double>, char> solveQuadEquation(const vector<double>& v);

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
