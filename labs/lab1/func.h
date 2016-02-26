#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

const static double eps = 1e-5;

enum TypeProblem { Gauss, Zeydel, SimpleIter };
enum TypeMatrix { Normal, Identity };

int strToInt(char* s);
