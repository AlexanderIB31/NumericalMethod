#pragma once
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

enum TypeSolve { Iter, Newtoon };

static double EPS = .01;

double F(double x);

double DF(double x);

double DDF(double x);

double fi(double x);

double Dfi(double x);
