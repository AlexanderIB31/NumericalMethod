#include "func.h"

double F(double x) {
    return sqrt(1.0 - x * x) - exp(x) + .1;
}

double DF(double x) {
    return x / sqrt(1.0 - x * x) - exp(x);
}

double DDF(double x) {
    return -1.0 / sqrt(pow(1.0 - x * x, 3.0)) - exp(x);
}

double fi(double x) {
    return x - F(x) / DF(x);
}

double Dfi(double x) {
    return x / (x * x - .1 * sqrt(1.0 - x * x) - 1.0);
}
