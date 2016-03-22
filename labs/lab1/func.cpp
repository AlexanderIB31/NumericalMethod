#include "func.h"

int strToInt(const char* s) {
    int res = 0;
    while (*s != '\0') {
        res = res * 10 + *s - '0';
        s++;
    }
    return res;
}

double EvklidNorm(const vector<double>& v) {
    double res = 0.0;
    for (size_t i = 0; i < v.size(); ++i) { res += v[i] * v[i];
    }
    return sqrt(res);
}

int sign(double a) {
    return (a >= 0 ? (a == 0 ? 0 : 1) : -1);
}

pair<pair<double, double>, char> solveQuadEquation(const vector<double>& v) {
    double D = v[1] * v[1] - 4.0 * v[0] * v[2];
    if (D > eps) {
        return make_pair(make_pair( (-v[1] + sqrt(D)) / (2 * v[0]), 
                (-v[1] - sqrt(D)) / (2 * v[0])), 0);    
    }
    else if (D < -eps) {
        return make_pair(make_pair(-v[1] / (2 * v[0]), 
                sqrt(-D) / (2 * v[0])), 1);    
    }
    else {
        return make_pair(make_pair(-v[1] / (2 * v[0]), 0.0), 2);
    }
}
