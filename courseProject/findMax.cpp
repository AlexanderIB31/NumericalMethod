#include <iostream>
#include <cmath>

double Dfi(double x1, double x2, double a) {
    return 1.0 - (2 * x1 * x2 + 2 * x2 - a) / 12.25 + 2 * x2 * x1 * (2 * x2 - a) / (12.25 * 12.25) - (2 * x1 * x1 * x1 + 2 * x1 * a * a - a * x1 * x1 - a * a * a) / (12.25 * 12.25);
}

struct TPair {
    double x1, x2;
};

int main() {
    TPair solution;
    double *max = NULL;
    double a = 4;
    double a1 = 5.0;//-2.0;
    double b1 = 6.0;//-1.0;
    double a2 = 1.0;//3.0;
    double b2 = 2.0;//4.0;

    double delta = 0.001;
    double eps = 1e-5;

    for (double stepX1 = a1; std::abs(stepX1 - b1) >= eps; stepX1 += delta) {
        for (double stepX2 = a2; std::abs(stepX2 - b2) >= eps; stepX2 += delta) {
            if (max == NULL) {
                max = new double; 
                *max = Dfi(stepX1, stepX2, a);
                solution.x1 = stepX1;
                solution.x2 = stepX2;
            }
            else {
                double res = Dfi(stepX1, stepX2, a);
                if (*max < res) {
                    *max = res;
                    solution.x1 = stepX1;
                    solution.x2 = stepX2;
                }
            }
        }
    }
    std::cout << solution.x1 << " | " << solution.x2 << std::endl; 
    return 0;
}
