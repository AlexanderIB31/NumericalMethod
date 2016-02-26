#include "TSolve.h"

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cout << "Error!!![2]" << endl; // Amount of params is incorrect
        return 0;
    }
    else {
        string pFrom = argv[1];
        string pTo = argv[2];
        TSolve solve(pFrom, pTo);
        int type = strToInt(argv[3]);
        /*
            0 --- LU method
            1 --- TripleDiagMatrix
            2 --- Simple Iterations
            3 --- Zeydel`s method
        */
        switch (type) {
            case 0:
                solve.ToSolveByGauss();
                break;
            case 1:
                solve.ToSolveByTripleDiagMatrix();
                break;
            case 2:
                solve.ToSolveBySimpleIterations();
                break;
            case 3:
                solve.ToSolveByZeydel();
                break;
            default:
                cout << "Error!!![6]" << endl;
                break;
        }
    }
    return 0;
}
