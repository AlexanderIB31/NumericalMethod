#include "TSolve.h"

int main(int argc, char* argv[]) {
    int type;
    string pFrom, pTo;
    try {
        pFrom = argv[1];
        pTo = argv[2];
        type = strToInt(argv[3]);
    }
    catch (const out_of_range& e) {
        cerr << "Out of range: " << e.what() << endl;
    }
    TSolve solve(pFrom, pTo);
    /*
            0 --- LU method
            1 --- TripleDiagMatrix
            2 --- Simple Iterations
            3 --- Zeydel`s method
    */
    switch (type) {
        case 0:
            cout << (!solve.ToSolveByGauss() ? "Ok" : "Error...") << endl;
            break;
        case 1:
            cout << (!solve.ToSolveByTripleDiagMatrix() ? "Ok" : "Error...") << endl;
            break;
        case 2:
            cout << (!solve.ToSolveBySimpleIterations() ? "Ok" : "Error...") << endl;
            break;
        case 3:
            cout << (!solve.ToSolveByZeydel() ? "Ok" : "Error...") << endl;
            break;
        case 4:
            cout << (!solve.ToSolveByRotateMethod() ? "Ok" : "Error...") << endl;
            break;
    }
    return 0;
}
