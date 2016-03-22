#include "TSolve.h"
int main(int argc, char* argv[]) {
    int type;
    string pFrom, pTo;
    try {
        if (argc == 4) {
            pFrom = argv[1];
            pTo = argv[2];
            type = strToInt(argv[3]);
        }
        else
            throw out_of_range("Number of params isn`t correct");   
    }
    catch (const out_of_range& e) {
        cerr << e.what() << endl;
        return -1;
    }

    TSolve solve(pFrom, pTo);
    /*
            0 --- LU method
            1 --- TripleDiagMatrix
            2 --- Simple Iterations
            3 --- Zeydel`s method
            4 --- Rotate method
            5 --- QR
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
        case 5:
            cout << (!solve.ToSolveByQR() ? "Ok" : "Error...") << endl;
            break;
    }
    return 0;
}
