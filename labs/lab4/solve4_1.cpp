#include "dependences/TMethodEuler.h"
#include "dependences/TMethodRungeKutta.h"
#include "dependences/TMethodAdams.h"

int main(int argc, char* argv[]) {
	if (argc != 2) {
		cerr << "Error! Incorrect number of arguments" << endl;
		exit(-1);
	}
	string outFileName = argv[1];
	/* Euler method */
	{
		TMethodEuler solution(outFileName + "Euler");
		solution.ToSolve();
		solution.RungeRomberg();
		cout << "Solution was completed![Euler]" << endl;
	}
	
	/* Runge-Kutta method */
	{
		TMethodRungeKutta solution(outFileName + "Runge-Kutta");
		solution.ToSolve();
		cout << "Solution was completed![Runge-Kutta]" << endl;
	}

	/* Adams method */
	{
		TMethodAdams solution(outFileName + "Adams");
		solution.ToSolve();
		solution.RungeRomberg();
		cout << "Solution was completed![Adams]" << endl;
	}
	
	return 0;
}