#include "TMethodEuler.h"

TMethodEuler::TMethodEuler(const string& outFile) {
	out.open(outFile, ios::out);
	log.open("MethodEuler.log", ios::out);
}

void TMethodEuler::ToSolve() {
	double x0 = this->_a, y0 = this->_funcY0;
	size_t N = (this->_b - this->_a) / this->_h;
	log << "h = " << this->_h << "; N = " << N << endl;
	double y_k = y0;
	double x_k = x0;
	double delta = this->_h * this->_FuncExpression(x_k, y_k);
	out.precision(5);
	log.precision(5);
	out << fixed << "x0\t" << x0 << ";\ty0\t" << y0 << endl;
	log << 0 << fixed << "\t" << "x0\t" << x0 << ";\ty0\t" << y0 << ";\td(y0)\t" << delta 
		<< ";\ty(true)\t" << this->_FuncY(x_k) << ";\teps(0)\t" << fabs(this->_FuncY(x_k) - y_k) << endl;
	for (size_t k = 1; k < N; ++k) {
		y_k = y_k + delta;
		x_k += this->_h;
		delta = this->_h * this->_FuncExpression(x_k, y_k);
		out << fixed << "x" << k << "\t" << x_k << ";\ty" << k << "\t" << y_k << endl;
		log << k << fixed << "\t" << "x" << k << "\t" << x_k << ";\ty" << k << "\t" << y_k 
			<< ";\td(y" << k << ")\t" << delta << ";\ty(true)\t" << this->_FuncY(x_k) 
			<< ";\teps(" << k << ")\t" << fabs(this->_FuncY(x_k) - y_k) << endl;
	}
	y_k = y_k + delta;
	x_k += this->_h;
	delta = this->_h * this->_FuncExpression(x_k, y_k);
	out << fixed << "x" << N << "\t" << x_k << ";\ty" << N << "\t" << y_k << endl;
	log << N << fixed << "\t" << "x" << N << "\t" << x_k << ";\ty" << N << "\t" << y_k 
			<< ";\t\t\t\ty(true)\t" << this->_FuncY(x_k) 
			<< ";\teps(" << N << ")\t" << fabs(this->_FuncY(x_k) - y_k) << endl;
}

TMethodEuler::~TMethodEuler() {
	out.close();
	log.close();
}