#include "TMethodRungeKutta.h"

TMethodRungeKutta::TMethodRungeKutta(const string& outFile) {
	out.open(outFile, ios::out);
	log.open("MethodRungeKutta.log", ios::out);
}

void TMethodRungeKutta::ToSolve() {
	string separator = "*************************************************************************";
	double x0 = this->_a, y0 = this->_funcY0;
	size_t N = (this->_b - this->_a) / this->_h;
	double K1, K2, K3, K4;
	log << "h = " << this->_h << "; N = " << N << endl;
	double y_k = y0;
	double x_k = x0;
	double delta;
	out.precision(5);
	log.precision(5);
	for (size_t k = 0; k < N; ++k) {
		K1 = this->_h * this->_FuncExpression(x_k, y_k);
		log << fixed << k << "/" << 1 << "\t" << x_k << "\t" << y_k << "\t" << K1 << "\t\t\t\t" 
			<< this->_FuncY(x_k) << "\t" << fabs(this->_FuncY(x_k) - y_k) << endl;		
		K2 = this->_h * this->_FuncExpression(x_k + 0.5 * this->_h, y_k + 0.5 * K1);
		log << fixed << k << "/" << 2 << "\t" << x_k + 0.5 * this->_h << "\t" << y_k + 0.5 * K1 << "\t" << K2 << endl;
		K3 = this->_h * this->_FuncExpression(x_k + 0.5 * this->_h, y_k + 0.5 * K2);
		log << fixed << k << "/" << 3 << "\t" << x_k + 0.5 * this->_h << "\t" << y_k + 0.5 * K2 << "\t" << K3 << endl;
		K4 = this->_h * this->_FuncExpression(x_k + this->_h, y_k + K3);
		delta = 1.0/6 * (K1 + 2.0 * K2 + 2.0 * K3 + K4);
		log << fixed << k << "/" << 4 << "\t" << x_k + this->_h << "\t" << y_k + K3 << "\t" << K4 
			<< "\t" << delta << "\t" << fabs((K2 - K3) / (K1 - K2)) << endl;
		out << fixed << "x" << k << "\t" << x_k << "\ty" << k << "\t" << y_k << endl;
		x_k += this->_h;
		y_k += delta;		
		log << separator << endl;
	}	
	out << fixed << "x" << N << "\t" << x_k << "\ty" << N << "\t" << y_k << endl;
	log << fixed << N << "\t" << x_k << "\t" << y_k << "\t\t\t\t\t" << this->_FuncY(x_k) << "\t" << fabs(this->_FuncY(x_k) - y_k) << endl;
}

TMethodRungeKutta::~TMethodRungeKutta() {
	out.close();
	log.close();
}