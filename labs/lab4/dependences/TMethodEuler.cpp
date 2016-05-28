#include "TMethodEuler.h"
#include <vector>

TMethodEuler::TMethodEuler(const string& outFile) {
	out.open(outFile, ios::out);
	log.open("MethodEuler.log", ios::out);
}

void TMethodEuler::ToSolve() {
	double x0 = this->_a, y0 = this->_funcY0, z0 = this->_funcZ0;
	size_t N = (this->_b - this->_a) / this->_h;
	log << "h = " << this->_h << "; N = " << N << endl;
	double y_k = y0;
	double x_k = x0;
	double z_k = z0;
	double deltaZ = this->_h * this->_FuncExpression(x_k, y_k);
	double deltaY = this->_h * z_k;
	out.precision(5);
	log.precision(5);
	out << "\tx\ty\tz" << endl;
	out << 0 << fixed << "\t" << x0 << "\t" << y0 << "\t" << z0 << endl;
	log << "\tx\ty\tz\t\td(z)\t\td(y)\t\ty(true)\t\teps" << endl;
	log << 0 << fixed << "\t" << x0 << "\t" << y0 << "\t" << z0 << "\t\t" << deltaZ 
		<< "\t\t" << deltaY << "\t\t" << this->_FuncY(x_k) << "\t\t" << fabs(this->_FuncY(x_k) - y_k) << endl;
	for (size_t k = 1; k < N; ++k) {		
		deltaY = this->_h * z_k;
		z_k = z_k + deltaZ;
		y_k = y_k + deltaY;
		x_k += this->_h;
		deltaZ = this->_h * this->_FuncExpression(x_k, y_k);
		out << fixed << k << "\t" << x_k << "\t" << y_k << "\t" << z_k << endl;
		log << k << fixed << "\t" << x_k << "\t" << y_k 
			<< "\t" << z_k << "\t\t" << deltaZ << "\t\t" << deltaY << "\t\t" << this->_FuncY(x_k) 
			<< "\t\t" << fabs(this->_FuncY(x_k) - y_k) << endl;
	}
	deltaY = this->_h * z_k;
	z_k = z_k + deltaZ;	
	y_k = y_k + deltaY;
	x_k += this->_h;
	deltaZ = this->_h * this->_FuncExpression(x_k, y_k);
	out << fixed << N << "\t" << x_k << "\t" << y_k << "\t" << z_k << endl;
	log << N << fixed << "\t" << x_k << "\t" << y_k 
			<< "\t" << z_k << "\t\t\t\t\t\t" << this->_FuncY(x_k) 
			<< "\t\t" << fabs(this->_FuncY(x_k) - y_k) << endl;
}

void TMethodEuler::RungeRomberg() {
	double x0 = this->_a, y0 = this->_funcY0, z0 = this->_funcZ0;
	double h1 = this->_h;
	double h2 = h1 / 2;
	size_t N1 = (this->_b - this->_a) / h1;
	size_t N2 = (this->_b - this->_a) / h2;
	log << "h1 = " << h1 << "; N1 = " << N1 << endl;
	log << "h2 = " << h2 << "; N2 = " << N2 << endl;
	double y_k = y0;
	double x_k = x0;
	double z_k = z0;
	double deltaZ = h1 * this->_FuncExpression(x_k, y_k);
	double deltaY = h1 * z_k;
	out.precision(5);
	log.precision(5);
	vector<double> X, Y1, Y2;
	for (size_t k = 1; k < N1; ++k) {		
		X.push_back(x_k);
		Y1.push_back(y_k);		
		deltaY = h1 * z_k;
		z_k = z_k + deltaZ;
		y_k = y_k + deltaY;
		x_k += h1;
		deltaZ = h1 * this->_FuncExpression(x_k, y_k);
	}
	deltaY = h1 * z_k;
	z_k = z_k + deltaZ;	
	y_k = y_k + deltaY;
	x_k += h1;
	X.push_back(x_k);
	Y1.push_back(y_k);			
	deltaZ = h1 * this->_FuncExpression(x_k, y_k);
	y_k = y0;
	x_k = x0;
	z_k = z0;
	deltaZ = h2 * this->_FuncExpression(x_k, y_k);
	deltaY = h2 * z_k;
	for (size_t k = 1; k < N2; ++k) {		
		if (!(k & 1)) {
			Y2.push_back(y_k);		
		}
		deltaY = h2 * z_k;
		z_k = z_k + deltaZ;
		y_k = y_k + deltaY;
		x_k += h2;
		deltaZ = h2 * this->_FuncExpression(x_k, y_k);
	}
	deltaY = h2 * z_k;
	z_k = z_k + deltaZ;	
	y_k = y_k + deltaY;
	x_k += h2;
	Y2.push_back(y_k);			
	deltaZ = h1 * this->_FuncExpression(x_k, y_k);
	log << "Runge-Romberg" << endl;
	log << "\tx\t\ty\t\ty*\t\terr" << endl;
	y_k = y0;
	for (size_t i = 0; i < min(Y1.size(), Y2.size()); ++i) {
		y_k = Y1[i] + (Y1[i] - Y2[i]);
		log << i << "\t" << X[i] << "\t\t" << Y1[i] << "\t\t" 
			<< y_k << "\t\t" << fabs(Y1[i] - Y2[i]) << endl;
	}
}

TMethodEuler::~TMethodEuler() {
	out.close();
	log.close();
}