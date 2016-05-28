#include "TMethodAdams.h"
#include <vector>

TMethodAdams::TMethodAdams(const string& outFile) {
	out.open(outFile, ios::out);
	log.open("MethodAdams.log", ios::out);
}

void TMethodAdams::ToSolve() {
	double x0 = this->_a, y0 = this->_funcY0, z0 = this->_funcZ0;
	size_t N = (this->_b - this->_a) / this->_h;
	double K1, K2, K3, K4, L1, L2, L3, L4;
	log << "h = " << this->_h << "; N = " << N << endl;
	double y_k = y0;
	double x_k = x0;
	double z_k = z0;
	double deltaZ, deltaY;
	out.precision(5);
	log.precision(5);
	log.width(10);
	vector<double> X, Y, Z;
	out << "\tx\ty" << endl;
	log << "\tx\t\ty\t\tz\t\td(y)\t\td(z)\t\ty(true)\t\teps\t\terr" << endl;
	size_t sz = min((size_t)4, N);
	for (size_t k = 0; k < sz; ++k) {
		L1 = this->_h * this->_FuncExpression(x_k, y_k);
		K1 = this->_h * z_k;		
		//log << fixed << k << "/" << 1 << "\t" << x_k << "\t" << y_k << "\t" << K1 << "\t\t\t\t" 
		//	<< this->_FuncY(x_k) << "\t" << fabs(this->_FuncY(x_k) - y_k) << endl;		
		K2 = this->_h * (z_k + 0.5 * L1);
		L2 = this->_h * this->_FuncExpression(x_k + 0.5 * this->_h, y_k + 0.5 * K1);
		//log << fixed << k << "/" << 2 << "\t" << x_k + 0.5 * this->_h << "\t" << y_k + 0.5 * K1 << "\t" << K2 << endl;
		K3 = this->_h * (z_k + 0.5 * L2);
		L3 = this->_h * this->_FuncExpression(x_k + 0.5 * this->_h, y_k + 0.5 * K2);
		//log << fixed << k << "/" << 3 << "\t" << x_k + 0.5 * this->_h << "\t" << y_k + 0.5 * K2 << "\t" << K3 << endl;
		K4 = this->_h * (z_k + L3);
		L4 = this->_h * this->_FuncExpression(x_k + this->_h, y_k + K3);
		deltaY = 1.0 / 6.0 * (K1 + 2.0 * K2 + 2.0 * K3 + K4);
		deltaZ = 1.0 / 6.0 * (L1 + 2.0 * L2 + 2.0 * L3 + L4);
		//log << fixed << k << "/" << 4 << "\t" << x_k + this->_h << "\t" << y_k + K3 << "\t" << K4 
		//	<< "\t" << delta << "\t" << fabs((K2 - K3) / (K1 - K2)) << endl;
		out << fixed << k << "\t" << x_k << "\t" << y_k << endl;
		log << fixed << k << "\t" << x_k << "\t\t" << y_k << "\t\t" << z_k 
			<< "\t\t" << deltaY << "\t\t" << deltaZ << "\t\t" << this->_FuncY(x_k) 
			<< "\t\t" << fabs(this->_FuncY(x_k) - y_k) << "\t\t" << fabs((K2 - K3) / (K1 - K2)) << endl;
		Y.push_back(y_k);
		X.push_back(x_k);
		Z.push_back(z_k);			
		x_k += this->_h;
		z_k += deltaZ;
		y_k += deltaY;		
	}		
	x_k -= this->_h;
	y_k -= deltaY;
	z_k -= deltaZ;
	for (size_t k = sz - 1; k < N; ++k) {
		z_k = z_k + this->_h / 24.0 * (55.0 * this->_FuncExpression(X[k], Y[k]) 
				- 59.0 * this->_FuncExpression(X[k - 1], Y[k - 1]) + 37.0 * this->_FuncExpression(X[k - 2], Y[k - 2]) 
				- 9.0 * this->_FuncExpression(X[k - 3], Y[k - 3]));
		y_k = y_k + this->_h * z_k;
		x_k += this->_h;
		log << fixed << k + 1 << "\t" << x_k << "\t\t" << y_k << "\t\t" << z_k 
			<< "\t\t" << deltaY << "\t\t" << deltaZ << "\t\t" << this->_FuncY(x_k) 
			<< "\t\t" << fabs(this->_FuncY(x_k) - y_k) << endl;
		out << fixed << "\t" << x_k << "\t" << y_k << endl;			
		Y.push_back(y_k);
		X.push_back(x_k);
		Z.push_back(z_k);
	}	
}

void TMethodAdams::RungeRomberg() {
	double x0 = this->_a, y0 = this->_funcY0, z0 = this->_funcZ0;
	double h1 = this->_h;
	double h2 = h1 / 2;
	size_t N1 = (this->_b - this->_a) / h1;
	size_t N2 = (this->_b - this->_a) / h2;
	double K1, K2, K3, K4, L1, L2, L3, L4;
	log << "h1 = " << h1 << "; N1 = " << N1 << endl;
	log << "h2 = " << h2 << "; N2 = " << N2 << endl;
	double y_k = y0;
	double x_k = x0;
	double z_k = z0;
	double deltaZ, deltaY;
	out.precision(5);
	log.precision(5);
	log.width(10);
	vector<double> X, Y, Z;
	size_t sz = min((size_t)4, N1);
	for (size_t k = 0; k < sz; ++k) {
		L1 = h1 * this->_FuncExpression(x_k, y_k);
		K1 = h1 * z_k;		
		K2 = h1 * (z_k + 0.5 * L1);
		L2 = h1 * this->_FuncExpression(x_k + 0.5 * h1, y_k + 0.5 * K1);
		K3 = h1 * (z_k + 0.5 * L2);
		L3 = h1 * this->_FuncExpression(x_k + 0.5 * h1, y_k + 0.5 * K2);
		K4 = h1 * (z_k + L3);
		L4 = h1 * this->_FuncExpression(x_k + h1, y_k + K3);
		deltaY = 1.0 / 6.0 * (K1 + 2.0 * K2 + 2.0 * K3 + K4);
		deltaZ = 1.0 / 6.0 * (L1 + 2.0 * L2 + 2.0 * L3 + L4);
		Y.push_back(y_k);
		X.push_back(x_k);
		Z.push_back(z_k);			
		x_k += h1;
		z_k += deltaZ;
		y_k += deltaY;		
	}		
	x_k -= h1;
	y_k -= deltaY;
	z_k -= deltaZ;
	for (size_t k = sz - 1; k < N1; ++k) {
		z_k = z_k + h1 / 24.0 * (55.0 * this->_FuncExpression(X[k], Y[k]) 
				- 59.0 * this->_FuncExpression(X[k - 1], Y[k - 1]) + 37.0 * this->_FuncExpression(X[k - 2], Y[k - 2]) 
				- 9.0 * this->_FuncExpression(X[k - 3], Y[k - 3]));
		y_k = y_k + h1 * z_k;
		x_k += h1;
		Y.push_back(y_k);
		X.push_back(x_k);
		Z.push_back(z_k);
	}

	vector<double> X2, Y2;
	y_k = y0;
	x_k = x0;
	z_k = z0;
	sz = min((size_t)4, N2);
	for (size_t k = 0; k < sz; ++k) {
		L1 = h2 * this->_FuncExpression(x_k, y_k);
		K1 = h2 * z_k;		
		K2 = h2 * (z_k + 0.5 * L1);
		L2 = h2 * this->_FuncExpression(x_k + 0.5 * h2, y_k + 0.5 * K1);
		K3 = h2 * (z_k + 0.5 * L2);
		L3 = h2 * this->_FuncExpression(x_k + 0.5 * h2, y_k + 0.5 * K2);
		K4 = h2 * (z_k + L3);
		L4 = h2 * this->_FuncExpression(x_k + h2, y_k + K3);
		deltaY = 1.0 / 6.0 * (K1 + 2.0 * K2 + 2.0 * K3 + K4);
		deltaZ = 1.0 / 6.0 * (L1 + 2.0 * L2 + 2.0 * L3 + L4);
		Y2.push_back(y_k);
		X2.push_back(x_k);
		x_k += h2;
		z_k += deltaZ;
		y_k += deltaY;		
	}		
	x_k -= h2;
	y_k -= deltaY;
	z_k -= deltaZ;
	for (size_t k = sz - 1; k < N2; ++k) {
		z_k = z_k + h2 / 24.0 * (55.0 * this->_FuncExpression(X2[k], Y2[k]) 
				- 59.0 * this->_FuncExpression(X2[k - 1], Y2[k - 1]) + 37.0 * this->_FuncExpression(X2[k - 2], Y2[k - 2]) 
				- 9.0 * this->_FuncExpression(X2[k - 3], Y2[k - 3]));
		y_k = y_k + h2 * z_k;
		x_k += h2;
		Y2.push_back(y_k);
		X2.push_back(x_k);
	}
	log << "Runge-Romberg" << endl;
	log << "\tx\t\ty\t\ty*\t\terr" << endl;
	for (size_t i = 0, j = 0; i < Y.size() && j < Y2.size(); ++i, j += 2) {
		y_k = Y[i] + (Y[i] - Y2[j]) / 15;
		log << i << "\t" << X[i] << "\t\t" << Y[i] << "\t\t" 
			<< y_k << "\t\t" << fabs(Y[i] - Y2[j]) / 15 << endl;

	}

}


TMethodAdams::~TMethodAdams() {
	out.close();
	log.close();
}