#include <iostream>
#include <cmath>
#include "gnuplot.h"
#include "..\lab1\TMatrix.cpp"

using namespace std;

double f(double x1, double x2) {
	return 2 * pow(x1,2) - x1 + pow(x2,2) - 1;
}

double g(double x1, double x2) {
	return x2 - tan(x1);
}

vector<double> phi(vector<double> x, TMatrix J1) {
	vector<double> F(2);
	F[0] = f(x[0], x[1]);
	F[1] = g(x[0], x[1]);
	vector<double> sub = J1 * F;
	for (int i = 0; i < x.size(); i++) {
		x[i] -= sub[i];
	}
	return x;
}


double f11(double x1, double x2) {
	return 4*x1 - 1;
}

double f12(double x1, double x2) {
	return 2*x2;
}

double g11(double x1, double x2) {
	return -1/(1 + pow(x1,2));
}

double g12(double x1, double x2) {
	return 1.0;
}

double phi1Max(double x1, double x2, TMatrix J1) {
	vector<double> x(2);
	vector<double> phi1(2);
	vector<double> phi2(2);
	x[0] = 1;
	x[1] = x2;
	vector<double> F(2);
	F[0] = f11(x1, x2);
	F[1] = g11(x1, x2);
	vector<double> sub = J1 * F;
	for (int i = 0; i < x.size(); i++) {
			phi1[i] = x[i] - sub[i];
		}

	x[0] = x1;
	x[1] = 1;
	F[0] = f12(x1, x2);
	F[1] = g12(x1, x2);
	sub = J1 * F;
	for (int i = 0; i < x.size(); i++) {
			phi2[i] = x[i] - sub[i];
		}
	double m1 = max(abs(phi1[0]) + abs(phi2[0]), abs(phi1[1]) + abs(phi2[1]));

	return m1;
}

TMatrix GetJ (double x1, double x2) {
	TMatrix J(2);
	J[0][0] = f11(x1, x2);
	J[0][1] = f12(x1, x2);
	J[1][0] = g11(x1, x2);
	J[1][1] = g12(x1, x2);
	return J;
}


vector<double> Newton(double a1, double b1, double a2, double b2, double eps, int& iter) {
	double epsk = eps + 1;

	vector<double> X;

	vector<double> b(2);
	vector<double> delta;
	X.push_back((a1 + b1) / 2);
	X.push_back((a2 + b2) / 2);

	while (epsk > eps){
		//cout << X[0] << " " << X[1] << '\n';
		TMatrix J = GetJ(X[0], X[1]);
		b[0] = -f(X[0], X[1]);
		b[1] = -g(X[0], X[1]);

		TMatrix U(J.Size()), L(J.Size());
    	LU(J, L, U);
    	delta = solveOfSystem(L, U, b);

    	for (int i = 0; i < delta.size(); ++i) {
    		X[i] += delta[i]; 
    	}

		abs(delta[0]) > abs(delta[1]) ? epsk = abs(delta[0]) : epsk = abs(delta[1]);
		iter++;
	}
	return X;
}

vector<double> SimpleIt(double a1, double a2, double eps, int& iter) {

	double epsk = eps + 1;
	TMatrix J = GetJ(a1, a2);
	TMatrix U(J.Size()), L(J.Size());
    TMatrix J1 = reverse2x2(J);
    vector<double> v(2);
    v[0] = a1;
    v[1] = a2;
    double q = phi1Max(a1,a2,J1);
    if (q > 1) {
    	cerr << "q > 1\n";
    }
    while (epsk > eps) {
    	vector<double> p = phi(v, J1);
    	epsk = q/(1-q) *max(abs(p[0] - v[0]), abs(p[1] - v[1]));
    	v = p;
    	iter++;
	}
	return v;
}

int main() {

	double eps;
	int iter = 0;
	cin >> eps;
	cout << "eps = " << eps << endl << endl;
	vector<double> x = SimpleIt(0.75, 0.75, eps, iter);
	cout << "X0 = [0.75 0.75]\n";
	cout << "Simple Itarations:\nx = " << x[0] << "  " << x[1] << "\nnumber of iterations: " << iter << "\n\n";
	iter = 0;
	cout << "Newton:\n";
	x = Newton(0.5, 1, 0.5, 1, eps, iter);
	cout << "x = " << x[0] << " " << x[1] << "\nnumber of iterations: " << iter << "\n\n";
	


	Gnuplot plot;
	plot("set xrange [0:+5]");
	plot("set yrange [0:+5]");
    plot("plot sqrt(-2*x**2 + x + 1), tan(x)");
    std::cin.get();
	return 0;
}