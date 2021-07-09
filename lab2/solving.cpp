#include <iostream>
#include <cmath>
#include "gnuplot.h"

using namespace std;

int sign(double x) {
	return x < 0 ? -1 : 1;
}

double f(double x) {
	return log(x + 2) - pow(x, 4) + 0.5;
}

double phi(double x, double lambda) {
	return x - lambda * f(x);
}

double f1(double x) {
	return 1/(x + 2) - 4 * pow(x, 3);
}

double phi1(double x, double lambda) {
	return 1 - lambda * f1(x);
}

double f2(double x) {
	return 1/pow(x + 2, 2) - 12 * pow(x, 2);
}

double Newton(double a, double b, double eps, int& iter) {
	double epsk = eps + 1;

	double x = f(a)* f2(a) > 0 ? a : b;

	if (f(x) * f2(x) <= 0) {
		cerr << "f(x) * f''(x) <= 0\n";
		exit(1);
	}

	cout << "x0 = " << x << endl;
	double x1;
	while (epsk > eps) {
		x1 = x - f(x) / f1(x);
		epsk = abs(x1 - x);
		x = x1;
		iter++;
	}
	return x;
}

double SimpleIt(double a, double b, double eps, int& iter) {
	
	double x1, x = (a + b) / 2;
	double ex = sign(f1(a));
	for (double i = a; i <= b; i += 0.1) {
		if (sign(f1(i)) != ex) {
			cerr << "sign f' != const\n";
			exit(2);
		}
	}
	double ex2 = sign(f2(a));
	for (double i = a; i <= b; i += 0.1) {
		if (sign(f2(i)) != ex2) {
			cerr << "sign f'' != const\n";
			exit(3);
		}
	}
	double maxF2 = max(abs(f1(a)), abs(f1(b)));
	double lambda = ex / maxF2;
	double q = max(phi1(a, lambda), phi1(b, lambda));
	double epsk = eps + 1;
	while (epsk > eps) {
		x1 = phi(x, lambda);
		epsk = q / (1 - q) * abs(x1 - x);
		x = x1;
		iter++;
	}
	return x;
}

int main() {

	double eps;
	int iter = 0;
	cin >> eps;
	cout << "eps = " << eps << endl << endl;
	double x1 = SimpleIt(1.0, 1.5, eps, iter);
	cout << "start segment [1.0 1.5]\n";
	cout << "Simple Itarations:\nx0 = 1.25\nx = " << x1 << "\nnumber of iterations: " << iter << "\n\n";
	iter = 0;
	cout << "Newton:\n";
	double x = Newton(1.0, 1.5, eps, iter);
	cout << "x = " << x << "\nnumber of iterations: " << iter << "\n\n";
	


	Gnuplot plot;
	plot("set xrange [0:+5]");
	plot("set yrange [-5:+5]");
    plot("plot log(x + 2) - x**4 + 0.5, 0");
    std::cin.get();
	return 0;
}