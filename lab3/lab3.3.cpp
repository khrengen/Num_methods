#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "..\lab2\gnuplot.h"
#include "..\lab1\TMatrix.cpp" 

using namespace std;

void ReadFromFile(vector<double>& x, vector<double>& y) {
	int c;
	double val;
	cin >> c;
	x.resize(c);
	y.resize(c);
	for (int i = 0; i < c; ++i) {
		cin >> val;
		x[i] = val;
	}
	for (int i = 0; i < c; ++i) {
		cin >> val;
		y[i] = val;
	}
	return; 
}

double F1(vector<double>& a, double x) {
	return a[0] + a[1] * x;
}

double F2(vector<double>& a, double x) {
	return a[0] + a[1] * x + a[2] * pow(x, 2);	
}

void MNK1 (vector<double>& x, vector<double>& y) {
	vector<double> b(2,0);
	double xSum = 0, xSqSum = 0;
	for (int i = 0; i < x.size(); i++) {
		b[0] += y[i];
		b[1] += x[i] * y[i];
		xSqSum += pow(x[i], 2);
		xSum += x[i];
	}

	TMatrix A(2), L(2), U(2);
	A[0][0] = x.size();
	A[0][1] = xSum;
	A[1][0] = xSum;
	A[1][1] = xSqSum;

	LU(A, L, U);
	for (int i = 0; i < U.swaps.size(); i++) {
        double tmp = b[U.swaps[i].first];
        b[U.swaps[i].first] = b[U.swaps[i].second];
        b[U.swaps[i].second] = tmp;
    }
	vector<double> a = solveOfSystem(L, U, b);

	cout << "polynomial of the first degree:\nF1(x) = " << a[0] << " + x * " << a[1] << '\n';
	printf("xi\t\tF1(xi)\n");
	for (int i = 0; i < x.size(); ++i) {
		printf("%f\t%f\n", x[i], F1(a, x[i]));
	}

	double sse = 0;
	for (int i = 0; i < y.size(); i++) {
		sse += pow(F1(a, x[i]) - y[i], 2);
	}
	cout << "sum of squared errors: " << sse << "\n\n";
}

void MNK2 (vector<double>& x, vector<double>& y) {
	vector<double> b(3,0);
	TMatrix A(3), L(3), U(3);
	A[0][0] = x.size();
	double xSum = 0, xSqSum = 0, xTrSum = 0, xQSum = 0;
	for (int i = 0; i < x.size(); i++) {
		xSum += x[i];
		xSqSum += pow(x[i], 2);
		xTrSum += pow(x[i],3);
		xQSum += pow(x[i],4);
		b[0] += y[i];
		b[1] += y[i]*x[i];
		b[2] += y[i]* pow(x[i],2); 
	}
	A[0][1] = xSum; 
	A[0][2] = xSqSum; 
	A[1][0] = xSum; 
	A[1][1] = xSqSum;
	A[1][2] = xTrSum; 
	A[2][0] = xSqSum; 
	A[2][1] = xTrSum;
	A[2][2] = xQSum;
	LU(A, L, U);
	for (int i = 0; i < U.swaps.size(); i++) {
        double tmp = b[U.swaps[i].first];
        b[U.swaps[i].first] = b[U.swaps[i].second];
        b[U.swaps[i].second] = tmp;
    }
	vector<double> a = solveOfSystem(L, U, b);

	cout << "polynomial of the second degree:\nF2(x) = " << a[0] << " + x * " << a[1] << " + x^2 * " << a[2] << '\n';
	printf("xi\t\tF2(xi)\n");
	for (int i = 0; i < x.size(); ++i) {
		printf("%f\t%f\n", x[i], F2(a, x[i]));
	}

	double sse = 0;
	for (int i = 0; i < y.size(); i++) {
		sse += pow(F2(a, x[i]) - y[i], 2);
	}
	cout << "sum of squared errors: " << sse << '\n';
}

int main() {
	vector<double> x, y;
	ReadFromFile(x, y);
	MNK1(x, y);
	MNK2(x, y);

	Gnuplot plot;
	plot("set xrange [0:+3]");
	plot("set yrange [-1:+11]");
	plot("set pointsize 2");
    plot("plot 6.591912 + x * -3.7283, 10.0988 + x* -14.1074 + x**2 *4.71779, 'test.dat' notitle with points");

	return 0;
}