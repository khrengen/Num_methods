#include <iostream>
#include <vector>
#include <cmath>


using namespace std;

vector<double> ReadX () {
	int num;
	cin >> num;
	vector<double> v(num);
	for (int i = 0; i < num; i++) {
		cin >> v[i];
	}
	return v;
}

double f(double x) {
	return 1 / x;
}

double f(double x1, double x2) {
	return (f(x1) - f(x2)) / (x1 - x2);
}

double f(double x1, double x2, double x3) {
	return (f(x1, x2) - f(x2, x3)) / (x1 - x3);
}

double f(double x1, double x2, double x3, double x4) {
	return (f(x1, x2, x3) - f(x2, x3, x4)) / (x1 - x4);
}

double NewtPol (double X, vector<vector<double>>& vals, vector<double>& x) {
	return vals[0][0] + (X - x[0]) * vals[1][0] + (X - x[0])*(X - x[1]) * vals[2][0] + 
		(X - x[0])*(X - x[1])*(X - x[2]) * vals[3][0];
}

double LagPol (double X, vector<double>& x, vector<double>& wi) {
	return f(x[0]) / wi[0] * (X - x[1]) * (X - x[2])*(X - x[3]) + 
			f(x[1]) / wi[1] * (X - x[0])*(X - x[2])*(X - x[3]) + 
			f(x[2]) / wi[2] * (X - x[0])*(X - x[1])*(X - x[3]) + 
			f(x[3]) / wi[3] * (X - x[0])*(X - x[1])*(X - x[2]);
}

void Lagrange(vector<double> x) {
	vector<double> fi(x.size());
	vector<double> wi(x.size());

	for (int i = 0; i < x.size(); ++i) {
		fi[i] = f(x[i]);
		double count = 1;
		for (int j = 0; j < x.size(); ++j) {
			if (i != j) {
				count *= x[i] - x[j];
			}
		}
		wi[i] = count;
		count = 1;
	}

	cout << "i" << "\t" << "xi" << "\t" << "fi" << "\t" << "w(xi)"<< "\t" << "fi/w(xi)\n";
	for (int i = 0; i < x.size(); i++) {
		cout << i << "\t" << x[i] << "\t" << f(x[i]) << "\t" << wi[i] << "\t" << f(x[i]) / wi[i] << '\n';
	}
	cout << "\n lagrange polynomial:\n L" << x.size()-1 <<"(x) = " 
	<< f(x[0]) / wi[0] << "(x - " << x[1] << ")(x - " << x[2] << ")( x - " << x[3] << ") \n" 
	<< "\t+ " << f(x[1]) / wi[1] << "(x - " << x[0] << ")(x - " << x[2] << ")(x - " << x[3] << ") \n" 
	<<'\t' <<f(x[2]) / wi[2] << "(x - " << x[0] << ")(x - " << x[1] << ")(x - " << x[3] << ") \n" 
	<< "\t+ " << f(x[3]) / wi[3] << "(x - " << x[0] << ")(x - " << x[1] << ")(x - " << x[2] << ")\n\n";

	double xStar = 0.8;
	cout << "X* = " << xStar << '\n';
	double lg = LagPol(xStar, x, wi);
	double trueVal = f(xStar);

	cout << "L3(" << xStar << ") = " << lg << '\n';
	cout << "F(" << xStar << ") = " << trueVal << '\n';
	cout << "absolute error: " << abs(lg - trueVal) << "\n\n";
}

void Newton (vector<double> x) {
	vector<vector<double>> vals(4);
	for (int i = 0; i < x.size(); i++) {
		vals[0].push_back(f(x[i]));
	}
	for (int i = 1; i < x.size(); i++) {
		vals[1].push_back(f(x[i - 1], x[i]));
	}
	for (int i = 2; i < x.size(); i++) {
		vals[2].push_back(f(x[i - 2], x[i - 1], x[i]));
	}
	for (int i = 3; i < x.size(); i++) {
		vals[3].push_back(f(x[i - 3], x[i - 2], x[i - 1], x[i]));
	}


	for (int i = 0; i < vals.size(); i++) {
		cout << "f" << i+1 << ": ";
		for (int j = 0; j < vals[i].size(); j++) {
			cout << vals[i][j] << "  ";
		}
		cout << '\n';
	}
	cout << endl;

	cout << "\n newton polynomial:\n";
	cout << "N3(x) = " << vals[0][0] << " + (x - " << x[0] << ") * " << vals[1][0] << '\n'
	<<'\t' << " + (x - " << x[0] << ")(x - " << x[1] << ") * " << vals[2][0] << '\n'
	<<'\t' << " + (x - " << x[0] << ")(x - " << x[1] << ")(x - " << x[2] << ") * " << vals[3][0] << "\n\n";

	double xStar = 0.8;
	cout << "X* = " << xStar << '\n';
	double np = NewtPol(xStar, vals, x);
	double trueVal = f(xStar);
	cout << "N3(" << xStar << ") = " << np << '\n';
	cout << "F(" << xStar << ") = " << trueVal << '\n';
	cout << "absolute error: " << abs(np - trueVal) << "\n";


}
int main() {
	vector<double> x = ReadX();
	cout << "Lagrange method:\n";
	Lagrange(x);
	cout << "Newton method:\n";
	Newton(x);

	return 0;
}