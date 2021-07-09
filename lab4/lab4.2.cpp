#include <iostream>
#include <vector>
#include <cmath>
#include "..\lab2\gnuplot.h"

using namespace std;

void ReadFromFile(double& z0, double& kr2, double& x0, double& x1, double& h) {
	cin >> z0;
	cin >> kr2;
	cin >> x0;
	cin >> x1;
	cin >> h;
	return; 
}

double NewNu (double nu0, double nu1, double F0, double F1) {
	return nu1 - (nu1-nu0)/(F1-F0)*F1;
}

double f(double x, double y, double z) {
	return (4*(pow(x,2)+3)*z - 6*x*y) / x / (pow(x,2)+6);
}

double RRR(double fh, double fkh, double h, double kh, double p) {
	double k = kh/h;
	return (fkh- fh)/(pow(k,p)-1);
}

pair<vector<double>,vector<double>> RungeKutta(vector<double>& x, double y0, double z0, double h) {
	vector<double> y(x.size()), z(x.size());
	y[0] = y0; z[0] = z0;
	vector<double> k1(x.size()), k2(x.size()), k3(x.size()), k4(x.size());
	vector<double> l1(x.size()), l2(x.size()), l3(x.size()), l4(x.size());

	for (int i = 0; i < x.size(); i++) {
		if (i > 0) {

			z[i] = z[i-1] + (l1[i-1] + 2*l2[i-1] + 2*l3[i-1] + l4[i-1])/6.0;
			y[i] = y[i-1] + (k1[i-1] + 2*k2[i-1] + 2*k3[i-1] + k4[i-1])/6.0;
		}
		l1[i] = h*f(x[i], y[i], z[i]);
		k1[i] = h*z[i];

		l2[i] = h*f(x[i] + h/2.0, y[i] + k1[i]/2.0, z[i] + l1[i]/2.0);
		k2[i] = h*(z[i] + l1[i]/2.0);

		l3[i] = h*f(x[i] + h/2.0, y[i] + k2[i]/2.0, z[i] + l2[i]/2.0);
		k3[i] = h*(z[i] + l2[i]/2.0);

		l4[i] = h*f(x[i] + h, y[i] + k3[i], z[i] + l3[i]);
		k4[i] = h*(z[i] + l3[i]);

	}
	return make_pair(y,z);
}

double Fi(double nu, vector<double> y, vector<double> z, double kr) {
	return abs(y[y.size()-1] - z[z.size()-1] - kr);
}

double Shooting(vector<double> x, double z0, double h, double kr) {
	double oldNu = 0;
	double nNu = 5;
	double nu;
	double res = 100.0;
	while(res > 0.0001) {
		pair<vector<double>,vector<double>> yz1 = RungeKutta(x, oldNu, z0, h);
		pair<vector<double>,vector<double>> yz2  = RungeKutta(x, nNu, z0, h);
		double F1 = Fi(oldNu, yz1.first, yz1.second, kr);
		double F2 = Fi(nNu, yz2.first, yz2.second, kr);
		nu = NewNu(oldNu, nNu, F1, F2);
		pair<vector<double>,vector<double>> yz3  = RungeKutta(x, nu, z0, h);
		res = Fi(nu, yz3.first, yz3.second, kr);
		oldNu = nNu;
		nNu=nu;

	}
	return nu;
}

double q(double x) {
	return 6/(pow(x,2)+6);
}
double p(double x) {
	return -4*(pow(x,2)+3)/x/(pow(x,2)+6);
}

double f(double x) {
	return 0.0;
}

vector<double> FiniteDifference(vector<double> x, double z0, double h, double kr) {

	vector<double> a(x.size()), b(x.size()), c(x.size()), d(x.size());
	a[0] = 0;
	b[0] = - 1 + pow(h,2)*p(x[0]) - p(x[0])*h/2;
	c[0] = 1 + p(x[0])*h/2;
	d[0] = h*z0*(1-p(x[0])*h/2);
	for (int k = 1; k < x.size()-1; ++k) {
		a[k] = (1-p(x[k])*h/2);
		b[k] = (-2 + pow(h,2)*q(x[k]));
		c[k] = (1 + p(x[k])*h/2);
		d[k] = f(x[k])*pow(h,2);
	}
	a[x.size()-1] = 1;
	b[x.size()-1] = h-1;
	c[x.size()-1] = 0;
	d[x.size()-1] = kr*h;


	vector<double> p(d.size());
    vector<double> q(d.size());

    p[0] = -c[0] / b[0];
    q[0] = d[0] / b[0];

    for (int i = 1; i < p.size(); i++) {
        p[i] = -c[i] / (b[i] + a[i] * p[i - 1]);
        q[i] = (d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]);
    }

    vector<double> y(d.size());
    y[x.size() - 1] = q[q.size() - 1];

    for (int i = y.size() - 2; i >= 0; i--) {
        y[i] = y[i + 1] * p[i] + q[i];
    }
    return y;
}

int main () {
	vector<double> x1,x2;
	double h,z0,kr2,x0,xk;
	ReadFromFile(z0, kr2, x0, xk, h);
	for (double i = x0; i <= xk+h; i += h) {
		x1.push_back(i);
	}
	double y0 = Shooting(x1, z0, h, kr2);
	vector<double> answS1 = RungeKutta(x1, y0, z0, h).first;
	vector<double> answFD1 = FiniteDifference(x1, z0, h, kr2);
	
	for (double i = x0; i <= xk+h; i += h*2) {
		x2.push_back(i);
	}

	vector<double> answS2 = RungeKutta(x2, y0, z0, h*2).first;
	vector<double> answFD2 = FiniteDifference(x2, z0, h*2, kr2);

	cout << "\nRunge-Romberg: inaccuracy\n";
	cout << "Shooting method: " << RRR(answS1[answS1.size()-1], answS2[answS2.size()-1], h, h*2, 4) << '\n';
	cout << "Finite Difference method: " << RRR(answFD1[answS1.size()-1], answFD2[answS2.size()-1], h, h*2, 1) << '\n';
	cout << "\nabsolute error:\n";
	cout << "Shooting method: " << abs(answS1[answS1.size()-1] - 82.0) / 82.0 << '\n';
	cout << "Finite Difference method: " << abs(answFD1[answS1.size()-1] - 82.0) / 82.0;

	Gnuplot plot;
	plot("set xrange [0:+4]");
    plot("plot x**3 + x**2 +2, 'RK1.dat'  with lines, 'FD1.dat'  with lines");
    //plot("plot x**3 + x**2 +2, 'RK2.dat'  with lines, 'FD2.dat'  with lines");

	return 0;

}