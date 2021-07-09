#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip> 

using namespace std;

void ReadFromFile(vector<double>& x, vector<double>& f, double& xStar) {
	int c;
	double val;
	cin >> c;
	x.resize(c);
	f.resize(c);
	for (int i = 0; i < c; ++i) {
		cin >> val;
		x[i] = val;
	}
	for (int i = 0; i < c; ++i) {
		cin >> val;
		f[i] = val;
	}
	cin >> xStar;
	return; 
}


double h(int i, vector<double> x) {
	return x[i] - x[i - 1];
}


vector<double> getC(vector<double>& x, vector<double>& f) {
	vector<double> a, b, c, d;
	a.push_back(0);
	b.push_back(2 * (h(1,x) + h(2,x)));
	c.push_back(h(2,x));

	a.push_back(h(2,x));
	b.push_back(2 * (h(2,x) + h(3,x)));
	c.push_back(h(4,x));

	a.push_back(h(3,x));
	b.push_back(2 * (h(3,x) + h(4,x)));
	c.push_back(0);

	for (int i = 2; i < x.size(); ++i) {
		d.push_back(3 * ((f[i] - f[i-1]) / h(i,x) - (f[i-1] - f[i-2]) / h(i-1,x)));
	}

	vector<double> p(d.size());
    vector<double> q(d.size());

    p[0] = -c[0] / b[0];
    q[0] = d[0] / b[0];

    for (int i = 1; i < p.size(); i++) {
        p[i] = -c[i] / (b[i] + a[i] * p[i - 1]);
        q[i] = (d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]);
    }

    vector<double> answ(d.size());
    answ[answ.size() - 1] = q[q.size() - 1];

    for (int i = answ.size() - 2; i >= 0; i--) {
        answ[i] = answ[i + 1] * p[i] + q[i];
    }

    answ.push_back(0);
    for (int i = answ.size()-2; i >= 0; i--) {
    	answ[i+1] = answ[i];
    }
    answ[0] = 0;
	return answ;
}

double fi(double xStar, int i, int xi, vector<vector<double>>& coefs, vector<double>& x) {
	return  coefs[0][i] + coefs[1][i]*(xStar-x[xi]) + coefs[2][i]*pow((xStar-x[xi]),2) + 
			coefs[3][i]*pow((xStar-x[xi]),3);
}

double fi1(double xStar, int i, int xi, vector<vector<double>>& coefs, vector<double>& x) {
	return  coefs[1][i] + coefs[2][i]*2*(xStar-x[xi]) + 
			coefs[3][i]*3*pow((xStar-x[xi]),2);
}

double fi2(double xStar, int i, int xi, vector<vector<double>>& coefs, vector<double>& x) {
	return coefs[2][i] * 2 + coefs[3][i]*6*(xStar-x[xi]);
}

void Spline(vector<double>& x, vector<double>& f, double xStar) {
	vector<double> c = getC(x, f);
	vector<double> a(c.size()), b(c.size()), d(c.size());
	for (int i = 1; i < c.size(); ++i) {
		a[i-1] = f[i-1];
		b[i-1] = (f[i]-f[i-1])/ h(i, x) - 1.0/3.0 * h(i, x) * (c[i] + 2*c[i-1]);
		d[i-1] = (c[i] - c[i-1]) / (3 * h(i,x));
	}
	int n = b.size() - 1;
	b[n] = (f[n+1] - f[n])/h(n,x) - 2.0/3.0 * h(n,x)*c[n];
	d[n] = -c[n]/ (3*h(n,x));
	a[n] = f[n];

	vector<vector<double>> coefs;
	coefs.push_back(a);
	coefs.push_back(b);
	coefs.push_back(c);
	coefs.push_back(d);

	printf("x\t\ta\t\tb\t\tc\t\td\n");
	for (int i = 0; i < a.size(); i++) {
		printf("[%.2f-%.2f]\t%f\t%f\t%f\t%f\n",x[i], x[i+1], a[i], b[i], c[i], d[i]);
	}

	cout << "\n X* = " << xStar << '\n'; 
	cout << "f(X*) = " << a[1] << " + " << b[1] << "(x - " << x[1] << ")"
		<< " + " << c[1] << "(x - " << x[1] << ")^2"
		<< " + " << d[1] << "(x - " << x[1] << ")^3\n";

	double res = fi(xStar, 1, 1, coefs, x);
	cout << "f(" << xStar << ") = " << res << '\n';

//check

	printf("\nCHECK:\n\tf(x)\t\t\tf'(x)\t\t\tf''(x)\n");
	for (int i = 1; i < 4; i++) {
		double res = fi(x[i], i, i, coefs, x);
		double res1 = fi(x[i], i-1, i-1, coefs, x);
		printf("x%d: %f %f\t", i ,res, res1);
	
		res = fi1(x[i], i, i, coefs, x);
		res1 = fi1(x[i], i-1, i-1, coefs, x);
		printf("%f %f\t",res, res1);

		
		res = fi2(x[i], i, i, coefs, x);
		res1 = fi2(x[i], i-1, i-1, coefs, x);
		printf("%f %f\n",res, res1);
	}
}

int main() {

	vector<double> x, f;
	double xStar;
	ReadFromFile(x, f, xStar);
	Spline(x,f, xStar);

	return 0;
}