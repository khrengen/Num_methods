#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

void ReadFromFile(double& y0, double& yx0, double& x0, double& x1, double& h) {
	cin >> y0;
	cin >> yx0;
	cin >> x0;
	cin >> x1;
	cin >> h;
	return; 
}

double f(double x, double y, double z) {
	return (-x*z + y + 3*pow(x,2)) / pow(x, 2);
}

double RRR(double fh, double fkh, double h, double kh, double p) {
	double k = kh/h;
	return (fkh- fh)/(pow(k,p)-1);
}

vector<double> Euler(vector<double>& x, double y0, double z0, double h) {
	vector<double> y(x.size()), z(x.size());
	y[0] = y0;
	z[0] = z0;

	for (int i = 1; i < y.size(); i++) {
		z[i] = z[i-1] + h*f(x[i-1], y[i-1], z[i-1]);
		y[i] = y[i-1] + h*z[i-1];
	}
	return y;
}

vector<double> BoostedEuler (vector<double>& x, double y0, double z0, double h) {

	vector<double> y(x.size()), z(x.size());
	y[0] = y0;
	z[0] = z0;

	for (int i = 1; i < y.size(); i++) {

		double xhalh = x[i-1] + h/2;
		double yhalf = y[i-1] + h/2*z[i-1];
		double zhalf = z[i-1] + h/2*f(x[i-1], y[i-1], z[i-1]);

		z[i] = z[i-1] + h*f(xhalh, yhalf, zhalf);
		y[i] = y[i-1] + h*zhalf;
	}
	return y;

}

vector<double> RungeKutta(vector<double>& x, double y0, double z0, double h) {
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

		double theta1 = abs((k2[i] - k3[i])/(k1[i]-k2[i]));
		double theta2 = abs((l2[i] - l3[i])/(l1[i]-l2[i]));
		//cout << theta1 << " " << theta2 << '\n';
		//if (theta1 > 0.1 || theta2 > 0.1) {
		//	h /= 2.0;
		//} else if (theta1 < 0.01 || theta2 < 0.01) {
	//		h *= 2.0;
	//	}
	}

	return y;
}

vector<double> Adams(vector<double>& x, vector<double>& yStart, double h) {
	vector<double> y(x.size());
	for (int i = 0; i < 4; i++) {
		y[i] = yStart[i];
	}

	for (int i = 4; i < y.size(); i++) {
		y[i] = y[i-1] + h/24*(55*y[i-1] - 59*y[i-2] + 37*y[i-3] - 9*y[i-4]);
	}

	return y;

}


int main() {
	vector<double> x1,x2;
	double h,y0,yx0,x0,xk;
	ReadFromFile(y0, yx0, x0, xk, h);
	for (double i = x0; i <= xk+h; i += h) {
		x1.push_back(i);
	}

	cout << "x^2*y'' + xy' - y - 3x^2 = 0\n replacement: z = y'\n";
	cout << "x^2*z'  + xz  - y - 3x^2 = 0\n\n";

	vector<double> answE = Euler(x1, y0, yx0, h);
	cout << "Euler:\n";
	cout << "Y = [ ";
	for (int i = 0; i < answE.size(); i++) {
		cout << answE[i] << " ";
	}
	cout << "]\n\n";

	vector<double> answBoost = BoostedEuler(x1, y0, yx0, h);
	cout << "Boosted Euler:\n";
	cout << "Y = [ ";
	for (int i = 0; i < answBoost.size(); i++) {
		cout << answBoost[i] << " ";
	}
	cout << "]\n\n";

	vector<double> answR = RungeKutta(x1, y0, yx0, h);
	cout << "RungeKutta:\n";
	cout << "Y = [ ";
	for (int i = 0; i < answR.size(); i++) {
		cout << answR[i] << " ";
	}
	cout << "]\n\n";

	vector<double> answA = Adams(x1, answR, h);
	cout << "Adams:\n";
	cout << "Y = [ ";
	for (int i = 0; i < answA.size(); i++) {
		cout << answA[i] << " ";
	}
	cout << "]\n";	


	for (double i = x0; i <= xk+h; i+= h*2) {
		x2.push_back(i);
	}
	vector<double> answE2 =  Euler(x2, y0, yx0, 2*h);
	vector<double> answBoost2 =  Euler(x2, y0, yx0, 2*h);
	vector<double> answR2 =  RungeKutta(x2, y0, yx0, 2*h);
	vector<double> answA2 =  Adams(x2,answR, 2*h);

	cout << "\nRunge-Romberg: inaccuracy\n";
	cout <<"Euler:\t"<<RRR(answE[10], answE2[10], h, h*2, 1) << '\n';
	cout <<"Boosted Euler:\t"<<RRR(answBoost[10], answBoost2[10], h, h*2, 2) << '\n';
	cout <<"RungeKutta:\t"<<RRR(answR[10], answR2[10], h, h*2, 4) << '\n'; 
	cout <<"Adams:\t"<<RRR(answA[10], answA2[10], h, h*2, 4) << '\n'; 

	cout << "\nabsolute error:\n";
	cout << "Euler:\t" << abs(answE[10] - 6.5) / 6.5 << '\n';
	cout << "Boosted Euler:\t" << abs(answBoost[10] - 6.5) / 6.5 << '\n';
	cout << "RungeKutta:\t" << abs(answR[10] - 6.5) / 6.5 << '\n'; 
	cout << "Adams:\t" << abs(answA[10] - 6.5) / 6.5 << '\n'; 

	return 0;
}