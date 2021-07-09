#include <iostream>
#include <vector>
#include <cmath>


using namespace std;

void ReadFromFile(double& x0, double& xk, double& h1, double& h2) {
	cin >> x0;
	cin >> xk;
	cin >> h1;
	cin >> h2;
	return; 
}

double f(double x) {
	return 1 / sqrt((2*x+7)*(3*x+4));
}

double rect(vector<double>& x, double h1) {
	double res = 0;
	for (int i = 0; i < x.size() - 1; i++) {
		res += f((x[i] + x[i+1]) / 2);
	}
	return res*h1;
}

double trape(vector<double>& x, double h1) {
	double res = 0;
	for (int i = 0; i < x.size(); i++) {
		if (i == 0 || i == x.size() - 1) {
			res += f(x[i]) / 2;
		} else {
			res += f(x[i]);
		}
	}
	return res*h1;
}

double simpson(vector<double>& x, double h1) {
	double res = 0;
	for (int i = 0; i < x.size(); i++) {
		if (i == 0 || i == x.size() - 1) {
			res += f(x[i]);
		} else if (i % 2){
			res += f(x[i])*4;
		} else {
			res += f(x[i])*2;
		}
	}
	return res*h1/3;
}

double RRR(double fh, double fkh, double h, double kh, double p) {
	double k = kh/h;
	return (fkh- fh)/(pow(k,p)-1);
}

int main() {
	vector<double> x1,x2;
	double h1,h2,x0,xk;
	ReadFromFile(x0, xk, h1, h2);
	for (double i = x0; i <= xk; i+= h1) {
		x1.push_back(i);
	}
	for (double i = x0; i <= xk; i+= h2) {
		x2.push_back(i);
	}
	cout << "step h1 = " << h1 << '\n';
	cout << "rect " << rect(x1, h1) << '\n';
	cout << "trape " << trape(x1, h1) << '\n';
	cout << "simpson " << simpson(x1, h1) << '\n';

	cout << "\nstep h2 = " << h2 << '\n';
	cout << "rect " << rect(x2, h2) << '\n';
	cout << "trape " << trape(x2, h2) << '\n';
	cout << "simpson " << simpson(x2, h2) << '\n';

	cout << "\nRunge-Romberg method: inaccuracy\n";
	cout << "\t    h2\n";
	cout <<"rect\t"<<RRR(rect(x2,h2), rect(x1,h1), h2, h1, 2) << '\n'; 
	cout <<"trape\t"<< RRR(trape(x2,h2), trape(x1,h1), h2, h1, 2) << '\n';
	cout <<"simpson\t"<< RRR(simpson(x2,h2), simpson(x1,h1), h2, h1, 4) << '\n';
	return 0;
}