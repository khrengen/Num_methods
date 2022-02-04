#include <iostream>
#include <vector>
#include <cmath>

#define MAX_ITER 1000

using namespace std;

double TrueSolution(double x, double y) {
	return cos(x)*cos(y); 
}

double MSE(vector<vector<double>>& u, vector<double>& x, vector<double>& y) {
	double err = 0;
	for (int i = 0; i < u.size(); i++) {
		for (int j = 0; j < u[0].size(); j++) {
			err+= pow(u[i][j]- TrueSolution(x[i], y[j]),2);
		}
	}
	return err/(u.size()*u[0].size());
}

double Norm(vector<vector<double>>& u1, vector<vector<double>>& u2) {
	double norm = 0;
	for (int i = 1; i < u1.size()-1; i++) {
		for (int j = 1; j < u1[0].size()-1; j++) {
			if (abs(u1[i][j] - u2[i][j]) > norm) {
				norm = abs(u1[i][j] - u2[i][j]);
			}
		}
	}
	return norm;
}

double xLeft (double y) { 
	return cos(y);
}

double xRight (double y) {
	return 0;
}

double yDown(double x) {
	return cos(x);
}

double yUp(double x) {
	return 0;
}

int main() {
	int n, m;
	double x0, xn, y0, ym;
	char chm, chg;
	cout << "Enter x0 xn y0 ym n m\n";
	cin >> x0 >> xn >> y0 >> ym >> n >> m;
	//init
	double h1 = (xn-x0)/n;
	double h2 = (ym-y0)/m;

	cout <<"h1 = " << h1 << "  h2 = " << h2 << '\n';

	vector<double> x, y;
	for (int i = 0; i <= n; i++) {
		x.push_back(x0+i*h1);
	}
	for (int i = 0; i <= m; i++) {
		y.push_back(y0+i*h2);
	}

	vector<vector<double>> u(x.size(), vector<double>(y.size(), 0));

	//граничные условия первого рода
	for (int i = 0; i < x.size(); i++) {
		u[i][0] = yDown(x[i]);
		u[i][y.size()-1] = yUp(x[i]);
	}

	for (int i = 0; i < y.size(); i++) {
		u[0][i] = xLeft(y[i]);
		u[x.size()-1][i] = xRight(y[i]);
	}

	//интерполяция
	for (int i = 1; i < x.size()-1; i++) {
		for (int j = 1; j < y.size()-1; j++) {
			u[i][j] = (ym-y[j])/(ym-y0)*u[i][0] + (y[j]-y0)/(ym-y0)*u[i][y.size()-1];
		}
	}

	for (int j = 1; j < y.size()-1; j++) {
		for (int i = 1; i < x.size()-1; i++) {
			double intx = (xn-x[i])/(xn-x0)*u[0][j] + (x[i]-x0)/(xn-x0)*u[x.size()-1][j];
			u[i][j] = (u[i][j]+intx)/2;
		}
	}

	double c = 1;
	double eps = 0.0001;
	cout << "eps = "<< eps <<"\n\n";
	cout << "select method:\n\tl -  Libman\n\tz - Zeidel\n";
	cin >> chm;
	cout << "need relaxation? (y/n)\n";
	cin >> chg;
	if (chg == 'y') {
		cout << "C = ";
		cin >> c;
	}

	int iter = 0;
	switch(chm) { 	
		case 'l':	//Либман

			for (int t = 0; t < MAX_ITER; t++) {
				vector<vector<double>> uNext = u;
				for (int i = 1; i < x.size()-1; i++) {
					for (int j = 1; j < y.size()-1; j++) {
						uNext[i][j] = c*(h2*h2*u[i+1][j] + h2*h2*u[i-1][j] + h1*h1*u[i][j+1] + h1*h1*u[i][j-1] + h1*h1*h2*h2*2*u[i][j])/(2*(h1*h1+h2*h2)) + u[i][j]*(1-c);
					}
				}
				iter++;
				if (Norm(uNext,u) < eps) {
					break;
				}
				u = uNext;
			}
			break;	
		default:	//Зейдель
			for (int t = 0; t < MAX_ITER; t++) {
				vector<vector<double>> uOld = u;
				for (int i = 1; i < x.size()-1; i++) {
					for (int j = 1; j < y.size()-1; j++) {
						//u[i][j] = c*(h2*h2*u[i+1][j] + h2*h2*u[i-1][j] + h1*h1*u[i][j+1] + h1*h1*u[i][j-1] + h1*h1*h2*h2*2*u[i][j])/(2*(h1*h1+h2*h2));
						u[i][j] = c*(h2*h2*u[i+1][j] + h2*h2*u[i-1][j] + h1*h1*u[i][j+1] + h1*h1*u[i][j-1])/(2*(h1*h1+h2*h2-h1*h1*h2*h2)) + u[i][j]*(1-c);
					}
				}
				iter++;
				if (Norm(uOld,u) < eps) {
					break;
				}
			}
			break;	
	}

	
	cout << fixed;
	//cout.precision(3);
	cout << "\nnumber of itaration: " << iter << '\n';
	
	cout.precision(10);
	cout << "\n MSE = " << MSE(u,x,y) << '\n';
	return 0;
}