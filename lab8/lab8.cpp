#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double TrueSolution(double x, double y, double t) {
	return x*y*cos(t);
}

double MSE(vector<vector<vector<double>>>& u, vector<double>& x, vector<double>& y, vector<double>& t) {
	double err = 0;
	for (int i = 1; i < u.size()-1; i++) {
		for (int j = 1; j < u[0].size()-1; j++) {
			for (int p = 1; p < u[0][0].size(); p++) {
				err+= pow(u[i][j][p] - TrueSolution(x[i], y[j], t[p]),2);
			}
		}
	}
	return err/((u.size()-2)*(u[0].size()-2)*(u[0][0].size()-1));
}

void Copy_level(vector<vector<double>> &u_half, vector<vector<vector<double>>> &u, int k) {
	for (int i = 0; i < u_half.size(); i++) {
		for (int j = 0; j < u_half[0].size(); j++) {
			u_half[i][j] = u[i][j][k];
		}
	}
	return;
}

vector<double> TomasRun (vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d) {
	vector<double> p(d.size());
    vector<double> q(d.size());

    p[0] = -c[0] / b[0];
    q[0] = d[0] / b[0];

    for (int i = 1; i < p.size(); i++) {
        p[i] = -c[i] / (b[i] + a[i] * p[i - 1]);
        q[i] = (d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]);
    }
    p[p.size()-1] = 0;

    vector<double> x(d.size());
    x[x.size() - 1] = q[q.size() - 1];

    for (int i = x.size() - 2; i >= 0; i--) {
        x[i] = x[i + 1] * p[i] + q[i];
    }

    return x;
}

double X0face (double y, double t) { //first kind
	return 0;
}

double X1face (double y, double t) { //first kind
	return y*cos(t);
}

double Y0face (double x, double t) { //first kind
	return 0;
}

double Y1face (double x, double t) { //first kind
	return x*cos(t);
}

double Z0face (double x, double y) { //first kind
	return x*y;
}

int main() {
	int a, b, c;
	double x0, xn, y0, ym, tp;
	char chm;
	cout << "Enter x0 x1 y0 y1 t1 nx ny nz\n";
	cin >> x0 >> xn >> y0 >> ym >> tp >> a >> b >> c;
	//init
	double hx = (xn-x0)/a;
	double hy = (ym-y0)/b;
	double ht = (tp)/c;
	cout << "hx = " << hx << "\thy = " << hy << "\tht = " << ht << '\n';

	vector<double> x, y, t;
	for (int i = 0; i <= a; i++) {
		x.push_back(x0+i*hx);
	}
	for (int i = 0; i <= b; i++) {
		y.push_back(y0+i*hy);
	}
	for (int i = 0; i <= c; i++) {
		t.push_back(i*ht);
	}

	vector<vector<vector<double>>> u(x.size(), vector<vector<double>>(y.size(), vector<double>(t.size(), 0)));
	
	//ну
	for (int i = 0; i < x.size(); i++) {
		for (int j = 0; j < y.size(); j++) {
			u[i][j][0] = Z0face(x[i], y[j]);
		}
	}

	//граничные условия первого рода
	for (int j = 0; j < y.size(); j++) {
		for (int p = 0; p < t.size(); p++) {
			u[0][j][p] = X0face(y[j], t[p]);
			u[x.size()-1][j][p] = X1face(y[j], t[p]);
		}
	}

	for (int i = 0; i < x.size(); i++) {
		for (int p = 0; p < t.size(); p++) {
			u[i][0][p] = Y0face(x[i], t[p]);
			u[i][y.size()-1][p] = Y1face(x[i], t[p]);
		}
	}

	cout << "select scheme:\n\tv - variable directions\n\tf - fractional steps\n";
	cin >> chm;
	double cf = min(hx, min(hy,ht));

	switch(chm) { 	
		case 'v':	//переменные направления
		cout << "v\n";
			for (int k = 0; k < t.size()-1; k++) {

				vector<vector<double>> u_half(x.size(), vector<double>(y.size(), 0));
				Copy_level(u_half, u, k);
				for (int j = 1; j < y.size()-1; j++) {
					vector<double> A(x.size(),0), B(x.size(),0), C(x.size(),0), D(x.size(),0);

					for (int i = 0; i < x.size(); i++) {
						A[i] = 0.5*ht/pow(hx,2);
						B[i] = -ht/pow(hx,2) - 1;
						C[i] = 0.5*ht/pow(hx,2);
						D[i] = -u[i][j][k] - cf*0.5*ht/pow(hy,2)*(u[i][j+1][k] - 2.0*u[i][j][k] + u[i][j-1][k]) + 0.5*ht*x[i]*y[j]*sin(t[k]+ht/2);
					}
					A[0] = 0;
					C[x.size()-1] = 0;
					D[0]-= 0.5*ht/pow(hx,2) * u[0][j][k];
					D[x.size()-1]-= 0.5*ht/pow(hx,2) * u[x.size()-1][j][k];
					vector<double> solve = TomasRun(A,B,C,D);
					for (int q = 1; q < x.size()-1; q++) {
						u_half[q][j] = solve[q];
					}
				}

				for (int i = 1; i < x.size()-1; i++) {
					vector<double> A(y.size(),0), B(y.size(),0), C(y.size(),0), D(y.size(),0);
					for (int j = 0; j < y.size(); j++) {
						A[j] = 0.5*ht/pow(hy,2);
						B[j] = -ht/pow(hy,2) - 1;
						C[j] = 0.5*ht/pow(hy,2);
						D[j] = -u_half[i][j] - cf*0.5*ht/pow(hx,2)*(u_half[i+1][j] -2*u_half[i][j] + u_half[i-1][j]) + 0.5*ht*x[i]*y[j]*sin(t[k]+ht/2);
					}
					A[0] = 0;
					C[y.size()-1] = 0;
					D[0]-= 0.5*ht/pow(hy,2) * u_half[i][0];
					D[y.size()-1]-= 0.5*ht/pow(hy,2) * u_half[i][y.size()-1];
					vector<double> solve = TomasRun(A,B,C,D);
					for (int q = 1; q < y.size()-1; q++) {
						u[i][q][k+1] = solve[q];
					}	
				}
			}
			break;	
			
		default:  //дробные шаги
			for (int k = 0; k < t.size()-1; k++) {
				vector<vector<double>> u_half(x.size(), vector<double>(y.size(), 0));
				Copy_level(u_half, u, k);
				for (int j = 1; j < y.size()-1; j++) {
					vector<double> A(x.size(),0), B(x.size(),0), C(x.size(),0), D(x.size(),0);
					for (int i = 0; i < x.size(); i++) {
						A[i] = ht;
						B[i] = -2*ht-pow(hx,2);
						C[i] = ht;
						D[i] = -u[i][j][k]*pow(hx,2)+ht*pow(hx,2)/2*x[i]*y[j]*sin(t[k]);
					}
					A[0] = 0;
					C[x.size()-1] = 0;
					D[0] -= X0face(y[j], t[k]+ht/2)*ht;
					D[x.size()-1] -= X1face(y[j], t[k]+ht/2)*ht;
					vector<double> solve = TomasRun(A,B,C,D);
					for (int q = 1; q < x.size()-1; q++) {
						u_half[q][j] = solve[q];
					}
				}

				for (int i = 1; i < x.size()-1; i++) {
					vector<double> A(y.size(),0), B(y.size(),0), C(y.size(),0), D(y.size(),0);
					for (int j = 0; j < y.size(); j++) {
						A[j] = ht;
						B[j] = -2*ht-pow(hy,2);
						C[j] = ht;
						D[j] = -pow(hy,2)*u_half[i][j] + ht*pow(hy,2)/2*x[i]*y[j]*sin(t[k]+ht);
					}
					A[0] = 0;
					C[y.size()-1] = 0;
					D[0] -= Y0face(x[i], t[k+1])*ht;
					D[x.size()-1] -= Y1face(x[i], t[k+1])*ht;
					vector<double> solve = TomasRun(A,B,C,D);
					for (int q = 1; q < y.size()-1; q++) {
						u[i][q][k+1] = solve[q];
					}	
				}
			}
			break;
	}

	
	cout << fixed;
	cout.precision(3);
	//for (int i = 0; i < x.size(); i++) { 
	//	for (int j = 0; j < y.size(); j++) {
	//		cout << u[i][j][t.size()-1] << ' ';
	//		cout << TrueSolution(x[i], y[j], t[t.size()-1]) << '\n';
	//	}
//		cout << '\n';
	//}

	cout.precision(10);
	cout << "\n MNE = " << MSE(u,x,y,t) << '\n';
	return 0;
}