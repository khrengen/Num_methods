#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double TrueSolution(double x, double t, double a) {
	return exp(-t*a)*sin(x);
}

double MNE(vector<vector<double>>& u, vector<double>& x, vector<double>& t, double a) {
	double err = 0;
	for (int i = 0; i < u.size(); i++) {
		for (int k = 0; k < u[0].size(); k++) {
			err+= pow(u[i][k]- TrueSolution(x[i], t[k], a),2);
		}
	}
	return err/(u.size()*u[0].size());
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

    vector<double> x(d.size());
    x[x.size() - 1] = q[q.size() - 1];

    for (int i = x.size() - 2; i >= 0; i--) {
        x[i] = x[i + 1] * p[i] + q[i];
    }

    return x;
}

double LeftBC (double t, double a) {
	return exp(-a*t);
}

double RightBC (double t, double a) {
	return -exp(-a*t);
}

void RenVectors(vector<double>& A, vector<double>& B, vector<double>& C, double h, double tau, double a, char chg) {
	if (chg == '1') {
		double sigma = a*tau/pow(h,2);
		B[0] = -1/h;
		C[0] = 1/h;
		A[A.size()-1] = -1/h;
		B[B.size()-1] = 1/h;
		for (int i = 1; i < A.size()-1; i++) {
			A[i] = sigma;
			C[i] = sigma;
			B[i] = -1 - 2*sigma;
		}
	} else {
		B[0] = 2*a/h + h/tau;
		C[0] = -2*a/h;
		A[A.size()-1] = -2*a/h;
		B[B.size()-1] = 2*a/h + h/tau;
		for (int i = 1; i < A.size()-1; i++) {
			A[i] = -a/pow(h,2);
			B[i] = 2*a/pow(h,2)+ 1/tau;
			C[i] = -a/pow(h,2);
		}
	}
}


int main() {
	int n, K;
	double L, R, T;
	double a;
	char chm, chg;
	cout << "Enter LRnTKa\n";
	cin >> L >> R >> n >> T >> K >> a;
	//init
	double h = (R-L)/n;
	double tau = T / K;
	double sigma = a*tau/pow(h,2);
	cout << "sigma = " << sigma << '\n';
	cout <<"h = " << h << "  tau = " << tau << '\n';
	vector<double> x, t;
	for (int i = 0; i <= n; i++) {
		x.push_back(L+i*h);
	}
	for (int i = 0; i <= K; i++) {
		t.push_back(i*tau);
	}

	vector<vector<double>> u(x.size(), vector<double>(t.size(), 0));

	//ну
	for (int i = 0; i < x.size(); i++) {
		u[i][0] = sin(x[i]);
	}

	//для прогонки
	vector<double> A(x.size(),0), B(x.size(),0), C(x.size(),0), D(x.size(),0);

	cout << "select scheme:\n\te - Explicit\n\ti - implicit\n\tc - Crank-Nicholson\n";
	cin >> chm;
	cout << "select approximation of the boundary conditions:\n\tfirst-order - 1\n\tsecond-order - 2\n";
	cin >> chg;

	switch(chm) { 	//явная схема
		case 'e':
			for (int q = 0; q < t.size()-1; q++) {
				for (int i = 1; i < x.size()-1; i++) {
					u[i][q+1] = sigma*u[i+1][q] + (1-2*sigma)*u[i][q] + sigma*u[i-1][q];
				}
				if (chg == '1'){
					u[0][q+1] = u[1][q+1] - h*LeftBC(t[q+1],a);
					u[x.size()-1][q+1] = u[x.size()-2][q+1] + h*RightBC(t[q+1],a);
				} else {
					u[0][q+1] = (-2*h*LeftBC(t[q+1],a) +4*u[1][q+1]-u[2][q+1])/3;
					u[x.size()-1][q+1] = (2*h*RightBC(t[q+1],a) +4*u[x.size()-2][q+1]-u[x.size()-3][q+1])/3;	
				}
			}
			break;

		case 'i':	//неявная схема
			RenVectors(A,B,C,h,tau,a,chg);
			for (int q = 0; q < t.size()-1; q++) {
				if (chg == '1') {
					D[0] = -LeftBC(t[q+1],a)*h;
					D[x.size()-1] = RightBC(t[q+1],a)*h;
					for (int i = 1; i < x.size()-1; i++) {
						D[i] = -u[i][q];
					}	
				} else {
					D[0] = h/tau*u[0][q]- LeftBC(t[q+1],a)*2*a;
					D[x.size()-1] = h/tau*u[x.size()-1][q]+RightBC(t[q+1],a)*2*a;
					for (int i = 1; i < x.size()-1; i++) {
						D[i] = u[i][q]/tau;
					}
				}

				vector<double> answ = TomasRun(A,B,C,D);
				for (int i = 0; i < x.size(); i++) {
					u[i][q+1] = answ[i];
				}
			}
			break;

		default:	//Кранка-Николсона
			RenVectors(A,B,C,h,tau,a/2,chg);

			for (int q = 0; q < t.size()-1; q++) {
				if (chg == '1'){
					D[0] = -LeftBC(t[q+1],a)*h ;
					D[x.size()-1] = RightBC(t[q+1],a)*h;
					for (int i = 1; i < x.size()-1; i++) {
						D[i] = -u[i][q] - sigma/2*(u[i+1][q]-2*u[i][q]+u[i-1][q]);
					}	
				} else {
					
					D[0] = h/tau*u[0][q]- LeftBC(t[q+1],a);
					D[x.size()-1] = h/tau*u[x.size()-1][q]+RightBC(t[q+1],a);
					for (int i = 1; i < x.size()-1; i++) {
						D[i] = (u[i][q]+ 0.5*sigma*(u[i+1][q]-2*u[i][q]+u[i-1][q]))/tau;
					}
				}

				vector<double> ImU = TomasRun(A,B,C,D);

				for (int i = 0; i < x.size(); i++) {
					u[i][q+1] = ImU[i];
				}

			}
			break;
	}

	
	cout << fixed;
	cout.precision(10);
	cout << "\n MNE = " << MNE(u,x,t,a) << '\n';
	return 0;
}