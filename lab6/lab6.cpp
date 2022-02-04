#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double TrueSolution(double x, double t) {
	return exp(-t-x)*cos(x)*cos(2*t);
}

double MNE(vector<vector<double>>& u, vector<double>& x, vector<double>& t) {
	double err = 0;
	for (int i = 0; i < u.size(); i++) {
		for (int k = 0; k < u[0].size(); k++) {
			err+= pow(u[i][k]- TrueSolution(x[i], t[k]),2);
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

double LeftBC (double t) { //first kind
	return exp(-t)*cos(2*t);
}

double RightBC (double t) { //first kind
	return 0;
}

double FirstInCond(double x) {
	return exp(-x)*cos(x);
}

double SecInCond(double x) {
	return -exp(-x)*cos(x);
}

int main() {
	int n, K;
	double L, R, T;
	double a;
	char chm, chg;
	cout << "Enter LRnTK\n";
	cin >> L >> R >> n >> T >> K;
	//init
	double h = (R-L)/n;
	double tau = T / K;
	double sigma = pow(tau,2)/pow(h,2);
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
		u[i][0] = FirstInCond(x[i]);
	}

	//граничные условия первого рода
	for (int i = 0; i < t.size(); i++) {
		u[0][i] = LeftBC(t[i]);
		u[x.size()-1][i] = RightBC(t[i]);
	}

	cout << "select scheme:\n\te - explicit\n\ti - implicit\n";
	cin >> chm;
	cout << "select approximation of the second initial condition:\n\tfirst-order - 1\n\tsecond-order - 2\n";
	cin >> chg;

	// второй слой
	for (int i = 1; i < x.size()-1; i++) {
		if (chg == '1') {
			u[i][1] = SecInCond(x[i])*tau + u[i][0];
		} else {
			u[i][1] = u[i][0] + SecInCond(x[i])*tau*(1-tau) + (2*exp(-x[i])*sin(x[i]) - 2*(exp(-x[i])*cos(x[i])+exp(-x[i])*sin(x[i])) -3*u[i][0])*pow(tau,2)/2;
		}
	}

	switch(chm) { 	
		case 'e':	//явная схема
			for (int k = 1; k < t.size()-1; k++) {
				for (int i = 1; i < x.size()-1; i++) {
					u[i][k+1] = (u[i][k-1]*pow(h,2)*(tau-1) + u[i-1][k]*pow(tau,2)*(1-h) + u[i][k]*(2*pow(h,2)-2*pow(tau,2)-3*pow(h,2)*pow(tau,2)) + u[i+1][k]*pow(tau,2)*(1+h))/(pow(h,2)*(tau+1));
				}
			}	
			break;	
			
		default:  //неявная схема
			//для прогонки
			vector<double> A(x.size()-2,0), B(x.size()-2,0), C(x.size()-2,0), D(x.size()-2,0);
			for (int i = 0; i < x.size()-2; i++) {
				A[i] = pow(tau,2)*(h-1);
				B[i] = pow(h,2) + pow(h,2)*tau + 2*pow(tau,2)+ 3*pow(h,2)*pow(tau,2);
				C[i] = -pow(tau,2)*(h+1);
			}

			for (int k = 1; k < t.size()-1; k++) {
				for (int i = 0; i < x.size()-2; i++) {
					D[i] = 2*pow(h,2)*u[i+1][k] + pow(h,2)*(tau-1)*u[i+1][k-1];
				}
				D[0] -= pow(tau,2)*(h-1)*LeftBC(t[k+1]);
				D[x.size()-3] += pow(tau,2)*(h+1)*RightBC(t[k+1]);
				vector<double> answ = TomasRun(A,B,C,D);
				for(int i = 1; i < x.size()-1; i++) {
					u[i][k+1] = answ[i-1];
				}
			}
			break;
	}

	cout << "\n[x , t]\t\tresult\ttrueValue\n";
	cout << fixed;
	cout.precision(3);
	for (int i = 0; i < min(x.size(), t.size()); i+=5) {
		cout << '['<<x[i]<< "," <<t[i]<<"]\t" << u[i][i]  << '\t' << TrueSolution(x[i], t[i]) << '\n'; 
	}
	cout << '['<<x[0]<< "," <<t[2]<<"]\t" << u[0][2]  << '\t' << TrueSolution(x[0], t[2]) << '\n';
	cout << '['<<x[x.size()-1]<< "," <<t[2]<<"]\t" << u[x.size()-1][2]  << '\t' << TrueSolution(x[x.size()-1], t[2]) << '\n';
	cout.precision(10);
	cout << "\n MNE = " << MNE(u,x,t) << '\n';
	return 0;
}