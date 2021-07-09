#include <iostream>
#include <vector>
#include "TMatrix.cpp"

using namespace std;

void ReadFromFile(vector<vector<double>>& vec, vector<double>& b, double& eps) {

	double c;
	int n;
	cin >> n;
	vec.resize(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cin >> c;
			vec[i].push_back(c);
		}		
	}
	for (int i = 0; i < n; i++) {
		cin >> c;
		b.push_back(c);
	}
	cin >> eps;
	return;
}

void Jacobi(TMatrix& a, vector<double>& b, TMatrix& alf, vector<double>& bet) {

	for (int i = 0; i < bet.size(); i++) {
		bet[i] = b[i] / a[i][i];
	}

	for (int i = 0; i < bet.size(); i++) {
		for (int j = 0; j < bet.size(); j++) {
			if (i == j) {
				alf[i][j] = 0;
			} else {
				alf[i][j] = -a[i][j] / a[i][i]; 
			}
		}
	}
	return;
}

double VecNorm(vector<double>& a) {
	double norm = 0;
	for (int i = 0; i < a.size(); i++) {
		norm += pow(a[i], 2);
	}
	return sqrt(norm);
}

vector<double> SimpleItSolution(TMatrix& alf, vector<double>& bet, double eps, int& iter)  {

	vector<double> x = bet;
	double epsK = eps + 1;
	double matNorm = alf.SqNorm();

	while (epsK > eps) {
		vector<double> x1 = alf * x;
		for (int i = 0; i < x1.size(); i++) {
			x1[i] += bet[i];
		}

		vector<double> dif = x1;
		for (int i = 0; i < x.size(); i++) {
			dif[i] -= x[i];
		}
		epsK = VecNorm(dif) * matNorm / (1 - matNorm);
		x = x1;
		iter++;
	}
	return x;
}

TMatrix FindUpTriangle (TMatrix& alf) {
	TMatrix C(alf.Size());
	for (int i = 0; i < C.Size(); i++) {
		for (int j = 0; j < C.Size(); j++) {
			if (i > j) {
				C[i][j] = 0;
			} else {
				C[i][j] = alf[i][j];
			}
		}
	}
	return C;
}

vector<double> ZeidelSolution(TMatrix& alf, vector<double>& bet, double eps, int& iter) {

	vector<double> x = bet;
	double epsK = eps + 1;
	double matNorm = alf.SqNorm();
	TMatrix C = FindUpTriangle(alf);
	double cNorm = C.SqNorm();
	
	while (epsK > eps) {
		vector<double> x1(x.size());

		for (int i = 0; i < x1.size(); i++) {
			for (int j = 0; j < x1.size(); j++) {
				if (i > j) {
					x1[i] += alf[i][j] * x1[j];
				} else {
					x1[i] += alf[i][j] * x[j]; 
				}
			}
			x1[i] += bet[i];
		}

		vector<double> dif = x1;
		for (int i = 0; i < x.size(); i++) {
			dif[i] -= x[i];
		}

		epsK = VecNorm(dif) * cNorm / (1 - matNorm); 
		x = x1;
		iter++;

	}

	return x;
}

int main() {

	double eps;
	int iter = 0;
	vector<vector<double>> vec;
	vector<double> b;
	ReadFromFile(vec, b, eps);
	TMatrix a(vec);
	TMatrix alf (a.Size());
	
	vector<double> bet(b.size());
	cout << "Rewriting the system to equivalent form x = beta + alfa * x :\nalfa = \n";
	Jacobi(a, b, alf, bet);
	Print(alf);
	cout << "beta = ( ";
	for (int i = 0; i < bet.size(); i++) {
		cout << bet[i] << ' ';
	}
	cout << ")\n\n";


	vector<double> answer = SimpleItSolution(alf, bet, eps, iter);
	cout << "Simple iterations (eps = " << eps << ") :\n\n";

	cout << "X = ( ";
	for (int i = 0; i < answer.size(); i++) {
		cout << answer[i] << ' ';
	}
	cout << ")\n number of iterations = " << iter << "\n";
	iter = 0;
	cout << "-----------------------------\n";

	answer = ZeidelSolution(alf, bet, eps, iter);
	cout << "Zeidel method (eps = " << eps << ") :\n\n";

	cout << "X = ( ";
	for (int i = 0; i < answer.size(); i++) {
		cout << answer[i] << ' ';
	}
	cout << ")\nnumber of iterations = " << iter << "\n";

	return 0;
}


