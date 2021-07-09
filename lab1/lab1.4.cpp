#include <iostream>
#include <vector>
#include "TMatrix.cpp"

using namespace std;

void ReadFromFile(vector<vector<double>>& vec, double& eps) {

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
	cin >> eps;
	return;
}

pair<int, int> MaxInd(TMatrix &a) {
	double value = a[0][1];
	int l = 0, k = 1;
	for (int i = 0; i < a.Size(); i++) {
		for (int j = i + 1; j < a.Size(); j++) {
			if (abs(a[i][j]) > value) {
				value = abs(a[i][j]);
				l = i;
				k = j;				
			}
		}
	}
	return make_pair(l, k);
}

TMatrix FindU(pair<int, int> ind, TMatrix& a) {
	TMatrix U(a.Size());
	for (int i = 0; i < a.Size(); i++) {
		U[i][i] = 1;
	}

	int i = ind.first;
	int j = ind.second;

	double phi = 0;
	if (a[i][i] != a[j][j]) {
		phi = 0.5 * atan( 2 * a[i][j] / (a[i][i] - a[j][j]));
	} else {
		phi = M_PI / 4;
	}

	U[i][i] = cos(phi);
	U[i][j] = -sin(phi);
	U[j][i] = sin(phi);
	U[j][j] = cos(phi);
	return U;
}

double SqSum(TMatrix& a) {
    double norm = 0;
    for (int i = 0; i < a.Size(); i++) {
        for (int j = i + 1; j < a.Size(); j++) {
           	norm += pow(a[i][j], 2);
        }
    }
    return sqrt(norm);
}

TMatrix RotationMeth(TMatrix &a, double eps, vector<double>& v, int& iter) {
	double epsk = eps + 1;
	TMatrix A = a;

	TMatrix Ucount(a.Size());
	for (int i = 0; i < Ucount.Size(); i++) {
		Ucount[i][i] = 1;
	}

	while (epsk > eps) {
		pair<int, int> indxs = MaxInd(A);
		TMatrix U = FindU(indxs, A);
		Ucount = Ucount * U;
		A = U.Transp() * A * U;
		epsk = SqSum(A);
		iter++;
	}
	v.resize(A.Size());
	for (int i = 0; i < A.Size(); i++) {
		v[i] = A[i][i];
	}
	return Ucount;
}

int main() {

	double eps;
	int iter = 0;
	vector<vector<double>> vec;
	ReadFromFile(vec, eps);
	TMatrix A(vec);
	vector<double> eigenvalues;
	TMatrix res = RotationMeth(A, eps, eigenvalues, iter);

	cout << "eigenvalues:\n\n";
	for (int i = 0; i < eigenvalues.size(); i++) {
		cout << "l" << i + 1 << " = " << eigenvalues[i] << '\n';
	}
	cout << "\neigenvectors:\n\n";
	cout << "x1\tx2\tx3\n";
	Print(res);

	cout << "number of iterations: " << iter << " (eps = " << eps << ")\n\n";

	cout << "Check:\t\t\tAx\t==\tlx\n";
	vector<vector<double>> v(res.Size());
	for (int i = 0; i < res.Size(); i++) {
		v[0].push_back(res[i][0]);
		v[1].push_back(res[i][1]);
		v[2].push_back(res[i][2]);
	}
	
	for (int k = 0; k < res.Size(); k ++) {
		vector<double> pr = A * v[k];
		cout << "( ";
		for (int i = 0; i < res.Size(); i++) {
			cout << pr[i] << ' ';
		}
		cout << " )  (";
		for (int i = 0; i < res.Size(); i++) {
			cout << v[k][i] * eigenvalues[k] << ' ';
		
		}
		cout << " )\n";
	}
	return 0;
}
