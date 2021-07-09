#include <iostream>
#include <vector>
#include "TMatrix.cpp"

using namespace std;

const int ITER_MAX = 100;

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

int sign(double a) {
	return a >= 0 ? 1 : -1;
}
/*
vector<double> VectorMult(const vector<double> &a, const vector<double> &b) {
	vector<double> res(a.size());
	for (int i = 0; i < res.size(); ++i) {
		res[i] = a[i] * b[i];
	}
	return res;
}

vector<double> VectorDiv(const vector<double> &a, const vector<double> &b) {
	vector<double> res(a.size());
	for (int i = 0; i < res.size(); ++i) {
		if (b[i] != 0) {
			res[i] = a[i] / b[i];
		} else {
			cerr << "0 element in division!\n";
		}
	}
	return res;
}
*/
double VecNorm(const vector<double>& a) {
	double norm = 0;
	for (int i = 0; i < a.size(); i++) {
		norm += pow(a[i], 2);
	}
	return sqrt(norm);
}

bool CheckEnd (TMatrix& A, double eps) {
	double sq1 = 0, lastsq1 = 10000, lastsq2 = 0,
	 sq2 = 0;
	for (int i = 0; i < A.Size() - 1; i++) {
		sq1 = 0;
		sq2 = 0;
		for (int j = i + 1; j < A.Size(); ++j) {
			sq1 += pow(A[j][i], 2);
			if (j < A.Size() - 1) {
				sq2 += pow(A[j + 1][i], 2);
			}
		}
		if (sqrt(sq1) > eps && (sqrt(sq2) > eps || sqrt(lastsq1) > eps)) {
			return false;
		}
		lastsq1 = sq1;
		lastsq2 = sq2;
	}
	return true;
}



vector<pair<double, double>> solveQv(double a, double b , double c) {
	vector<pair<double, double>> answ(2);
	double d = pow(b, 2) - 4 * a * c;
	if (d < 0) {
		answ[0].first = -b/(2*a);
		answ[1].first = -b/(2*a);
		answ[0].second = sqrt(abs(d))/(2*a);
		answ[1].second = -sqrt(abs(d))/(2*a);
	} else {
		cerr << "no complex solution\n";
	}
	return answ;
}


TMatrix GetH(vector<double>& v) {
	TMatrix E;
	E.GetE(v.size());
	TMatrix mat(v.size());
	for (int i = 0; i < v.size(); ++i) {
		for (int j = 0; j < v.size(); ++j) {
			mat[i][j] = v[i] * v[j];
		}
	}

	double scal = 0;
	for (int i = 0; i < v.size(); ++i) {
		scal += v[i] * v[i];
	}
	return E -  mat / scal * 2.0;
}


TMatrix HouseHolder(TMatrix& A, TMatrix& Q) {
	TMatrix H(A.Size());
	TMatrix A0 = A;
	vector<double> v(A.Size());

	for (int i = 0; i < v.size() - 1; ++i) {
		for(int k = 0; k < v.size(); ++k) {
			if (k < i) {
				v[k] = 0;
			} else if (k == i) {
				double sq = 0;
				for (int j = k; j < v.size(); ++j) {
					sq += pow(A0[j][k], 2);
				}
				sq = sqrt(sq);
				v[k] = A0[k][k]+ sign(A0[k][k]) * sq;
			} else {
				v[k] = A0[k][i];
			}
		}
		H = GetH(v);
		A0 = H * A0;
		Q = Q * H;
	}	
	return A0;
}


void QRMethod (TMatrix &A0, double& eps) {
	int iter = 0;
	TMatrix Q;
	TMatrix A = A0;
	while (iter < ITER_MAX) {
		Q.GetE(A0.Size());
		A = HouseHolder(A, Q);
		A = A * Q;
	//	Print(A);
		iter++;
		if (CheckEnd(A, eps)) {
			break;
		}

	}
	
	//Print(A);
	for (int i = 0; i < A.Size(); ++i) {
		double sq1 = 0;
		for (int j = i + 1; j < A.Size(); ++j) {
			sq1 += pow(A[j][i], 2);
		}
		if (sqrt(sq1) <= eps) {
			cout << "l" << i + 1 << " = " << A[i][i] << '\n';
		} else {
			vector<pair<double, double>> ans = solveQv(1.0, -A[i][i]*A[i+1][i+1], -A[i][i+1]*A[i+1][i]);
			cout << "l" << i+1 << " = " << ans[0].first << " + " << ans[0].second << "i\n";
			cout << "l" << i+2 << " = " << ans[1].first << " - " << abs(ans[1].second) << "i\n";
			++i;
		}
	}
	return;
}

int main() {

	double eps;
	vector<vector<double>> vec;
	ReadFromFile(vec, eps);
	TMatrix A(vec);
	cout << "input matrix: \n";
	Print(A);
	cout << "example of decomposition:\n";
	TMatrix Q;
	Q.GetE(A.Size());
	TMatrix R = HouseHolder(A, Q);
	cout << " Q = \n";
	Print(Q);
	cout << "R =\n";
	Print(R);
	cout << "Q * R = \n";
	R = Q * R;
	Print(R);
	cout << "eigenvalues by QR (eps = " << eps << "):\n";
	QRMethod(A, eps);
}