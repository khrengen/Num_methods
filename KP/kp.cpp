#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

using namespace std;


class DisMatrix {
public:
	DisMatrix(vector<double>& d,vector<double>& e,vector<int>& j,vector<int>& i):
		diag(d), elem(e), jptr(j), iptr(i) {

	}
	DisMatrix(vector<double>& d):
		diag(d){
		iptr.resize(d.size()+1,0);
	}

	DisMatrix(){
	}

	vector<double> multV(const vector<double>& x) {
		vector<double> z(x.size());
		for (int i = 0; i < x.size(); ++i) {
			z[i] = x[i] * diag[i];
		}
		for (int i = 0; i < x.size(); ++i) {
			for (int j = iptr[i]; j < iptr[i+1]; ++j) {
				z[i] += x[jptr[j]] * elem[j];
				z[jptr[j]] += x[i] * elem[j];
			}
		}
		return z;
	}

	void PrintDiag() {
		for (int i = 0; i < diag.size(); i++) {
			cout << diag[i] << '\n';
		}
	}
private:
	vector<double> diag, elem;
	vector<int> jptr, iptr;
};

double scalMult(vector<double> v1, vector<double> v2){
	double res;
	for (int i = 0; i < v1.size(); ++i) {
		res += v1[i]*v2[i];
	}
	return res;
}

vector<double> vecMult(vector<double> v1, double a){
	for (int i = 0; i < v1.size(); ++i) {
		v1[i] *= a;
	}
	return v1;
}

vector<double> vecMult(vector<double> v1, vector<double> v2){
	for (int i = 0; i < v1.size(); ++i) {
		v1[i] *= v2[i];
	}
	return v1;
}

vector<double> vecAdd(vector<double> v1, vector<double> v2){
	vector<double> res(v1.size());
	for (int i = 0; i < v1.size(); ++i) {
		res[i] = v1[i] + v2[i];
	}
	return res;
}

vector<double> vecSub(vector<double> v1, vector<double> v2){
	vector<double> res(v1.size());
	for (int i = 0; i < v1.size(); ++i) {
		res[i] = v1[i] - v2[i];
	}
	return res;
}

vector<double> vecMinus(vector<double> v1){
	for (int i = 0; i < v1.size(); ++i) {
		v1[i] =  -v1[i];
	}
	return v1;
}

double vecNorm(vector<double>& x) {
	double res = 0;
	for (int i = 0; i < x.size(); ++i) {
		res += pow(x[i],2);
	}
	return sqrt(res);
}

void Print(vector<double> answ) {
	for (int i = 0; i < answ.size(); i++) {
		cout << answ[i] << '\n';
	}
	cout <<'\n';
}

DisMatrix ReadFromFile(vector<double> &b, DisMatrix& M1) {
	int n;
	double v;
	vector<double> elem, diag;
	vector<int> jptr,iptr(1,0);
	cin >> n;
	vector<vector<double>> fullvec(n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			cin >> v;
			fullvec[i].push_back(v);
			if (j == i) {
				diag.push_back(v);
			} else if (j < i) {
				if (v) {
					elem.push_back(v);
					jptr.push_back(j);
				}
			}
		}
		iptr.push_back(elem.size());
	}

	for (int i = 0; i < n; ++i) {
		cin >> v;
		b.push_back(v);
	}
	DisMatrix mat(diag, elem, jptr, iptr);
	for (int i = 0; i < diag.size(); i++) {
		diag[i] = 1/diag[i];
	}
	M1 = DisMatrix(diag); 
	return mat;
}

vector<double> ConjugateGradient(DisMatrix& A, vector<double>& b) {
	
	vector<double> xOld(b.size(), 0);
	vector<double> rOld = vecSub(b, A.multV(xOld));
	vector<double> zOld(rOld);
	vector<double> x,r,z;
	double end = 0.001;
	double eps = 100;
	int n = 0;
	while (eps > end) {
		double alfa = scalMult(rOld,rOld) / scalMult(A.multV(zOld), zOld);
		x = vecAdd(xOld, vecMult(zOld, alfa));
		r = vecSub(rOld, vecMult(A.multV(zOld),alfa));

		double beta = scalMult(r, r) / scalMult(rOld, rOld);
		z = vecAdd(r, vecMult(zOld, beta));
		xOld = x;
		rOld = r;
		zOld = z;
		eps = vecNorm(r);
		n++;
	}
	cout << n<<"\n";
	return x;
} 


vector<double> PrecConjugateGradient(DisMatrix& A, DisMatrix&Minv, vector<double>& b) {
	vector<double> xOld(b.size(), 0);
	vector<double> rOld = vecSub(b, A.multV(xOld));
	vector<double> zOld = Minv.multV(rOld);
	vector<double> pOld(zOld);
	vector<double> x,r,z,p;
	double end = 0.001;
	double eps = 100;
	int n = 0;
	while (eps > end) {
		double alfa = scalMult(rOld, zOld) / scalMult(A.multV(pOld),pOld);
		x = vecAdd(xOld, vecMult(pOld, alfa));
		r = vecSub(rOld, vecMult(A.multV(pOld),alfa));
		z = Minv.multV(r);
		double beta = scalMult(r, z) / scalMult(rOld, zOld);
		p = vecAdd(z, vecMult(pOld, beta));
		xOld = x;
		rOld = r;
		zOld = z;
		pOld = p;
		eps = vecNorm(r);
		n++;
	}
	cout << n<< "\n";
	return x;

}
int main () {
	vector<double> v;
	DisMatrix M1;
	DisMatrix A = ReadFromFile(v, M1);
	unsigned int start_time =  clock();
	vector<double> answ = ConjugateGradient(A, v);
	unsigned int end_time =  clock();
	cout << double(end_time - start_time) / 1000 << '\n';
	//cout << "X = \n";
	//Print(answ);
	unsigned int start_time2 =  clock();
	vector<double> answ2 = PrecConjugateGradient(A, M1, v);
	unsigned int end_time2 =  clock();
	cout << double(end_time2 - start_time2) / 1000;
	//Print(answ2);
	return 0;
}