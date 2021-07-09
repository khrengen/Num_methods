#include <iostream>
#include <vector>
#include <cmath>


using namespace std;

void ReadFromFile(vector<double>& x, vector<double>& y, double& xStar) {
	int c;
	double val;
	cin >> c;
	x.resize(c);
	y.resize(c);
	for (int i = 0; i < c; ++i) {
		cin >> val;
		x[i] = val;
	}
	for (int i = 0; i < c; ++i) {
		cin >> val;
		y[i] = val;
	}

	cin >> val;
	xStar = val;
	return; 
}

int findInd (const vector<double> x, double a) {
	for (int i = 0; i < x.size(); i++) {
		if (a <= x[i]){
			return i;
		}
	}
	cerr << "x out of range\n";
}

double derivativeA1 (const vector<double>& x, const vector<double>&y, double xStar) {
	int i = findInd(x, xStar);
	if (i == 0) {
		i++;
	}
	return (y[i] - y[i-1]) / (x[i] - x[i-1]);
}

double derivativeA1Right (const vector<double>& x, const vector<double>&y, double xStar) {
	int i = findInd(x, xStar);
	if (i == x.size()-1) {
		i--;
	}
	return (y[i+1] - y[i]) / (x[i+1] - x[i]);
}


double derivativeA2 (const vector<double>& x, const vector<double>&y, double xStar) {
	int i = findInd(x, xStar);
	if (i == 0) {
		i++;
	}

	if (i == x.size()){
		cerr << "x out of range\n";
	}


	return (y[i] - y[i-1]) / (x[i] - x[i-1]) + 
		((y[i+1] - y[i]) / (x[i+1] - x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1])) / 
			(x[i+1]-x[i-1]) * (2*xStar - x[i-1] - x[i]);
}

double derivative2 (const vector<double>& x, const vector<double>&y, double xStar) {
	int i = findInd(x, xStar);
	if (i == 0) {
		i++;
	}

	if (i == x.size()){
		cerr << "x out of range\n";
	}


	return ((y[i+1] - y[i]) / (x[i+1] - x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1])) / 
			(x[i+1]-x[i-1]) * 2;
}


int main(){

	vector<double> x, y;
	double xStar;
	ReadFromFile(x, y, xStar);

	cout << "X* = " << xStar << '\n';
	cout << "first-order precision first derivative: " << derivativeA1(x, y, xStar) << "\n\n";
	cout << "second-order precision first derivative: " << derivativeA2(x, y, xStar) << "\n\n";
	cout << "second derivative: " << derivative2(x, y, xStar) << "\n\n";

	return 0;
}