#include <iostream>
#include <vector>
#include "TMatrix.cpp"

using namespace std;

void ReadFromFile(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& q) {
	
	double c1, c2, c3, size;
	cin >> size;
	for (int i = 0; i < size; i++) {
		if (i == 0) {
            cin >> c2 >> c3;
            c1 = 0;
        } else if (i == size - 1){
            cin >> c1 >> c2;
            c3 = 0;
        } else {
            cin >> c1 >> c2 >> c3;
        }
		a.push_back(c1);
        b.push_back(c2);
        c.push_back(c3);
	}

    for (int i = 0; i < size; i++) {
        cin >> c1;
        q.push_back(c1);
    }
}

int main() {

	vector<double> d, a, b, c;
    ReadFromFile(a, b, c, d);

    vector<double> p(d.size());
    vector<double> q(d.size());

    p[0] = -c[0] / b[0];
    q[0] = d[0] / b[0];

    for (int i = 1; i < p.size(); i++) {
        p[i] = -c[i] / (b[i] + a[i] * p[i - 1]);
        q[i] = (d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]);
    }

    cout << "run factors :\n P = ( ";
    for(int i = 0; i < p.size(); i++) {
        cout << p[i] << " ";
    }
    cout << ")\n\n Q = ( ";
    for(int i = 0; i < q.size(); i++) {
        cout << q[i] << " ";
    }
    cout << ")\n\n";

    vector<double> x(d.size());
    x[x.size() - 1] = q[q.size() - 1];

    for (int i = x.size() - 2; i >= 0; i--) {
        x[i] = x[i + 1] * p[i] + q[i];
    }

    cout << "answer\nX = ( ";
    for (int i = 0; i < x.size(); i++) {
        cout << x[i] << ' ';
    }
    cout << ")\n";

	return 0;
}


