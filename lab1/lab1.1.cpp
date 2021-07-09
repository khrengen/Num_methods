#include <iostream>
#include <vector>
#include "TMatrix.cpp"

using namespace std;



double Udet(TMatrix &U) {
    double det = 1;
    for (int i = 0; i < U.Size(); i++) {
        det *= U[i][i];
    }
    return pow(-1, U.swaps.size()) * det;
}


void ReadFromFile(vector<vector<double>>& vec, vector<double>& b) {

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
    return;
}

int main() {

    vector<vector<double>> vec;
    vector<double> b;
    ReadFromFile(vec, b);

    TMatrix A(vec);
    cout << "input matrix:\n";
    Print(A);
    TMatrix U(A.Size()), L(A.Size());
    LU(A, L, U);
    
    cout << "LU-decomposition:\n L =\n";
    Print(L);
    cout <<" U =\n";
    Print(U);

    cout << "L * U =\n";
    TMatrix res = L * U;
    Print(res);

    for (int i = 0; i < U.swaps.size(); i++) {
        double tmp = b[U.swaps[i].first];
        b[U.swaps[i].first] = b[U.swaps[i].second];
        b[U.swaps[i].second] = tmp;
    }
   vector<double> res2 = solveOfSystem(L, U, b);

    cout << "system solution:\n\n x = (  ";
    for (int i = 0; i < res2.size(); i++) {
        cout  << res2[i] << "  ";
    }
    cout << ")\n\n";


    double det = Udet(U);
    cout << "determinant of matrix:\ndet = " << det << "\n\n";

    TMatrix Arev = reverse(L,U);
    cout << "reverse matrix:\n Arev= \n";
    Print(Arev);

   //TMatrix res3 = Arev * A;
    //cout << "Arev * A\n";
    //Print(res3);
    
    
    return 0;
}