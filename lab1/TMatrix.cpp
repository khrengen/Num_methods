#include <iostream>
#include <cmath>
#include <vector>


class TMatrix {

public: 
    TMatrix(): n(0), m(0) {
    }

    TMatrix(int n1, int m1) {
        n = n1;
        m = m1;
        stor.resize(n);
        for (int i = 0; i < n; i++) {
            stor[i].resize(m);
        }
    }

    TMatrix(int n1) {
        n = n1;
        m = n1;
        stor.resize(n);
        for (int i = 0; i < n; i++) {
            stor[i].resize(n);
        }
    }

    TMatrix(std::vector<std::vector<double>> &vec) {
        n = vec.size();
        m = vec[0].size();
        stor = vec; 
    }

    TMatrix(int n1, int m1, std::vector<double> &vec) {
        n = n1;
        m = m1;

        stor.resize(n);
        for (int i = 0; i < n; i++) {
            stor[i].resize(m);
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                stor[i][j] = vec[m*i + j];
            }
        }   
    }

    TMatrix& operator= (const TMatrix& b) {
        n = b.n;
        m = b.m;
        stor = b.stor;
        return *this;
    }


    TMatrix operator+ (TMatrix& b) {
        if (IsSameSize(*this, b)) {
            TMatrix c(n, m);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    c.stor[i][j] = stor[i][j] + b.stor[i][j];
                }
            }
            return c;
        }

    }

    TMatrix operator- (const TMatrix& b) {
        if (IsSameSize(*this, b)) {
            TMatrix c(n, m);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    c.stor[i][j] = stor[i][j] - b.stor[i][j];
                }
            }
            return c;
        }

    }

    TMatrix operator* (double c) {
        TMatrix mtrx(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                mtrx.stor[i][j] = stor[i][j] * c;
            }
        }
        return mtrx;
    }

    TMatrix operator* (TMatrix& b) {
        if (IsCompatible(*this, b)) {
            TMatrix c(n, b.m);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < b.m; j++) {
                    c.stor[i][j] = 0;
                    for (int k = 0; k < m; k++) {
                        c.stor[i][j] += stor[i][k] * b.stor[k][j];
                    }
                }
            }
            return c;
        }

    }

    std::vector<double> operator* (std::vector<double> &v) {
        
        std::vector<double> res(n);
        for (int i = 0; i < res.size(); i++) {
            for (int j = 0; j < v.size(); j++) {
                res[i] += stor[i][j] * v[j];
            }
        }
        return res;

    }

    TMatrix operator/ (double c) {
        TMatrix mtrx(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                mtrx.stor[i][j] = stor[i][j] / c;
            }
        }
        return mtrx;
    }

    TMatrix& SwapRow(int l, int k) {
        for (int i = 0; i < m; i++) {
            double tmp = stor[l][i];
            stor[l][i] = stor[k][i];
            stor[k][i] = tmp;
        }
        return *this;
    }

        TMatrix& SwapColumn(int l, int k) {
        for (int i = 0; i < n; i++) {
            double tmp = stor[i][l];
            stor[i][l] = stor[i][k];
            stor[i][k] = tmp;
        }
        return *this;
    }


    friend void LU(TMatrix &a, TMatrix &L, TMatrix &U) {

        if (a.n != a.m) {
            std::cerr << "matrix isn`t square\n";
            return;
        }
        U = a;
        for (int k = 0; k < a.n; k++) {

            if (U.CheckSwaps(k)) {
                int l = U.swaps.back().first;
                int m = U.swaps.back().second;
                U.SwapRow(l, m);
                L.SwapRow(l,m);
                L.SwapColumn(l,m);
            }
            
            for (int i = k; i < a.n; i++) {
                L.stor[i][k] = U.stor[i][k] / U.stor[k][k];
            } //Print(L);
        
            for (int i = k + 1; i < a.n; i++) {
                for (int j = k; j < a.n; j++) {
                    U.stor[i][j] -= L.stor[i][k] * U.stor[k][j];
                }
            } //Print(U);
        }
        return;
    }

    friend void Print(TMatrix& c) {
        for (int i = 0; i < c.n; i++) {
            for (int j = 0; j < c.m; j++) {
                 printf ("%.3f\t", c.stor[i][j]);
            }
            std::cout << '\n';
        }
        std::cout << '\n';
        return;
    }

    TMatrix& GetE(int n1) {
        n = n1;
        m = n1;
        stor.resize(n);
        for (int i = 0; i < n; i++) {
            stor[i].resize(n);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    stor[i][j] = 1;
                } else {
                    stor[i][j] = 0;
                }
            }
        }
        return *this;
    }

    TMatrix Transp() {
        TMatrix c(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                c.stor[j][i] = stor[i][j];
            }
        }
        return c;
    }

    double SqNorm() {
        double norm = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                norm += pow(stor[i][j], 2);
            }
        }
        return sqrt(norm);
    }

    double InfNorm() {
        double norm = 0;
        double max = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                norm += abs(stor[i][j]);
            }
            if (norm > max) {
                max = norm;
            }
            norm = 0;
        }
        return max;
    }

    size_t Size() {
        return n;
    }

    std::vector<double>& operator[] (int i) {
        return this->stor[i];
    }

    std::vector<std::pair<int, int>> swaps;

private:

    int n, m;
    std::vector<std::vector<double>> stor;
    

    bool IsCompatible (TMatrix& a, TMatrix& b) {

        if (a.m = b.n) {
            return true;
        }
        std::cerr << "matrices arn`t compatible\n";
        return false;
    }

    bool IsSameSize (const TMatrix& a, const TMatrix& b) {

        if (a.n == b.n && a.m == b.m) {
            return true;
        }
        std::cerr << "matrices arn`t same size\n";
        return false;
    }

    bool CheckSwaps(int k) {
        double max = abs(stor[k][k]);
        int index;
        bool sw = false;
        for (int i = k + 1; i < n; i++) {
            if (abs(stor[i][k]) > max) {
                max = abs(stor[i][k]);
                index = i;
                sw = true;
            }
        }
        if (sw) {
            swaps.push_back(std::make_pair(k, index));
        }
        return sw;
    }

};


std::vector<double> solveLyb(TMatrix& L, std::vector<double>& b) {

    std::vector<double> y(b.size());
    for (int i = 0; i < b.size(); i++) {
        double sum = 0;
        for (int p = 0; p < i; p++) {
            sum += y[p] * L[i][p];
        }
        y[i] = b[i] - sum;
    }

    return y;
}

std::vector<double> solveUxy(TMatrix& U, std::vector<double>& y) {

    std::vector<double> x(y.size());
    for (int i = y.size() - 1; i >= 0; i--) {
        double sum = 0;
        for (int p = y.size() - 1; p > i; p--) {
            sum += U[i][p] * x[p];
        }
        x[i] = (y[i] - sum) / U[i][i]; 
    }

    return x;
}

std::vector<double> solveOfSystem(TMatrix &L, TMatrix &U, std::vector<double> &b) {

    std::vector<double> x(b.size()), y(b.size());
    y = solveLyb(L, b);
    x = solveUxy(U, y);
    return x;
}

TMatrix reverse(TMatrix& L, TMatrix &U) {

    TMatrix X(L.Size());

    TMatrix E;
    E.GetE(L.Size());

    for (int i = 0; i < E.Size(); i++) {
        X[i] = solveOfSystem(L, U, E[i]);
    }
    X.Transp();

    for (int i = 0; i < U.swaps.size(); i++) {
       E.SwapRow(U.swaps[i].first, U.swaps[i].second);
    }
    TMatrix res = X * E;
    return res;
}

TMatrix reverse2x2(TMatrix A) {

   double det = A[0][0]*A[1][1]-A[0][1]*A[1][0];
   double tmp = A[0][0];
   A[0][0] = A[1][1] / det;
   A[1][1] = tmp / det;
   A[0][1] = -A[0][1] / det;
   A[1][0] = -A[1][0] / det;
   return A;
}
