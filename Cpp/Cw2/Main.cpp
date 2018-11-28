#include "../Cw1/importy.h"
#include "../Cw1/Normy.cpp"

using namespace std;


template<typename T, typename T1, typename T2>	//np: T = float		T1 = float*		T2 = float**
void metodaJakobiego(T2 A, T1 b, T1 x, int n, int iter, T _) {
	T1 x2 = new T[n];

    // x == poprzednia iteracja (x^k) | x2 = obecna iteracja (x^k+1)
    for (int k = 0; k < iter; ++k) {
        for (int i = 0; i < n; ++i) {
            T sum = 0;
            for (int j = 0; j < n; ++j) {
                if (i != j)
                    sum += A[i][j] * x[j];
            }

            x2[i] = (1 / A[i][i]) * (b[i] - sum);
        }
        for (int l = 0; l < n; ++l) {
            x[l] = x2[l];
        }
    }
}

//===================================================
template<typename T2, typename T>
void matrixA(T2 A, int n, T m, T k){
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j)
                A[i][i] = k;
            else
                A[i][j] = 1 / (abs(i - j) + m);
        }
    }
}

void macierzX(int* X, int n) {
    int pom[] = { 1, -1 };

    for (int i = 0; i < n; i++)
        X[i] = pom[rand() % 2];
}
template<typename T, typename T1, typename T2>
void matrixMul(T A, int* X, T1 B, int n, T2 _) {
    T2 sum;

    for (int i = 0; i < n; i++) {
        sum = 0;
        for (int j = 0; j < n; j++) {
            sum += (A[i][j] * X[j]);
        }
        B[i] = sum;
    }
}
//===================================================

int main() {
    freopen(".\\Cw2\\output.txt","w",stdout);
    cout << fixed;

    int n = 100;
    double** A;
    double *b;
    double *x_po;
    int *x;
    double pom = 0, k = 5, m = 9, ro = 0.001;

    A = new double*[n];
    for (int i = 0; i < n; ++i) {
        A[i] = new double[n];
    }
    b = new double[n];

    x_po = new double[n];
    x = new int[n];
    macierzX(x, n);

    matrixA(A, n, m, k);

    matrixMul(A, x, b, n, pom);
    cout << "Prawidlowy euklides: " << euklidesNorm(x, n, pom) << endl << endl;

    for (int i = 0; i < 10; ++i) {
        for (int i = 0; i < n; ++i) {
            x[i] = 0;
        }
        metodaJakobiego(A, b, x_po, n, i, pom);

        cout << "Ile iteracji: " << i << endl;
        cout << "Euklides: " << euklidesNorm(x_po, n, pom) << endl;
        cout << "Max: " << maxNorm(x_po, n, pom) << endl << "------------------------------" << endl;
    }


}