#include "importy.h"

// wypelnienie macierzy A
template<typename T>
void matrixA1(T A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == 0)
                A[i][j] = 1;
            else
                A[i][j] = ((1.0 / (i + j + 2.0 - 1.0)));
        }
    }
}

template<typename T>
void matrixA2(T A, int n) {
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (j >= i)
                A[i][j] = ((2 * i + 1.0) / (j + 1.0));
            else
                A[i][j] = A[j][i];
        }
    }
}

template<typename T>
void matrixA3a(T A, int n, int k, int m) {
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j) A[i][j] = k;
            else if (j == i + 1) A[i][j] = (1.0 / (i + m));
            else if (j == i - 1) A[i][j] = (k / (i + m + 1.0));
            else A[i][j] = 0;
        }
    }
}
