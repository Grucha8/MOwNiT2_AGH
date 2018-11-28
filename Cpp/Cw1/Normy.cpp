#include "importy.h"

template<typename T, typename T1>
T1 euklidesNorm(T X, int n, T1 _) {
    T1 sum = 0;

    for (int i = 0; i < n; i++){
        sum += (pow(X[i], 2));
    }

    return (sqrt(sum));
}

template<typename T, typename T1>
T1 maxNorm(T x, int n, T1 _) {
    T1 max = 0;

    for (int i = 0; i < n; i++){
        if (fabs(x[i]) > max) max = fabs(x[i]);
    }

    return max;
}

template<typename T, typename T1>
double uwarunkowanie(T Pom, T A, int n, T1 _) {
    double maxA = 0;
    double maxPom = 0;
    double sumA, sumPom;

    for (int i = 0; i < n; i++){
        sumA = 0;
        sumPom = 0;
        for (int j = 0; j < n; j++){
            sumA += abs(A[i][j]);
            sumPom += abs(Pom[i][j]);
        }
        if (sumA > maxA) maxA = sumA;
        if (sumPom > maxPom) maxPom = sumA;
    }

    return maxA * maxPom;
}