import numpy as np
from decimal import *

#params
RO = 0.00001
Flag = 0        #1 = x^i+1 - x^i    0 = Ax - b
K = 11
M = 3
NUMTYPE = np.float64


# macierz A
def matrixA(n):
    global K, M

    A = np.zeros((n, n), dtype=NUMTYPE)

    for i in range(n):
         for j in range(n):
              if i == j:
                  A[i][j] = K
              else:
                  A[i][j] = 1 / (np.abs(i - j) + M)
    return A


def metodaJakobiego(A, x, b, n, iter_num):
    x2 = np.zeros(n, dtype=NUMTYPE)
    #x2 = np.random.ranf(n)
    x = np.copy(x2)

    for k in range(iter_num):
        for i in range(n):
            s = 0.
            for j in range(n):
                if i != j:
                    s += A[i][j] * x[j]
            x2[i] = (1 / A[i][i]) * (b[i] - s)
        # wylicz norme
        if k > 0:
            if Flag:
                pom = x2 - x
            else:
                pom = np.dot(A, x) - b

            curNorm = np.linalg.norm(pom)
            if curNorm < RO:
                print(f"k = {k} | Norma Euklidesowaska: {curNorm}")
                return x

        for l in range(n):
            x[l] = x2[l]

    print("Brak chcianego wyniku")
    return x


def spectral_radius(A: np.array):
    Pom = np.copy(A)
    D = np.zeros(A.shape)

    for i in range(A[0].size):
        D[i] = Pom[i][i]
        Pom[i][i] = 0


    w, _ = np.linalg.eig(np.dot(-(1 / D), Pom))
    return np.linalg.norm(w, np.inf)


def main():
    np.set_printoptions(suppress=True)
    getcontext().prec = 15

    sizes = (2, 4, 8, 20, 100)
    ros = (0.01, 0.001, 0.0001, 0.00001)

    for n in sizes:
        print(f"=========================================="
              f"\nSize: {n}")

        A = matrixA(n)
        x = np.random.choice([-1., 1.], n)

        print(f"promien spektralny: {spectral_radius(A)}")

        b = np.dot(A, x)
        for RO in ros:
            print(f"\nRO = {RO}")

            FLAG = 1
            print("Pierwsze kryterium stopu")
            metodaJakobiego(A, x, b, n, 25)

            FLAG = 0
            print("Drugie kryterium stopu")
            metodaJakobiego(A, x, b, n, 25)


if __name__ == '__main__':
    main()