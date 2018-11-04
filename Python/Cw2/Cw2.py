import numpy as np
from decimal import *
import time

#params
K = 11
M = 3
NUMTYPE = np.float64
HIGHEST_MATRIX_SZIE = 200


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


# x - wektor wejsciowy
def metodaJakobiego(A, x, b, n, iter_num, ro, kryteriumStopu):
    x2 = np.zeros(n, dtype=NUMTYPE)

    for k in range(iter_num):
        for i in range(n):
            s = 0.
            for j in range(n):
                if i != j:
                    s += A[i][j] * x[j]
            x2[i] = (1 / A[i][i]) * (b[i] - s)
        # wylicz norme
        if k > 0:
            if kryteriumStopu == 0:
                pom = x2 - x
            else:
                pom = np.dot(A, x2) - b

            curNorm = np.linalg.norm(pom, ord=np.inf)
            if curNorm < ro:
                print(f"Iteracja nr = {k}\tNorma Maksimum: {curNorm}\t", end='')
                return x

        x = np.copy(x2)

    print(f"Brak chcianego wyniku\t Norma Maksimum: {curNorm}\t", end='')
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

    xh = np.random.choice([-1., 1.], HIGHEST_MATRIX_SZIE)
    sizes = (2, 4, 8, 20, 100, HIGHEST_MATRIX_SZIE)
    ros = (0.01, 0.001, 0.0001, 0.00001, 1e-10)
    xs = (
        np.zeros(HIGHEST_MATRIX_SZIE, dtype=NUMTYPE),              # same zera
        np.ones(HIGHEST_MATRIX_SZIE, dtype=NUMTYPE),               # same jedynki
        np.full(HIGHEST_MATRIX_SZIE, -1, dtype=NUMTYPE),           # same -1
        np.full(HIGHEST_MATRIX_SZIE, 5, dtype=NUMTYPE),            # same piatki
        2 * np.random.ranf(HIGHEST_MATRIX_SZIE) - 1                # losowe liczby zmiennoprzecinkowe od -1 do 1
    )
    print(f"Wektor x:\n {xh}\n======================================")

    print(f"Losowy wektor wejsciowy x: {xs[4]}")

    for n in sizes:
        print(f"===================================================="
              f"\nSize: {n}")

        A = matrixA(n)
        x = np.resize(xh, n)

        print(f"promien spektralny: {spectral_radius(A)}")

        b = np.dot(A, x)

        for ro in ros:
            print(f"\n----------------\nRO = {ro}")

            for xd, i in zip(xs, range(5)):

                print("\nWektor poczatkowy: ", end='')
                if i == 0: print("Same zera")
                elif i == 1: print("Same 1")
                elif i == 2: print("Same -1")
                elif i == 3: print("Same 5")
                else: print("Losowe")

                x_p = np.resize(xd, n)
                print("Pierwsze kryterium stopu")
                start = time.perf_counter()
                print(np.linalg.norm(metodaJakobiego(A, x_p, b, n, 25, ro, 0), np.inf))
                end = time.perf_counter()
                print(f"\tCzas = {end-start}")

                x_p = np.resize(xd, n)
                print("Drugie kryterium stopu")
                start = time.perf_counter()
                print(np.linalg.norm(metodaJakobiego(A, x_p, b, n, 25, ro, 1), np.inf))
                end = time.perf_counter()
                print(f"\tCzas = {end-start}\n")


if __name__ == '__main__':
    main()