import os

import numpy as np
import matplotlib.pyplot as plt


def f(x):
    return np.sin(x*2) * np.sin((2 * x**2) / np.pi)


# =====================================================
def list_of_coords(a, b, n):
    points = []
    d = (b - a) / (n - 1)
    x = a
    for _ in range(n):
        points.append((x, f(x)))
        x += d

    return points


def list_of_chebyshev_coords(a, b, k):
    points = []
    for j in range(1, k + 1):
        x = np.cos(((2 * j - 1) * np.pi) / float(2 * k)) * (b-a)/2 + (a+b)/2
        points.append((x, f(x)))
    return points


def list_of_random_points(a, b, n):
    return (np.random.ranf(n) * b) - a


# =====================================================
def least_squares_approximation(points, m):
    n = len(points)

    x = [d[0] for d in points]
    y = [d[1] for d in points]

    def s(k):
        sum = 0
        for i in range(n):
            sum += x[i] ** k
        return sum

    def t(k):
        sum = 0
        for i in range(n):
            sum += y[i] * x[i] ** k
        return sum

    A = np.zeros((m + 1, m + 1))
    B = np.zeros(m + 1)
    for i in range(m + 1):
        for j in range(m + 1):
            A[i][j] = s(i+j)
        B[i] = t(i)

    return np.linalg.solve(A, B)


def least_squares_approximation_t(points, m, M, LX):
    n = len(points)

    dx = [d[0] for d in points]
    dy = [d[1] for d in points]

    def s(x):
        sum = 0
        for k in range(1, m):
            sum += b(k) * np.sin(k * x) + a(k) * np.cos(k * x)
        return sum + a(0) / 2

    def b(k):
        sum = 0
        for i in range(n):
            sum += dy[i] * np.sin(k * dx[i])
        return sum * 2 / (n - 1)

    def a(k):
        sum = 0
        for i in range(n):
            sum += dy[i] * np.cos(k * dx[i])
        return sum * 2 / (n - 1)

    LY = np.zeros(M)
    for i, x in zip(range(M), LX):
        LY[i] = s(x)
    return LY


# ======================================
def get_points(a, LX):
    LY = np.zeros(len(LX))

    for i in range(len(LX)):
        tmp = 0
        for j in range(len(a)):
            tmp += a[j] * LX[i] ** j
        LY[i] = tmp
    return LY


# ======================================
def main1():
    start, end = -np.pi, np.pi
    M = 1000
    ns = [10, 20, 30, 50]
    ms = [1, 2, 3, 4, 5, 6, 7, 8, 10, 15]

    for n in ns:
        for m in ms:
            x = lambda n: np.linspace(start, end, n)
            LX = x(M)

            points = list_of_coords(start, end, n)

            a = least_squares_approximation(points, m)
            LY = get_points(a, LX)

            a = np.linalg.norm(LY - f(x(M)))
            print(f'EVELNY:\t\tn = {n}; m = {m}; norm = {a}')
            plot(LX, LY, M, n, points, x, m, "Evenly")

            points = list_of_chebyshev_coords(start, end, n)

            a = least_squares_approximation(points, m)
            LY = get_points(a, LX)

            a = np.linalg.norm(LY - f(x(M)))
            print(f'CHEBYSHEV:\tn = {n}; m = {m}; norm = {a}')
            plot(LX, LY, M, n, points, x, m, "Chebyshev")
        print()


def main2():
    start, end = -np.pi, np.pi
    M = 1000
    ns = [10, 20, 30, 50]
    ms = [1, 2, 3, 4, 5, 6, 7, 8, 10, 15]

    for n in ns:
        for m in ms:

            x = lambda n: np.linspace(start, end, n)
            LX = x(M)

            points = list_of_coords(start, end, n)
            LY = least_squares_approximation_t(points, m, M, LX)
            a = np.linalg.norm(LY - f(x(M)))
            print(f'Trygonometria:\t\tn = {n}; m = {m}; norm = {a}')

            plot(LX, LY, M, n, points, x, m, "Trygonometryczne")
        print()


def plot(LX, LY, M, n, points, x, m, type_):
    plt.ion()
    plt.plot(x(M), f(x(M)), 'k')

    dx = [d[0] for d in points]
    dy = [d[1] for d in points]
    l1, _ = plt.plot(LX, LY, 'r', dx, dy, 'ro')
    plt.grid()
    plt.title('Liczba węzłów = ' + str(n) + "; Stopień wielomiany = " + str(m) + " " + type_)

    name = "plot" + str(n) + "_" + str(m) + "_" + type_
    if type_ == 'Trygonometryczne':
        tmp = 'm_' + str(m) + '/'
        if not os.path.exists('./plotsT/' + tmp):
            os.mkdir('./plotsT/' + tmp)
        plt.savefig("./plotsT/" + tmp + name)
        tmp = 'n_' + str(n) + '/'
        if not os.path.exists('./plotsT/' + tmp):
            os.mkdir('./plotsT/' + tmp)
        plt.savefig("./plotsT/" + tmp + name)

    else:
        tmp = 'm_' + str(m)
        if type_ == 'Evenly':
            tmp += './evenly/'
        elif type_ == 'Chebyshev':
            tmp += './chebyshev/'
        plt.savefig("./plots/" + tmp + name)

    #plt.show()
    plt.close()


if __name__ == '__main__':
    if not os.path.exists('./plots'):
        os.mkdir('./plots')
    main1()
    print('\n==========================\n========================\n')
    main2()
