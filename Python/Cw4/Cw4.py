import numpy as np
import matplotlib.pyplot as plt
import os

from LagrangeInterpolation import lagrange_interpolation
from NewtonInterpolation import newton_interpolation
from HermiteInterpolation import hermite_interpolation


def f(x):
    return np.sin(x*2) * np.sin((2 * x**2) / np.pi)


def df(x):
    a = 4 * x * np.sin(2 * x) * np.cos(2 * x ** 2 / np.pi) / np.pi
    b = 2 * np.sin(2 * x ** 2 / np.pi) * np.cos(2 * x)
    return a + b


# ==============================================
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


# ==============================================
def find_max_diff(LX, LY, M):
    max = -np.inf
    x = 0
    for i in range(M):
        d = np.abs(LY[i] - f(LX[i]))
        if d > max:
            max = d
            x = LX[i]
    return (x, max)


# ==============================================
def main():
    start, end = -np.pi, np.pi
    M = 1000
    x = lambda n: np.linspace(start, end, n)
    LX = x(M)

    plt.plot(x(M), f(x(M)), 'k')
    plt.grid()
    plt.xlabel("x")
    plt.ylabel('y')
    plt.savefig("./plots/original")
    plt.close()

    lagr = [[], []]     # 0 - chebyshev, 1 - evenly
    ns = [2, 3, 4, 5, 7, 10, 13, 15, 20, 30, 35, 40, 50]
    print("\n=====================\nLAGRANGEA\n=====================")
    for n in ns:
        points1 = list_of_chebyshev_coords(start, end, n)
        LY1 = lagrange_interpolation(points1, LX)
        p = find_max_diff(LX, LY1, M)
        print(n, " - Lagrangea-Chebyshev:\t\t", p)
        lagr[0].append(p)

        points2 = list_of_coords(start, end, n)
        LY2 = lagrange_interpolation(points2, LX)
        p = find_max_diff(LX, LY2, M)
        print(n, " - Lagrangea-Evenly:\t\t", p)
        lagr[1].append(p)

        plot(LX, LY1, LY2, M, n, points1, points2, x, 'L')

    newt = [[], []]
    print("\n=====================\nNEWTON\n=====================")
    for n in ns:
        points1 = list_of_chebyshev_coords(start, end, n)
        LY1 = newton_interpolation(points1, LX)
        p = find_max_diff(LX, LY1, M)
        print(n, " - Newton-Chebyshev:\t\t", p)
        newt[0].append(p)

        points2 = list_of_coords(start, end, n)
        LY2 = newton_interpolation(points2, LX)
        p = find_max_diff(LX, LY2, M)
        print(n, " - Newton-Evenly:\t\t", p)
        newt[1].append(p)

        plot(LX, LY1, LY2, M, n, points1, points2, x, 'N')
    return lagr, newt


def hermit():
    start, end = -np.pi, np.pi
    M = 1000
    x = lambda n: np.linspace(start, end, n)
    LX = x(M)

    ns = [2, 3, 4, 5, 7, 10, 13, 15, 20, 30, 50]
    herm = [[], []]
    print("\n=====================\nHERMIT\n=====================")
    for n in ns:
        points1 = list_of_chebyshev_coords(start, end, n)
        points1 = [(point[0], point[1], df(point[0])) for point in points1]
        LY1 = hermite_interpolation(points1, LX)
        p = find_max_diff(LX, LY1, M)
        print(n, " - Hermite-Chebyshev:\t\t", p)
        herm[0].append(p)

        points2 = list_of_coords(start, end, n)
        points2 = [(point[0], point[1], df(point[0])) for point in points2]
        LY2 = hermite_interpolation(points2, LX)
        p = find_max_diff(LX, LY2, M)
        print(n, " - Hermite-Evenly:\t\t", p)
        herm[1].append(p)

        plot(LX, LY1, LY2, M, n, points1, points2, x, 'H')

    return herm


# ============================================================
def plot(LX, LY1, LY2, M, n, points1, points2, x, inter):
    # black plot == actual function
    plt.ion()
    plt.plot(x(M), f(x(M)), 'k')

    dx = [d[0] for d in points1]
    dy = [d[1] for d in points1]
    # red plot == interpolated function chebyshev
    l1, _ = plt.plot(LX, LY1, 'r', dx, dy, 'ro')

    # blue plot == interpolated function evenly
    dx = [d[0] for d in points2]
    dy = [d[1] for d in points2]

    l2, _ = plt.plot(LX, LY2, 'b', dx, dy, 'bo')

    plt.legend((l1, l2), ('Zera Czebyszewa', 'Równoodlegle'))
    plt.grid()
    plt.xlabel("x")
    plt.ylabel('y')
    title = 'Interpolacja '
    if inter == 'L':
        title += "Lagrange\'a"
    elif inter == 'N':
        title += 'Newton\'a'
    elif inter == 'H':
        title += 'Hermit\'a'
    plt.title(title + ', Liczba węzłów = ' + str(n))

    name = "plot" + str(n) + inter + ".png"
    plt.savefig("./plots/" + name)
    plt.show()
    plt.close()


if __name__ == "__main__":
    if not os.path.exists('./plots'):
        os.mkdir('./plots')
    main()
    hermit()

