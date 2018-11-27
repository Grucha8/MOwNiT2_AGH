import numpy as np
import matplotlib.pyplot as plt

from LagrangeInterpolation import lagrange_interpolation
from NewtonInterpolation import newton_interpolation
from HermiteInterpolation import hermite_interpolation


def f(x):
    return np.sin(x*2) * np.sin((2 * x**2) / np.pi)


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
def main():
    start, end = -np.pi, np.pi
    M = 1000

    ns = [2, 3, 4, 5, 7, 10, 13, 15, 20, 30]
    for n in ns:
        print("=====================\nLiczba wezlow: " + str(n))
        x = lambda n: np.linspace(start, end, n)
        LX = x(M)

        points1 = list_of_chebyshev_coords(start, end, n)
        print(points1)
        LY1 = lagrange_interpolation(points1, LX)

        points2 = list_of_coords(start, end, n)
        print(points2)
        LY2 = lagrange_interpolation(points2, LX)

        plot(LX, LY1, LY2, M, n, points1, points2, x, 'L')

    print("\n=====================\nNEWTON\n=====================")
    for n in ns:
        print("===============\nLiczba wezlow: " + str(n))
        x = lambda n: np.linspace(start, end, n)
        LX = x(M)

        points1 = list_of_chebyshev_coords(start, end, n)
        print(points1)
        LY1 = newton_interpolation(points1, LX)

        points2 = list_of_coords(start, end, n)
        print(points2)
        LY2 = newton_interpolation(points2, LX)

        plot(LX, LY1, LY2, M, n, points1, points2, x, 'N')


def xd():
    start, end = -np.pi, np.pi
    M = 1000
    x = lambda n: np.linspace(start, end, n)
    LX = x(M)
    n = 7

    points1 = list_of_chebyshev_coords(start, end, n)
    print(points1)
    LY1 = hermite_interpolation(points1, LX)

    points2 = list_of_coords(start, end, n)
    print(points2)
    LY2 = hermite_interpolation(points2, LX)

    plot(LX, LY1, LY2, M, n, points1, points2, x, 'H')


# ============================================================
def plot(LX, LY1, LY2, M, n, points1, points2, x, inter):
    # putting vertical lines which indicates nodes for noramal <blue>
    for i in range(n):
        plt.axvline(x=points2[i][0], color='b')

    # chebushev nodes <red>
    for i in range(n):
        plt.axvline(x=points1[i][0], color='r')

    # black plot == actual function
    plt.ion()
    plt.plot(x(M), f(x(M)), 'k')

    # red plot == interpolated function chebyshev
    plt.plot(LX, LY1, 'r')

    # _ plot == interpolated function
    plt.plot(LX, LY2, 'b')

    name = "plot" + str(n) + inter + ".png"
    plt.savefig("./plots/" + name)
    plt.show()


if __name__ == "__main__":
   # main()
    xd()