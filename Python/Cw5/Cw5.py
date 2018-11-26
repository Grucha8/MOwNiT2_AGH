import numpy as np
import matplotlib.pyplot as plt

from ThirdDegreeSpline import third_degree_spline
from SecondDegreeSpline import second_degree_spline


def f(x):
    # return x**14 + x**3 + 12 * x - 1
    return np.exp(4 * np.cos(5 * x))
    # return np.sin(x * 2) * np.sin((2 * x ** 2) / np.pi)


# =========================================
def list_of_coords(a, b, n):
    points = []
    d = (b - a) / (n - 1)
    x = a
    for _ in range(n):
        points.append((x, f(x)))
        x += d

    return points


# =========================================
def main():
    start, end = -np.pi, np.pi
    n = 9
    M = 1000

    x = lambda n: np.linspace(start, end, n)
    LX = x(M)
    p = list_of_coords(start, end, n)
    # LY = second_degree_spline(p, LX, 0)
    LY = third_degree_spline(p, LX)
    plot(LX, LY, M, n, p, x)


def plot(LX, LY, M, n, points, x):
    for i in range(n):
        plt.axvline(x=points[i][0], color='b')

    plt.ion()
    plt.plot(x(M), f(x(M)), 'k')

    plt.plot(LX, LY, 'b')

    plt.show()


if __name__ == "__main__":
    main()
