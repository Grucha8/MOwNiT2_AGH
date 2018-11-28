import numpy as np
import matplotlib.pyplot as plt

from ThirdDegreeSpline import third_degree_spline
from SecondDegreeSpline import second_degree_spline


def f(x):
    # return x**14 + x**3 + 12 * x - 1
    # return np.exp(4 * np.cos(5 * x))
    return np.sin(x * 2) * np.sin((2 * x ** 2) / np.pi)


def df(x):
    a = 4 * x * np.sin(2*x) * np.cos(2*x**2/np.pi) / np.pi
    b = 2 * np.sin(2*x**2/np.pi) * np.cos(2*x)
    return a + b


# =========================================
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


# =========================================
def main():
    start, end = -np.pi, np.pi
    n = 2
    M = 1000

    ns = [4, 5, 7, 10, 13, 15, 20, 30]
    for n in ns:
        x = lambda n: np.linspace(start, end, n)
        LX = x(M)
       # p = list_of_chebyshev_coords(start, end, n)
        p = list_of_coords(start, end, n)
        #LY = second_degree_spline(p, LX, 5)
        LY = third_degree_spline(p, LX)
        plot(LX, LY, M, n, p, x, 3, 'natural')


def plot(LX, LY, M, n, points, x, degree, type):
    for i in range(n):
        plt.axvline(x=points[i][0], color='b')

    plt.ion()
    plt.plot(x(M), f(x(M)), 'k')

    plt.plot(LX, LY, 'b')

    name = "plot" + str(n) + "_" + str(degree) + "_" + type
    plt.savefig("./plots/" + name)
    #plt.show()


if __name__ == "__main__":
    main()
