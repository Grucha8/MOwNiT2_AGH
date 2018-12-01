import numpy as np
import matplotlib.pyplot as plt
import os

from ThirdDegreeSpline import third_degree_spline
from SecondDegreeSpline import second_degree_spline


def f(x):
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


# =========================================
def solo_second_and_third():
    start, end = -np.pi, np.pi
    M = 1000

    ns = [4, 5, 7, 10, 13, 15, 20, 30]
    for n in ns:
        x = lambda n: np.linspace(start, end, n)
        LX = x(M)
        p = list_of_coords(start, end, n)

        LY = third_degree_spline(p, LX)
        print(n, " - thirdN: ", find_max_diff(LX, LY, M))
        plot(LX, LY, M, n, p, x, 'Natural', 3)
        plt.close()

        LY = third_degree_spline(p, LX, df(LX[0]), df(LX[-1]))
        print(n, " - thirdC: ", find_max_diff(LX, LY, M))
        plot(LX, LY, M, n, p, x, 'Clamped', 3)
        plt.close()

        LY = second_degree_spline(p, LX, True)
        print(n, " - secondN: ", find_max_diff(LX, LY, M))
        plot(LX, LY, M, n, p, x, 'Natural', 2)
        plt.close()

        LY = second_degree_spline(p, LX, False)
        print(n, " - secondNon: ", find_max_diff(LX, LY, M))
        plot(LX, LY, M, n, p, x, 'NoNNatural', 2)
        plt.close()
        print()


def mixed_second_and_third():
    start, end = -np.pi, np.pi
    M = 1000

    ns = [4, 5, 7, 10, 13, 15, 20, 30]
    for n in ns:
        x = lambda n: np.linspace(start, end, n)
        LX = x(M)
        p = list_of_coords(start, end, n)

        LY = third_degree_spline(p, LX)
        LY2 = second_degree_spline(p, LX, True)
        plot(LX, LY, M, n, p, x, 'Natural', LY2=LY2)

        plt.close()

        LY = third_degree_spline(p, LX, df(LX[0]), df(LX[-1]))
        LY2 = second_degree_spline(p, LX, False)
        plot(LX, LY, M, n, p, x, 'NoNNatural', LY2=LY2)

        plt.close()


# =====================================
def find_max_diff(LX, LY, M):
    max = -np.inf
    for i in range(M):
        d = np.abs(LY[i] - f(LX[i]))
        if d > max:
            max = d
    return max


# =====================================
def plot(LX, LY, M, n, points, x, type, degree=None, LY2=None):
    for i in range(n):
        plt.axvline(x=points[i][0], color='b')

    plt.ion()
    plt.plot(x(M), f(x(M)), 'k')

    if degree is None:
        l1, l2 = plt.plot(LX, LY, 'b', LX, LY2, 'r')
        plt.legend((l1, l2), ('third degree', 'second degree'))
    else:
        plt.plot(LX, LY, 'b')

    name = "plot" + str(n) + "_" + type
    if degree == 2:
        plot_dir = './plots_2nd_degree/'
    elif degree == 3:
        plot_dir = './plots_3rd_degree/'
    else:
        plot_dir = './mixed/'

    plt.savefig(plot_dir + name)
    plt.show()


if __name__ == "__main__":
    # create dir
    if not os.path.exists('./plots_2nd_degree'):
        os.mkdir('./plots_2nd_degree')
    if not os.path.exists('./plots_3rd_degree'):
        os.mkdir('./plots_3rd_degree')
    if not os.path.exists('./mixed'):
        os.mkdir('./mixed')

    mixed_second_and_third()
    solo_second_and_third()