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
def main():
    start, end = -10, 10
    M = 1000
    n = 50
    m = 6

    n += 1
    x = lambda n: np.linspace(start, end, n)
    LX = x(M)
    points = list_of_coords(start, end, n)

    a = least_squares_approximation(points, m)
    LY = get_points(a, LX)

    plot(LX, LY, M, n, points, x)


def plot(LX, LY, M, n, points, x):
    plt.ion()
    plt.plot(x(M), f(x(M)), 'k')

    dx = [d[0] for d in points]
    dy = [d[1] for d in points]
    l1, _ = plt.plot(LX, LY, 'r', dx, dy, 'ro')
    plt.grid()
    plt.savefig("xd")
    plt.show()


if __name__ == '__main__':
    main()
