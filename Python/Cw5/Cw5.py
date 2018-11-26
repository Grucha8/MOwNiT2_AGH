import numpy as np
import matplotlib.pyplot as plt


def f(x):
    return np.sin(x * 2) * np.sin((2 * x ** 2) / np.pi)


# =======================================================
# points = [(x, f(x))]
def second_degree_spline(points, LX, z0):
    n = len(points) - 1
    nx = len(LX)

    dt = [d[0] for d in points]
    dy = [d[1] for d in points]

    z = [0.0] * (n + 1)
    z[0] = z0
    for i in range(1, n+1):
        z[i] = -z[i-1] + 2 * ((dy[i] - dy[i-1]) / (dt[i] - dt[i-1]))

    Q = (lambda x, c:
         ((z[c+1] - z[c]) / (2 * (dt[c+1] - dt[c]))) * (x - dt[c])**2 + z[c] * (x - dt[c]) + dy[c])

    L = [0.0] * nx

    for i, xi in enumerate(LX):
        for j, ti in enumerate(points):
            if xi < ti[0]:
                L[i] = Q(xi, j-1)
                break
    L[0] = dy[0]

    return L


# =========================================
def list_of_coords(a, b, n):
    points = []
    d = (b - a) / (n - 1)
    x = a
    for _ in range(n):
        points.append((x, f(x)))
        x += d

    return points


def main():
    start, end = -np.pi, np.pi
    n = 6
    M = 1000

    x = lambda n: np.linspace(start, end, n)
    LX = x(M)
    p = list_of_coords(start, end, n)
    LY = second_degree_spline(p, LX, 0)

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