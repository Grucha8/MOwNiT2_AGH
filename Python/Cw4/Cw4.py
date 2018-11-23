import numpy as np
import matplotlib.pyplot as plt


def f(x):
    return np.sin(x*2) * np.sin((2 * x**2) / np.pi)


# ==============================================
def list_of_coords(a, b, n):
    points = []
    d = (b - a) / n
    x = a
    n += 1
    for _ in range(n):
        points.append((x, f(x)))
        x += d

    return points


def list_of_chebyshev_coords(a, b, k):
    k += 1
    points = []
    for j in range(1, k+1):
        x = np.cos(((2 * j - 1) * np.pi) / float(2 * k)) * (b-a)/2 + (a+b)/2
        points.append((x, f(x)))
    return points


# ================================================
def lagrange_interpolation(points, x):
    n = len(points)
    nx = len(x)

    dx = [d[0] for d in points]
    dy = [d[1] for d in points]

    L = [0.0] * (nx)

    def b(j, xi):
        v = 1.0
        for k in range(n):
            if k != j:
                v *= (xi - dx[k]) / (dx[j] - dx[k])
        return v

    for i, xi in enumerate(x):
        for j in range(n):
            L[i] += dy[j] * b(j, xi)

    return L


# ==============================================
def newton_interpolation(points, x):
    pts = len(points)
    nx = len(x)

    dx = [d[0] for d in points]
    dy = [d[1] for d in points]

    L = [0.0] * (nx)

    def a(j0, j1=None):
        if j1 is None:
            j1, j0 = j0, 0

        if j0 == j1:
            return dy[j0]
        elif j1 - j0 == 1:
            return (dy[j1] - dy[j0]) / (dx[j1] - dx[j0])
        else:
            return (a(j0 + 1, j1) - a(j0, j1 - 1)) / (dx[j1] - dx[j0])

    def n(j, x_):
        v = 1.0
        for i in range(0, j):
            v *= float(x_ - dx[i])
        return v

    for i in range(nx):
        for j in range(pts):
            L[i] += a(j) * n(j, x[i])

    return L


# ==============================================

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


# ============================================================
def plot(LX, LY1, LY2, M, n, points1, points2, x, inter):
    # putting vertical lines which indicates nodes for noramal <blue>
    for i in range(n+1):
        plt.axvline(x=points2[i][0], color='b')

    # chebushev nodes <red>
    for i in range(n+1):
        plt.axvline(x=points1[i][0], color='r')

    # black plot == actual function
    plt.ion()
    plt.plot(x(M), f(x(M)), 'k')

    # red plot == interpolated function chebyshev
    plt.plot(LX, LY1, 'r')

    # _ plot == interpolated function
    plt.plot(LX, LY2, 'b')

    name = "plot" + str(n) + inter + ".png"
    plt.savefig(name)
    plt.show()


if __name__ == "__main__":
    main()