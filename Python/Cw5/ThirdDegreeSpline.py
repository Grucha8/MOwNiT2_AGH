import numpy as np

def third_degree_spline(points, LX):
    n = len(points) - 1
    nx = len(LX)

    t = [d[0] for d in points]
    y = [d[1] for d in points]
    h = [(t[i+1] - t[i]) for i in range(n)]
    z = tridiagonal_solver(h, y, n + 1)

    def S(i, x):
        p1 = z[i+1] * (x - t[i])**3 / (6 * h[i]) + z[i] * (t[i+1] - x)**3 / (6 * h[i])
        p2 = ((y[i+1] / h[i]) - (h[i] * z[i+1] / 6)) * (x - t[i])
        p3 = ((y[i] / h[i]) - h[i] * z[i] / 6) * (t[i+1] - x)
        return p1 + p2 + p3

    L = [0.0] * nx

    for i, xi in enumerate(LX):
        for j, ti in enumerate(points):
            if xi < ti[0]:
                L[i] = S(j - 1, xi)
                break
    L[0] = y[0]

    return L


def tridiagonal_solver(h, y, n):
    u = np.zeros(n)
    v = np.zeros(n)

    for i in range(1, n-1):
        u[i] = 2 * (h[i-1] + h[i])
        v[i] = 6 * ((y[i+1] - y[i]) / h[i]) - ((y[i] - y[i-1]) / h[i-1])

    for i in range(2, n-1):
        u[i] -= h[i-1]**2 / u[i-1]
        v[i] -= h[i-1] * v[i-1] / u[i-1]

    z = np.zeros(n)
    z[n-2] = v[n-2] / u[n-2]
    for i in range(n-3, 0, -1):
        z[i] = (v[i] - h[i] * z[i+1]) / u[i]

    return z