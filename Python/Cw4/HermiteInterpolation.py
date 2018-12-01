import numpy as np
import math


def hermite_interpolation(points, x):
    pts = len(points)
    nx = len(x)

    dx = np.array([d[0] for d in points])
    dy = np.array([(d[1], d[2]) for d in points])
    m = 2

    L = np.zeros(nx)

    dd = divided_diference(dx, dy)
    w = np.ones(nx)

    for i in range(pts):
        p = x - dx[i]
        for j in range(m):
            l = i * m + j
            L += dd[l][l] * w
            w *= p
    return L


# dy[0][i] = f(xi), dy[1][i] = f'(xi)
def divided_diference(x, dy):
    n = len(dy)
    m = 2
    dd = np.zeros((n*m, n*m))
    z = np.zeros(n*m)

    for i in range(n):
        for j in range(m):
            k = i * m + j
            z[k] = x[i]
            dd[k][0] = dy[i][0]
            for l in range(1, k+1):
                if dd[k][l-1] == dd[k-1][l-1]:
                    dd[k][l] = dy[i][l] / math.factorial(l-1)
                else:
                    dd[k][l] = (dd[k][l-1] - dd[k-1][l-1]) / (z[k] - z[k-l])
    return dd
