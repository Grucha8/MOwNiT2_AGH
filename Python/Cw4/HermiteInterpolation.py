import numpy as np


def hermite_interpolation(points, x):
    pts = len(points)
    nx = len(x)
    m = pts

    dx = [d[0] for d in points]
    dy = [d[1] for d in points]

    L = [0.0] * nx

    def divided_diference(d):
        if len(d) == 1:
            for i in range(pts):
                if dx[i] == d[0]:
                    return dy[i]


        l = len(d) // 2
        return (divided_diference(d[l:]) - divided_diference(d[:l])) / (d[-1] - d[0])

    # A = np.zeros((pts, m))
    A = [0.0] * pts
    for i in range(pts):
        A[i] = [0.0] * m
    A[0][0] = dy[0]
    for i in range(1, pts):
        for j in range(pts + 1):
            divided_diference(dy[(i - j):(i+1)])

    def n(k, x_):
        v = 1.0
        for l in range(0, k):
            v *= float(x_ - dx[l])
        return v

    for i in range(nx):
        for j in range(pts):
            L[i] += A[j][j] * n(j, x[i])

    return L

