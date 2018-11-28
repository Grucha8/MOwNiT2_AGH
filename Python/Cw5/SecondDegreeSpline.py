import numpy as np


def second_degree_spline(points, LX, is_natural):
    n = len(points) - 1
    nx = len(LX)

    dt = [d[0] for d in points]
    dy = [d[1] for d in points]

    z = np.zeros(n + 1)

    if is_natural:
        z[0] = 0
        for i in range(1, n+1):
            z[i] = -z[i-1] + 2 * ((dy[i] - dy[i-1]) / (dt[i] - dt[i-1]))
    else:
        # Q[0] - linear; Q''[0] = 0; z[0] == z[1]
        z[0] = z[1] = (dy[1] - dy[0]) / (dt[1] - dt[0])
        for i in range(2, n+1):
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