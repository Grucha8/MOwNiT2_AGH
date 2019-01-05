import numpy as np


def euler_method(x0, y0, n, xn, f):
    h = (xn - x0) / (n - 1)

    x = np.linspace(x0, xn, n)
    y = np.zeros([n])
    y[0] = y0

    for i in range(1, n):
        y[i] = h * f(x[i-1], y[i-1]) + y[i-1]

    return x, y