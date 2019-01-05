import numpy as np


def runge_kutty_metohd(x0, y0, n, xn, degree, f):
    h = (xn - x0) / (n - 1)

    x = np.linspace(x0, xn, n)
    y = np.zeros([n])
    y[0] = y0

    def dy(i):
        k1 =  h * f(x[i], y[i])
        if degree == 1:
            return k1

        k2 = h * f(x[i] + h/2, y[i] + k1/2)
        if degree == 2:
            return k2

        k3 = h * f(x[i] + h/2, y[i] + k2/2)
        if degree == 3:
            return k3   # probably wrong

        k4 = h * f(x[i] + h, y[i] + k3)
        if degree == 4:
            return (k1 + 2*k2 + 2*k3 + k4) / 6
        exit(-54)

    for i in range(1, n):
        y[i] = y[i-1] + dy(i-1)

    return x, y