import numpy as np


# ddf_eq = {'u': u(x), 'v': v(x), 'w': w(x)} \ - lambda
def finite_difference_method(x0, y0, xn, yn, n, ddf_eq):
    h = (xn - x0) / n

    x = np.linspace(x0, xn, n)
    a, b, d, c = np.zeros([n]), np.zeros([n]), np.zeros([n]), np.zeros([n])
    for i in range(1, n-1):
        a[i] = -(1 + ddf_eq['w'](x[i]) * h / 2)
        b[i] = -ddf_eq['u'](x[i]) * h**2
        c[i] = -(1 - ddf_eq['w'](x[i]) * h / 2)
        d[i] = (2 + ddf_eq['v'](x[i]) * h**2)

    A = np.zeros([n, n])
    for i in range(2, n-2):
        A[i][i] = d[i]
        A[i][i-1] = a[i]
        A[i][i+1] = c[i]
    A[1][1], A[n-2][n-2] = d[1], d[n-2]
    A[1][2] = c[1]
    A[n-2][n-3] = a[n-2]

    b[1] = b[1] - a[1] * y0
    b[n-2] = b[n-2] - c[n-2] * yn
    b[0], b[n-1] = y0, yn
    A[0][0] = 1
    A[n-1][n-1] = 1

    y = np.linalg.solve(A, b)   # Ay = b

    return x, y



