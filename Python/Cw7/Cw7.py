from matplotlib import pyplot as plt
from Euler_method import *
from Runge_Kutty_method import *
from aaa import *


# y'' = u + vy + wy'
#      [u, v, w]
def dd_fun():
    return {
        'u': lambda x: np.exp(x) - 3*np.sin(x),
        'v': lambda x: -1,
        'w': lambda x: 1
    }


def main():
    f = lambda x, y: -2 * x * y

    x0, y0 = 1, 1.09737491
    n = 100
    xn, yn = 2, 8.63749661

    #x, y = runge_kutty_metohd(x0, y0, n, xn, 4, f)
    x, y = finite_difference_method(x0, y0, xn, yn, n, dd_fun())
    for i in range(n):
        print(x[i], y[i])

    plot(x, y)


def plot(x, y):
    plt.plot(x, y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()


if __name__ == '__main__':
    main()