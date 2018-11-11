import numpy as np


def f(x):
    return (x+1)**2 - 10


def fprim(x, dx=1e-5):
    return (f(x + dx) - f(x)) / dx


def secant(a, b, nmax, eps):
    fa = f(a)
    fb = f(b)

    if np.abs(fa) > np.abs(fb):
        a, b = b, a
        fa, fb = fb, fa

    for n in range(2, nmax):
        if np.abs(fa) > np.abs(fb):
            a, b = b, a
            fa, fb = fb, fa

        d = (b - a) / (fb - fa)
        b, fb = a, fa
        d = d * fa
        if np.abs(d) < eps:
            print("spoko")
            return
        a = a - d
        fa = f(a)
        print(n, a, fa, sep="\t")


def newton(x, nmax, eps, delta):
    fx = f(x)

    for n in range(nmax):
        fp = fprim(x)
        if np.abs(fp) < delta:
            print("mala pochodna")
            return

        d = fx / fp
        x -= d
        fx = f(x)

        print(n, x, fx, sep="\t")
        if np.abs(d) < eps:
            print("spoko")
            return


def main():
    secant(-10, 10, 30, 1e-5)
    print("===============")
    newton(-10, 30, 1e-5, 1e-5)


if __name__ == '__main__':
    main()