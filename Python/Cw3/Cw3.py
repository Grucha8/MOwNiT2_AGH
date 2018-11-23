import numpy as np
from math import floor
import csv

#[-2.5, 0.5]
NMAX = 500  # maksymalna liczba iteracji

def f(x):
    return x**14 + x**16


def fprim_(x):
    return 14 * x**13 + 16 * x**15


def secant(a, b, nmax, eps, stop):
    a_prev = 1
    fa = f(a)
    fb = f(b)
    if fa == fb:
        print("skiped")
        return

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

        if stop and n > 2:
            if np.abs(a - a_prev) < eps:
                print(n, a, fa, sep="\t\t")
                return
        else:
            if np.abs(fa) < eps:
                print(n, a, fa, sep="\t\t")
                return
        a_prev = a
        a = a - d
        fa = f(a)


def newton(x, nmax, eps, stop):
    fx = f(x)

    for n in range(nmax):

        fp = fprim_(x)

        try:
            d = fx / fp
        except ZeroDivisionError:
            print("Division by zer. Aborting this x0")
            return
        x_prev = x
        x -= d
        fx = f(x)

        # warunek stopu
        if stop == 0:
            if np.abs(x - x_prev) < eps:
                print(n, x, fx, sep="\t\t")
                return
        else:
            if np.abs(fx) < eps:
                print(n, x, fx, sep="\t\t")
                return

        if n == nmax-1:
            print("Brak")


##############################################################################3
def f1(x1, x2, x3):
    return np.array(
        [x1**2 + x2**2 + x3 - 1,
         2*x1**2 - x2**2 - 4*x3**2 + 3,
         x1**2 + x2 + x3 - 1]
    )


def df1(x1, x2, x3):
    return np.array(
        [[2*x1, 2*x2, 1],
         [4*x1, -2*x2, -8*x3],
         [2*x1, 1, 1]]
    )


def non_linear_solver(x0, max_iter, ro, stop):
    x = x0
    k = 0
    while k < max_iter:
        f_val = f1(x[0], x[1], x[2])
        df_val = df1(x[0], x[1], x[2])

        try:
            s = np.linalg.solve(df_val, -f_val)
        except np.linalg.LinAlgError:
            print("Singluar matrix, skiping this x0")
            return

        if stop:
            l = np.linalg.norm(np.abs(s), ord=2)
            if l < ro:
                print(k, x+s, l, sep='\t\t\t')
                return
        else:
            l = np.linalg.norm(f1((x+s)[0], (x+s)[1], (x+s)[2]), np.Inf)
            if l < ro:
                print(k, x + s, l, sep='\t\t\t')
                return

        x = x + s
        #print(f1(x[0],x[1],x[2]))
        k += 1

    print("Nie znaleziono wyniku")
################################################################################


def part1():
    epss = [0.1, 0.01, 0.001, 0.0001, 1e-5, 1e-6, 1e-7, 1e-10]

    print("METODA NEWTONA")
    for eps in epss:
        print("---------------------------\neps =", eps)
        x = -2.5
        while (x < 0.6):
            print("\nx =", x)
            print("Pierwszy warunek stopu".capitalize())
            newton(x, NMAX, eps, 0)
            print("Drugi warunek stopu".capitalize())
            newton(x, NMAX, eps, 1)
            x = (int)((x + 0.1) * 10) / 10
    print("=========================================================================")
    print("METODA SIECZNYCH")
    for eps in epss:
        print("--------------------------------------------------\neps =", eps)
        x = -2.4
        while (x < 0.6):
            print("\nx =", x)
            print("Pierwszy warunek stopu".capitalize())
            secant(-2.5, x, NMAX, eps, 0)
            print("Drugi warunek stopu".capitalize())
            secant(-2.5, x, NMAX, eps, 1)
            x = x + 0.1
        print("######################")
        x = 0.4
        while (x > -2.5):
            print("\nx =", x)
            print("Pierwszy warunek stopu".capitalize())
            secant(x, 0.5, NMAX, eps, 0)
            print("Drugi warunek stopu".capitalize())
            secant(x, 0.5, NMAX, eps, 1)
            x = x - 0.1


def part2():
    epss = [0.01, 0.0001, 1e-6, 1e-10]
    xs = [np.array([i, j, k])
          for i in np.arange(-2.0, 2.0, 0.1)
          for j in np.arange(-2.0, 2.0, 0.1)
          for k in np.arange(-2.0, 2.0, 0.1)
          ]

    for eps in epss:
        print("---------------------------\neps =", eps)
        for xd, i in zip(xs, range(len(xs))):
            print("\nWektor poczatkowy: ", end='')
            print("Wektor: ", xd)

            print("Pierwsze kryterium stopu")
            non_linear_solver(xd, 20, eps, 1)
            print("Drugie kryterium stopu")
            non_linear_solver(xd, 20, eps, 0)


if __name__ == '__main__':
    np.set_printoptions(precision=5, suppress=True)
    part1()

    #part2()
    print("===========")