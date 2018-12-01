def newton_interpolation(points, x):
    pts = len(points)
    nx = len(x)

    dx = [d[0] for d in points]
    dy = [d[1] for d in points]

    L = [0.0] * nx

    A = [None] * (pts + 1)
    for i in range(pts + 1):
        A[i] = [None] * (pts + 1)

    def a(j0, j1=None):
        if j1 is None:
            j1, j0 = j0, 0

        if j0 == j1:
            A[j0][j1] = dy[j0]
            return dy[j0]
        elif j1 - j0 == 1:
            A[j0][j1] = (dy[j1] - dy[j0]) / (dx[j1] - dx[j0])
            return (dy[j1] - dy[j0]) / (dx[j1] - dx[j0])
        else:
            if A[j0 + 1][j1] is None:
                A[j0 + 1][j1] = a(j0 + 1, j1)

            if A[j0][j1 - 1] is None:
                A[j0][j1 - 1] = a(j0, j1 - 1)

            return (A[j0 + 1][j1] - A[j0][j1 - 1]) / (dx[j1] - dx[j0])

    def n(k, x_):
        v = 1.0
        for l in range(0, k):
            v *= float(x_ - dx[l])
        return v

    for i in range(nx):
        for j in range(pts):
            L[i] += a(j) * n(j, x[i])

    return L
