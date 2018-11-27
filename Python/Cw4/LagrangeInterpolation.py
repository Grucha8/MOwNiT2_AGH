def lagrange_interpolation(points, x):
    n = len(points)
    nx = len(x)

    dx = [d[0] for d in points]
    dy = [d[1] for d in points]

    L = [0.0] * (nx)

    def b(j, xi):
        v = 1.0
        for k in range(n):
            if k != j:
                v *= (xi - dx[k]) / (dx[j] - dx[k])
        return v

    for i, xi in enumerate(x):
        for j in range(n):
            L[i] += dy[j] * b(j, xi)

    return L