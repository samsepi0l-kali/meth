import math

#Matrix
c = [[-1, 0, 6],
     [1, -1, 5],
     [8, 4, -3]]

#Known eigenvalues from QR decomposition method (Program8.py)
#Replace these with the output from your QR run
eigenvalues = [6.497306221842258, -10.147158601184731, -1.350147620657524]

#Given Accuracy
e = 0.001


def check_eigenpair(n, c, lam, x):
    #Residual vector cx - lambdax (should be close to zero)
    r = [0.0] * n
    for i in range(n):
        r[i] = 0.0
        for j in range(n):
            r[i] = r[i] + c[i][j] * x[j]
        r[i] = r[i] - lam * x[i]
    return r


def mat_inverse(n, m):
    #Gauss-Jordan inversion of a real n x n matrix
    aug = [[float(m[i][j]) for j in range(n)] + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for col in range(n):
        #Pivot
        piv = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[piv][col]) < 1e-15:
            raise ValueError("Matrix is singular")
        aug[col], aug[piv] = aug[piv], aug[col]
        pv = aug[col][col]
        for j in range(2 * n):
            aug[col][j] /= pv
        for r in range(n):
            if r != col:
                fac = aug[r][col]
                for j in range(2 * n):
                    aug[r][j] -= fac * aug[col][j]
    return [[aug[i][n + j] for j in range(n)] for i in range(n)]


def mat_inverse_complex(n, m):
    #Gauss-Jordan inversion of a complex n x n matrix
    aug = [[complex(m[i][j]) for j in range(n)] + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[piv][col]) < 1e-15:
            raise ValueError("Matrix is singular")
        aug[col], aug[piv] = aug[piv], aug[col]
        pv = aug[col][col]
        for j in range(2 * n):
            aug[col][j] /= pv
        for r in range(n):
            if r != col:
                fac = aug[r][col]
                for j in range(2 * n):
                    aug[r][j] -= fac * aug[col][j]
    return [[aug[i][n + j] for j in range(n)] for i in range(n)]


def shift_matrix(n, a, mu):
    #Return (A - mu*E)
    return [[a[i][j] - (mu if i == j else 0) for j in range(n)] for i in range(n)]


def inverse_iteration_method(n, a, b, e, mu, label=""):
    #Inverse iteration for a given shift mu (b = (A - mu*E)^(-1))
    y = [1.0] * n
    f = 2 * e
    k = 0
    w = 0.0

    while f > e:
        #z = B * y
        z = [0.0] * n
        r = 0.0
        for i in range(n):
            z[i] = 0.0
            for j in range(n):
                z[i] = z[i] + b[i][j] * y[j]
            if z[i] > r:
                r = z[i]
            if -z[i] > r:
                r = -z[i]

        if k > 0:
            f = abs(w)

        #Normalization of the iterated vector
        for i in range(n):
            y[i] = z[i] / r

        #Refining the eigenvalue: lambda = (y^T A y) / (y^T y)
        z = [0.0] * n
        for i in range(n):
            z[i] = 0.0
            for j in range(n):
                z[i] = z[i] + a[i][j] * y[j]

        w_new = 0.0
        r = 0.0
        for j in range(n):
            w_new = w_new + y[j] * z[j]
            r = r + y[j] * y[j]
        w_new = w_new / r

        if k > 0:
            f = w_new - w
        w = w_new

        if f < 0:
            f = -f

        #Iteration Counter
        k = k + 1

    print(f"Inverse Iteration Method ({label}):")
    print("Number of iterations:", k)
    print("Given accuracy:", e)
    print("Refined eigenvalue:", f"{w:.7f}")
    print("Corresponding eigenvector y (normalized):")
    print("  y =", "\t".join(f"{yi:.7f}" for yi in y))

    #Check: Av - lambda*v (should be close to zero)
    print("Check A*y - lambda*y:")
    r = check_eigenpair(n, a, w, y)
    print("  ", "\t".join(f"{t:.7f}" for t in r))
    print()

    return w, y, k


def inverse_iteration_complex(n, a, mu, e, label=""):
    #Complex version of inverse iteration (for complex shifts)
    m = shift_matrix(n, a, mu)
    b = mat_inverse_complex(n, m)

    y = [1.0 + 0j] * n
    f = 2 * e
    k = 0
    w = 0.0 + 0j

    while f > e:
        z = [0j] * n
        r = 0.0
        for i in range(n):
            z[i] = 0j
            for j in range(n):
                z[i] = z[i] + b[i][j] * y[j]
            if abs(z[i]) > r:
                r = abs(z[i])

        if k > 0:
            f = abs(w)

        for i in range(n):
            y[i] = z[i] / r

        z = [0j] * n
        for i in range(n):
            z[i] = 0j
            for j in range(n):
                z[i] = z[i] + a[i][j] * y[j]

        w_new = 0j
        r = 0.0
        for j in range(n):
            w_new = w_new + y[j].conjugate() * z[j]
            r = r + abs(y[j]) ** 2
        w_new = w_new / r

        if k > 0:
            f = abs(w_new - w)
        w = w_new

        #Iteration Counter
        k = k + 1

    print(f"Inverse Iteration Method ({label}):")
    print("Number of iterations:", k)
    print("Given accuracy:", e)
    print("Refined eigenvalue:", w)
    print("Corresponding eigenvector y (normalized):")
    print("  y =", "\t".join(f"{yi:.7f}" for yi in y))
    print()

    return w, y, k



n = len(c)

print("Inverse Iteration Method:")
print("Matrix:")
for row in c:
    print("  ", row)
print("Given accuracy:", e)
print()

for idx, lam in enumerate(eigenvalues):
    if isinstance(lam, complex):
        #Complex eigenvalue -> use complex inverse iteration
        inverse_iteration_complex(n, c, lam, e, label=f"lambda_{idx + 1} = {lam}")
    else:
        #Real eigenvalue -> use real inverse iteration
        m = shift_matrix(n, c, lam)
        b = mat_inverse(n, m)
        inverse_iteration_method(n, c, b, e, lam, label=f"lambda_{idx + 1} = {lam:.7f}")