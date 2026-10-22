import math

#Matrix
c = [[-8, 5, -7],
     [5, 1, 4],
     [-7, 4, 4]]

#Given Accuracy
e = 0.01


def check_eigenpair(n, c, lam, x):
    #Residual vector cx - lambdax (should be close to zero)
    r = [0.0] * n
    for i in range(n):
        r[i] = 0.0
        for j in range(n):
            r[i] = r[i] + c[i][j] * x[j]
        r[i] = r[i] - lam * x[i]
    return r


def jacobi_rotation_method(n, c, e):
    a = [[float(c[i][j]) for j in range(n)] for i in range(n)]
    u = [[0.0] * n for _ in range(n)]
    b = [[0.0] * n for _ in range(n)]
    v = [[0.0] * n for _ in range(n)]
    w = [0.0] * n

    #Filling the eigenvectors matrix with the unity matrix
    for i in range(n):
        for j in range(n):
            if i == j:
                v[i][j] = 1.0
            else:
                v[i][j] = 0.0

    #Critertion check
    f = 0.0
    for i in range(n - 1):
        for j in range(i + 1, n):
            f = f + a[i][j] * a[i][j]
    f = math.sqrt(f)
    k = 0

    #Iteration starts here
    while f > e:
        #Search for the maximum in modulus non-diagonal element
        g = 0.0
        l = 0
        m = 1

        for i in range(n - 1):
            for j in range(i + 1, n):
                if a[i][j] > g:
                    g = a[i][j]
                    l = i
                    m = j

                if -a[i][j] > g:
                    g = -a[i][j]
                    l = i
                    m = j

        #Rotation matrix U filled with the unity matrix
        for i in range(n):
            for j in range(n):
                if i == j:
                    u[i][j] = 1.0
                else:
                    u[i][j] = 0.0

        #Rotation angle
        if a[l][l] == a[m][m]:
            h = math.pi / 4
        else:
            h = math.atan(2 * a[l][m] / (a[l][l] - a[m][m])) / 2

        u[l][l] = math.cos(h)
        u[l][m] = -math.sin(h)
        u[m][l] = math.sin(h)
        u[m][m] = math.cos(h)

        #b = U^T * A (only rows and columns l, m change)
        for i in range(n):
            for j in range(n):
                if i == l or i == m or j == l or j == m:
                    b[i][j] = 0.0

                    for p in range(n):
                        b[i][j] = b[i][j] + u[p][i] * a[p][j]

        #A = b * U (only rows and columns l, m change)
        for i in range(n):
            for j in range(n):
                if i == l or i == m or j == l or j == m:
                    a[i][j] = 0.0

                    for p in range(n):
                        a[i][j] = a[i][j] + b[i][p] * u[p][j]

        #b = V * U (only rows and columns l, m change)
        for i in range(n):
            for j in range(n):
                if i == l or i == m or j == l or j == m:
                    b[i][j] = 0.0

                    for p in range(n):
                        b[i][j] = b[i][j] + v[i][p] * u[p][j]

        #V = b
        for i in range(n):
            for j in range(n):
                if i == l or i == m or j == l or j == m:
                    v[i][j] = b[i][j]

        #Criterion for iterations termination
        f = 0.0
        for i in range(n - 1):
            for j in range(i + 1, n):
                f = f + a[i][j] * a[i][j]
        f = math.sqrt(f)
        #Iteration Counter
        k += 1

    #Eigenvalues are the diagonal elements
    for i in range(n):
        w[i] = a[i][i]

    print("Eigenvalues w:")
    print("  ", "\t".join(f"{x:.7f}" for x in w))

    print("Eigenvectors (columns of V):")
    for j in range(n):
        print("  v" + str(j + 1) + " =", "\t".join(f"{v[i][j]:.7f}" for i in range(n)))

    #Check: Av - lambda*v for each eigenpair
    print("Check A*v - lambda*v:")
    for j in range(n):
        x = [v[i][j] for i in range(n)]
        r = check_eigenpair(n, c, w[j], x)
        print("  v" + str(j + 1) + ":", "\t".join(f"{t:.7f}" for t in r))

    print("Number of iterations:", k)
    print("Given accuracy:", e)
    print("Final criterion value:", f)

    return w, v, k


print("Jacobi Rotation Method:")
result = jacobi_rotation_method(len(c), c, e)