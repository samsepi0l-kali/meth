import math

#Matrix
c = [[-1, 0, 6],
     [1, -1, 5],
     [8, 4, -3]]

#Given Accuracy
e = 0.001


def qr_decomposition_method(n, c, e):
    #Working copy of the matrix
    a = [[float(c[i][j]) for j in range(n)] for i in range(n)]

    #Initial criterion value f
    f = 0.0
    g = a[1][0] * a[1][0]
    for i in range(1, n):
        for j in range(i + 1, n):
            f = f + a[j][i - 1] * a[j][i - 1]
        if a[i][i - 1] * a[i][i - 1] < g:
            g = a[i][i - 1] * a[i][i - 1]
    f = math.sqrt(f + g)

    k = 0

    #Iteration starts here
    while f > e:
        #Filling the initial Q matrix with the unity matrix
        q = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j:
                    q[i][j] = 1.0
                else:
                    q[i][j] = 0.0

        #QR-decomposition
        for p in range(1, n):
            v = [0.0] * n
            g = 0.0
            for i in range(1, p):
                v[i - 1] = 0.0
            for i in range(p, n + 1):
                v[i - 1] = a[i - 1][p - 1]
                g = g + a[i - 1][p - 1] * a[i - 1][p - 1]
            if a[p - 1][p - 1] > 0:
                v[p - 1] = v[p - 1] + math.sqrt(g)
            if a[p - 1][p - 1] < 0:
                v[p - 1] = v[p - 1] - math.sqrt(g)
            g = g - a[p - 1][p - 1] * a[p - 1][p - 1] + v[p - 1] * v[p - 1]

            #Householder matrix H
            h = [[0.0] * n for _ in range(n)]
            for i in range(n):
                for j in range(n):
                    if i == j:
                        h[i][j] = 1.0
                    else:
                        h[i][j] = 0.0
                    h[i][j] = h[i][j] - 2 * v[i] * v[j] / g

            #b = H * a
            b = [[0.0] * n for _ in range(n)]
            for i in range(n):
                for j in range(n):
                    b[i][j] = 0.0
                    for s in range(n):
                        b[i][j] = b[i][j] + h[i][s] * a[s][j]
            for i in range(n):
                for j in range(n):
                    a[i][j] = b[i][j]

            #q = q * h
            b = [[0.0] * n for _ in range(n)]
            for i in range(n):
                for j in range(n):
                    b[i][j] = 0.0
                    for s in range(n):
                        b[i][j] = b[i][j] + q[i][s] * h[s][j]
            for i in range(n):
                for j in range(n):
                    q[i][j] = b[i][j]

        #Matrices Q and R are multiplied in reverse order: a = R * Q
        b = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                b[i][j] = 0.0
                for s in range(n):
                    b[i][j] = b[i][j] + a[i][s] * q[s][j]
        for i in range(n):
            for j in range(n):
                a[i][j] = b[i][j]

        #Recalculate criterion value f
        f = 0.0
        g = a[1][0] * a[1][0]
        for i in range(1, n):
            for j in range(i + 1, n):
                f = f + a[j][i - 1] * a[j][i - 1]
            if a[i][i - 1] * a[i][i - 1] < g:
                g = a[i][i - 1] * a[i][i - 1]
        f = math.sqrt(f + g)

        #Iteration Counter
        k = k + 1

    #Extract eigenvalues from the quasi-triangular matrix
    w = [0.0] * n
    j = 0
    for i in range(1, n):
        if j == 0:
            #Checking for the presence of a complex-conjugate pair
            if abs(a[i][i - 1]) <= 10 * f:
                w[i - 1] = a[i - 1][i - 1]
            else:
                #Complex-conjugate pair from 2x2 block
                tr = a[i - 1][i - 1] + a[i][i]
                det = a[i - 1][i - 1] * a[i][i] - a[i - 1][i] * a[i][i - 1]
                disc = tr * tr - 4 * det
                
                re = tr / 2
                if disc >= 0:
                    d = math.sqrt(disc) / 2
                    w[i - 1] = re + d
                    w[i] = re - d
                else:
                    im = math.sqrt(-disc) / 2
                    w[i - 1] = complex(re, im)
                    w[i] = complex(re, -im)
                j = 1
        else:
            j = 0
    if j == 0:
        w[n - 1] = a[n - 1][n - 1]

    print("QR Decomposition Method:")
    print("Number of iterations:", k)
    print("Given accuracy:", e)
    print("Final criterion value:", f)
    print("Eigenvalues:")
    for i in range(n):
        print(f"  lambda_{i + 1} = {w[i]}")
    print()

    return w, k


print("QR Decomposition Method:")
print("Matrix:")
for row in c:
    print("  ", row)
print("Given accuracy:", e)
print()

eigenvalues, iters = qr_decomposition_method(len(c), c, e)