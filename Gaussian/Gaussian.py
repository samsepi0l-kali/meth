#Matrix
a = [[7, 8, 4, -6],
     [-1, 6, -2, -6],
     [2, 9, 6, -4],
     [5, 9, 1, 1]]

#Vector
b = [-126, -42, -115, -67]

def gaussian_pivot(a, b, l):
    n = len(a)
    z = [[float(i == j) for j in range(n)] for i in range(n)]
    p = 0
    x = [0.0] * n

    #Forward path
    for k in range(n - 1):
        #Pivot selection
        if l != 0:
            g = abs(a[k][k])
            m = k
            for i in range(k + 1, n):
                f = abs(a[i][k])
                if f > g:
                    g, m = f, i
            if m != k:
                a[k], a[m] = a[m], a[k]
                b[k], b[m] = b[m], b[k]
                z[k], z[m] = z[m], z[k]
                p += 1

        if a[k][k] == 0:
            print("Singular matrix!")
            return

        #Elimination
        for i in range(k + 1, n):
            h = -a[i][k] / a[k][k]
            for j in range(k, n):
                a[i][j] += h * a[k][j]
            b[i] += h * b[k]
            for j in range(n):
                z[i][j] += h * z[k][j]


    #Determinant
    y = 1 - 2 * (p % 2)

    #Backward path 
    for k in range(n-1, -1 ,-1):
        y *= a[k][k]

        h = sum(a[k][i] * x[i] for i in range(k + 1, n))
        x[k] = (b[k] - h) / a[k][k]

        for j in range(n):
            h = sum(a[k][i] * z[i][j] for i in range(k + 1, n))
            z[k][j] = (z[k][j] - h) / a[k][k]

    
    print("x =", " ".join(f"{v:.7f}" for v in x))
    print("y =", f"{y:.7f}")
    print("z:")
    for row in z:
        print("  ", "\t".join(f"{v:.7f}" for v in row))

#Function call
gaussian_pivot(a, b, 1)