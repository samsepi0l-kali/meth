#Matrix
a = [[7, 8, 4, -6],
     [-1, 6, -2, -6],
     [2, 9, 6, -4],
     [5, 9, 1, 1]]

#Vector
b = [-126, -42, -115, -67]


a_original = [row[:] for row in a]


def lup_decomposition(a, b, l):
    n = len(a)
    q = list(range(n))
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
                q[k], q[m] = q[m], q[k]
                p += 1

        if a[k][k] == 0:
            print("Singular matrix!")
            return

        #Elimination
        for i in range(k + 1, n):
            a[i][k] = a[i][k] / a[k][k]
            for j in range(k + 1, n):
                a[i][j] -= a[i][k] * a[k][j]

    #Check PA = LU
    print("\nCheck PA = LU:")
    for i in range(n):
        row_check = []
        for j in range(n):
            pa = a_original[q[i]][j]
            lu = 0.0
            for k in range(n):
                l_ik = a[i][k] if i > k else (1.0 if i == k else 0.0)
                u_kj = a[k][j] if k <= j else 0.0
                lu += l_ik * u_kj
            if abs(pa - lu) < 1e-9:
                row_check.append(0.0)
            else:
                row_check.append(pa - lu)
        print("  ", "\t".join(f"{v:.7f}" for v in row_check))

    #Determinant
    y = 1 - 2 * (p % 2)
    for k in range(n):
        y *= a[k][k]

    #Backward path
    d = [b[q[i]] for i in range(n)]

    #Solving for LZ = d
    z = [0.0] * n
    for i in range(n):
        s = d[i]
        for j in range(i):
            s -= a[i][j] * z[j]
        z[i] = s
    #Solving for UX = z
    for i in range(n - 1, -1, -1):
        s = z[i]
        for j in range(i + 1, n):
            s -= a[i][j] * x[j]
        x[i] = s / a[i][i]

    z_matrix = [[0.0] * n for _ in range(n)]

    #Inverse Matrix
    for col in range(n):
        d = [1.0 if q[i] == col else 0.0 for i in range(n)]

        #Solving LZ = d for each column after each iteration
        z = [0.0] * n
        for i in range(n):
            s = d[i]
            for j in range(i):
                s -= a[i][j] * z[j]
            z[i] = s

        xcol = [0.0] * n
        for i in range(n - 1, -1, -1):
            s = z[i]
            for j in range(i + 1, n):
                s -= a[i][j] * xcol[j]
            xcol[i] = s / a[i][i]

        for i in range(n):
            z_matrix[i][col] = xcol[i]

    #Print
    print("x =", " ".join(f"{v:.7f}" for v in x))
    print("y =", f"{y:.7f}")
    print("z:")

    for row in z_matrix:
        print("  ", "\t".join(f"{v:.7f}" for v in row))

    return x, y, z_matrix


#Function call
x, y, z = lup_decomposition(a, b, 1)


#Verify A · A⁻¹ = I
print("\nCheck A · A⁻¹:")
for i in range(4):
    row_check = []
    for j in range(4):
        total = 0.0
        for k in range(4):
            total += a_original[i][k] * z[k][j]
        if abs(total) < 1e-9:
            total = 0.0
        elif abs(total - 1) < 1e-9:
            total = 1.0
        row_check.append(total)
    print("  ", "\t".join(f"{v:.7f}" for v in row_check))