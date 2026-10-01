#Matrix
a = [[7, 8, 4, -6],
     [-1, 6, -2, -6],
     [2, 9, 6, -4],
     [5, 9, 1, 1]]

#Vector
b = [-126, -42, -115, -67]


a_original = [row[:] for row in a]


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
    for k in range(n - 1, -1, -1):
        y *= a[k][k]

        h = sum(a[k][i] * x[i] for i in range(k + 1, n))
        x[k] = (b[k] - h) / a[k][k]

        for j in range(n):
            h = sum(a[k][i] * z[i][j] for i in range(k + 1, n))
            z[k][j] = (z[k][j] - h) / a[k][k]

    #Print
    print("x =", " ".join(f"{v:.7f}" for v in x))
    print("y =", f"{y:.7f}")
    print("z:")
    for row in z:
        print("  ", "\t".join(f"{v:.7f}" for v in row))

    return x, y, z


#Function call
x, y, z = gaussian_pivot(a, b, 1)


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