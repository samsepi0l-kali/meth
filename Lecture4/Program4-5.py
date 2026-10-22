#Matrix 
c = [[10, -1, -2, 5],
     [4, 28, 7, 9],
     [6, 5, -23, 4],
     [1, 4, 5, -15]]
     

#Vector
d = [-99, 0, 67, 58]

#Given Accuracy 
e = 0.01


def fixed_point_iteration_method_seidel_method(n, c, d, e, l):
    a = [[0.0] * n for _ in range(n)]
    b = [0.0] * n
    x = [0.0] * n
    y = [0.0] * n

    #Reduction to equivalent form x = b + a*x 
    for i in range(n):
        if c[i][i] == 0:
            print("Zero diagonal element found, exiting.")
            return

        for j in range(n):
            if i != j:
                a[i][j] = -c[i][j] / c[i][i]
            else:
                a[i][j] = 0.0

        b[i] = d[i] / c[i][i]

    #Norm of matrix A
    g = 0.0

    for i in range(n):
        h = 0.0

        for j in range(n):
            if a[i][j] > 0:
                h += a[i][j]
            else:
                h -= a[i][j]

        if h > g:
            g = h

    print("Matrix norm g =", f"{g:.7f}")

    #Convergence condition
    if g >= 1:
        print("The convergence condition is not fulfilled!")
        return

    #Initial approximation x = b and initial error estimate
    f = 0.0

    for i in range(n):
        x[i] = b[i]

        if x[i] > f:
            f = x[i]

        if -x[i] > f:
            f = -x[i]

    f = f * g / (1 - g)

    k = 0

    #Iteration process
    while f > e:
        for i in range(n):
            y[i] = x[i]

        for i in range(n):
            x[i] = b[i]

            if l == 0:
                for j in range(n):
                    x[i] = x[i] + a[i][j] * y[j]

            else:
                for j in range(i):
                    x[i] = x[i] + a[i][j] * x[j]

                for j in range(i, n):
                    x[i] = x[i] + a[i][j] * y[j]

        # Error estimate
        f = 0.0

        for i in range(n):
            if x[i] - y[i] > f:
                f = x[i] - y[i]

            if y[i] - x[i] > f:
                f = y[i] - x[i]

        f = f * g / (1 - g)
        
        #Iteration Counter
        k += 1

    print("Solution vector x:")
    print("  ", "\t".join(f"{v:.7f}" for v in x))
    print("Number of iterations:", k)
    print("Given accuracy:", e)
    print("Final error estimate:", f)

    return x, k


print("Fixed-point Iteration Method:")
result = fixed_point_iteration_method_seidel_method(
    len(c), c, d, e, 0
)

print("\n-----------------------------\n")

print("Seidel Method:")
result = fixed_point_iteration_method_seidel_method(
    len(c), c, d, e, 1
)