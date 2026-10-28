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

def power_iteration_method(n, c, e):
    a = [[float(c[i][j]) for j in range(n)] for i in range(n)]

    #Filling the initial eigenvector with ones
    y = [1.0] * n

    #Initial criterion value
    f = 2 * e
    k = 0

    #Previous eigenvalue estimate
    lam_old = 0.0

    #Iteration starts here
    while f > e:
        #z = A * y
        z = [0.0] * n
        for i in range(n):
            z[i] = 0.0
            for j in range(n):
                z[i] = z[i] + a[i][j] * y[j]

        #Eigenvalue estimate (use j = 0, basically the first component)
        lam = z[0] / y[0]

        #Normalization of the iterated vector 
        norm = max(abs(zi) for zi in z)
        y_new = [zi / norm for zi in z]

        #Criterion for iterations termination
        f = abs(lam - lam_old)

        y = y_new
        lam_old = lam
        #Iteration Counter
        k += 1

    #The last eigenvalue estimate is the maximum in modulus eigenvalue
    w = lam_old

    print("Maximum in modulus eigenvalue w:")
    print("  ", f"{w:.7f}")

    print("Corresponding eigenvector y (normalized):")
    print("  y =", "\t".join(f"{yi:.7f}" for yi in y))

    #Check: Av - lambda*v (should be close to zero)
    print("Check A*y - lambda*y:")
    r = check_eigenpair(n, c, w, y)
    print("  ", "\t".join(f"{t:.7f}" for t in r))

    print("Number of iterations:", k)
    print("Given accuracy:", e)
    print("Final criterion value:", f)

    return w, y, k


print("Power Iteration Method:")
result = power_iteration_method(len(c), c, e)