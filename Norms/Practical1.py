# Matrix
A = [[7, 4, 2],[-8, -6, 5], [9, 0, -8]]

#Vector
b = [4, -5, 9]

# N.O.V
v1 = sum(abs(x) for x in b)
v2 = (sum(x**2 for x in b)) ** (1/2)
v3 = sum(abs(x)**3 for x in b) ** (1/3)
vc = max(abs(x) for x in b)

# N.O.M
m1 = max(sum(abs(A[i][j]) for i in range(3)) for j in range(3))
m2 = (sum(A[i][j]**2 for i in range(3) for j in range(3))) ** (1/2)
mc = max(sum(abs(A[i][j]) for j in range(3)) for i in range(3))

# A*b
Ab = [sum(A[i][j] * b[j] for j in range(3)) for i in range(3)]

# Norms of Ab
ab1 = sum(abs(x) for x in Ab)
ab2 = (sum(x**2 for x in Ab)) ** (1/2)
abc = max(abs(x) for x in Ab)

print("Vector norms:", round(v1, 4), round(v2, 4), round(v3, 4), round(vc, 4))
print("Matrix norms:", round(m1, 4), round(m2, 4), round(mc, 4))
print("Ab =", Ab)
print("Norms of Ab:", round(ab1, 4), round(ab2, 4), round(abc, 4))

print("Consistency:")
print("||Ab||_1 <= ||A||_1 * ||b||_1 :", ab1 <= m1*v1)
print("||Ab||_2 <= ||A||_2 * ||b||_2 :", ab2 <= m2*v2)
print("||Ab||_c <= ||A||_c * ||b||_c :", abc <= mc*vc)