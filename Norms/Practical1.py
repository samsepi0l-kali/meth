import math

n = 3
p = 3

# Variant 1 matrix
A = [
    [3, -7, 1],
    [9, 8, -7],
    [5, -2, 2]
]

# Variant 1 vector
b = [1, 5, -3]

# Vector p-norms
x = []

for j in range(1, p + 1):
    s = 0

    for i in range(n):
        h = abs(b[i])
        s += h ** j

    x.append(s ** (1.0 / j))

# Vector infinity norm
y = 0

for i in range(n):
    y = max(y, abs(b[i]))

# Matrix 1-norm:
# maximum absolute column sum
u = 0

for j in range(n):
    h = 0

    for i in range(n):
        h += abs(A[i][j])

    u = max(u, h)

# Matrix 2-norm according to the lecture program:
# square root of the sum of squares
v = 0

for i in range(n):
    for j in range(n):
        v += A[i][j] ** 2

v = math.sqrt(v)

# Matrix infinity norm:
# maximum absolute row sum
w = 0

for i in range(n):
    h = 0

    for j in range(n):
        h += abs(A[i][j])

    w = max(w, h)

# Output
print("Vector norms:")
print(f"||b||1   = {x[0]:.6f}")
print(f"||b||2   = {x[1]:.6f}")
print(f"||b||3   = {x[2]:.6f}")
print(f"||b||inf = {y:.6f}")

print()

print("Matrix norms:")
print(f"||A||1   = {u:.6f}")
print(f"||A||2   = {v:.6f}")
print(f"||A||inf = {w:.6f}")