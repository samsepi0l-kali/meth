n = 5 #Matrix Dimension

#Coefficients Here
a = [0, 2, -8, 6, 9]
b = [-6, 10, 18, -17, 14]
c = [6, -7, 9, -6, 0]
d = [30, -31, 108, -114, 124]

p = [0] * n
q = [0] * n
det = b[0] #Determinant Accumulator here 

p[0] = -c[0] / b[0]
q[0] = d[0] / b[0]

#Forward Path
for i in range(1, n):
    y = b[i] + a[i] * p[i - 1] #Denominator
    p[i] = -c[i] / y
    q[i] = (d[i] - a[i] * q[i - 1]) / y
    det *= y 

x = [0] * n
x[-1] = q[-1]

#Back Substitution
for i in range(n-2, -1, -1):
    x[i] = p[i] * x[i + 1] + q[i]

print("x:", x)
print("Determinant:", det)