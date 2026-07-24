from sage.all import *

R3 = PolynomialRing(QQ, ['x', 'y', 'z'])
x, y, z = R3.gens()
#g basis for abaBAbaBabAB
basis_strings = ["x*y*z^4 + x^2*y^2*z - x^2*z^3 - y^2*z^3 - z^5 - x^3*y - x*y^3 - 4*x*y*z^2 + 3*x^2*z + 3*y^2*z + 6*z^3 + 3*x*y - 9*z", "x^3*y*z + x*y*z^3 - x^4 - x^2*y^2 - 2*x^2*z^2 - y^2*z^2 - z^4 - 3*x*y*z + 6*x^2 + 3*y^2 + 6*z^2 - 9", "x^2*y*z^2 - x^3*z - 2*x*y^2*z - x*z^3 + x^2*y + y^3 + y*z^2 + 3*x*z - 3*y"]

basis = [R3(s.replace("^", "**")) for s in basis_strings]

for f in basis: 
    print(factor(f))

g = basis[0]
for f in basis[1:]:
    g = gcd(g,f)

print('gcd ', factor(g))
I = ideal(basis)
print(I.dimension())

Q = x**2 + y**2 + z**2 - x*y*z - 3

quotients = []
remainders = []

for f in basis:
    q, r = f.quo_rem(Q)
    quotients.append(q)
    remainders.append(r)

print("Remainders:")
print(remainders)

print("Factored quotients:")
for q in quotients:
    print(factor(q))