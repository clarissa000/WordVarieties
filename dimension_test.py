from sage.all import *
from time import perf_counter
import random
R3 = PolynomialRing(QQ, ['x', 'y', 'z'])
x, y, z = R3.gens()

# Bigger ring for the ideal with u
R10 = PolynomialRing(QQ, ['x', 'y', 'z', 'u', 'v11', 'v12', 'v21', 'v22', 'v31', 'v32'])
x4, y4, z4, u, v11, v12, v21, v22, v31, v32 = R10.gens()


memo = {}

#A = a^-1
def invert_letter(c):
    return{'a':'A', 'b':'B', 'A':'a', 'B':'b'}[c]

#This removes any neighbouring letters
def free_reduce(word):
    stack = []
    for c in word:
        if stack and invert_letter(c) == stack[-1]:
            stack.pop()
        else:
            stack.append(c)
    return tuple(stack)

#This removes inverse letters from ends - as trace is conjugate invariant
def cyclic_reduce(word):
    word = free_reduce(word)
    while len(word) >= 2 and invert_letter(word[0]) == word[-1]:
        word = word[1:-1]
    return word

def invert_word(word):
    return tuple(invert_letter(c) for c in reversed(word))

#GOLDMAN METHOD
def choose_split(word):
    word = cyclic_reduce(word)
    n = len(word)

    positions = {}
    for i, c in enumerate(word):
        positions.setdefault(c, []).append(i)

    #words are cyclically invariant
    for s, pos in positions.items():
        if len(pos) >= 2:
            i, j = pos[0], pos[1]
            # rotate so the word starts just after j
            # and ends at j, giving p s q s
            p = word[j+1:] + word[:i]
            q = word[i+1:j]

            u1 = p + (s,)
            u2 = q + (s,)
            return u1, u2
    return None

def trace_poly(word):
    word = cyclic_reduce(word)
    if word in memo:
        return memo[word]

    #base cases
    n = len(word)
    if n == 0:
        return R3(2)
    if n == 1:
        return x if word[0] in ('a', 'A') else y
    if n == 2:
        c, d = word
        # same generator twice
        if c in ('a', 'A') and d in ('a', 'A'):
            return x**2 - 2
        if c in ('b', 'B') and d in ('b', 'B'):
            return y**2 - 2
        # mixed generator cases
        if (c, d) in [('a', 'b'), ('b', 'a'), ('A', 'B'), ('B', 'A')]:
            return z
        if (c, d) in [('a', 'B'), ('B', 'a'), ('A', 'b'), ('b', 'A')]:
            return x*y - z
    if n == 4 and len(set(word)) == 4:
        return x**2 + y**2 + z**2 - x*y*z - 2
    
    #otherwise there is a repeated letter
    split = choose_split(word)
    u1, u2 = split
    u3 = free_reduce(u1 + invert_word(u2))

    #using the identity that tr(AB) = tr(A)tr(B) - tr(AB^-1)
    ans = trace_poly(u1) * trace_poly(u2) - trace_poly(u3)
    memo[word] = ans
    return ans

def dimension_singular(word): 
    pw = trace_poly(word)
    paw = trace_poly('a' + word)
    pbw = trace_poly('b' + word)

    F = x**2 + y**2 + z**2 - x*y*z - 4
    I = ideal([pw - 2, paw - x, pbw - y, u*F - 1])
    Iwd = ideal([pw - 2, paw - x, pbw - y])

    Is = I._singular_()
    Gs = Is.groebner()

    Iwds = Iwd._singular_()
    Gwds = Iwds.groebner()

    print('Groebner Basis', Gs)
    print('Dimension' , singular.eval('dim(%s)' % Gs.name()))

    #Not sure what to do here as still has u in the ring - I think just subtract 1? 
    print('w D Groebner Basis', Gwds)
    print('w D Dimension' , singular.eval('dim(%s)' % Gwds.name()))

def dimension_sage(word): 
    pw = trace_poly(word)
    paw = trace_poly('a' + word)
    pbw = trace_poly('b' + word)

    F = x**2 + y**2 + z**2 - x*y*z - 4
    I = ideal([pw - 2, paw - x, pbw - y, u*F - 1])
    Iwd = ideal([pw - 2, paw - x, pbw - y])

    Gs = I.groebner_basis()

    print('Groebner Basis', I.groebner_basis())
    print('Dimension' , I.dimension())

    #Not sure what to do here as still has u in the ring - I think just subtract 1? 
    print('w D Groebner Basis', Iwd.groebner_basis())
    print('w D Dimension' , Iwd.dimension() - 1)

def dimension_singular(word): 
    pw = trace_poly(word)
    paw = trace_poly('a' + word)
    pbw = trace_poly('b' + word)

    F = x**2 + y**2 + z**2 - x*y*z - 4
    I = ideal([pw - 2, paw - x, pbw - y, u*F - 1])
    Iwd = ideal([pw - 2, paw - x, pbw - y])

    Is = I._singular_()
    Gs = Is.groebner()

    Iwds = Iwd._singular_()
    Gwds = Iwds.groebner()

    print('Groebner Basis', Gs)
    print('Dimension' , singular.eval('dim(%s)' % Gs.name()))

    #Not sure what to do here as still has u in the ring - I think just subtract 1? 
    print('w D Groebner Basis', Gwds)
    print('w D Dimension' , int(singular.eval('dim(%s)' % Gwds.name()))-1 )

def additional_variables(word):
    pw  = trace_poly(word)
    paw = trace_poly('a'+ word)
    pbw = trace_poly('b' + word)

    F = x**2 + y**2 + z**2 - x*y*z - 4

    # Ideal in R10 = R3[u, v11, ...]
    I10 = ideal([R10(pw) - 2, R10(paw) - x4, R10(pbw) - y4, u*R10(F) - 1, v11*x4 + v12*z4 - 1, v21*x4 + v22*y4 - 1, v31*y4 + v32*z4 - 1])
    return I10.groebner_basis(), I10.dimension()

def saturation_method(word):
    pw  = trace_poly(word)
    paw = trace_poly('a'+ word)
    pbw = trace_poly('b' + word)

    F = x**2 + y**2 + z**2 - x*y*z - 4

    # Ideal in R3
    I3 = ideal([pw - 2, paw - x, pbw - y])
    #saturation method?
    ideal_F = ideal([F])
    ideal_axes = ideal([x*y, x*z, y*z])

    I3_cleaned = I3.saturation(ideal_F)[0]
    I3_cleaned = I3_cleaned.saturation(ideal_axes)[0]
    return I3_cleaned.groebner_basis(), I3_cleaned.dimension()

def compare_one_word(word, f1, f2, repeats=20):
    #to avoid calculation of trace polynomial contributing
    pw = trace_poly(word)
    paw = trace_poly('a' + word)
    pbw = trace_poly('b' + word)

    start = perf_counter()
    for _ in range(repeats):
        f1(word)
    val1 = perf_counter() - start

    start = perf_counter()
    for _ in range(repeats):
        f2(word)
    val2 = perf_counter() - start

    print(f"Length {len(word)} Word {word}")
    print(f"Function 1: {val1:.6f}s")
    print(f"Function 2: {val2:.6f}s")
    return val1 > val2

def random_word(length):
    alphabet = ['a', 'b', 'A', 'B']
    return ''.join(random.choice(alphabet) for _ in range(length))

sum = 0
count = 0
for i in range(1, 8):
    for _ in range(2):
        count += 1
        sum = sum + 1 if compare_one_word(random_word(5*i), additional_variables, additional_variables) else sum
print(sum / count)