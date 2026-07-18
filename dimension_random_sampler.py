SAMPLE_SIZE = 20
LENGTH = 20

from sage import *

# Base ring for trace polynomials
R3 = PolynomialRing(QQ, ['x', 'y', 'z'])
x, y, z = R3.gens()

# Bigger ring for the ideal with u
#R10 = PolynomialRing(QQ, ['x', 'y', 'z', 'u', 'v11', 'v12', 'v21', 'v22', 'v31', 'v32'])
#x4, y4, z4, u, v11, v12, v21, v22, v31, v32 = R10.gens()

memo = {}

#A = a^-1
def invert_letter(c):
    return{'a':'A', 'b':'B', 'A':'a', 'B':'b'}[c]

def swap_ab_letter(c):
    return{'a':'b', 'b':'a', 'A':'B', 'B':'A'}[c]

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

def swap_ab_word(word):
    return tuple(swap_ab_letter(c) for c in word)

def swap_inverse_word(word):
    return tuple(invert_letter(c) for c in word)

def rotations(word):
    word = tuple(word)
    if len(word) == 0:
        return [word]
    return [word[i:] + word[:i] for i in range(len(word))]

def word_key(word):
    order = {'a': 0, 'b': 1, 'A': 2, 'B': 3}
    return tuple(order[c] for c in word)

def canonical_word(word, use_swap=False):
    """
    Return a 'canonical' representative for the word up to:
      - free reduction
      - cyclic permutation
      - inversion
      - optionally swapping a <-> b

    Set use_swap=True if you want to identify words up to renaming a and b. - Different variety but only with naming of coordinates
    """
    word = cyclic_reduce(word)

    candidates = []

    variants = [word, invert_word(word)]
    if use_swap:
        variants += [swap_ab_word(word), invert_word(swap_ab_word(word))]
        variants += [swap_inverse_word(word), invert_word(swap_inverse_word(word))]

    for v in variants:
        for r in rotations(v):
            candidates.append(r)

    return min(candidates, key=word_key)

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

def dimension_two_rings(word):
    pw  = trace_poly(word)
    paw = trace_poly('a'+ word)
    pbw = trace_poly('b' + word)

    F = x**2 + y**2 + z**2 - x*y*z - 4

    # Ideal in R3
    I3 = ideal([pw - 2, paw - x, pbw - y])

    #print('\n\nWord', word)
    #print('Groebner basis for I3:')
    #print(I3.groebner_basis())
    #print('Dimension of I3:')
    #print(I3.dimension())

    #saturation method?
    ideal_F = ideal([F])
    ideal_axes = ideal([x*y, x*z, y*z])

    I3_cleaned = I3.saturation(ideal_F)[0]
    I3_cleaned = I3_cleaned.saturation(ideal_axes)[0]
    basis = I3_cleaned.groebner_basis()
    dim = I3_cleaned.dimension()
    if dim == 0:
        dimv = I3_cleaned.vector_space_dimension()
        dimr = I3_cleaned.radical().vector_space_dimension()
    else:
        dimv, dimr = -1, -1
    #print('\nGroebner basis for I3 cleaned:')
    #print(I3_cleaned.groebner_basis())
    #print('Dimension of I3 cleaned:')
    #print(I3_cleaned.dimension())
    return I3.dimension(), dim, basis, dimv, dimr

alphabet = ['a', 'b', 'A', 'B']
inverse = {'a': 'A', 'A': 'a', 'b': 'B', 'B': 'b'}

def word_generator(n):
    #Generates a freely reduced word of length n
    pass

new_rows = []
for _ in range(SAMPLE_SIZE):
    w = word_generator(LENGTH)
    dim, dim2, basis, dimv, dimr = dimension_two_rings(w)
    new_rows.append({"length":LENGTH, "word":w, "dimension":dim, "dimensionZDense":dim2, "gb basis":basis, "cardinality":dimv, "radical_cardinality":dimr})
        
if new_rows <= 20:
    print(new_rows)