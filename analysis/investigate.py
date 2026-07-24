from sage.all import *

memo = {}

R3 = PolynomialRing(QQ, ['x', 'y', 'z'])
x, y, z = R3.gens()

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

alphabet = ['a', 'b', 'A', 'B']
inverse = {'a': 'A', 'A': 'a', 'b': 'B', 'B': 'b'}

def reduced_words(n):
    if n == 0:
        yield ""
        return

    def extend(prefix):
        if len(prefix) == n:
            yield prefix
            return
        for letter in alphabet:
            if not prefix or letter != inverse[prefix[-1]]:
                yield from extend(prefix + letter)

    yield from extend("")

def canonical_reduced_words(n):
    if n == 0:
        yield ""
        return

    def extend(prefix):
        if len(prefix) == n:
            w = tuple(prefix)
            canon = canonical_word(w, use_swap=True)
            if w == canon:
                yield ''.join(canon)
            return
        for letter in alphabet:
            if not prefix or letter != inverse[prefix[-1]]:
                prefix.append(letter)
                yield from extend(prefix)
                prefix.pop()

    yield from extend([])

def uru(n):
    if n%2 != 0:
        return []
    for w in reduced_words(n//2):
        neww = w + ''.join(inverse[c] for c in w)
        if tuple(neww) == canonical_word(neww, use_swap=True):
            yield neww
    return

F = x**2 + y**2 + z**2 - x*y*z - 4
axes = ideal([x*y, x*z, y*z])

def bar_letterwise(u):
    return ''.join(inverse[c] for c in u)


def saturated_ideal_for_word(w):
    Pw  = trace_poly(w)
    Paw = trace_poly('a' + w)
    Pbw = trace_poly('b' + w)

    I = ideal([Pw - 2, Paw - x, Pbw - y])
    I = I.saturation(ideal([F]))[0]
    I = I.saturation(axes)[0]
    return I

def check_family(max_len=14):
    bad = []
    for n in range(1, max_len + 1):
        for u in uru(n):
            bu = bar_letterwise(u)
            w = u + bu
            I = saturated_ideal_for_word(w)
            d = I.dimension()

            if d != 2:
                bad.append((u, w, d, I.groebner_basis()))
                print("BAD:", u, "->", w, "dimension =", d)
    return bad

bad = check_family()
print("Number of bad words:", len(bad))