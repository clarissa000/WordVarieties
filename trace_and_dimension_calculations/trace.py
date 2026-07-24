from sage.all import *
from time import perf_counter
import random

R = PolynomialRing(QQ, ['x', 'y', 'z'])
x,y,z = R.gens()

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

#----------------------
#GOLDMAN METHOD
#----------------------
def choose_split(word):
    word = cyclic_reduce(word)
    n = len(word)

    positions = {}
    for i, c in enumerate(word):
        positions.setdefault(c, []).append(i)

    #words are cyclically invariant
    for s, pos in positions.items():
        if len(pos) >= 2:
            i, j = pos[0], pos[len(pos)//2]

            # rotate so the word starts just after j
            # and ends at j, giving p s q s
            p = word[j+1:] + word[:i]
            q = word[i+1:j]

            u1 = p + (s,)
            u2 = q + (s,)
            return u1, u2
    return None

def trace_poly_goldman(word, memo = None):
    word = cyclic_reduce(word)
    if word in memo:
        return memo[word]

    #base cases
    n = len(word)
    if n == 0:
        return R(2)
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
    ans = trace_poly_goldman(u1, memo) * trace_poly_goldman(u2, memo) - trace_poly_goldman(u3, memo)
    memo[word] = ans
    return ans

def trace_poly_goldman_fast(word):
    return trace_poly_goldman(word, {})
#----------------------
#TRAINA METHOD
#----------------------
# Chebyshev recurrences
# T_0 = 2, T_1 = t, T_n = t*T_{n-1} - T_{n-2}
def T(n, t):
    if not hasattr(T, "_cache"):
        T._cache = {
            x: {0: R(2), 1: x},
            y: {0: R(2), 1: y},
        }

    n = int(n)
    cache = T._cache[t]

    if n < 0:
        n = -n

    if n in cache:
        return cache[n]

    cache[n] = t * T(n - 1, t) - T(n - 2, t)
    return cache[n]

# P_{-1} = 0, P_0 = 1, P_n = t*P_{n-1} - P_{n-2}
def P(n, t):
    if not hasattr(P, "_cache"):
        P._cache = {
            x: {-1: R(0), 0: R(1)},
            y: {-1: R(0), 0: R(1)},
        }

    n = int(n)
    cache = P._cache[t]

    if n in cache:
        return cache[n]

    cache[n] = t * P(n - 1, t) - P(n - 2, t)
    return cache[n]

def trailing_block(word):
    """
    Return (prefix, last_letter, run_length)
    where last_letter is the final repeated letter.
    """
    word = cyclic_reduce(word)
    if not word:
        return (), None, 0

    c = word[-1]
    k = 1
    i = len(word) - 2

    while i >= 0 and word[i] == c:
        k += 1
        i -= 1

    prefix = word[:i + 1]
    return prefix, c, k

def trace_poly_syllable(word, memo=None):
    if memo is None:
        memo = {}

    word = cyclic_reduce(word)

    if word in memo:
        return memo[word]

    n = len(word)

    # base cases
    if n == 0:
        memo[word] = R(2)
        return memo[word]

    if n == 1:
        ans = x if word[0] in ('a', 'A') else y
        memo[word] = ans
        return ans

    # pure powers
    if len(set(word)) == 1:
        c = word[0]
        t = x if c in ('a', 'A') else y
        ans = T(n, t)
        memo[word] = ans
        return ans

    # small cases
    if n == 2:
        c, d = word

        if c in ('a', 'A') and d in ('a', 'A'):
            ans = x**2 - 2
            memo[word] = ans
            return ans
        if c in ('b', 'B') and d in ('b', 'B'):
            ans = y**2 - 2
            memo[word] = ans
            return ans

        if (c, d) in [('a', 'b'), ('b', 'a'), ('A', 'B'), ('B', 'A')]:
            ans = z
            memo[word] = ans
            return ans

        if (c, d) in [('a', 'B'), ('B', 'a'), ('A', 'b'), ('b', 'A')]:
            ans = x*y - z
            memo[word] = ans
            return ans

    if n == 4 and len(set(word)) == 4:
        ans = x**2 + y**2 + z**2 - x*y*z - 2
        memo[word] = ans
        return ans

    # strip a final block c^k if possible
    prefix, c, k = trailing_block(word)

    if k >= 2:
        t = x if c in ('a', 'A') else y
        ans = P(k - 1, t) * trace_poly_syllable(prefix + (c,), memo) - P(k - 2, t) * trace_poly_syllable(prefix, memo)
        memo[word] = ans
        return ans

    # otherwise fall back to Goldman split
    split = choose_split(word)

    u1, u2 = split
    u3 = free_reduce(u1 + invert_word(u2))

    ans = trace_poly_syllable(u1, memo) * trace_poly_syllable(u2, memo) - trace_poly_syllable(u3, memo)
    memo[word] = ans
    return ans

def trace_poly_syllable_fast(word):
    return trace_poly_syllable(word, {})

#----------------------
#CAYLEY-HAMILTON DIRECT METHOD
#----------------------
def letter_var(c):
    return x if c.lower() == 'a' else y

def find_ch_step(word):
    """
    Find the first place where a Cayley-Hamilton rewrite applies.
    """
    word = cyclic_reduce(word)

    # 1) adjacent equal letters
    for i in range(len(word) - 1):
        if word[i] == word[i + 1]:
            c = word[i]
            t = letter_var(c)

            # U c c V = t * (U c V) - (U V)
            keep_one = word[:i] + (c,) + word[i + 2:]
            drop_both = word[:i] + word[i + 2:]
            return ("pair", t, keep_one, drop_both)

    # 2) inverses
    for i, c in enumerate(word):
        if c in ('A', 'B'):
            lower = c.lower()
            t = letter_var(c)

            # U C V = t * (U V) - (U lower V)
            drop_capital = word[:i] + word[i + 1:]
            replace_lower = word[:i] + (lower,) + word[i + 1:]
            return ("capital", t, drop_capital, replace_lower)

    return None

def trace_poly_CH(word, memo=None):
    """
    Cayley-Hamilton style trace polynomial.

    Strategy:
      - free/cyclic reduce
      - use CH rewrites on repeated letters or capitals
      - if no CH step applies, fall back to Goldman split
    """
    if memo is None:
        memo = {}

    word = cyclic_reduce(word)

    if word in memo:
        return memo[word]

    n = len(word)

    # base cases
    if n == 0:
        memo[word] = R(2)
        return memo[word]

    if n == 1:
        ans = x if word[0] in ('a', 'A') else y
        memo[word] = ans
        return ans

    if n == 2:
        c, d = word

        if c in ('a', 'A') and d in ('a', 'A'):
            ans = x**2 - 2
            memo[word] = ans
            return ans
        if c in ('b', 'B') and d in ('b', 'B'):
            ans = y**2 - 2
            memo[word] = ans
            return ans

        if (c, d) in [('a', 'b'), ('b', 'a'), ('A', 'B'), ('B', 'A')]:
            ans = z
            memo[word] = ans
            return ans

        if (c, d) in [('a', 'B'), ('B', 'a'), ('A', 'b'), ('b', 'A')]:
            ans = x*y - z
            memo[word] = ans
            return ans

    if n == 4 and len(set(word)) == 4:
        ans = x**2 + y**2 + z**2 - x*y*z - 2
        memo[word] = ans
        return ans

    # pure power shortcut
    if len(set(word)) == 1:
        c = word[0]
        t = x if c.lower() == 'a' else y
        ans = T(n, t)
        memo[word] = ans
        return ans

    # CH rewrite step
    step = find_ch_step(word)
    if step is not None:
        kind, t, w1, w2 = step
        ans = t * trace_poly_CH(w1, memo) - trace_poly_CH(w2, memo)
        memo[word] = ans
        return ans

    # fallback to Goldman split
    split = choose_split(word)
    if split is None:
        raise ValueError(f"No split found for word {word}")

    u1, u2 = split
    u3 = free_reduce(u1 + invert_word(u2))

    ans = trace_poly_CH(u1, memo) * trace_poly_CH(u2, memo) - trace_poly_CH(u3, memo)
    memo[word] = ans
    return ans

def trace_poly_CH_fast(word):
    return trace_poly_CH(word, {})

#----------------------
#Method 4 Matrix method
#----------------------
Pxyz = PolynomialRing(ZZ, names=('x', 'y', 'z'))
x3, y3, z3 = Pxyz.gens()

# Quadratic extension ring:
#   R = Z[x,y,z,zeta]/(zeta^2 - z*zeta + 1)
S = PolynomialRing(ZZ, names=('x', 'y', 'z', 'u'))
xS, yS, zS, uS = S.gens()
Rquad = S.quotient(uS**2 - zS*uS + 1, names=('x', 'y', 'z', 'u'))
xq, yq, zq, uq = Rquad.gens()

# Matrices from the paper
Aq = Matrix(Rquad, [[xq, -1], [1, 0]])
Bq = Matrix(Rquad, [[0, uq], [uq - zq, yq]])

Aq_inv = Matrix(Rquad, [[0, 1], [-1, xq]])
Bq_inv = Matrix(Rquad, [[yq, -uq], [zq - uq, 0]])

GEN_MATS = {
    'a': Aq,
    'A': Aq_inv,
    'b': Bq,
    'B': Bq_inv,
}

def to_tuple_word(word):
    if isinstance(word, str):
        return tuple(word)
    return tuple(word)

def trace_poly_generic(word, reduce_first=True):
    word = to_tuple_word(word)

    if reduce_first:
        word = free_reduce(word)

    M = identity_matrix(Rquad, 2)
    for c in word:
        M = M * GEN_MATS[c]

    tr = M.trace()

    # lift from the quotient ring back to S
    tr_lift = tr.lift()

    # evaluate u = 0 to get the polynomial in x,y,z
    eval_u0 = S.hom([xS, yS, zS, 0], Pxyz)
    return eval_u0(tr_lift)

def trace_poly_generic_fast(word):
    return trace_poly_generic(word, reduce_first=True)

#----------------------
# COMPARISON / BENCHMARKING
#----------------------
def compare_one_word(word, repeats=100):

    start = perf_counter()
    for _ in range(repeats):
        trace_poly_goldman_fast(word)
    gold = perf_counter() - start

    start = perf_counter()
    for _ in range(repeats):
        trace_poly_syllable_fast(word)
    syl = perf_counter() - start

    '''start = perf_counter()
    for _ in range(repeats):
        trace_poly_CH_fast(word)
    ch = perf_counter() - start

    start = perf_counter()
    for _ in range(repeats):
        trace_poly_generic_fast(word)
    gen = perf_counter() - start'''

    print(f"Length {len(word)} Word {word}")
    print(f"Goldman : {gold:.6f}s")
    print(f"Syllable: {syl:.6f}s")
    #print(f"Cayley-Hamilton: {ch:.6f}s")
    #print(f"Matrix: {gen:.6f}s")
    #print(f"Speedup G/S: {gold/syl:.2f}x")
    #print(f"Speedup C/S: {ch/syl:.2f}x")
    #print('\n')
    return True if gold < syl else False

def random_word(length):
    alphabet = ['a', 'b', 'A', 'B']
    return ''.join(random.choice(alphabet) for _ in range(length))

sum = 0
count = 0
for i in range(1, 10):
    for _ in range(10):
        count += 1
        sum = sum + 1 if compare_one_word(random_word(5*i)) else sum
print(sum / count)