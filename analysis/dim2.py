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

for w in uru(14):
    print(w)