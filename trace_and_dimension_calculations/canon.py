alphabet = {'a', 'b', 'A','B'}
inverse = {'a': 'A', 'A': 'a', 'b': 'B', 'B': 'b'}
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

def reverse_word(word):
    return tuple(reversed(word))

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
      - optionally swapping a <-> b and reversing

    Set use_swap=True if you want to identify words up to renaming a and b. - Different variety but only with naming of coordinates
    """
def canonical_word(word, use_swap=False):
    word = cyclic_reduce(word)

    variants = {word, invert_word(word)}

    if use_swap:
        swapped = swap_ab_word(word)
        reversed_w = reverse_word(word)

        variants.update({
            swapped,
            invert_word(swapped),
            reversed_w,
            invert_word(reversed_w),
            swap_ab_word(reversed_w),
            invert_word(swap_ab_word(reversed_w)),
        })

    candidates = []
    for v in variants:
        candidates.extend(rotations(cyclic_reduce(v)))

    return min(candidates, key=word_key)

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

'''count = 0
for w in canonical_reduced_words(15):
    count += 1
print(count)'''

#print(''.join(canonical_word('aaabbabbbaBab', use_swap=True)))

print(''.join(canonical_word('abABabABABabABab')))