def invert_letter(c):
    return {'a': 'A', 'b': 'B', 'A': 'a', 'B': 'b'}[c]

def free_reduce(word):
    stack = []
    for c in word:
        if stack and invert_letter(c) == stack[-1]:
            stack.pop()
        else:
            stack.append(c)
    return tuple(stack)

def invert_word(word):
    return tuple(invert_letter(c) for c in reversed(word))

def rotations(word):
    word = tuple(word)
    if len(word) == 0:
        return [word]
    return [word[i:] + word[:i] for i in range(len(word))]

def power_decomposition(word):
    """
    If word = u^k with k >= 2, return (u, k).
    Otherwise return None.
    """
    word = tuple(word)
    n = len(word)
    if n < 2:
        return None

    pi = [0] * n
    for i in range(1, n):
        j = pi[i - 1]
        while j > 0 and word[i] != word[j]:
            j = pi[j - 1]
        if word[i] == word[j]:
            j += 1
        pi[i] = j

    period = n - pi[-1]
    if period < n and n % period == 0:
        k = n // period
        if k >= 2:
            return word[:period], k

    return None

def one_cancel_insertions(word):
    """
    Yield the word itself, and every word obtained by inserting one
    cancelling pair c c^{-1} at every possible position.
    """
    word = tuple(word)
    yield word
    for pos in range(len(word) + 1):
        for c in ('a', 'b', 'A', 'B'):
            yield word[:pos] + (c, invert_letter(c)) + word[pos:]

def is_power_or_power_conjugate_product(word):
    """
    Returns True if the reduced word is either:

      1) a proper power u^k with k >= 2, or
      2) (up to cyclic permutation and one optional inserted cancelling pair)
         of the form u^k v u^k v^{-1}, with k >= 2.

    Returns:
        (True, "power", rotation, u, k)
        (True, "power_conjugate", rotation, u, k, v)
    or
        (False, None)
    """
    w = free_reduce(word)

    # First check: plain power words
    for rot in rotations(w):
        info = power_decomposition(rot)
        if info is not None:
            u, k = info
            return True, ("power", rot, u, k)

    # Second check: power-conjugate product words
    for rot in rotations(w):
        rot = tuple(rot)

        for expanded in one_cancel_insertions(rot):
            m = len(expanded)

            for v_len in range(0, m // 2 + 1):
                rem = m - 2 * v_len
                if rem < 0 or rem % 2:
                    continue

                p_len = rem // 2
                P = expanded[:p_len]
                info = power_decomposition(P)
                if info is None:
                    continue

                u, k = info
                v = expanded[p_len:p_len + v_len]

                if expanded == P + v + P + invert_word(v):
                    return True, ("power_conjugate", rot, u, k, v)

    return False, None

'''print(is_power_or_power_conjugate_product('aaa'))
print(is_power_or_power_conjugate_product('ababab'))
print(is_power_or_power_conjugate_product('aababbab'))
print(is_power_or_power_conjugate_product('abAB'))'''