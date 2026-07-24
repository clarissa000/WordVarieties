from sage.all import FreeGroup
from tqdm import tqdm

# --- word utilities on strings like "ababAbABABaB" ---

def word_to_free_group(F, word):
    """
    Convert a string word in {a,b,A,B} to an element of the free group F.
    """
    a, b = F.gens()
    out = F.one()
    for c in word:
        if c == 'a':
            out *= a
        elif c == 'b':
            out *= b
        elif c == 'A':
            out *= a**-1
        elif c == 'B':
            out *= b**-1
        else:
            raise ValueError(f"Bad letter: {c}")
    return out

def free_reduce_str(word):
    inv = {'a':'A', 'A':'a', 'b':'B', 'B':'b'}
    stack = []
    for c in word:
        if stack and inv[c] == stack[-1]:
            stack.pop()
        else:
            stack.append(c)
    return ''.join(stack)

def is_power_word_str(word):
    """
    Check whether a reduced string word is a proper power.
    Returns (True, base, k) or (False, None, None).
    """
    w = free_reduce_str(word)
    n = len(w)
    for d in range(1, n):
        if n % d != 0:
            continue
        base = w[:d]
        if base * (n // d) == w:
            return True, base, n // d
    return False, None, None

def vanishes_in_normal_closure(target_word, relator_word):
    """
    Return True iff target_word becomes trivial in the quotient by the
    normal closure of relator_word.

    Both inputs are strings in {a,b,A,B}.
    """
    F = FreeGroup("a, b")
    a, b = F.gens()

    target = word_to_free_group(F, free_reduce_str(target_word))
    relator = word_to_free_group(F, free_reduce_str(relator_word))

    Q = F / [relator]   # quotient by the normal closure of relator
    qtarget = Q(target)

    return qtarget == Q.one()

alphabet = ['a', 'b', 'A', 'B']

def inverse_letter(c):
    return {'a':'A', 'A':'a', 'b':'B', 'B':'b'}[c]

def reduced_words_of_length(n):
    if n == 0:
        yield ""
        return

    def extend(prefix):
        if len(prefix) == n:
            yield ''.join(prefix)
            return
        for c in alphabet:
            if not prefix or inverse_letter(c) != prefix[-1]:
                yield from extend(prefix + [c])

    yield from extend([])

def proper_power_candidates(min_len, max_len):
    """
    Generate reduced proper power words up to length max_len.
    Deduplicate by cyclic rotation and inversion only crudely via a set.
    """
    seen = set()
    for n in range(min_len, max_len + 1):
        for w in reduced_words_of_length(n):
            ok, base, k = is_power_word_str(w)
            if ok:
                canon = min(w[i:] + w[:i] for i in range(len(w)))
                inv = ''.join(inverse_letter(c) for c in reversed(w))
                canon = min(canon, min(inv[i:] + inv[:i] for i in range(len(inv))))
                if canon not in seen:
                    seen.add(canon)
                    yield w, base, k

def find_power_relators_for_target(target_word, min_len = 2, max_len=8):
    hits = []
    for relator, base, k in proper_power_candidates(min_len, max_len):
        if vanishes_in_normal_closure(target_word, relator):
            hits.append((relator, base, k))
    return hits

#hits = find_power_relators_for_target("ababABAB", max_len=8)
#print(hits)

'''words = 'aabAAB aaabAAAB aabbAABB ababABAB abAbABaB aaaaaabaaB aaaaaabAAB aaaabAAAAB aaabbAAABB aababAABAB aabaBAAbAB aabaBAABAb aabAbAABaB aaBabAAbAB aaaaaabaaaB aaaaaabAAAB aaaaaabbaaBB aaaaaabbAABB aaaaabAAAAAB aaaabaabaaBB aaaabaabAABB aaaabaaBBaab aaaabaaBBAAb aaaababaaBaB aaaababAABaB aaaabaBaabaB aaaabaBAAbaB aaaabbaaBaaB aaaabbaaBAAB aaaabbAAAABB aaaabbAABaaB aaaabbAABAAB aaaabAbaaBAB aaaabAbAABAB aaaabAAbaaBB aaaabAAbAABB aaaabAABBAAb aaaabABaabAB aaaabABAAbAB aaaaBaabbaaB aaaaBaabbAAB aaaaBAAbbAAB aaabaaBabaaB aaabaaBabAAB aaabaaBaBaab aaabaaBaBAAb aaababAAABAB aaabaBAAAbAB aaabaBAAABAb aaabbbAAABBB aaabAbAAABaB aaabAABabaaB aaabAABabAAB aaabAABaBAAb aaaBaababaaB aaaBaababAAB aaaBabAAAbAB aaaBAAbabAAB aabaabaaBaaB aabaabAABaaB aabaabAABAAB aabaaBAAbAAB aababbAABABB aababbAABBAB aababbABAABB aababAbABABB aababABBABAb aabaBBAAbbAB aabaBBAABAbb aabbaBAAbABB aabbaBAABBAb aabbAbAABaBB aabbAbAABBaB aabbAbABABaB aabbAAbaBBAB aabbAAbABBaB aabbAABaBBAb aabAbbAABaBB aabAAbaaBAAB aabAAbAABaaB aabABabbaBAb aaBaBAbbAbAB aaBAbaBBabAB aaBAbbAAbaBB abababABABAB ababAbABABaB abaBabABAbAB abaBAbaBabAB abAbAbABaBaB'.split(' ')
upto8 = []
for word in words:
    if find_power_relators_for_target(word, max_len=8) == []:
        print(word)
        upto8.append(word)

print(upto8)
'''

'''upto8 = ['aababbAABABB', 'aabaBBAABAbb', 'aabbaBAABBAb', 'aabbAbAABBaB', 'aabAbbAABaBB', 'aabABabbaBAb', 'aaBAbaBBabAB', 'aaBAbbAAbaBB', 'ababAbABABaB', 'abaBabABAbAB', 'abaBAbaBabAB']
upto10 = [] 
for word in upto8:
    print('n')
    if find_power_relators_for_target(word, max_len=10) == []:
        print(word)
        upto10.append(word)
print(upto10)
print(len(upto8), len(upto10))
'''
#print(find_power_relators_for_target('aababbAABABB',min_len=1, max_len=10))

upto5 = ['aabaBAABAb', 'aaBabAAbAB', 'aaabaBAAABAb', 'aaaBabAAAbAB', 'aababbAABABB', 'aabaBBAABAbb', 'aabbaBAABBAb', 'aabbAbAABBaB', 'aabAbbAABaBB', 'aabABabbaBAb', 'aaBAbaBBabAB', 'aaBAbbAAbaBB', 'abababABABAB', 'ababAbABABaB', 'abaBabABAbAB', 'abaBAbaBabAB', 'abAbAbABaBaB', 'aaaabaBAAAABAb', 'aaaaBabAAAAbAB', 'aaababbAAABABB', 'aaabaBBAAABAbb', 'aaabbabAAABBAB', 'aaabbaBAAABBAb', 'aaabbAbAAABBaB', 'aaabAbbAAABaBB', 'aaabABababaBAb', 'aaaBabbAAAbABB', 'aaaBAbaBaBabAB', 'aaaBAbbAAAbaBB', 'aabaabbABAABAB', 'aababaBAABABAb', 'aabaBaBAABAbAb', 'aabAbaBAABaBAb', 'aabAbbAbABBaBB', 'aaBababAAbABAB', 'aaBaBabAAbAbAB', 'aaBAbabAAbaBAB', 'aaaaaaabAAAAAAAB', 'aaaaabaaababaaab', 'aaaaabaBAAAAABAb', 'aaaaaBaaaBaBaaaB', 'aaaaaBabAAAAAbAB', 'aaaababbAAAABABB', 'aaaabaBBAAAABAbb', 'aaaabbaabbbbaabb', 'aaaabbabAAAABBAB', 'aaaabbaBAAAABBAb', 'aaaabbAbAAAABBaB', 'aaaabAbbAAAABaBB', 'aaaabABabaabaBAb', 'aaaaBabbAAAAbABB', 'aaaaBAbaBaaBabAB', 'aaaaBAbbAAAAbaBB', 'aaaaBBaaBBBBaaBB', 'aaabaabbAAABAABB', 'aaabaabAbAABAABB', 'aaabaabAbABAABAB', 'aaabaabABABAABAb', 'aaabaabABBAABAAb', 'aaabaaBBAAABAAbb', 'aaababaBAAABABAb', 'aaababbbAAABABBB', 'aaababAbAAABABaB', 'aaababABAAABABab', 'aaabaBabaaaBabaB', 'aaabaBaBAAABAbAb', 'aaabaBAbAAABAbaB', 'aaabaBBBAAABAbbb', 'aaabbaabAAABBAAB', 'aaabbaaBAAABBAAb', 'aaabbabbAABABBAB', 'aaabbbaBAAABBBAb', 'aaabbbAbAAABBBaB', 'aaabbAAbAAABBaaB', 'aaabbAAbAABABaaB', 'aaabAbabAAABaBAB', 'aaabAbaBAAABaBAb', 'aaabAbbbAAABaBBB', 'aaabAbbAbAABBaBB', 'aaabAbAAbABABaaB', 'aaabAAbbAAABaaBB', 'aaabABabAAABabAB', 'aaabABAbaaaBAbAB', 'aaaBaabbAAAbAABB', 'aaaBaaBAbbAAbAAB', 'aaaBaaBAbAbAAbAB', 'aaaBababAAAbABAB', 'aaaBabAbAAAbABaB', 'aaaBaBabAAAbAbAB', 'aaaBAbabAAAbaBAB', 'aaaBAbbbAAAbaBBB', 'aaaBAAbbAAAbaaBB', 'aabaababAABAABAB', 'aabaabaBAABAABAb', 'aabaabAbAABAABaB', 'aabaaBabaaBaabaB', 'aabaaBabAABAAbAB', 'aabaaBAbaaBaabAB', 'aababaabAABABAAB', 'aabababbAABABABB', 'aabababAbABABABB', 'aabababABBABABAb', 'aababaBBAABABAbb', 'aababAAbAABABaaB', 'aabaBaabAABAbAAB', 'aabaBaaBAABAbAAb', 'aabaBabbAABAbABB', 'aabaBAbbAABAbaBB', 'aabaBAAbaaBabAAB', 'aabaBAAbAABAbaaB', 'aabaBABBAABAbabb', 'aabaBBaBAABAbbAb', 'aabbaabABBAABBAb', 'aabbaBabAABBAbAB', 'aabbaBaBAABBAbAb', 'aabbaBAbabaBabAB', 'aabbaBAbAABBAbaB', 'aabbAbabAABBaBAB', 'aabbAbAbAABBaBaB', 'aabbAbAbABABaBaB', 'aabbABabAbabaBAb', 'aabbABabAABBabAB', 'aabAbaabAABaBAAB', 'aabAbaBAbbaBabAB', 'aabAbAbbAABaBaBB', 'aabAbAAbAABaBaaB', 'aabAAbabAABaaBAB', 'aabAAbaBAABaaBAb', 'aabAAbAbAABaaBaB', 'aabAABabAABaabAB', 'aabABabaBaBAbaBB', 'aabABabaBBAbaBAB', 'aabABabbAABabABB', 'aabABabABBAbaBAb', 'aabABAbbAABabaBB', 'aaBaaBabAAbAAbAB', 'aaBabAbbAAbABaBB', 'aaBaBaBAbbAbAbAB', 'aaBAbaBAbbABabAB', 'aaBAbbAAbbABaaBB', 'aaBABAbbAAbabaBB', 'abababAbABABABaB', 'ababaBabaBaBabaB', 'ababaBabABABAbAB', 'ababAbAbABABaBaB', 'abaBabaBAbABAbAB', 'abaBaBabABAbAbAB', 'abaBAbabABAbaBAB', 'abaBAbaBAbABabAB', 'abAbaBabABaBAbAB']
upto8 = ['aababbAABABB', 'aabaBBAABAbb', 'aabbaBAABBAb', 'aabbAbAABBaB', 'aabAbbAABaBB', 'aabABabbaBAb', 'aaBAbaBBabAB', 'aaBAbbAAbaBB', 'ababAbABABaB', 'abaBabABAbAB', 'abaBAbaBabAB', 'aaaabaBAAAABAb', 'aaaaBabAAAAbAB', 'aaababbAAABABB', 'aaabaBBAAABAbb', 'aaabbabAAABBAB', 'aaabbaBAAABBAb', 'aaabbAbAAABBaB', 'aaabAbbAAABaBB', 'aaabABababaBAb', 'aaaBabbAAAbABB', 'aaaBAbaBaBabAB', 'aaaBAbbAAAbaBB', 'aababaBAABABAb', 'aabaBaBAABAbAb', 'aabAbaBAABaBAb', 'aaBababAAbABAB', 'aaBaBabAAbAbAB', 'aaBAbabAAbaBAB', 'aaaaabaBAAAAABAb', 'aaaaaBabAAAAAbAB', 'aaaababbAAAABABB', 'aaaabaBBAAAABAbb', 'aaaabbabAAAABBAB', 'aaaabbaBAAAABBAb', 'aaaabbAbAAAABBaB', 'aaaabAbbAAAABaBB', 'aaaabABabaabaBAb', 'aaaaBabbAAAAbABB', 'aaaaBAbaBaaBabAB', 'aaaaBAbbAAAAbaBB', 'aaabaabbAAABAABB', 'aaabaaBBAAABAAbb', 'aaababaBAAABABAb', 'aaababbbAAABABBB', 'aaababAbAAABABaB', 'aaababABAAABABab', 'aaabaBabaaaBabaB', 'aaabaBaBAAABAbAb', 'aaabaBAbAAABAbaB', 'aaabaBBBAAABAbbb', 'aaabbaabAAABBAAB', 'aaabbaaBAAABBAAb', 'aaabbbaBAAABBBAb', 'aaabbbAbAAABBBaB', 'aaabbAAbAAABBaaB', 'aaabAbabAAABaBAB', 'aaabAbaBAAABaBAb', 'aaabAbbbAAABaBBB', 'aaabAAbbAAABaaBB', 'aaabABabAAABabAB', 'aaabABAbaaaBAbAB', 'aaaBaabbAAAbAABB', 'aaaBababAAAbABAB', 'aaaBabAbAAAbABaB', 'aaaBaBabAAAbAbAB', 'aaaBAbabAAAbaBAB', 'aaaBAbbbAAAbaBBB', 'aaaBAAbbAAAbaaBB', 'aabaababAABAABAB', 'aabaabaBAABAABAb', 'aabaabAbAABAABaB', 'aabaaBabaaBaabaB', 'aabaaBabAABAAbAB', 'aabaaBAbaaBaabAB', 'aababaabAABABAAB', 'aabababbAABABABB', 'aababaBBAABABAbb', 'aababAAbAABABaaB', 'aabaBaabAABAbAAB', 'aabaBaaBAABAbAAb', 'aabaBabbAABAbABB', 'aabaBAbbAABAbaBB', 'aabaBAAbaaBabAAB', 'aabaBAAbAABAbaaB', 'aabaBABBAABAbabb', 'aabaBBaBAABAbbAb', 'aabbaBabAABBAbAB', 'aabbaBaBAABBAbAb', 'aabbaBAbabaBabAB', 'aabbaBAbAABBAbaB', 'aabbAbabAABBaBAB', 'aabbAbAbAABBaBaB', 'aabbABabAbabaBAb', 'aabbABabAABBabAB', 'aabAbaabAABaBAAB', 'aabAbaBAbbaBabAB', 'aabAbAbbAABaBaBB', 'aabAbAAbAABaBaaB', 'aabAAbabAABaaBAB', 'aabAAbaBAABaaBAb', 'aabAAbAbAABaaBaB', 'aabAABabAABaabAB', 'aabABabaBaBAbaBB', 'aabABabaBBAbaBAB', 'aabABabbAABabABB', 'aabABAbbAABabaBB', 'aaBaaBabAAbAAbAB', 'aaBabAbbAAbABaBB', 'aaBABAbbAAbabaBB', 'abababAbABABABaB', 'ababaBabABABAbAB', 'ababAbAbABABaBaB', 'abaBaBabABAbAbAB', 'abaBAbabABAbaBAB', 'abAbaBabABaBAbAB']