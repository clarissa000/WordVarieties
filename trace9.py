# -------------------------------------------------
# GENERIC MATRIX METHOD (Section 9 style)
# -------------------------------------------------

from sage.all import *
from time import perf_counter

# Base polynomial ring for the final answer
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
    """
    Accept a string like 'abbaAABaBA' or an iterable of letters.
    """
    if isinstance(word, str):
        return tuple(word)
    return tuple(word)

def trace_poly_generic(word, reduce_first=True):
    """
    Compute the Fricke polynomial using the generic matrix model.

    The result is returned as a polynomial in ZZ[x,y,z].
    """
    word = to_tuple_word(word)

    if reduce_first:
        word = free_reduce(word)

    M = identity_matrix(Rquad, 2)
    for c in word:
        M = M * GEN_MATS[c]

    tr = M.trace()

    # Proposition 9.3 says the trace has no u-part.
    tr_lift = tr.lift()          # polynomial in S = ZZ[x,y,z,u]
    tr_xyz = tr_lift.subs(uS=0)  # should remove nothing if the computation is correct

    # Convert to ZZ[x,y,z]
    return Pxyz(tr_xyz)

def trace_poly_generic_fast(word):
    return trace_poly_generic(word, reduce_first=True)

def benchmark_generic(words, repeats=1):
    start = perf_counter()
    out = None
    for _ in range(repeats):
        for w in words:
            out = trace_poly_generic_fast(w)
    return perf_counter() - start, out

words = [
    "abAB",
    "aaab",
    "aBab",
    "abbaAABaBA",
]

for w in words:
    print(w, trace_poly_generic_fast(w))