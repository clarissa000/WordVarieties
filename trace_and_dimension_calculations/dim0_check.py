from sage.all import *
import pandas as pd

# Base ring for trace polynomials
R3 = PolynomialRing(QQ, ['x', 'y', 'z'])
x, y, z = R3.gens()

# Bigger ring for the ideal with u
R10 = PolynomialRing(QQ, ['x', 'y', 'z', 'u', 'v11', 'v12', 'v21', 'v22', 'v31', 'v32'])
x4, y4, z4, u, v11, v12, v21, v22, v31, v32 = R10.gens()

sqrt2 = QQbar(sqrt(2))
phi = (1 + QQbar(sqrt(5))) / 2

E = {
    QQbar(0),
    QQbar(1), QQbar(-1),
    sqrt2, -sqrt2,
    phi, -phi,
    1 - phi, -(1 - phi)
}

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

def dimension_two_rings(word):
    pw  = trace_poly(word)
    paw = trace_poly('a'+ word)
    pbw = trace_poly('b' + word)

    F = x**2 + y**2 + z**2 - x*y*z - 4

    # Ideal in R3
    I3 = ideal([pw - 2, paw - x, pbw - y])

    print('\n\nWord', word)
    print('Groebner basis for I3:')
    print(I3.groebner_basis())
    print('Dimension of I3:')
    print(I3.dimension())

    # Ideal in R10 = R3[u, v11, ...]
    I10 = ideal([R10(pw) - 2, R10(paw) - x4, R10(pbw) - y4, u*R10(F) - 1, v11*x4 + v12*z4 - 1, v21*x4 + v22*y4 - 1, v31*y4 + v32*z4 - 1])
    print('\nGroebner basis for I10:')
    gb = I10.groebner_basis()
    newgb = []
    for term in gb:
        if term.subs({u:0, v11: 0, v12:0, v21:0, v22:0, v31:0, v32:0}) == term:
            newgb.append(term)
    print(newgb)
    print('Dimension of I10:')
    #3 bigger? 
    print(I10.dimension()-3)

    #saturation method?
    ideal_F = ideal([F])
    ideal_axes = ideal([x*y, x*z, y*z])

    I3_cleaned = I3.saturation(ideal_F)[0]
    I3_cleaned = I3_cleaned.saturation(ideal_axes)[0]
    print('\nGroebner basis for I3 cleaned:')
    print(I3_cleaned.groebner_basis())
    print('Dimension of I3 cleaned:')
    print(I3_cleaned.dimension())

def additional_variables(word):
    pw  = trace_poly(word)
    paw = trace_poly('a'+ word)
    pbw = trace_poly('b' + word)

    F = x**2 + y**2 + z**2 - x*y*z - 4

    # Ideal in R10 = R3[u, v11, ...]
    I10 = ideal([R10(pw) - 2, R10(paw) - x4, R10(pbw) - y4, u*R10(F) - 1, v11*x4 + v12*z4 - 1, v21*x4 + v22*y4 - 1, v31*y4 + v32*z4 - 1])
    print('\nGroebner basis for I10:')
    gb = I10.groebner_basis().subs({u:0, v11: 0, v12:0, v21:0, v22:0, v31:0, v32:0})
    print(gb)

    print('Dimension of I10, subtracting 3:')
    #3 bigger? 
    print(I10.dimension() - 3)

def analyse_groebner_basis(gb):
    # Compute common gcd
    g = gb[0]
    for f in gb[1:]:
        g = gcd(g, f)

    print("Common gcd:")
    print(factor(g))
    print()

    print("Generators divided by gcd:")
    for i, f in enumerate(gb, start=1):
        q = f // g
        print(f"g{i} =")
        print(factor(q))
        print()

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
    J = I3_cleaned
    print('\nGroebner basis for I3 cleaned:')
    print(I3_cleaned.groebner_basis())
    print('Dimension of I3 cleaned:')
    print(I3_cleaned.dimension())
    if I3_cleaned.dimension() == 0:
        print('Vector space dimension for ideal & radical')
        print(I3_cleaned.vector_space_dimension())
        print(I3_cleaned.variety(QQbar))
        print(J.radical() == J)
        print(J.radical().vector_space_dimension())
    print(analyse_groebner_basis(I3_cleaned.groebner_basis()))

def finite_strongly_irreducible_point(point):
    xv = QQbar(point[x])
    yv = QQbar(point[y])
    zv = QQbar(point[z])

    # kappa = tr([A,B])
    kappa = xv**2 + yv**2 + zv**2 - xv*yv*zv - 2

    # No two of x,y,z are zero
    if ((xv == 0 and yv == 0) or
        (xv == 0 and zv == 0) or
        (yv == 0 and zv == 0)):
        return False

    # kappa != 0
    if kappa == 0:
        return False

    # x,y,z,kappa all lie in E
    return (
        xv in E and
        yv in E and
        zv in E and
        kappa in E
    )

def analyse_dim0_word(word, verbose=True):

    pw = trace_poly(word)
    paw = trace_poly('a' + word)
    pbw = trace_poly('b' + word)

    F = x**2 + y**2 + z**2 - x*y*z - 4

    I3 = ideal([
        pw - 2,
        paw - x,
        pbw - y
    ])

    ideal_F = ideal([F])
    ideal_axes = ideal([x*y, x*z, y*z])

    # Same cleaning you already use
    J = I3.saturation(ideal_F)[0]
    J = J.saturation(ideal_axes)[0]

    if J.dimension() != 0:
        raise ValueError(
            f"{word}: expected dimension 0, got {J.dimension()}"
        )

    solutions = J.variety(QQbar)

    finite_points = []
    zdense_points = []

    for point in solutions:
        if finite_strongly_irreducible_point(point):
            finite_points.append(point)
        else:
            zdense_points.append(point)

    if len(zdense_points) == 0:
        zdense_dimension = -1
    else:
        zdense_dimension = 0

    if verbose:
        print("\n" + "=" * 70)
        print("WORD:", word)
        print("=" * 70)
        print("Number of solutions:", len(solutions))
        print("Finite strongly irreducible:", len(finite_points))
        print("Non-finite:", len(zdense_points))
        print("ZD dimension:", zdense_dimension)

        if zdense_points:
            print("\nNon-finite points:")
            for p in zdense_points:
                xv = QQbar(p[x])
                yv = QQbar(p[y])
                zv = QQbar(p[z])
                kappa = xv**2 + yv**2 + zv**2 - xv*yv*zv - 2

                print("  ", p)
                print("     kappa =", kappa)

    return zdense_dimension

csv_in = "databases/word_dimensions_nogb.csv"
csv_out = "databases/word_dimensions_canon_cardinality_ZDense.csv"

df = pd.read_csv(csv_in)
df = df[df["word"].str.len() <= 14]

# Only the original dimension-0 cases
dim0 = df[df["dimensionZDense"] == 0]

print("Total words:", len(df))
print("Dimension-0 words:", len(dim0))
df["dimensionZDense"] = df["dimensionZDense"].astype(int)
df["dimensionZDense2"] = df["dimensionZDense"]

counter = 0

for idx, row in dim0.iterrows():

    counter += 1
    word = row["word"]

    #print(f"[{counter}/{len(dim0)}] {word}")

    try:
        zd_dim = analyse_dim0_word(word, verbose=False)
        df.loc[idx, "dimensionZDense2"] = zd_dim

    except Exception as e:
        print(f"ERROR on {word}: {e}")
        continue

    # Save progress every 1000 words
    if counter % 1000 == 0:
        df.to_csv(csv_out, index=False)
        print(f"--- Saved progress: {counter} words ---")

# Final save
df.to_csv(csv_out, index=False)

print("DONE")
print("Saved to:", csv_out)