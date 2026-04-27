def parse_poly(s):
    # Split a coefficient sequence given as "a_n,a_{n-1},...,a_0"
    parts = s.split(',')
    res = []
    for item in parts:
        # Convert to float
        res.append(float(item))
    return res

# Read polynomial coefficients (from highest degree to constant term)
A = parse_poly(input("A: "))

# Reverse list so a[i] is the coefficient at x^i
a = []
i = len(A) - 1
while i >= 0:
    a.append(A[i])
    i -= 1

# Polynomial degree
n = len(a) - 1

# Read point x
x = float(input("x: "))

# Read parameter p
p = int(input("p (best value: 1): "))
# Check whether p divides n+1
if (n + 1) % p != 0:
    print("Invalid parameter p (it must be a divisor of n+1)")
    raise SystemExit

q = (n + 1) // p

# Helper functions s(j) and r(j)
def s(j):
    return (n - j) % q

def r(j):
    if j % q == 0:
        return q
    else:
        return 0

# Prepare tables T_tbl ((n+1)x(n+1)) and T_neg1 (for j = -1)
T_tbl = []
i = 0
while i <= n:
    # for each i create a row with n+1 zeros
    row = []
    j = 0
    while j <= n:
        row.append(0.0)
        j += 1
    T_tbl.append(row)
    i += 1

T_neg1 = [0.0] * (n + 1)

# Print s(j)
print("\ns(j) =")
for j in range(n + 1):
    print("{} for j = {}".format(s(j), j))

# Print r(j)
print("\nr(j) =")
for j in range(n + 1):
    print("{} for j = {}".format(r(j), j))

# Fill table T[-1][j]
print("\nFilling T[-1][j]:")
i = 0
while i <= n:
    si1 = s(i + 1)
    T_neg1[i] = a[n - i - 1] * (x ** si1) if si1 >= 0 else 0.0
    print("T[-1][{}] = a[{}] * x^s({}) = {} * {}^{} = {}".format(i, n - i - 1, i + 1, a[n - i - 1], x, si1, T_neg1[i]))
    i += 1


# Fill tables
if x == 0.0:
    # If x == 0, then for each i and j>=0: T(i,j) = a[j]
    # (according to the original code, function T returns a[j] when x == 0)
    i = 0
    while i <= n:
        si1 = s(i + 1)
        T_neg1[i] = a[n - i - 1] * (0.0 ** si1) if si1 >= 0 else 0.0

        # For j from 0 to i: T(i,j) = a[j]
        j = 0
        while j <= i:
            T_tbl[i][j] = a[j]
            j += 1

        i += 1
else:
    # x != 0: compute iteratively
    # First compute s0 = s(0)
    s0 = s(0)

    # i = 0: base row
    T_tbl[0][0] = a[n] * (x ** s0)
    T_neg1[0] = a[n - 1] * (x ** s(1))

    # Print step 0
    print("\ns(0) = {}, r(0) = {}".format(s0, r(0)))
    print("T[0][0] = a[n] * x^{} = {} * {}^{} = {}".format(s0, a[n], x, s0, T_tbl[0][0]))

    print("T[0][0]", end="")
    for i in range(1, n + 1):
        print(" = T[{}][{}]".format(i, i), end="")
    print("\n")

    # For i = 1..n:
    i = 1
    while i <= n:
        # Compute T_neg1[i] first (used in the next step)
        si1 = s(i + 1)
        T_neg1[i] = a[n - i - 1] * (x ** si1)

        # Now fill T_tbl[i][j] for j = 0..i
        j = 0
        while j <= i:
            if j == i:
                T_tbl[i][j] = a[n] * (x ** s0)
            else:
                if j == 0:
                    prev_neg = T_neg1[i - 1]
                    prev_zero = T_tbl[i - 1][0]
                    T_tbl[i][0] = prev_neg + prev_zero * (x ** r(i - 0))
                    print("T[{}][0] = T[-1][{}] + T[{}][0] * x^{} = {} + {} * {}^{} = {}".format(i, i - 1, i - 1, 0, prev_neg, prev_zero, x, r(i - 0), T_tbl[i][0]))
                    if j == i:
                        break
                else:
                    prev_left = T_tbl[i - 1][j - 1]
                    prev_same = T_tbl[i - 1][j]
                    T_tbl[i][j] = prev_left + prev_same * (x ** r(i - j))
                    print("T[{}][{}] = T[{}][{}] + T[{}][{}] * x^{} = {} + {} * {}^{} = {}".format(i, j, i - 1, j - 1, i - 1, j, r(i - j), prev_left, prev_same, x, r(i - j), T_tbl[i][j]))

            j += 1

        i += 1



# Now print full table T[i][j] for i,j = 0..n
print("\nTable T[j][i] (j = -1..", n, ", i = 0..", n, "):")
i = 0

row_str = "j=-1:"
j = 0
while j < n:
    row_str += " " + str(T_neg1[j])
    j += 1
print(row_str)

while i <= n:
    row_str = "j=" + str(i) + ": "
    j = 0
    while j <= n:
        row_str += " " + str(T_tbl[j][i])
        j += 1
    print(row_str)
    i += 1

# After filling T_tbl and T_neg1, compute normalized derivatives for all degrees k = 0..n
print("\nNormalized derivatives for degrees 0..n at point", x, ":")
k = 0
while k <= n:
    if x == 0.0:
        der = T_tbl[n][k]
    else:
        der = T_tbl[n][k] / (x ** (k % q))
    print("f({})_Z = T[{}][{}] / ({}^({}%{}) = {}".format(k, n, k, x, k, q, der))
    k += 1