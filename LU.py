# ------------------------------------------------------------------
# Euclid: greatest common divisor
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# ------------------------------------------------------------------
# Fraction as tuple (num, den), always in reduced form
def make_frac(num, den):
    if den == 0:
        raise ZeroDivisionError
    if den < 0:
        num, den = -num, -den
    g = gcd(abs(num), den)
    return (num // g, den // g)

def add_frac(a, b):
    return make_frac(a[0]*b[1] + b[0]*a[1], a[1]*b[1])

def sub_frac(a, b):
    return make_frac(a[0]*b[1] - b[0]*a[1], a[1]*b[1])

def mul_frac(a, b):
    return make_frac(a[0]*b[0], a[1]*b[1])

def div_frac(a, b):
    if b[0] == 0:
        raise ZeroDivisionError
    return make_frac(a[0]*b[1], a[1]*b[0])

# ------------------------------------------------------------------
# Parse text "1/2, 0.5, -3, .25" -> list of tuples (num, den)
def parse_line(s):
    parts = s.split(',')
    coeffs = []
    for t in parts:
        t = t.strip()
        if not t:
            continue
        if '/' in t:
            i = t.find('/')
            p, q = t[:i], t[i+1:]
            coeffs.append(make_frac(int(p), int(q)))
        elif '.' in t:
            sign = -1 if t[0] == '-' else 1
            raw = t.lstrip('+-')
            if raw[0] == '.':
                raw = '0' + raw
            if raw[-1] == '.':
                raw = raw[:-1]
            d = raw.find('.')
            ip, fp = raw[:d], raw[d+1:]
            if fp == '':
                coeffs.append(make_frac(sign*int(ip), 1))
            else:
                den = 10 ** len(fp)
                num = int(ip)*den + int(fp)
                coeffs.append(make_frac(sign*num, den))
        else:
            coeffs.append(make_frac(int(t), 1))
    return coeffs

# ------------------------------------------------------------------
# Format fraction as text
def format_frac(fr):
    if fr[1] == 1:
        return str(fr[0])
    else:
        return str(fr[0]) + "/" + str(fr[1])

# ------------------------------------------------------------------
def abs_frac(fr):
    return (abs(fr[0]), fr[1])

def is_zero(fr):
    return fr[0] == 0

def is_greater_abs(a, b):
    return abs(a[0]) * b[1] > abs(b[0]) * a[1]

def print_matrix(M, name):
    print(name + ":")
    for row in M:
        line = []
        for fr in row:
            line.append(format_frac(fr))
        print("  " + " ".join(line))
    print("")

# ------------------------------------------------------------------
def lu_factor(A, deg):
    pivot = [i+1 for i in range(deg)]
    swap_count = 0
    U = [row.copy() for row in A]
    L = [[ (1,1) if i==j else (0,1) for j in range(deg)] for i in range(deg)]

    for k in range(deg):
        # partial pivoting
        max_row = k
        for i in range(k+1, deg):
            if is_greater_abs(U[i][k], U[max_row][k]):
                max_row = i
        if is_zero(U[max_row][k]):
            print("Zero pivot at (" + str(k+1) + "," + str(k+1) + "), skipping further eliminations.\n")
            continue
        if max_row != k:
            U[k], U[max_row] = U[max_row], U[k]
            for j in range(k):
                L[k][j], L[max_row][j] = L[max_row][j], L[k][j]
            pivot[k], pivot[max_row] = pivot[max_row], pivot[k]
            swap_count += 1
            print("Swap: " + str(k+1) + " with " + str(max_row+1))
        print("Pivoting: " + str(pivot) + "\n")

        # elimination
        for i in range(k+1, deg):
            l = div_frac(U[i][k], U[k][k])
            print("l[" + str(i+1) + "][" + str(k+1) + "] = " +
                  format_frac(U[i][k]) + " / " + format_frac(U[k][k]) +
                  " = " + format_frac(l))
            L[i][k] = l
            for j in range(k, deg):
                U[i][j] = sub_frac(U[i][j], mul_frac(l, U[k][j]))
            print_matrix(U, "U after elimination of row " + str(i+1))

    print_matrix(L, "Matrix L")
    print_matrix(U, "matrix U")
    # determinant
    det = (1,1)
    for i in range(deg):
        det = mul_frac(det, U[i][i])
    if swap_count % 2 == 1:
        det = (-det[0], det[1])
    print("Determinant for A = " + format_frac(det) + "\n")

    return pivot, L, U
# ------------------------------------------------------------------
def invert_via_lu(pivot, L, U, deg):
    # prepare inverse matrix as a list of columns
    inv = [[(0,1) for _ in range(deg)] for _ in range(deg)]

    # for each unit column
    for k in range(deg):
        # 1) P*e_k
        b = []
        for i in range(deg):
            if pivot[i] == k+1:
                b.append((1,1))
            else:
                b.append((0,1))
        print("L*U*X(" + str(k+1) + ") = " +
              "[" + ", ".join("1" if i==k else "0" for i in range(deg)) +
              "]^T")
        print("P*e" + str(k+1) + " = b = [" +
              ", ".join(format_frac(x) for x in b) + "]^T\n")

        # 2) Solve L*y = b
        y = [(0,1)] * deg
        print("Solve for L * y = b:")
        for i in range(deg):
            # sum L[i][j] * y[j]
            s = (0,1)
            terms = []
            for j in range(i):
                t = mul_frac(L[i][j], y[j])
                s = add_frac(s, t)
                terms.append(format_frac(L[i][j]) + "*y" + str(j+1))
            rhs = sub_frac(b[i], s)
            y[i] = rhs  # because L[i][i] = 1
            if i == 0:
                print("y1 = b1 = " + format_frac(b[0]))
            else:
                print(
                  "y" + str(i+1) + " = (b" + str(i+1) + " - " +
                  " + ".join(terms) + ") = " + format_frac(b[i]) +
                  " - " + format_frac(s) + " = " + format_frac(y[i])
                )
        print("y = [" + ", ".join(format_frac(v) for v in y) + "]^T\n")

        # 3) Solve U*x = y
        x = [(0,1)] * deg
        print("Solve for U * x = y:")
        for i in range(deg-1, -1, -1):
            s = (0,1)
            terms = []
            for j in range(i+1, deg):
                t = mul_frac(U[i][j], x[j])
                s = add_frac(s, t)
                terms.append(format_frac(U[i][j]) + "*x" + str(j+1))
            num = sub_frac(y[i], s)
            x[i] = div_frac(num, U[i][i])
            if i == deg-1:
                print(
                  "x" + str(deg) + " = y" + str(deg) +
                  " / U[" + str(deg) + "][" + str(deg) + "] = " +
                  format_frac(y[i]) + " / " + format_frac(U[i][i]) +
                  " = " + format_frac(x[i])
                )
            else:
                print(
                  "x" + str(i+1) + " = (y" + str(i+1) + " - " +
                  " - ".join(terms) + ") / " +
                  format_frac(U[i][i]) + " = (" +
                  format_frac(y[i]) + " - " +
                  format_frac(s) + ") / " + format_frac(U[i][i]) +
                  " = " + format_frac(x[i])
                )
        print("X(" + str(k+1) + ") = [" +
              ", ".join(format_frac(v) for v in x) + "]^T\n")

        # store computed X(k) as a column of the inverse matrix
        for i in range(deg):
            inv[i][k] = x[i]

    # finally print the full inverse matrix
    print_matrix(inv, "A^-1:")

def main():
    deg = int(input("S: "))
    A = []
    for i in range(1, deg+1):
        A.append(parse_line(input(str(i) + ": ")))
    print_matrix(A, "Matrix A")
    pivot, L, U = lu_factor(A, deg)
    k = input("Reverse matrix? ")
    if k == "0":
        return
    invert_via_lu(pivot, L, U, deg)

main()
