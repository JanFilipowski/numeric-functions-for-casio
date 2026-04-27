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

def print_matrix_frac(M, name):
    print(name + ":")
    for row in M:
        line = "  "
        first = True
        for x in row:
            if not first:
                line += "  "
            line += format_frac(x)
            first = False
        print(line)
    print("")

def print_matrix_float(M, name):
    print(name + ":")
    for row in M:
        line = "  "
        first = True
        for x in row:
            if not first:
                line += "  "
            line += "{:.6g}".format(x)
            first = False
        print(line)
    print("")

# ------------------------------------------------------------------

def chol_decomp(A_frac, deg):
    # 1) first convert A_frac to A as float
    A = [[a[0]/a[1] for a in row] for row in A_frac]
    # check symmetry
    for i in range(deg):
        for j in range(i+1, deg):
            if abs(A[i][j] - A[j][i]) > 1e-9:
                print("Error: matrix is not symmetric.")
                return
    # prepare L as a zero float matrix
    L = [[0.0]*deg for _ in range(deg)]

    print_matrix_frac(A_frac, "Matrix A (fraction)")

    # iterate through subsequent rows/columns
    for k in range(deg):
        # sum of squares l[k][0..k-1]
        sum_sq = 0.0
        for j in range(k):
            sum_sq += L[k][j]**2
        # diagonal step
        rad = A[k][k] - sum_sq
        if rad <= 0:
            print("Error: radicand at step {0} is not positive ({1:.6g})."
                  .format(k+1, rad))
            return

        # build expression without f-strings
        expr = "a[{0}][{0}]".format(k+1)
        if k > 0:
            parts = []
            for j in range(k):
                parts.append("({0:.6g})^2".format(L[k][j]))
            subs = " + ".join(parts)
            expr += " - (" + subs + ")"

        lkk = rad**0.5
        # print using .format
        print("l[{0}][{0}] = Sqr({1}) = Sqr({2:.6g}) = {3:.6g}"
              .format(k+1, expr, rad, lkk))
        L[k][k] = lkk

        # off-diagonal elements
        for i in range(k+1, deg):
            sum_pr = 0.0
            for j in range(k):
                sum_pr += L[i][j] * L[k][j]
            num = A[i][k] - sum_pr
            lik = num / lkk
            if k > 0:
                terms = []
                for j in range(k):
                    terms.append("({0:.6g})*({1:.6g})".format(L[i][j], L[k][j]))
                term_str = " + ".join(terms)
                print("l[{0}][{1}] = (a[{0}][{1}] - ({2})) / {3:.6g} = "
                      "({4:.6g} - {5:.6g}) / {3:.6g} = {6:.6g}"
                      .format(i+1, k+1, term_str, lkk, A[i][k], sum_pr, lik))
            else:
                print("l[{0}][{1}] = a[{0}][{1}] / {2:.6g} = "
                      "{3:.6g} / {2:.6g} = {4:.6g}"
                      .format(i+1, k+1, lkk, A[i][k], lik))
            L[i][k] = lik

        print("")

    # print results
    print_matrix_float(L, "Matrix L")
    # transposed matrix
    LT = [[L[j][i] for j in range(deg)] for i in range(deg)]
    print_matrix_float(LT, "Matrix L^T")

    if input("Show L.^2? ") != "0":
        L2 = [[L[i][j]**2 for j in range(deg)] for i in range(deg)]
        print_matrix_float(L2, "Matrix L.^2")
    # --- new section: forward and back substitution ---
    answer = input("Add right-hand-side vector b and show solving steps? "
                   "(yes/no): ").strip().lower()
    if answer in ("yes", "y", "1"):
        raw = input("Enter {0} coefficients of vector b, separated by commas: "
                    .format(deg))
        b_frac = parse_line(raw)
        b = []
        for num, den in b_frac:
            b.append(num/den)

        # forward: Ly = b
        y = [0.0] * deg
        print("\n=== Forward substitution: Ly = b ===")
        for i in range(deg):
            # 1) symbolic form
            terms = []
            for j in range(i):
                terms.append("l[{0}][{1}]*y{1}".format(i+1, j+1))
            sym = " + ".join(terms) or "0"
            print("y{0} = (b{0} - ({1})) / l[{0}][{0}]"
                  .format(i+1, sym))
            # 2) numeric substitution form
            terms = []
            for j in range(i):
                terms.append("{0:.6g}*{1:.6g}".format(L[i][j], y[j]))
            num_sym = " + ".join(terms) or "0"
            print("      = ({0:.6g} - ({1})) / {2:.6g}"
                  .format(b[i], num_sym, L[i][i]))
            # 3) result
            yi = (b[i] - sum(L[i][j] * y[j] for j in range(i))) / L[i][i]
            print("      = {0:.6g}\n".format(yi))
            y[i] = yi

        # back: L^T x = y
        x = [0.0] * deg
        print("\n=== Back substitution: L^T x = y ===")
        for idx in range(deg-1, -1, -1):
            # 1) symbolic form
            terms = []
            for j in range(idx+1, deg):
                terms.append("l[{0}][{1}]*x{0}".format(j+1, idx+1))
            sym = " + ".join(terms) or "0"
            print("x{0} = (y{0} - ({1})) / l[{0}][{0}]"
                  .format(idx+1, sym))
            # 2) numeric substitution form
            terms = []
            for j in range(idx+1, deg):
                terms.append("{0:.6g}*{1:.6g}".format(L[j][idx], x[j]))
            num_sym = " + ".join(terms) or "0"
            print("       = ({0:.6g} - ({1})) / {2:.6g}"
                  .format(y[idx], num_sym, L[idx][idx]))
            # 3) result
            xi = (y[idx] - sum(L[j][idx] * x[j] for j in range(idx+1, deg))) / L[idx][idx]
            print("       = {0:.6g}\n".format(xi))
            x[idx] = xi

        print("Solution x:")
        for idx in range(len(x)):
            i = idx + 1
            xi = x[idx]
            print("  x{0} = {1:.6g}".format(i, xi))


# ------------------------------------------------------------------
def main():
    deg = int(input("S (matrix size): "))
    A_frac = []
    for i in range(deg):
        row_in = input("Row {0}: ".format(i+1))
        A_frac.append(parse_line(row_in))
    chol_decomp(A_frac, deg)

main()
