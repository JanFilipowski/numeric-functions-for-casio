def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def make_frac(num, den):
    if den == 0:
        raise ZeroDivisionError
    if den < 0:
        num, den = -num, -den
    g = gcd(abs(num), den)
    return (num // g, den // g)

def add_frac(a, b):
    return make_frac(a[0] * b[1] + b[0] * a[1], a[1] * b[1])

def sub_frac(a, b):
    return make_frac(a[0] * b[1] - b[0] * a[1], a[1] * b[1])

def mul_frac(a, b):
    return make_frac(a[0] * b[0], a[1] * b[1])

def div_frac(a, b):
    if b[0] == 0:
        raise ZeroDivisionError
    return make_frac(a[0] * b[1], a[1] * b[0])

def parse_poly(s):
    parts = s.split(",")
    coeffs = []
    for t in parts:
        t = t.strip()
        if not t:
            continue
        if "/" in t:
            i = t.find("/")
            p = t[:i]
            q = t[i + 1:]
            coeffs.append(make_frac(int(p), int(q)))
        elif "." in t:
            sign = -1 if t[0] == "-" else 1
            raw = t.lstrip("+-")
            if raw[0] == ".":
                raw = "0" + raw
            if raw[-1] == ".":
                raw = raw[:-1]
            d = raw.find(".")
            ip = raw[:d]
            fp = raw[d + 1:]
            if fp == "":
                coeffs.append(make_frac(sign * int(ip), 1))
            else:
                den = 10 ** len(fp)
                num = int(ip) * den + int(fp)
                coeffs.append(make_frac(sign * num, den))
        else:
            coeffs.append(make_frac(int(t), 1))
    return coeffs

def format_frac(fr):
    return str(fr[0]) if fr[1] == 1 else "{}/{}".format(fr[0], fr[1])

# Neville algorithm with indices P[i][k]: i - start, k - sub-interval length (degree)
def neville_table(x, y, x0):
    n = len(x)
    P = [[None for k in range(n)] for i in range(n)]

    # P[i][0] = f[x_i]
    for i in range(n):
        P[i][0] = y[i]
        print("P[0][{}] = f[{}] = {}".format(i, i, format_frac(y[i])))
    print("\nComputation steps (Neville table):\n")

    for k in range(1, n):
        for i in range(n - k):
            # (x0 - x_{i+k}) * P[i][k-1] + (x_{i} - x0) * P[i+1][k-1])
            num1 = mul_frac(sub_frac(x0, x[i + k]), P[i][k-1])
            num2 = mul_frac(sub_frac(x[i], x0), P[i+1][k-1])
            num = add_frac(num1, num2)
            den = sub_frac(x[i], x[i + k])
            P[i][k] = div_frac(num, den)
            print("P[{}][{}] = ((x0 - x[{}]) * P[{}][{}] + (x[{}] - x0) * P[{}][{}]) / (x[{}] - x[{}])".format(
                k+i, k, i + k, i, k-1, i, i+1, k-1, i, i + k))
            print("         = (({} - {}) * {} + ({} - {}) * {}) / ({} - {})".format(
                format_frac(x0), format_frac(x[i + k]), format_frac(P[i][k-1]),
                format_frac(x[i]), format_frac(x0), format_frac(P[i+1][k-1]),
                format_frac(x[i]), format_frac(x[i + k])))
            print("         = ({} * {} + {} * {}) / {}".format(
                format_frac(sub_frac(x0, x[i + k])), format_frac(P[i][k-1]),
                format_frac(sub_frac(x[i], x0)), format_frac(P[i+1][k-1]),
                format_frac(den)))
            print("         = ({} + {}) / {} = {} / {} = {}\n".format(
                format_frac(num1), format_frac(num2), format_frac(den),
                format_frac(num), format_frac(den), format_frac(P[i][k])))
    print("Table P (Neville):")
    for i in range(n):
        row = []
        for k in range(n):
            if P[i][k] is not None and i + k < n:
                row.append(format_frac(P[i][k]))
            else:
                row.append("-")
        print("\t".join(row))
    print("\nInterpolation result: P({}) = {}".format(format_frac(x0), format_frac(P[0][n-1])))
    return P[0][n-1]

# ------------------------- Main program -------------------------

x_input = input("x: ")
y_input = input("y: ")
x0_input = input("x0 (point): ")

x_vals = parse_poly(x_input)
y_vals = parse_poly(y_input)
x0 = parse_poly(x0_input.strip())[0]  # single value

if len(x_vals) != len(y_vals):
    raise ValueError("x and y must have the same length")

neville_table(x_vals, y_vals, x0)