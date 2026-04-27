# ------------------------------------------------------------------

def enumerate(seq):
    i = 0
    result = []
    for item in seq:
        result.append((i, item))
        i = i + 1
    return result


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
    """Return a fraction in num/den format (without f-strings)."""
    return str(fr[0]) if fr[1] == 1 else "{}/{}".format(fr[0], fr[1])


def newton_coeffs(x, y):
    n = len(x)
    coef = y[:]
    table = [coef[:]]
    print("\ndivided differences:")
    for j in range(1, n):
        row = coef[:]
        print("\nOrder {}:".format(j))
        for i in range(n - 1, j - 1, -1):
            a = coef[i]
            b = coef[i - 1]
            xi = x[i]
            xij = x[i - j]
            num = sub_frac(a, b)
            den = sub_frac(xi, xij)
            coef[i] = div_frac(num, den)
            row[i] = coef[i]
            print(
                "f[x{},...,x{}] = (f[x{},...,x{}] - f[x{},...,x{}]) / (x{} - x{})".format(
                    i - j,
                    i,
                    i - j + 1,
                    i,
                    i - j,
                    i - 1,
                    i,
                    i - j,
                )
            )
            print(
                "               = ({} - {}) / ({} - {})".format(
                    format_frac(a),
                    format_frac(b),
                    format_frac(xi),
                    format_frac(xij),
                )
            )
            print(
                "               = ({}) / ({}) = {}".format(
                    format_frac(num),
                    format_frac(den),
                    format_frac(coef[i]),
                )
            )
        table.append(row[:])
    print("\nDivided differences table:")
    for r in range(n):
        row = [format_frac(val) if i >= r else "" for i, val in enumerate(table[r])]
        print("\t".join(row))
    return coef


def newton_polynomial(x, coef):
    terms = [format_frac(coef[0])]
    for i in range(1, len(coef)):
        part = format_frac(coef[i])
        for j in range(i):
            part += "*(x - {})".format(format_frac(x[j]))
        terms.append(part)
    return " + ".join(terms)


# ------------------------- Main program -------------------------

x_input = input("x: ")
y_input = input("y: ")

x_vals = parse_poly(x_input)
y_vals = parse_poly(y_input)

if len(x_vals) != len(y_vals):
    raise ValueError("x != y")

coef = newton_coeffs(x_vals, y_vals)
poly = newton_polynomial(x_vals, coef)

print("\nNewton polynomial:")
print("P(x) =", poly)

print("\nDivided differences coefficients:")
for i, c in enumerate(coef):
    print("a[{}] = {}".format(i, format_frac(c)))
