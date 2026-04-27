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


def lagrange_polynomial(x, y):
    n = len(x)
    terms = []
    print("\nComputing Lagrange interpolation polynomial:")
    for k in range(n):
        print("\nStep {}: for node x_{} = {}, y_{} = {}".format(k+1, k, format_frac(x[k]), k, format_frac(y[k])))
        denom = (1, 1)
        print("  Computing denominator D_{}:".format(k))
        for j in range(n):
            if j == k:
                continue
            diff = sub_frac(x[k], x[j])
            print("    x_{} - x_{} = {} - {} = {}".format(k, j, format_frac(x[k]), format_frac(x[j]),format_frac(diff)))
            denom = mul_frac(denom, diff)
            print("    Multiplication: new denominator = {}".format(format_frac(denom)))
        print("  Total denominator D_{} = {}".format(k, format_frac(denom)))
        # Build symbolic numerator (sequence of factors (x - x_j))
        num_str = ""
        for j in range(n):
            if j == k:
                continue
            term = "(x - {})".format(format_frac(x[j]))
            num_str += term
        print("  Basis polynomial L_{}(x) = ({}) / {}".format(k, num_str, format_frac(denom)))
        # Multiply by y_k
        coef_frac = div_frac(y[k], denom)
        print(
            "  Multiply by y_{}: coefficient = {} / {} = {}".format(k,format_frac(y[k]),format_frac(denom),format_frac(coef_frac)))
        term_str = "{}".format(format_frac(coef_frac))
        for j in range(n):
            if j == k:
                continue
            term_str += "*(x - {})".format(format_frac(x[j]))
        print("  Term in the sum: {}".format(term_str))
        terms.append(term_str)
    poly = " + ".join(terms)
    print("\nLagrange polynomial in summation form:")
    print("P(x) =", poly)
    return poly


# ------------------------- Main program -------------------------

x_input = input("x: ")
y_input = input("y: ")

x_vals = parse_poly(x_input)
y_vals = parse_poly(y_input)

if len(x_vals) != len(y_vals):
    raise ValueError("x != y")

# Instead of Newton form, run Lagrange interpolation with steps
lagrange_polynomial(x_vals, y_vals)
