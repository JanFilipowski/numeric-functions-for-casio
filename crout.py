# ------------------------------------------------------------------
# Euclid: greatest common divisor
# ------------------------------------------------------------------
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# ------------------------------------------------------------------
# Fraction as tuple (num, den), always in reduced form
# ------------------------------------------------------------------
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
# ------------------------------------------------------------------
def parse_line(s):
    parts = s.split(',')
    coeffs = []
    for t in parts:
        t = t.strip()
        if not t:
            continue
        if '/' in t:
            p, q = t.split('/')
            coeffs.append(make_frac(int(p), int(q)))
        elif '.' in t:
            sign = -1 if t.startswith('-') else 1
            raw = t.lstrip('+-')
            if raw.startswith('.'):
                raw = '0' + raw
            if raw.endswith('.'):
                raw = raw[:-1]
            ip, fp = raw.split('.')
            den = 10 ** len(fp)
            num = int(ip) * den + int(fp)
            coeffs.append(make_frac(sign * num, den))
        else:
            coeffs.append(make_frac(int(t), 1))
    return coeffs

# ------------------------------------------------------------------
# Format fraction as text
# ------------------------------------------------------------------
def format_frac(fr):
    if fr[1] == 1:
        return str(fr[0])
    return "{}/{}".format(fr[0], fr[1])

# ------------------------------------------------------------------
# Thomas algorithm for a tridiagonal matrix
# ------------------------------------------------------------------
def thomas(A, d_list, deg):
    # Allocate arrays with 1-based indexing
    a = [None] * (deg+1)
    b = [None] * (deg+1)
    c = [None] * (deg+1)
    d = [None] * (deg+1)
    # Load coefficients from A and D
    for i in range(1, deg+1):
        a[i] = A[i-1][i-1]
        d[i] = d_list[i-1]
        b[i] = A[i-1][i] if i < deg else make_frac(0,1)
        c[i] = A[i-1][i-2] if i > 1 else make_frac(0,1)

    # Print c, a, b, d
    print("Coefficients:")
    for i in range(1, deg+1): print("c{} = {}".format(i, format_frac(c[i])))
    for i in range(1, deg+1): print("a{} = {}".format(i, format_frac(a[i])))
    for i in range(1, deg+1): print("b{} = {}".format(i, format_frac(b[i])))
    for i in range(1, deg+1): print("d{} = {}".format(i, format_frac(d[i])))
    print("")

    # Initialize y and u
    u = [None] * (deg+1)
    y = [None] * (deg+1)

    # Initial step
    y[1] = div_frac(d[1], a[1])
    print("y1 = d1/a1 = {} / {} = {}".format(format_frac(d[1]), format_frac(a[1]), format_frac(y[1])))
    u[1] = div_frac(b[1], a[1])
    print("u1 = b1/a1 = {} / {} = {}".format(format_frac(b[1]), format_frac(a[1]), format_frac(u[1])))

    # Forward elimination
    for i in range(2, deg+1):
        denom = sub_frac(a[i], mul_frac(u[i-1], c[i]))
        if i < deg:
            u[i] = div_frac(b[i], denom)
            print("u{} = b{} / (a{} - u{}*c{}) = {} / ({} - {}*{}) = {}".format(
                i, i, i, i-1, i,
                format_frac(b[i]), format_frac(a[i]), format_frac(u[i-1]), format_frac(c[i]), format_frac(u[i])
            ))
        y[i] = div_frac(sub_frac(d[i], mul_frac(c[i], y[i-1])), denom)
        print("y{} = (d{} - c{}*y{}) / (a{} - u{}*c{}) = ({} - {}*{}) / {} = {}".format(
            i, i, i, i-1, i, i-1, i,
            format_frac(d[i]), format_frac(c[i]), format_frac(y[i-1]), format_frac(denom), format_frac(y[i])
        ))

    print("")
    # Back substitution
    x = [None] * (deg+1)
    x[deg] = y[deg]
    print("x{} = y{} = {}".format(deg, deg, format_frac(x[deg])))
    for i in range(deg-1, 0, -1):
        x[i] = sub_frac(y[i], mul_frac(u[i], x[i+1]))
        print("x{} = y{} - u{}*x{} = {} - {}*{} = {}".format(
            i, i, i, i+1,
            format_frac(y[i]), format_frac(u[i]), format_frac(x[i+1]), format_frac(x[i])
        ))
    return x

# ------------------------------------------------------------------
# Main program
# ------------------------------------------------------------------
def main():
    deg = int(input("S: "))
    A = []
    for i in range(1, deg+1):
        A.append(parse_line(input("{}: ".format(i))))
    print("Solution vector:")
    D = parse_line(input("D: "))
    print("")
    thomas(A, D, deg)

main()
