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
def parse_poly(s):
    parts = s.split(',')
    coeffs = []
    for t in parts:
        t = t.strip()
        if not t:
            continue

        # common fraction
        if '/' in t:
            i = t.find('/')
            p = t[:i]
            q = t[i+1:]
            coeffs.append(make_frac(int(p), int(q)))

        # decimal fraction
        elif '.' in t:
            sign = -1 if t[0] == '-' else 1
            raw = t.lstrip('+-')
            if raw[0] == '.':
                raw = '0' + raw
            if raw[-1] == '.':
                raw = raw[:-1]
            d = raw.find('.')
            ip = raw[:d]
            fp = raw[d+1:]
            if fp == '':
                coeffs.append(make_frac(sign*int(ip), 1))
            else:
                den = 10 ** len(fp)
                num = int(ip)*den + int(fp)
                coeffs.append(make_frac(sign*num, den))

        # integer number
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
# Horner synthetic division:
# divide poly by (x - x0), return (quotient_coeffs, remainder)
def synthetic_division(poly, x0):
    n = len(poly)
    W = [None] * n
    W[0] = poly[0]

    # successive Horner steps
    for i in range(1, n):
        W[i] = add_frac( mul_frac(W[i-1], x0), poly[i] )

    # print w_n ... w_0
    deg = n - 1
    for idx in range(n):
        print("w%d= %s" % (deg - idx, format_frac(W[idx])))

    # quotient is W[0..n-2], remainder is W[n-1]
    return W[:-1], W[-1]

# ------------------------------------------------------------------
# Main function
def main():
    try:
        # read polynomial
        print("Polynomial coefficients: ", end="")
        Q = parse_poly(input())

        # read x0
        print("x0: ", end="")
        x0 = parse_poly(input())[0]

        print("\n=== Horner synthetic division results ===")
        while True:
            Q, R = synthetic_division(Q, x0)

            # print quotient
            print("\nQuotient: ", end="")
            L = len(Q)
            if L == 0:
                print("0", end="")
            else:
                for idx in range(L):
                    coeff = Q[idx]
                    deg = L - idx - 1
                    print(format_frac(coeff), end="")
                    if deg > 0:
                        print("*x^%d" % deg, end="")
                        if idx < L-1:
                            print(" + ", end="")

            # print remainder
            print("\nRemainder R = %s" % format_frac(R))

            # ask whether to continue
            print("Continue? (0 = no): ", end="")
            a = input()
            if a == "0":
                break

    except Exception as e:
        print("Error:", e)

# program entry point
main()
