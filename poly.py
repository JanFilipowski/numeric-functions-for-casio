# ------------------------------------------------------------------
# Euclid: greatest common divisor
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# ------------------------------------------------------------------
# Create and normalize a fraction as (num, den)
def make_frac(num, den):
    if den == 0:
        raise ZeroDivisionError
    if den < 0:
        num, den = -num, -den
    g = gcd(abs(num), den)
    return (num // g, den // g)

# Fraction operations (a and b are tuples)
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
# Parse input "1/2,0.5,-3,4." -> list of fractions (num,den)
def parse_poly(s):
    parts = s.split(',')
    coeffs = []
    for t in parts:
        t = t.strip()
        if not t:
            continue
        # common fraction
        if '/' in t:
            p, q = t.split('/', 1)
            coeffs.append(make_frac(int(p), int(q)))
        # decimal number
        elif '.' in t:
            sign = -1 if t[0]=='-' else 1
            raw = t.lstrip('+-')
            if raw.startswith('.'):
                raw = '0' + raw
            if raw.endswith('.'):
                raw = raw[:-1]
            ip, fp = raw.split('.',1)
            if fp == '':
                coeffs.append(make_frac(sign*int(ip), 1))
            else:
                den = 10 ** len(fp)
                num = int(ip)*den + int(fp)
                coeffs.append(make_frac(sign*num, den))
        # integer
        else:
            coeffs.append(make_frac(int(t), 1))
    return coeffs

# ------------------------------------------------------------------
# Polynomial formatting: [(2,1),(-3,2),...] -> "2x^4 -3/2x^3 + ..."
def format_poly(coefs):
    terms = []
    n = len(coefs)
    for i in range(n):
        num, den = coefs[i]
        if num == 0:
            continue
        deg = n -i -1
        # build coefficient string
        if den == 1:
            coeff = str(abs(num))
        else:
            coeff = str(abs(num)) + "/" + str(den)
        if coeff == "1" and deg>0:
            coeff = ""
        # monomial
        if deg == 0:
            mon = coeff
        elif deg == 1:
            mon = coeff + "x"
        else:
            mon = coeff + "x^" + str(deg)
        # sign
        if not terms:
            prefix = "-" if num<0 else ""
        else:
            prefix = " - " if num<0 else " + "
        terms.append(prefix + mon)
    if not terms:
        return "0"
    # join the list into one string
    res = ""
    for t in terms:
        res += t
    return res

# ------------------------------------------------------------------
# Format a single quotient term
def format_term(fr, deg):
    num, den = fr
    sign = "-" if num<0 else ""
    base = (str(abs(num)) if den==1 else str(abs(num))+"/"+str(den))
    if base=="1" and deg>0:
        base = ""
    if deg==0:
        return sign + base
    elif deg==1:
        return sign + base + "x"
    else:
        return sign + base + "x^" + str(deg)

# ------------------------------------------------------------------
# Long division of polynomials A, B - lists of fractions
def long_division(A, B):
    m, n = len(A), len(B)
    if n==0 or B[0][0]==0:
        raise ZeroDivisionError("Divisor must have a non-zero leading coefficient")
    # quotient = all zeros
    q = [make_frac(0,1) for _ in range(m-n+1)]
    cur = A[:]  # copy

    print(format_poly(cur), " -Dividend polynomial")
    print(format_poly(B), " -Divisor")
    print("")

    step = 1
    while len(cur) >= n:
        # check whether cur != 0
        all_zero = True
        for x in cur:
            if x[0]!=0:
                all_zero = False
                break
        if all_zero:
            break

        deg_diff = len(cur) -n
        # quotient term = cur[0]/B[0]
        lt = div_frac(cur[0], B[0])
        q[deg_diff] = lt

        # subtrahend = shifted B*lt
        sub = []
        for bj in B:
            sub.append(mul_frac(bj, lt))
        for _ in range(deg_diff):
            sub.append(make_frac(0,1))

        print("Step ", str(step), ":")
        print(format_poly(cur), " -current polynomial")
        print(format_term(lt, deg_diff), " -quotient term")
        print(format_poly(sub), " -polynomial to subtract")

        # subtraction cur - sub
        L = len(cur) if len(cur)>len(sub) else len(sub)
        new_cur = []
        for j in range(L):
            a = cur[j] if j<len(cur) else make_frac(0,1)
            b = sub[j] if j<len(sub) else make_frac(0,1)
            new_cur.append(sub_frac(a, b))
        # remove leading zeros
        while len(new_cur)>1 and new_cur[0][0]==0:
            del new_cur[0]
        cur = new_cur

        print(format_poly(cur), " -subtraction result")
        print("")
        step += 1

    print(format_poly(q[::-1]), " -Quotient")
    print(format_poly(cur), " -Remainder")

# ------------------------------------------------------------------
# Main function
def main():
    try:
        A = parse_poly(input("A: "))
        B = parse_poly(input("B: "))
        long_division(A, B)
    except Exception as e:
        print("Error:", e)

main()
