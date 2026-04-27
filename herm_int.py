# Hermite Interpolation (Micropython 1.9.3 Compatible)
# no external modules, no f-strings, ASCII-friendly

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
    return str(fr[0]) if fr[1] == 1 else "{}/{}".format(fr[0], fr[1])

def input_list(msg, n):
    s = input(msg)
    out = []
    parts = s.split(",")
    for i in range(n):
        if i >= len(parts) or parts[i].strip() == "":
            out.append("")
        else:
            out.append(parts[i].strip())
    return out

def frac_or_empty(s):
    s = s.strip()
    if s == "":
        return None
    return make_frac(int(s), 1)

def print_step(msg):
    print(msg)

# ------------------------- Main program -------------------------

print("Input data:")
x_input = input("x: ")
m_input = input("m: ")

x_vals = [int(z.strip()) for z in x_input.split(",") if z.strip() != ""]
m_vals = [int(z.strip()) for z in m_input.split(",") if z.strip() != ""]
n = len(x_vals)
M = sum(m_vals)

# read f, f', f'' (only up to 2nd derivative for this example)
f_in   = input_list("f(xi): ", n)
fp_in  = input_list("f'(xi): ", n)
fpp_in = input_list("f''(xi): ", n)

# expand Hermite nodes
z = []
for i in range(n):
    for j in range(m_vals[i]):
        z.append(x_vals[i])

# prepare initial table of values and derivatives
d = []
for i in range(n):
    d.append(frac_or_empty(f_in[i]))
    for j in range(1, m_vals[i]):
        if j == 1:
            d.append(frac_or_empty(fp_in[i]))
        elif j == 2:
            d.append(frac_or_empty(fpp_in[i]))
        else:
            d.append(None)

# compute s(i)
s = [0]*n
acc = 0
for i in range(n):
    s[i] = acc
    acc += m_vals[i]

print("\ns(i):")
for i in range(n):
    print("s({}) = {}".format(i, s[i]))

print("\nl (Hermite index):")
l = 0
for i in range(n):
    for j in range(m_vals[i]):
        print("l = {} = s({}) + {}".format(l, i, j))
        l += 1



# -----------------------------------------------------------
# 4. Derive base polynomials p_k(x)
print("\nPolynomials p_k(x):")
def pow_poly(x, m):
    if m == 0:
        return "1"
    out = []
    for _ in range(m):
        out.append("(x - {})".format(x))
    return "*".join(out)

for k in range(M):
    poly = []
    last = None
    cnt  = 0
    for j in range(k):
        if z[j] == last:
            cnt += 1
        else:
            if last is not None:
                poly.append(pow_poly(last, cnt))
            last = z[j]
            cnt  = 1
    if last is not None and cnt>0:
        poly.append(pow_poly(last, cnt))
    txt = "*".join([p for p in poly if p!="1"])
    print("p{}(x) = {}".format(k, txt if txt else "1"))


# === CHANGE #1: list of b_k definitions ===
print("\nList of b_k definitions:")
bracket_strs = []
for k in range(M):
    parts = []
    for i in range(n):
        if k >= s[i]:
            cnt = k - s[i] + 1
            if cnt > m_vals[i]:
                cnt = m_vals[i]
            parts.append("x{},{}".format(i, cnt))
    bs = "[" + "; ".join(parts) + "; f]"
    bracket_strs.append(bs)
    print("b{} = {}".format(k, bs))



print("\nCompute the remaining part manually")
