def get_precision(eps):
    if eps == 0:
        return 10
    s = "{:.0e}".format(eps)  # e.g. 1e-4 -> '1e-04'
    exponent = int(s.split('e')[-1])
    return abs(exponent) + 3

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
            coeffs.append(float(p) / float(q))
        else:
            coeffs.append(float(t))
    return coeffs

def poly_val(coeffs, x):
    res = 0.0
    n = len(coeffs)
    for i in range(n):
        res += coeffs[n-1-i] * x**i
    return res

def poly_deriv(coeffs):
    n = len(coeffs)
    result = []
    for i in range(n-1):
        result.append(coeffs[i]*(n-1-i))
    return result

def bisection(coeffs, a, b, eps, precision):
    print("\n== Bisection method ==")
    fmt = f"{{:.{precision}f}}"
    n = 1
    fa = poly_val(coeffs, a)
    fb = poly_val(coeffs, b)
    if fa*fb > 0:
        print("The function does not change sign on this interval!")
        return None
    while abs(b-a) > eps:
        c = (a+b)/2.0
        fc = poly_val(coeffs, c)
        print(f"Iteracja {n}: a={fmt.format(a)}, b={fmt.format(b)}, c={fmt.format(c)}, W(a)={fmt.format(fa)}, W(b)={fmt.format(fb)}, W(c)={fmt.format(fc)}")
        if abs(fc) < eps:
            print("Target accuracy reached at point c.")
            return c
        if fa*fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
        n += 1
    x = (a + b) / 2.0
    print(f"Done: x={fmt.format(x)}, W(x)={fmt.format(poly_val(coeffs, x))}")
    return x

def regula_falsi(coeffs, a, b, eps, precision):
    print("\n== Regula falsi method ==")
    fmt = f"{{:.{precision}f}}"
    n = 1
    fa = poly_val(coeffs, a)
    fb = poly_val(coeffs, b)
    if fa*fb > 0:
        print("The function does not change sign on this interval!")
        return None
    while True:
        print(f"-- Iteracja {n} --")
        print(f"a = {fmt.format(a)}, b = {fmt.format(b)}")
        print(f"W(a) = {fmt.format(fa)}, W(b) = {fmt.format(fb)}")
        numerator = a * fb - b * fa
        denominator = fb - fa
        if abs(denominator) < 1e-20:
            print("Denominator is near zero!")
            return None
        c = numerator / denominator
        fc = poly_val(coeffs, c)
        print(f"c = {fmt.format(numerator)} / {fmt.format(denominator)} = {fmt.format(c)}")
        print(f"W(c) = {fmt.format(fc)}")
        if abs(fc) < eps:
            print(f"Done: c = {fmt.format(c)}, W(c) = {fmt.format(fc)}")
            return c
        if fa*fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
        n += 1

def secant(coeffs, x0, x1, eps, precision):
    print("\n== Secant method ==")
    fmt = f"{{:.{precision}f}}"
    n = 1
    while True:
        fx0 = poly_val(coeffs, x0)
        fx1 = poly_val(coeffs, x1)
        delta_x = x1 - x0
        delta_fx = fx1 - fx0
        if abs(delta_fx) < 1e-20:
            print("Difference W(x1)-W(x0) is too small!")
            return None
        product = fx1 * delta_x
        x2 = x1 - product / delta_fx
        fx2 = poly_val(coeffs, x2)
        print(f"-- Iteracja {n} --")
        print(f"x0 = {fmt.format(x0)}, x1 = {fmt.format(x1)}")
        print(f"W(x0) = {fmt.format(fx0)}, W(x1) = {fmt.format(fx1)}")
        print(f"x2 = {fmt.format(x1)} - {fmt.format(product)} / {fmt.format(delta_fx)} = {fmt.format(x2)}")
        print(f"W(x2) = {fmt.format(fx2)}")
        if abs(fx2) < eps:
            print(f"Done: x2 = {fmt.format(x2)}, W(x2) = {fmt.format(fx2)}")
            return x2
        x0, x1 = x1, x2
        n += 1

def newton(coeffs, x0, eps, precision):
    deriv = poly_deriv(coeffs)
    fmt = f"{{:.{precision}f}}"
    print("\n== Newton method ==")
    n = 1
    while True:
        fx = poly_val(coeffs, x0)
        dfx = poly_val(deriv, x0)
        if abs(dfx) < 1e-20:
            print("Derivative is near zero!")
            return None
        x1 = x0 - fx / dfx
        print(f"Iteracja {n}: x[{n}]={fmt.format(x0)}, W(x[{n}])={fmt.format(fx)}, W'(x[{n}])={fmt.format(dfx)}, x[{n+1}]={fmt.format(x1)}")
        if abs(x1-x0) < eps or abs(fx) < eps:
            print(f"Done: x={fmt.format(x1)}, W(x)={fmt.format(poly_val(coeffs, x1))}")
            return x1
        x0 = x1
        n += 1

def menu():
    print("\nChoose method:")
    print("1 - bisection")
    print("2 - regula falsi")
    print("3 - secant")
    print("4 - Newton")
    return input("Method number: ").strip()

def main():
    print("Enter polynomial coefficients W(x) from highest degree, comma-separated (e.g. 2,0,-3 for 2x^2-3):")
    w_input = input("w(x): ").strip()
    coeffs = parse_poly(w_input)
    eps = float(input("Enter required accuracy (e.g. 0.0001): ").strip())
    precision = get_precision(eps)
    method = menu()

    if method == "1":
        a = float(input("Enter left interval endpoint a: ").strip())
        b = float(input("Enter right interval endpoint b: ").strip())
        bisection(coeffs, a, b, eps, precision)
    elif method == "2":
        a = float(input("Enter left interval endpoint a: ").strip())
        b = float(input("Enter right interval endpoint b: ").strip())
        regula_falsi(coeffs, a, b, eps, precision)
    elif method == "3":
        x0 = float(input("Enter first approximation x0: ").strip())
        x1 = float(input("Enter second approximation x1: ").strip())
        secant(coeffs, x0, x1, eps, precision)
    elif method == "4":
        x0 = float(input("Enter initial approximation x0: ").strip())
        newton(coeffs, x0, eps, precision)
    else:
        print("Invalid method")

main()
