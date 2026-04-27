# Quick helper: sum_powers(X, k): sum of x_j^k

def sum_powers(x, k):
    s = 0.0
    for xi in x:
        s += xi**k
    return s

def sum_powers_weighted(x, y, k):
    s = 0.0
    for j in range(len(x)):
        s += y[j]*x[j]**k
    return s

def input_list(prompt):
    # Reads comma-separated numbers
    lst = []
    t = input(prompt)
    items = t.split(",")
    for el in items:
        lst.append(float(el.strip()))
    return lst

def print_symbolic_equations(x, y, n):
    m = len(x)
    print("\nSYMBOLIC EQUATION SYSTEM:")
    for k in range(n+1):
        row = ""
        for i in range(n+1):
            row += "a{}(".format(i)
            ss = []
            for j in range(m):
                ss.append("{}^{}".format(int(x[j]), i+k))
            row += "+".join(ss)
            row += ")"
            if i < n:
                row += " + "
        row += " = "
        rhs = []
        for j in range(m):
            rhs.append("{}*{}^{}".format(int(y[j]), int(x[j]), k))
        row += " + ".join(rhs)
        print(row)

def print_numeric_equations(x, y, n):
    m = len(x)
    print("\nNUMERIC EQUATION SYSTEM:")
    for k in range(n+1):
        row = ""
        for i in range(n+1):
            total = sum_powers(x, i+k)
            row += "{}a{}".format(int(total), i)
            if i < n:
                row += " + "
        # computed right-hand side
        rhs = sum_powers_weighted(x, y, k)
        row += " = {}".format(int(rhs))
        print(row)

def gauss(a, b):
    # Simple Gaussian elimination with in-place updates
    n = len(b)
    for k in range(n):
        # find max pivot in the column
        max_row = k
        for i in range(k+1, n):
            if abs(a[i][k]) > abs(a[max_row][k]):
                max_row = i
        # swap rows
        a[k], a[max_row] = a[max_row], a[k]
        b[k], b[max_row] = b[max_row], b[k]
        # elimination
        for i in range(k+1, n):
            c = a[i][k]/a[k][k]
            for j in range(k, n):
                a[i][j] -= c*a[k][j]
            b[i] -= c*b[k]
    # back substitution
    x = [0.0]*n
    for i in range(n-1, -1, -1):
        s = b[i]
        for j in range(i+1, n):
            s -= a[i][j]*x[j]
        x[i] = s/a[i][i]
    return x

def solve_coefficients(x, y, n):
    # Build the system matrix and right-hand side
    mat = []
    rhs = []
    for k in range(n+1):
        row = []
        for i in range(n+1):
            row.append(sum_powers(x, i+k))
        mat.append(row)
        rhs.append(sum_powers_weighted(x, y, k))
    coeffs = gauss(mat, rhs)
    return coeffs

# -------------------- MAIN PROGRAM ---------------

x = input_list("x: ")
y = input_list("y: ")
n = int(input("degree: "))

print_symbolic_equations(x, y, n)
print_numeric_equations(x, y, n)

coeffs = solve_coefficients(x, y, n)
print("\nCOEFFICIENTS")
for i in range(n+1):
    print("a{} = {:.8g}".format(i, coeffs[i]))