import numpy as np

EPSILON = 1e-6

# http://th-www.if.uj.edu.pl/zfs/gora/metnum15/wyklad10.pdf
# strona 14
def Laggure(P: np.poly1d, z = complex(0, 0), loopbreak = 32):
    zk = complex(0, 0)
    n = P.order+1
    dx = P.deriv()
    ddx = dx.deriv()
    for i in range(loopbreak):
        a = dx(z) - np.sqrt((n-1)*((n-1) * dx(z)**2 - n * P(z)*ddx(z)))
        b = dx(z) + np.sqrt((n-1)*((n-1) * dx(z)**2 - n * P(z)*ddx(z)))
        if abs(a) > abs(b):
            zk = z - (n * P(z) / a)
        else:
            zk = z - (n * P(z) / b)

        if abs(zk - z) < EPSILON:
            return zk
        z = zk

# http://th-www.if.uj.edu.pl/zfs/gora/metnum15/wyklad10.pdf
# strona 18
def deflation(P: np.poly1d, x):
    new_P, _ = P / np.poly1d([1, -x]) 
    return new_P

# http://th-www.if.uj.edu.pl/zfs/gora/metnum15/wyklad10.pdf
# strona 20
def solve(P: np.poly1d):
    roots = np.zeros(P.order, dtype=complex)
    roots[0] = Laggure(P)
    temp_P = P

    for i in range(1, P.order):
        new_P = deflation(temp_P, roots[i-1])
        temp = Laggure(new_P)
        roots[i] = Laggure(P, temp)
        temp_P = new_P

    return roots

if __name__ == "__main__":

    p1 = np.poly1d([243, -486, 783, -990, 558, -28, -72, 16], variable='z')
    p2 = np.poly1d([1, 1, 3, 2, -1, -3, -11, -8, -12, -4, -4], variable='z')
    p3 = np.poly1d([1, 1j, -1, -1j, 1], variable='z')

    x1 = solve(p1)
    x2 = solve(p2)
    x3 = solve(p3)

    print("1:")
    print(f"Wielomian: \n{p1}\n")
    print("Miejsca zerowe:")
    for root in x1:
        print(f"{round(root, 8)}")

    print("\n2:")
    print(f"Wielomian: \n{p2}\n")
    print("Miejsca zerowe:")
    for root in x2:
        print(f"{round(root, 8)}")

    print("\n3:")
    print(f"Wielomian: \n{p3}\n")
    print("Miejsca zerowe:")
    for root in x3:
        print(f"{round(root, 8)}")