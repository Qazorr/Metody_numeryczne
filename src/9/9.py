import numpy as np

# http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad06.pdf
# strona 34
def A(xj, xj1):
    #xj1 - x
    licz = np.poly1d([-1, xj1])
    mian = xj1 - xj
    return licz/mian

def B(xj, xj1):
    licz = np.poly1d([1, -xj])
    mian = xj1 - xj
    return licz/mian

def C(xj, xj1):
    first = A(xj, xj1)**3 - A(xj, xj1)
    second = (xj1 - xj)**2
    return 1/6 * first * second

def D(xj, xj1):
    first = B(xj, xj1)**3 - B(xj, xj1)
    second = (xj1 - xj)**2
    return 1/6 * first * second

#funkcja z zadania 1
def thomas(a,b,c,d):
    beta, gamma = [], []
    n = len(d)

    beta.append(c[0]/b[0])
    gamma.append(d[0]/b[0])

    for i in range (1, n):
        beta.append(c[i] / (b[i] - a[i]*beta[i-1]))
        gamma.append((d[i] - a[i] * gamma[i-1]) / (b[i] - a[i] * beta[i-1]))
    
    x = np.zeros(n)

    x[-1] = gamma[-1]
    for i in range(n-2, -1, -1):
        x[i] =  gamma[i] - beta[i] * x[i+1]

    return x

# http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad06.pdf
# strona 38
def find_zeta(x, fx, h):
    d = np.zeros(x.shape[0] - 2)
    
    zeta = np.array([0])

    #tworze macierz A (diag = 4, diag ponizej = 1, diag powyzej = 1)
    a = [0, 1, 1, 1, 1, 1]
    b = [4 for i in range(6)]
    c = a[::-1]

    #wyliczam b
    mul = 6/(h**2)
    for i in range(x.shape[0]-2):
        val = fx[i] - 2 * fx[i+1] + fx[i+2] 
        d[i] = mul * val

    print(d)
    #macierz jest trojdiagonalna, korzystam z algorytmu z zadania (1)
    solved = thomas(a,b,c, d)


    zeta = np.append(zeta, solved)
    zeta = np.append(zeta, [0])

    return zeta

#yj(x) = A*f(j) + B*f(j+1) + C*ξ(j) + D*ξ(j+1) ,
def splajn(x, fx, zeta):
    polynomials = []

    for i in range(x.shape[0]-1):
        yj = A(x[i], x[i+1]) * fx[i] + B(x[i], x[i+1]) * fx[i+1] + C(x[i], x[i+1]) * zeta[i] + D(x[i], x[i+1]) * zeta[i+1]
        polynomials.append(yj)
    return polynomials

if __name__ == "__main__":
    x = np.arange(-7/8, 7/8+2/8, 2/8)
    fx = np.array([1/(1+5*(val**2)) for val in x])

    zeta = find_zeta(x, fx, (x[1] - x[0]))

    print(zeta)

    polynomials = splajn(x, fx, zeta)

    i = 1
    for p in polynomials:
        print(f"{i} wielomian:\n{p}", end = "\n\n")
        i+=1

    with open("9.data", "w") as out:
        for p in polynomials:
            for coof in p.c:
                out.write(str(coof) + " ")
            out.write("\n")