import numpy as np


"""
W zadaniu uzywam algorytmu Thomasa
czyli faktoryzacji LU dla macierzy trójdiagonalnej
ze względu na zlozonosc obliczeniowa faktoryzacji bedacej O(n)

W podpunkcie B korzystam z algorytmu Shermana Morissona
"""

#a -> diagonala pod glowna
#b -> glowna diagonala
#c -> diagonala nad glowna
# https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
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

#http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad03.pdf
#str 55
def sherman(a, b, c, d, u):
    y = thomas(a, b, c, d)
    q = thomas(a, b, c, u)
    x = y - (np.dot(u, y)/(1+np.dot(u,q))) * q
    return x

if __name__ == "__main__":
    A = np.array(
        [
        [4,1,0,0,0,0,0],
        [1,4,1,0,0,0,0],
        [0,1,4,1,0,0,0],
        [0,0,1,4,1,0,0],
        [0,0,0,1,4,1,0],
        [0,0,0,0,1,4,1],
        [0,0,0,0,0,1,4]
        ]
    )
    B = np.array(
        [
        [4,1,0,0,0,0,1],
        [1,4,1,0,0,0,0],
        [0,1,4,1,0,0,0],
        [0,0,1,4,1,0,0],
        [0,0,0,1,4,1,0],
        [0,0,0,0,1,4,1],
        [1,0,0,0,0,1,4]
        ]
    )

    print("==============A==============\n")

    print(f"Matrix:\n{A}\n")

    #zapamietuje tylko wazne wektory dla tej macierzy (A)
    a = [0, 1, 1, 1, 1, 1, 1]
    b = [4 for i in range(7)]
    c = a[::-1]
    d = [i+1 for i in range(7)]

    print(f"Vector b = {d}\n")

    x = thomas(a,b,c,d)

    print(f"Wynik dla rownania Ax = b")

    for i in range(7):
        print(f"x{i+1} = {x[i]}")

    print("\n==============B==============\n")

    u = np.array([1, 0, 0, 0, 0, 0, 1])
    v = np.array([[1], [0], [0], [0], [0], [0], [1]])

    print(f"Matrix:\n{B}\n")

    #powstaje nam trójdiagonalna macierz
    B = B - (u * v)

    print(f"Vector u = v.T = {u}\n")
    print(f"Matrix B po modyfikacji (B - (u * v)):\n{B}\n")

    #tworze wektory wazne dla tej macierzy ktore wykorzystam w algorytmie Shermana Morissona
    a_B = B.diagonal(-1)
    b_B = B.diagonal()
    c_B = B.diagonal(1)

    a_B = np.append(np.array([0]), a_B, axis = 0)
    c_B = np.append(c_B, np.array([0]))

    x = sherman(a_B, b_B, c_B, d, u)
    
    print(f"Wynik dla rownania Bx = b za pomoca algorytmu Shermana-Morissona:")

    for i in range(7):
        print(f"x{i+1} = {x[i]}")
