import numpy as np
from numpy.linalg.linalg import matrix_rank

def sum(A, k):
    s = 0.0
    for j in range(k+1, A.shape[0]):
        s += (A[j][k] ** 2)
    return s

def step11_sum(A, v, j, k):
    s = 0.0
    for i in range(k+1, A.shape[0]):
        s += (A[j][i] * v[i])
    return s

# https://johnfoster.pge.utexas.edu/numerical-methods-book/LinearAlgebra_EigenProblem2.html
def householder(A):
    matrix_size = A.shape[0]
    v,u,z = np.zeros(matrix_size), np.zeros(matrix_size), np.zeros(matrix_size)

    for k in range(matrix_size-2):
        if A[k+1][k] == 0:
            alfa = -np.sqrt(sum(A, k))
        else:
            alfa = -np.sign(A[k+1][k]) * np.sqrt(sum(A,k))

        r = alfa ** 2 - alfa * A[k+1][k]
        v[k] = 0.0
        v[k+1] = A[k+1][k] - alfa

        for j in range(k+2, matrix_size):
            v[j] = A[j][k]
        for j in range(k, matrix_size):
            u[j] = (1.0/r) * step11_sum(A, v, j, k)
        for j in range(k, matrix_size):
            z[j] = u[j] - ((np.dot(u,v))/(2.0 * r) * v[j])

        for l in range(k+1, matrix_size-1):
            for j in range(l+1, matrix_size):
                A[j][l] = A[j][l] - v[l]*z[j] - v[j]*z[l]
                A[l][j] = A[j][l]
            A[l][l] = A[l][l] - (2.0 * v[l]*z[l])

        A[-1][-1] = A[-1][-1] - 2.0 * v[-1] * z[-1]
        for j in range(k+2, matrix_size):
            A[k][j] = 0.0
            A[j][k] = 0.0
        
        A[k+1][k] -= v[k+1] * z[k]
        A[k][k+1] = A[k+1][k]

#http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad03.pdf
#str 45
def Givens_matrix(A, i):
    I = np.identity(A.shape[0])
    x = np.array(A[i])
    c = x[i] / np.sqrt(x[i] ** 2 + x[i+1] ** 2)
    s = x[i+1] / np.sqrt(x[i] ** 2 + x[i+1] ** 2)

    I[i][i] = c
    I[i][i+1] = s
    I[i+1][i+1] = c
    I[i+1][i] = -s

    return I

# http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad05.pdf
# str 26
def QR(A):
    A_prime = A.copy()
    for i in range(A.shape[0] - 1):
        G = Givens_matrix(A_prime, i)
        A_prime = G @ A_prime @ G.T

    return A_prime

def find_eigens(A, iterations = 100):
    A_prime = A.copy()
    for _ in range(iterations):
        A_prime = QR(A_prime)

    return A_prime

if __name__ == "__main__":
    
    A = np.array(
        [
            [19/12, 13/12, 5/6, 5/6, 13/12, -17/12],
            [13/12, 13/12, 5/6, 5/6, -11/12, 13/12],
            [5/6,   5/6,   5/6, -1/6, 5/6,   5/6],
            [5/6,   5/6,   -1/6, 5/6, 5/6,   5/6],
            [13/12, -11/12, 5/6, 5/6, 13/12, 13/12],
            [-17/12, 13/12, 5/6, 5/6, 13/12, 19/12]
        ]
    )

    print(f"Matrix z zadania:\n{A}\n")
    print()

    #korzystam z biblitecznej funkcji do pozniejszego sprawdzenia wyniku
    D, V = np.linalg.eig(A)

    #sprowadzam macierz do postaci trojdiagonalnej
    householder(A)

    for i in range(A.shape[0]):
        for j in range(A.shape[0]):
            A[i][j] = round(A[i][j], 10)

    print(f"Matrix z zadania w postaci trojdiagonalnej:\n{A}\n")

    A_prime = find_eigens(A)

    for i in range(A.shape[0]):
        for j in range(A.shape[0]):
            A_prime[i][j] = round(A_prime[i][j], 10)

    print(f"Matrix z eigenvalues policzona za pomoca obrotow Givensa i QR\n{A_prime}\n")

    eigenvalues = np.diagonal(A_prime)

    print(f"Eigenvalues z uzyciem numpy = {D}\n")
    print(f"Eigenvalues z uzyciem obrotow Givensa + QR = {eigenvalues}\n")