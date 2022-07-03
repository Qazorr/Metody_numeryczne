import numpy as np

# http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad05.pdf
# strona 13
def metoda_potegowa(A, iterations = 250):
    b = np.zeros((A.shape[0], 1))
    #wektor z norma 1
    b[0][0] = 1

    for i in range(iterations):        
        #mnozymy A * b
        new_b = A @ b
        norm = np.linalg.norm(new_b)
        new_b = new_b / norm
        b = new_b
    return b, norm

# https://math.stackexchange.com/questions/1114777/approximate-the-second-largest-eigenvalue-and-corresponding-eigenvector-given
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

    vec, eigenvalue = metoda_potegowa(A)
    eigenvalue = round(eigenvalue, 4)
    for i in range(vec.shape[0]):
        vec[i] = np.round(vec[i], 8)

    print(f"Najwiekszy eigenvalue: {eigenvalue}")
    print(f"Vector:\n{vec}")
    
    C = A - eigenvalue * vec * vec.T
    vec, eigenvalue = metoda_potegowa(C)
    eigenvalue = round(eigenvalue, 4)
    for i in range(vec.shape[0]):
        vec[i] = np.round(vec[i], 8)
    print(f"\nDrugi najwiekszy eigenvalue: {eigenvalue}")
    print(f"Vector:\n{vec}")

