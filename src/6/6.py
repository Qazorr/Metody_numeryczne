import numpy as np

# http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad05.pdf
# strona 43
# https://en.wikipedia.org/wiki/Inverse_iteration
def inverse_iteration(A, iterations = 1000):
    y = np.zeros(A.shape[0])
    y[-1] = 1

    for _ in range(iterations):
        z = np.linalg.solve(A, y)
        y = z / np.linalg.norm(z)
    return z/np.linalg.norm(z)

if __name__ == "__main__":
    A = np.array(
        [
            [2., -1., 0., 0., 1.],
            [-1., 2, 1., 0., 0.],
            [0., 1., 1., 1., 0.],
            [0., 0., 1., 2., -1.],
            [1., 0., 0., -1., 2.]
        ]
    )
    eigenvalue = 0.38197

    print(f"Matrix z zadania:\n{A}")
    
    D, V = np.linalg.eig(A)

    for i in range(V.shape[0]):
        for j in range(V.shape[0]):
            V[i][j] = round(V[i][j], 6)


    A -= np.identity(A.shape[0]) * eigenvalue
    eigenvector = inverse_iteration(A)

    for i in range(len(eigenvector)):
        eigenvector[i] = round(eigenvector[i], 6)
    
    print(f"\nEigenvector using numpy:\t     {V[:,1]}")
    print(f"Eigenvector using Inverse Iteration: {eigenvector}")