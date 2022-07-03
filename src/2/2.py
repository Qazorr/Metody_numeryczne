import numpy as np
import copy as cp

#przyblizenie
EPSILON = 1e-10

# http://th-www.if.uj.edu.pl/zfs/gora/metnum18/wyklad04.pdf
# strona 4
def gauseid(A, e):

    #wstepna estymacja
    x = np.zeros(128)
    k = 0
    with open("GAUSS-SEIDEL_diff.data", "w") as GSD:            
        before = np.linalg.norm(x)
        while True:

            x[0] = (e[0] - x[1]*A[0][1] - x[4]*A[0][4]) / A[0][0]

            for i in range(1,4):
                x[i] = (e[i] - x[i-1]*A[i][i-1] - x[i+1]*A[i][i+1] - x[i+4]*A[i][i+4]) / A[i][i]

            for i in range(4, 124):
                x[i] = (e[i] - x[i-4]*A[i][i-4] - x[i-1]*A[i][i-1] - x[i+1]*A[i][i+1] - x[i+4]*A[i][i+4]) / A[i][i]

            for i in range(124, 127):
                x[i] = (e[i] - x[i-4]*A[i][i-4] - x[i-1]*A[i][i-1] - x[i+1]*A[i][i+1])/A[i][i]

            x[-1] = (e[-1] - x[-1-4]*A[-1][-1-4] - x[-1-1]*A[-1][-1-1]) / A[-1][-1]

            after = np.linalg.norm(x)
            GSD.write(str(abs(before-after)) + "\n")
            k += 1

            if abs(before - after) < EPSILON:
                break
            before = after

    return x, k

"""
a.b = dot(a,b)

Algorytm:
r1 = b - A.x1, p1 = r1
while ||rk|| > epsilon 
{
1)    a(k) = r(k)^T.r(k) / p(k)^T.A.p(k)
2)    r(k+1) = r(k) - (a(k)* A).p(k)
3)    B(k) = r(k+1)^T.r(k+1) / r(k)^T.r(k)
4)    p(k+1) = r(k+1) + B(k).p(k)
5)    x(k+1) = x(k) + a(k).p(k)
}

http://th-www.if.uj.edu.pl/zfs/gora/metnum18/wyklad04.pdf
strona 15
"""
def gradients(A, e):
    #wektor zerowy, stojacy
    x = np.zeros((128, 1))
    k = 0

#    Wstepnie:
#    r(1) = b - A.x(1)
#    p(1) = r(1) 
    r = e - np.dot(A, x)
    p = r.copy()
    with open("GRADIENTS_diff.data", "w") as GDD:            
        before = np.linalg.norm(r)
        while True:
            # 1)
            a = np.dot(r.T, r) / (np.dot(np.dot(p.T, A), p))
            
            # 2)
            previous_r = r
            r = previous_r - np.dot(a * A, p)
            
            # 3)
            B = np.dot(r.T, r) / np.dot(previous_r.T, previous_r)

            # 4)
            previous_p = p
            p = r + B * previous_p

            # 5)
            x = x + a * previous_p

            after = np.linalg.norm(r)
            k += 1

            GDD.write(str(abs(before-after)) + "\n")

            if abs(after - before) < EPSILON:
                break
            before = after
    return x, k


if __name__ == "__main__":    
    A = np.zeros((128, 128))
    e = np.ones(128)

    np.fill_diagonal(A, 4)

    for i in range(-1, 2, 2):
        a = A.diagonal(i)
        a.setflags(write=True)
        a.fill(1)
        a.setflags(write=False)

    for i in range(-4, 5, 8):
        a = A.diagonal(i)
        a.setflags(write=True)
        a.fill(1)
        a.setflags(write=False)

    x_gs, k_gs = gauseid(A, e)
    x_g, k_g = gradients(A, e)

    print(f"xi = GAUSS-SEIDEL || GRADIENTS")
    for i in range(x_gs.shape[0]):
        print(f"x{i+1} = {round(x_gs[i], 8)} || {round(x_g[i][0], 8)}")

    print(f"Number of iterations:\nGauss-Seidel = {k_gs} || Gradients = {k_g}")