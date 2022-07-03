import numpy as np
import math

EPSILON = 1e-7

#? Wskazowka dotyczaca kranca przedzialu
def findA():
    for i in range(1,10000):
        if math.exp(-i) < EPSILON:
            return i
    return 0

#funckja z zadania
def f(x):
    inside_sin = math.pi * (1+math.sqrt(x))/(1+x**2)
    return math.sin(inside_sin) * math.exp(-x);


# http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad07.pdf
# strona 20
def trapezoidal(start,end, subintervals):
    h = (end - start) / subintervals
    
    integration = f(start) + f(end)
    
    for i in range(1, subintervals):
        k = start + i*h #kolejny punkt 
        integration += 2 * f(k)
    
    integration = integration * h/2
    
    return integration

# http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad07.pdf
# strona 28
def Romberg(start, end, romberg_lenght = 100):
    A = np.zeros((romberg_lenght,romberg_lenght))
    for k in range(1, romberg_lenght):
        A[k,0] = trapezoidal(start, end, 2**k)

        # http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad07.pdf <- strona 30
        for j in range(k):
            A[k][j+1] = (4**(j+1) * A[k, j] - A[k-1][j]) / (4**(j+1) - 1)
    
        if abs(A[k-1][k-1] - A[k][k]) < EPSILON and k != 0:
            return A[k][k]

        print(f"∫({k+1}) = {A[k-1][k-1]}", end = "\r")

    print()
    return A[-1][-1]

if __name__ == "__main__":
    start = 0
    end = findA()

    for i in range(1, 5001):
        integral = trapezoidal(start,end,i)
        print(f"∫({i}) = {integral}", end = "\r")

    print()

    romberg = Romberg(start, end)
    # print(f"∫ value using trapezoidal = {integral}")
    print(f"∫ value using trapezoidal + Romberg = {romberg}")