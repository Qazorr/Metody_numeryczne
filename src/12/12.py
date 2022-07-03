import numpy as np

EPSILON = 1e-08

def f(x):
    return np.cos((1+x)/(x**2 + 0.04)) * np.exp(-(x**2))

#patrzac na wykres funkcji widzimy ze mozemy ograniczyc sie do przedzialu ~ [-2;2], poniewaz reszta jest zaniedbywalnie mala
#najpierw policze starting i ending point dla calki
def start():
    interval = np.arange(-2, -10, -0.001)
    for x in interval:
        if(f(x) < EPSILON):
            return x

def end():
    interval = np.arange(2, 10, 0.001)
    for x in interval:
        if(f(x) < EPSILON):
            return x

# http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad07.pdf
# strona 20
def trapezoidal(start,end, subintervals = 1000):
    h = (end - start) / subintervals
    
    integration = f(start) + f(end)
    
    for i in range(1, subintervals):
        k = start + i*h #kolejny punkt 
        integration += 2 * f(k)
    
    integration = integration * h/2
    
    return integration

if __name__ == "__main__":
    start = start()
    end = end()

    print(f"Punkty wybrane w taki sposob aby f(x) <= epsilon\nStart = {start}\nEnd = {end}")

    I = trapezoidal(start, end, 1000)

    print(f"∫ value using trapezoidal = {I}")