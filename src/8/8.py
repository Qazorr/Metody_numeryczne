import numpy as np

def PI(x, to_miss):
    val = 1
    for i in range(x.shape[0]):
        if i != to_miss:
            val *= (x[to_miss] - x[i])
    return val

def PI_polynomial(x, to_miss):
    polyn = np.poly1d([1])
    for i in range(x.shape[0]):
        if i != to_miss:
            #mnoze razy x - x[i]
            polyn *= np.poly1d([1, -x[i]])
    to_ret = np.array([wsp for wsp in polyn])
    return to_ret

# http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad06.pdf
# strona 16
def Lagrange(x, fx):
    m = x.shape[0]
    wsp = np.zeros(m)
    for i in range(m):
        wsp += (fx[i] * PI_polynomial(x, i)) / PI(x,i)
    return wsp



if __name__ == "__main__":
    #poprzednia wersja zadania
    prev_x = np.arange(-1, 1 + 1/32, 1/32)
    prev_y = np.array([1/(1+5*(val**2)) for val in prev_x])

    x = np.arange(-7/8, 7/8+2/8, 2/8)
    fx = np.array([1/(1+5*(val**2)) for val in x])

    coefficients = Lagrange(x, fx)

    print("Wezel\t\tWartosc funkcji w wezle\n")
    for i in range(x.shape[0]):
        print(f"{x[i]}\t||\t{fx[i]}")
    print(f"\nWspolczynniki wielomianu interpolowanego metoda Lagrange'a:")
    with open("8.data", "w") as out:   
        for i in range(len(coefficients)):
            print(f"a{i} = {coefficients[i]}")
            out.write(str(coefficients[i]) + "\n")

    polynomial = np.poly1d(coefficients)
    print(f"\nWielomian ma postac:\n\n{polynomial}")