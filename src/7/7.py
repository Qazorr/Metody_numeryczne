import numpy as np

# http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad06.pdf
# strona 14
def solve(x, fx):
    N = x.shape[0]
    #macierz Vandermonde'a
    A = np.array(
        [ [x[j] ** i for i in range(N-1, -1, -1)] for j in range(N) ]
    )    

    coefficients = np.linalg.solve(A, fx)
    return coefficients

if __name__ == "__main__":
    x = np.array([0.062500, 0.187500, 0.312500, 0.437500, 0.562500, 0.687500, 0.812500, 0.937500])
    fx = np.array([0.687959, 0.073443, -0.517558, -1.077264, -1.600455, -2.080815, -2.507266, -2.860307])

    coefficients = solve(x,fx)

    print("Wezel\t\tWartosc funkcji w wezle\n")
    for i in range(x.shape[0]):
        print(f"{x[i]}\t||\t{fx[i]}")
    print()

    with open("7.data", "w") as out:
        for i in range(coefficients.shape[0]):
            out.write(str(coefficients[i]) + "\n")

    print("Wspolczynniki wielomianu interpolowanego za pomoca macierzy Vandermonde'a")
    for i in range(coefficients.shape[0]):
        print(f"a{i} = {coefficients[i]}")

    polynomial = np.poly1d(coefficients)
    print(f"\nWielomian ma postac:\n\n{polynomial}")