import numpy as np
from numpy.fft import fft, ifft
import cmath

# calculate rho_H(n)
# n : integer in 0..N-1
def rho(n, H):
    assert(type(n) is int and n >= 0)
    n = np.cast["float64"](n)
    if n == 0:
        return 1.0
    return (-2*n**(2*H) + (n+1)**(2*H) + (n-1)**(2*H))/2

# calculate c_0, ..., c_{M-1} where M = 2(N-1)
# N : integer in R+
def circulant_coefficients(N, M, H):
    assert(N > 0 and type(M) is int)
    lower = np.array([rho(n, H) for n in range(N)], dtype=np.float64)
    upper = np.array([rho(M-n, H) for n in range(N, M)], dtype=np.float64)
    return np.concatenate((lower, upper))

# given a complex array, take the real part
def get_reals(arr):
    return np.array([arr[i].real for i in range(len(arr))])

# driver for all the simulation steps
def simulate(N, M, H):
    circ_coeffs = circulant_coefficients(N, M, H)
    Lambda = get_reals(fft(circ_coeffs))
    zeta = np.random.normal(size=(M, ))
    q_star = get_reals(ifft(zeta))
    fGn = get_reals(fft(np.multiply(Lambda, q_star))[:N])

    increments = fGn * N**(-H)

    fBm = np.cumsum(increments)
    return fBm

########## extra (debugging) ##########

def print_matrix(m):
    for j in range(min(len(m[0]), 8)):
        print("---------------------", end='')
    print("\n")
    for i in range(min(len(m[0]), 8)):
        for j in range(len(m[0])):
            print(m[i][j], ",\t", sep='', end='')
        print("\n")
    for j in range(min(len(m[0]), 8)):
        print("---------------------", end='')
    print("\n\n")

def print_array(a):
    for i in range(min(len(a), 8)):
        print("---------------------", end='')
    print("\n")
    for i in range(min(len(a), 8)):
        print(a[i], ",\t", sep='', end='')
    print("\n")
    for i in range(min(len(a), 8)):
        print("---------------------", end='')
    print("\n\n")

def main():
    # H : float in (0, 1)
    H = 0.99
    # q : int in (0, \inf)
    q = 10
    
    # N : int in (0, \inf) which is one more than a power of 2
    N = 2**q + 1
    M = 2**(q+1)
    assert(0 < H and H < 1 and type(N) is int and 0 < N)
    
    print_array(simulate(N, M, H))    

if __name__ == "__main__":
    main()