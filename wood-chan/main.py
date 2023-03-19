import numpy as np
from numpy.fft import fft, ifft
import cmath
import matplotlib.pyplot as plt
import pandas as pd

# calculate rho_H(n)
# n : integer in 0..N-1
def rho(n, H):
    assert(type(n) is int and n >= 0)
    n = np.cast["float64"](n)
    if n == 0:
        return 1.0
    return (-2*n**(2*H) + (n+1)**(2*H) + (n-1)**(2*H))/2

# calculate c_0, ..., c_{M-1} where M = 2(N-1)
# N : integer in R^+
def circulant_coefficients(N, M, H):
    assert(N > 0 and type(M) is int)
    lower = np.array([rho(n, H) for n in range(N)], dtype=np.float64) # N of these
    upper = np.array([rho(M-n, H) for n in range(N, M)], dtype=np.float64) # N - 2 of these
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

    fBm = np.cumsum(increments) # xi_1, ..., xi_N
    return np.concatenate(([0.0], fBm))

# plot all the fBms next to each other
def plot(Hs, Ns, Ts, fBms):
    assert(len(Hs) == len(Ns) and len(Ns) == len(Ts) and len(Ts) == len(fBms))
    numrows = len(Ns)
    fig, axs = plt.subplots(numrows, 1, figsize=(14, 6*numrows), sharex=True)
    
    for n in range(numrows):
        axs[n].plot(Ts[n], fBms[n], label="H={:.2f}, N={}".format(Hs[n], Ns[n]))
        axs[n].legend(fontsize = 18)

    fig.suptitle("Paths of fBm for different values of H", fontsize=24)
    fig.subplots_adjust(wspace=0.1, hspace=0.1, top=0.95, bottom=0.05)

    for ax in axs.flat:
        ax.set_xlabel("Time (over the unit interval)", fontsize=18)
        ax.set_ylabel("Particle position", fontsize=18)
    for ax in axs.flat:
        ax.label_outer()

    plt.savefig("fBm.png")

def main():
    epsilon = 0.2 # distance between Hs
    # q : int in (0, infty)
    q = 10

    N = 2**q + 1 # 1/delta where delta is the granularity of the discretized time interval
    M = 2**(q+1) #2(N-1)
    delta = 1/float(N)

    # H : float in (0, 1)
    Hs = np.arange(epsilon, 1, epsilon)
    n = len(Hs) # number of trials
    Ns = np.full((n, ), N) # for now we won't vary this parameter

    Ts = np.array([np.arange(0, 1+delta, delta) for k in range(n)])

    fBms = np.empty((n, N+1, ))

    for k in range(n):
        H = Hs[k]
        assert(0 < H and H < 1 and type(N) is int and 0 < N)
        
        fBms[k] = simulate(N, M, H)

    plot(Hs, Ns, Ts, fBms)

if __name__ == "__main__":
    main()