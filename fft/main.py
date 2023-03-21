# FFT from EECS 477
# i.e. Cooley-Tukey algorithm for powers of 2, where N_1 is always 2
# i.e. "the radix-2 case" https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#The_radix-2_DIT_case

from math import ceil, log2, cos, sin, pi
import numpy as np
import cmath


# outermost caller (wrapper for fft_recurse with pre & post processing)
# x : sequence of data points (any numeric type)
def fft(x):
    pad(x)
    return np.array(fft_recurse(x))


# helper for recursive calls
def fft_recurse(x):
    if len(x) == 1:
        x[0] = complex(x[0], 0)
    else:
        # python list splicing is cool
        x0, x1 = x[::2], x[1::2]
        d, e = fft_recurse(x0), fft_recurse(x1)

        for i in range(len(x) // 2):
            # primary nth root of unity
            omega = complex(cos(-2 * pi / len(x)), sin(-2 * pi / len(x)))
            tmp = omega**i * e[i] # omega**i * e[i]
            x[i] = d[i] + tmp
            x[i + len(x)//2] = d[i] - tmp
    
    return x


# pad input to a power of 2, plus an extra copy of all zeros
def pad(x):
    n = 2 ** ceil(log2(len(x)))
    if n != len(x):
        print("** warning: your input length is not optimal as it is not a power of 2 **")
    x += np.zeros(2*n - len(x))


# test runs
if __name__ == "__main__":
    m = 10 # number of trials per length
    p = 12 # number of array lengths to test (all powers of 2)
    tol = 1e-10 # tolerance for floating point comparisons
    for i in range(p):
        n = 2**i
        for j in range(m):
            x = (np.random.rand(n) * 100).tolist()
            assert(all(np.isclose(np.fft.fft(x), fft(x), rtol=tol, atol=tol)))
        print("{} trials passed for arrays of length {}".format(m, n))
