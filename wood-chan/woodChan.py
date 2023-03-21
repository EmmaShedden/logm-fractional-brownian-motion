from scipy.fft import fft, ifft
import numpy as np

#covariance function
def rho(n, h):
    return 0.5*(pow(n+1, 2*h) + pow(n-1, 2*h) - 2*pow(n, 2*h))

#wood chan method as a function
def woodChan(q, h, time):
    n = pow(2, q) + 1
    m = 2*(q-1)
    c = list()
    for k in range(m):
        if k == 0:
            c.append(1)
        if k < q:
            c.append(rho(k, h))
        else:
            c.append(rho(m - k, h))
    lambdas = pow(np.real(fft(c)), 0.5)
    zetas = np.random.normal(0, 1, q)
    qstars = np.real(ifft(zetas))
    for i in range(len(qstars)):
        qstars[i] = qstars[i]*lambdas[i]
    noise = np.real(fft(qstars))
    increments = noise
    for x in increments:
        x = x*pow(time/n, h)
    values = list()
    for x in range(len(increments)):
        csum = 0
        for i in range(x):
            csum = csum + increments[i]
        values.append(csum)
    return values



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    v = woodChan(10, 0.2, 1)
    print(v)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
