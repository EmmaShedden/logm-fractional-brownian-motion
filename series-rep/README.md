# Series Representation of Fractional Brownian Motion

## Mathematical Fundation
Source: [_A mini-course on Wavelets and Fractional Processes_, Antoine Ayache](http://math.univ-lille1.fr/~ayache/COURSE-WavFrac.pdf)

The series representation of fractional Brownian motion is the generalization of that of stantion Brownian motion [here](https://gitlab.eecs.umich.edu/logm/wn23/fractional-brownian-motion/standard-brownian-motion/-/tree/main/series-representation) 

We use the Meyer wavelet for simulation.
The Meyer mother wavelet is defined as follows:

Let $$\psi_1(t) = \frac{\frac{4}{3\pi}(t-1/2)\cos{(\frac{2\pi}{3}(t-1/2))} - \frac{1}{\pi}\sin{(\frac{4\pi}{3}(t-1/2))}}{(t-1/2) - \frac{16}{9}(t-1/2)^3}$$
and 
$$\psi_2(t) = \frac{\frac{8}{3\pi}(t-1/2)\cos{(\frac{8\pi}{3}(t-1/2))} + \frac{1}{\pi}\sin{(\frac{4\pi}{3}(t-1/2))}}{(t-1/2) - \frac{64}{9}(t-1/2)^3}.$$
Then, the Meyer mother wavelet is
$$\psi(t) = \psi_1(t) + \psi_2(t).$$

We have the following theorem:
### Theorem 3.1
Let $\{\varepsilon_{j, k} : (j, k) \in \mathbb{Z} \times \mathbb{Z}\}$ be the sequence of i.i.d. real-valued standard Gaussian random variables defined as
$$\varepsilon_{j,k} = \int_{\mathbb R} 2^{j/2}\psi(-2^{j}t-k) ~ dW(\xi).$$
Then $\{B^H(t) : t \in \mathbb {R}\}$ can be represented as 
$$B^H(t) = \sum_{j, k \in \mathbb {Z}} 2^{-jH}\varepsilon_{j,k}(\Psi_H(2^jt-k) - \Psi_H(-k)).$$
where the series is almost surely uniformly convergent in $t$ on each compact subset of $\mathbb{R}$.

Here, $\Psi_H$ is the left-sided fractional primitive of order $H + \frac{1}{2}$ of $\psi$:
$$\Psi_H (x) = \int_{\mathbb{R}} e^{ix \xi} \frac{\hat{\psi}(\xi)}{(i \xi)^{H + 1/2}} d \xi.$$

## Simulation

### Remark
- We only simulate $B^H(t)$ for finite summands.

### Algorithm
- Representing $\varepsilon_{j,k}$ as a standard Gaussian random variable.
- Use the `scipy.integrate` module in Python to evaluate the left-sided fractional primitive and simulate fractional Brownian motion.
