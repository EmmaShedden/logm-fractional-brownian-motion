# Wood-Chan (Circulant Method) Algorithm

## Mathematical foundation
Source: [_Fractional Brownian motion in a Nutshell_, Georgiy Shevchenko](https://drive.google.com/file/d/1BEjP1AHJWwW1HtJDZcKPLzWJ1wXDoxcW/view)

### Problem setup
- Let $T := \{ 0, \frac{1}{N}, \frac{2}{N}, \cdots, \frac{N-1}{N}, 1 \}$ be the discretized unit time interval with $N$ increments.

- Let $t_k = t_k^N := \frac{k}{N} \in T$ for $k \in I_0 := \{0, \cdots, N\}$ denote the timestamps.

- Let $B_t^H$ for $t \in T$ be the value of an fBm with Hurst parameter $H$ at time $t$.
    - Since $B^H$ is self-stationary, $B_{t_k}^H = B_{k/N}^H = N^{-2H} B_k^H$ where $k \in I_0$.

- We know the covariance function of $B^H$ from its definition:
$$ \text{Cov}(B^H_s, B^H_t) = \mathbb{E}[B^H_s B^H_t] = \frac{1}{2}(t^{2H} + s^{2H} - |t - s|^{2H}) $$

- Let $\xi_k = \xi_k^H := B_k^H - B_{k-1}^H = N^{-2H} \left ( B_{t_k}^H - B_{t_{k-1}}^H \right )$ for $k \in I_1 := \{1, \cdots, N\}$.
    - $\xi_0 = \xi_0^H := B_0^H = 0$ with probability $1$.

- Now, $\xi := (\xi_1, \cdots, \xi_N)$ is a centered Gaussian vector:
$$ \mathbb{E}(\xi_k)
 = \mathbb{E}[B_k^H - B_{k-1}^H]
 = \mathbb{E}[B_k^H] - \mathbb{E}[B_{k-1}^H]
 = N^{-2H} \left ( \mathbb{E}[B_{t_k}^H] - \mathbb{E}[B_{t_{k-1}}^H] \right )
 = N^{-2H} \left ( 0 - 0 \right )
 = 0 ~~ \forall k \in I_1 $$

- $\xi$ represents fGn, or fractional Gaussian noise.

- It is a fact, which we do not prove here, that a Gaussian vector $x$ with mean (vector) $\mu$ and variance-covariance (matrix) $C = S S^T$ can be written as $x = \mu + S \zeta$, where $\zeta$ is a standard Gaussian vector.

    Here, $C$ will be the _circulant matrix_, and $\mu$ will be zero.

### Covariance matrix of $\xi$
- We can calculate the covariance function (matrix) of $\xi$ using its relationship to $B^H$:
$$ \begin{align*}
\text{Cov}(\xi_i, ~ \xi_j) 

 &= \text{Cov}(B_i^H - B_{i-1}^H, ~ B_j^H - B_{j-1}^H) & \text{definition of } \xi \\

 &= \mathbb{E} \left[ \left(B_i^H - B_{i-1}^H \right) \left(B_j^H - B_{j-1}^H \right) \right] - \mathbb{E} \left[B_i^H - B_{i-1}^H \right] \cdot \mathbb{E} \left[B_j^H - B_{j-1}^H \right] 
    & \text{definition of covariance} \\

 &= \mathbb{E} \left[ \left(B_i^H - B_{i-1}^H \right) \left(B_j^H - B_{j-1}^H \right) \right] - \left(\mathbb{E} \left[B_i^H \right] - \mathbb{E} \left[B_{i-1}^H \right] \right) \cdot \left(\mathbb{E} \left[B_j^H \right] - \mathbb{E} \left[B_{j-1}^H \right] \right) 
    & \text{linearity of expectation} \\
 
 &= \mathbb{E} \left[ \left(B_i^H - B_{i-1}^H \right) \left(B_j^H - B_{j-1}^H \right) \right] - (0 - 0) (0 - 0) 
    & B^H \text{ is centered by definition} \\
 
 &= \mathbb{E} \left[B_i^H B_j^H - B_i^H B_{j-1}^H - B_{i-1}^H B_j^H + B_{i-1}^H B_{j-1}^H \right] 
    & \\
 
 &= \mathbb{E} \left[B_i^H B_j^H \right]
  - \mathbb{E} \left[B_i^H B_{j-1}^H \right]
  - \mathbb{E} \left[B_{i-1}^H B_j^H \right]
  + \mathbb{E} \left[B_{i-1}^H B_{j-1}^H \right] 
    & \text{linearity of expectation} \\
 
 &=                 \frac{1}{2} \left(i^{2H} + j^{2H} - \left| i - j \right|^{2H} \right) & \\
  &\hspace{0.4cm} - \frac{1}{2} \left(i^{2H} + (j-1)^{2H} - \left| i - (j-1) \right|^{2H} \right) & \\
  &\hspace{0.4cm} - \frac{1}{2} \left((i-1)^{2H} + j^{2H} - \left| (i-1) - j \right|^{2H} \right)  & \\
  &\hspace{0.4cm} + \frac{1}{2} \left((i-1)^{2H} + (j-1)^{2H} - \left| (i-1) - (j-1) \right|^{2H} \right)
    & \text{definition of fBm} \\
 
 &= \frac{1}{2} \left(- | i - j |^{2H}
                      + | i - (j-1) |^{2H}
                      + | (i-1) - j |^{2H}
                      - | (i-1) - (j-1) |^{2H}
                \right) 
    & \\
 
 &= \frac{1}{2} \left(-2 | i - j |^{2H}
                       + | i - (j-1) |^{2H}
                       + | (i-1) - j |^{2H}
                \right) 
    & \\
 
 &= \frac{1}{2} \left(-2 | i-j |^{2H}
                       + | i-j+1 |^{2H}
                       + | i-j-1 |^{2H}
                \right) 
    & \\
 
 &= \frac{1}{2} \left(-2 | i-j |^{2H}
                       + | i-j+1 |^{2H}
                       + | i-j-1 |^{2H}
                \right)
    & (*) \\
 
 &= \frac{1}{2} \left(-2 n^{2H}
                       + (n+1)^{2H}
                       + (n-1)^{2H}
                \right)
    & \text{where } n := |i-j| \text{, } i \neq j \in I_1 \\ \\

\text{and ~ Cov}(\xi_i, ~ \xi_i)
 
 &= \frac{1}{2} \left(-2 | i-i |^{2H}
                       + | i-i+1 |^{2H}
                       + | i-i-1 |^{2H}
                \right)
    & \text{same as above, up to } (*); ~ j := i \\

 &= \frac{1}{2} \left(1^{2H} + 1^{2H} \right)
    & \\

 &= 1 & \text{where } i \in I_1
\end{align*}$$

- Observe that the covariance depends only on the difference between the timestamps, not their values; in other words, the sequence $\xi_1, \cdots, \xi_N$ is a stationary Gaussian process.
    - This is essentially a special case of the proof that fBm has stationary increments; the general version is not much different.
    - As a result, the elements of each diagonal of the covariance matrix for $\xi$ are equal.
    - Like all variance-covariance matrices, this matrix is also symmetric.

- Let $\rho_H(n) := \text{cov} \left(\xi_1, ~ \xi_{n+1} \right) = \text{cov} \left(\xi_i, ~ \xi_{j} \right) ~~ \forall i, j \in I_1$ such that $|i-j| = n$ as calculated above.
    - So all the principle diagonal elements are $\rho_H(0) = 1$, all the super-diagonal and sub-diagonal elements are $\rho_H(1)$, the next two diagonals have $\rho_H(2)$, and so on:

$$ \text{cov}(\xi) = \begin{pmatrix}
1           & \rho_H(1)   & \rho_H(2)   & \cdots & \rho_H(N-2)  & \rho_H(N-1)   \\
\rho_H(1)   & 1           & \rho_H(1)   & \cdots & \rho_H(N-3)  & \rho_H(N-2)   \\
\rho_H(2)   & \rho_H(1)   & 1           & \cdots & \rho_H(N-4)  & \rho_H(N-3)   \\
\vdots      & \vdots      & \vdots      & \ddots & \vdots       & \vdots        \\
\rho_H(N-2) & \rho_H(N-3) & \rho_H(N-4) & \cdots & 1            & \rho_H(1)     \\
\rho_H(N-1) & \rho_H(N-2) & \rho_H(N-3) & \cdots & \rho_H(1)    & 1             \\
\end{pmatrix}$$

- It is a fact, which we do not prove here, that the matrix $\text{cov}(\xi)$ is positive definite, so it has positive, real eigenvalues.

### Circulant matrix
- Let $M := 2(N-1)$.

- Let the coefficients $c_k$ for $k \in \{ 0, \cdots, M-1 \}$ be defined as:
$$\begin{align*}
c_0 &= 1, & \\
c_k &= \rho_H(k), & 1 \leq k \leq N-1  \\
c_k &= \rho_H(M-k), & N \leq k \leq M-1 \\
\text{which imply } ~~~ c_{M-k} &= \rho_H(k), & 1 \leq k \leq N-1
\end{align*}$$

- Let $C$ denote the circulant matrix on $c_0, ~ \cdots, ~ c_{M-1}$:

$$ C := \text{circ}(c_0, ~ \cdots, ~ c_{M-1}) = 
\begin{pmatrix}
 c_0        & c_1       & c_2       & \cdots    & c_{M-2}   & c_{M-1}   \\
 c_{M-1}    & c_0       & c_1       & \cdots    & c_{M-3}   & c_{M-2}   \\
 c_{M-2}    & c_{M-1}   & c_0       & \cdots    & c_{M-4}   & c_{M-3}   \\
 \vdots     & \vdots    & \vdots    & \ddots    & \vdots    & \vdots    \\
 c_2        & c_3       & c_4       & \cdots    & c_0       & c_1       \\
 c_1        & c_2       & c_3       & \cdots    & c_{M-1}   & c_0       \\
\end{pmatrix} = 
\begin{pmatrix}
 1          & \rho_H(1) & \rho_H(2) & \cdots    & \rho_H(2) & \rho_H(1) \\
 \rho_H(1)  & 1         & \rho_H(1) & \cdots    & \rho_H(3) & \rho_H(2) \\
 \rho_H(2)  & \rho_H(1) & 1         & \cdots    & \rho_H(4) & \rho_H(3) \\
 \vdots     & \vdots    & \vdots    & \ddots    & \vdots    & \vdots    \\
 \rho_H(2)  & \rho_H(3) & \rho_H(4) & \cdots    & 1         & \rho_H(1) \\
 \rho_H(1)  & \rho_H(2) & \rho_H(3) & \cdots    & \rho_H(1) & 1         \\
\end{pmatrix}
$$

- Observe that $C_{jk} = c_{k-j} = c_{|k-j|}$ for $j \leq k$ (i.e. on/above the main diagonal), and $C_{jk} = c_{M-j+k} = c_{M-|k-j|}$ for $j > k$ (i.e. below the main diagonal). But in fact, $c_r = c_{M-r} = \rho_H(r)$ for all $r \in \{ 1, ~ \cdots, ~ M-1 \}$. So $C_{jk} = c_{|k-j|} = c_{M-|k-j|}$ for all $j, ~ k \in \{ 0, ~ \cdots, ~ M-1 \}$ such that $|k-j| \not \in \{ 0, ~ M \}$ i.e. $M\not| ~~ j-k$.

- Let $Q := ( q_{jk} )_{j, k \in \{ 0, ~ \cdots, ~ M-1 \}}$ with coefficients proportional to $M$-th roots of unity:

$$ q_{jk} = \frac{1}{\sqrt{M}} \exp {\left( -2\pi i \frac{jk}{M} \right)} $$

- Observe that $Q^* =: ( q^*_{jk} ) = ( \overline{ q_{kj}} ) = ( \overline{ q_{jk}} )$ so that $\forall j, ~ k \in \{0, ~ \cdots, ~ M-1 \}$:
$$ \begin{align*}
(QQ^*)_{jk}
 &= \sum_{r=0}^{M-1} q_{jr} q^*_{rk} & \\
 &= \sum_{r=0}^{M-1} q_{jr} \overline{q_{rk}} & \\
 &= \sum_{r=0}^{M-1} \frac{1}{\sqrt{M}} \exp {\left( -2\pi i \frac{jr}{M} \right)} \frac{1}{\sqrt{M}} \exp {\left( 2\pi i \frac{kr}{M} \right)} & \\
 &= \sum_{r=0}^{M-1} \frac{1}{M} \exp {\left( 2\pi i \frac{(k-j)r}{M} \right)} & \\
 &= \begin{cases} 
        \displaystyle{\sum_{r=0}^{M-1} \frac{1}{M} \cdot 1}, & k = j, \\
         &  \\
        \displaystyle{ \frac{1}{M} \cdot \frac{1 - \exp \{2 \pi i (k-j)\} }{1-\exp \{ 2 \pi i \frac{k-j}{M} \} } }, & k \neq j \\
    \end{cases} & \text{by the geometric series formula} \\
 &= \begin{cases} 
        \displaystyle{\frac{1}{M} \cdot M}, & k = j, \\
         &  \\
        \displaystyle{ \frac{1}{M} \cdot \frac{1 - 1^{k-j} }{1-\exp \{ 2 \pi i \frac{k-j}{M} \} } }, & k \neq j \\
    \end{cases} & \text{note the denominator is nonzero by assumption} \\
 &= \delta_{jk}
\end{align*} $$

- Therefore by symmetry, $Q Q^* = Q^* Q = I_M$, and $Q$ is unitary.

- Right-multiplying $Q$ and $Q^*$, respectively, by a column vector $a = (a_0, ~ \cdots, ~ a_{M-1})^T \in \mathbb{C}^M$ gives:
$$\begin{align*}
(Qa)_j &= \sum_{k=0}^{M-1} q_{jk} a_k
 = \frac{1}{\sqrt{M}} \sum_{k=0}^{M-1} a_k \exp {\left( -2\pi i \frac{jk}{M} \right)} & \text{and} \\
(Q^* a)_j &= \sum_{k=0}^{M-1} \overline{q_{jk}} a_k
 = \frac{1}{\sqrt{M}} \sum_{k=0}^{M-1} a_k \exp {\left( 2\pi i \frac{jk}{M} \right)} & \\
\end{align*}$$

### Applying DFT
- The discrete Fourier transform takes as input the coefficient representation of a polynomial and outputs a point-value representation of the polynomial, where the points lie on the unit circle. In particular:
    - Let $f(x) = \displaystyle{ \sum_{k=0}^{n-1} a_k x^k }$ be a polynomial in $\mathbb{C}$ of degree $n-1$.
    - Then the input is $(a_0, ~ \cdots, ~ a_{n-1}) \in \mathbb{C}^n$.
    - DFT returns the point-value representation where the points are the $n$-th roots of unity, i.e. $(X_0, ~ \cdots, ~ X_{n-1})$ where $X_k = f(x_k)$ for $\displaystyle{ x_k = \exp(-2 \pi i \frac{k}{n}) }$, $ k \in \{ 0, ~ \cdots, ~ n-1 \} $.

- In other words, for input vector $a = (a_0, ~ \cdots, ~ a_{n-1}) \in \mathbb{C}^n$, DFT outputs $\text{DFT}(a) = X = (X_0, ~ \cdots, ~ X_{n-1}) \in \mathbb{C}^n$ such that:
$$ X_j
 = \sum_{k=0}^{n-1} a_k \exp \left (-2 \pi i \frac{j}{n} \right )^k
 = \sum_{k=0}^{n-1} a_k \exp \left (-2 \pi i \frac{jk}{n} \right ) ~~~~ \forall j \in \{ 0, ~ \cdots, ~ n-1 \} $$

- Therefore $\displaystyle{ Qa = \frac{1}{\sqrt{M}} \text{DFT}(a) }$ for any vector $a \in \mathbb{C}^{M}$.

- The inverse DFT simply takes the roots of unity in the opposite direction around the unit circle, i.e. the conjugates of the sequence of points given above; it also scales by a factor of $\frac{1}{n}$.

- Consider $a$, $X$, and $n$ from above.
$$\begin{align*}
(\text{DFT}^{-1} \circ \text{DFT})(a) &= \text{DFT}^{-1}(X) & \\

 &= \left ( \frac{1}{n} \sum_{j=0}^{n-1} X_j \exp \left (2 \pi i \frac{r}{n} \right )^j \right )_{r \in \{ 0, ~ \cdots, ~ n-1 \}} & \\
 
 &= \frac{1}{n} \left ( \sum_{j=0}^{n-1} X_j \exp \left (2 \pi i \frac{jr}{n} \right ) \right )_{r \in \{ 0, ~ \cdots, ~ n-1 \}} & \\
 
 &= \frac{1}{n} \left ( \sum_{j=0}^{n-1} \sum_{k=0}^{n-1} a_k \exp \left (-2 \pi i \frac{jk}{n} \right ) \exp \left (2 \pi i \frac{jr}{n} \right ) \right )_{r \in \{ 0, ~ \cdots, ~ n-1 \}} & \\
 
 &= \frac{1}{n} \left ( \sum_{k=0}^{n-1} a_k \sum_{j=0}^{n-1} \exp \left (2 \pi i \frac{j(r-k)}{n} \right ) \right )_{r \in \{ 0, ~ \cdots, ~ n-1 \}} & \\
 
 &= \frac{1}{n} \left ( \sum_{k=0}^{n-1} a_k n \delta_{kr} \right )_{r \in \{ 0, ~ \cdots, ~ n-1 \}} & \text{by } (\dag) \text{ with } m = r - k \\
 
 &= \left ( \sum_{k=0}^{n-1} a_k \delta_{kr} \right )_{r \in \{ 0, ~ \cdots, ~ n-1 \}} & \\
 
 &=  ( a_r )_{r \in \{ 0, ~ \cdots, ~ n-1 \}} & \\
 
 &= a & \\
\end{align*}$$

- Proof of $(\dag)$:

$$\begin{align*}
\sum_{j=0}^{n-1} \exp \left (2 \pi i \frac{jm}{n} \right )
 &= \begin{cases} 
        \displaystyle{\sum_{j=0}^{n-1} 1^{jm / n}}, & n ~ | ~ m, \text{ i.e. } \frac{m}{n} \in \mathbb{Z} \\
         &  \\
        \displaystyle{ \frac{1 - \exp \{2 \pi i m\} }{1-\exp \{ 2 \pi i \frac{m}{n} \} } }, & n ~ \not | ~ m \\
    \end{cases} & \text{by the geometric series formula} \\
 &= \begin{cases} 
        n, & n ~ | ~ m, \\
         &  \\
        \displaystyle{ \frac{1 - 1^{m} }{1-\exp \{ 2 \pi i \frac{m}{n} \} } }, & n ~ \not | ~ m \\
    \end{cases} & \text{note the denominator is nonzero by assumption} \\
 &= n \delta_{n ~ | ~ m} & \text{for any } n \neq 0, m \in \mathbb{Z}
\end{align*}$$

- So $\displaystyle{ Q^* a = \sqrt{M} \cdot \text{DFT}^{-1}(a) }$ for any vector $a \in \mathbb{C}^{M}$.

- So $\text{DFT} \circ \text{DFT}^{-1} = \text{DFT}^{-1} \circ \text{DFT} = I_n$. This can also be seen using the relationship with $Q$ and $Q^*$ given above.

### Theorem 6.1
- Let $\displaystyle{ \lambda_k := \sum_{j=0}^{M-1}c_j \exp(2 \pi i \frac{jk}{M}) }$ for $k \in \{ 0, ~ \cdots, ~ M-1 \}$.

- Let $\Lambda := \text{diag}(\lambda_0, ~ \cdots, ~ \lambda_{M-1})$ be the diagonal matrix with $\lambda_k$'s along the diagonal.

$$ (Q \Lambda)_{jk} = q_{jk} \lambda_k = \frac{1}{\sqrt{M}} \exp {\left( -2\pi i \frac{jk}{M} \right)} \sum_{r=0}^{M-1}c_r \exp(2 \pi i \frac{kr}{M}) $$

$$\begin{align*}
(Q \Lambda Q^*)_{js} &= \sum_{k=0}^{M-1} (Q \Lambda)_{jk} Q^*_{ks} & \\

 &= \sum_{k=0}^{M-1} (Q \Lambda)_{jk} \overline{q_{ks}} & \\

 &= \sum_{k=0}^{M-1} q_{jk} \lambda_k \overline{q_{ks}} & \\

 &= \sum_{k=0}^{M-1} \frac{1}{\sqrt{M}} \exp {\left(-2\pi i \frac{jk}{M} \right)} \lambda_k \overline{\frac{1}{\sqrt{M}} \exp {\left(-2\pi i \frac{ks}{M} \right)}} & \\

 &= \sum_{k=0}^{M-1} \frac{1}{\sqrt{M}} \exp {\left(-2\pi i \frac{jk}{M} \right)} \lambda_k \frac{1}{\sqrt{M}} \exp {\left(2\pi i \frac{ks}{M} \right)} & \\

 &= \frac{1}{M} \sum_{k=0}^{M-1} \exp {\left(2\pi i \frac{k(s-j)}{M} \right)} \lambda_k & \\

 &= \frac{1}{M} \sum_{k=0}^{M-1} \exp {\left(2\pi i \frac{k(s-j)}{M} \right)} \sum_{r=0}^{M-1} c_r \exp \left ( 2 \pi i \frac{kr}{M} \right ) & \\

 &= \frac{1}{M} \sum_{r=0}^{M-1} c_r \sum_{k=0}^{M-1} \exp {\left(2\pi i \frac{k(s-j + r)}{M} \right)} & \\

 &= \frac{1}{M} \sum_{r=0}^{M-1} c_r \cdot M \delta_{M ~ | ~ s-j+r} & \text{by } (\dag) \\

 &= \sum_{r=0}^{M-1} c_r \cdot \delta_{(s-j+r) \in \{ 0, ~ M \}} & -(M-1) \leq s-j+r \leq 2(M-1) \\

 &= \sum_{r=0}^{M-1} c_r \cdot \delta_{r \in \{ j-s, ~ M+j-s \}} & \\

 &= \begin{cases} 
        c_{M-|j-s|}, & j < s, \\
         &  \\
        c_{|j-s|}, & j \geq s \\
    \end{cases} & \\

 &= C_{js} & \forall j, ~ s \in \{ 0, ~ \cdots, ~ M-1 \} \\
 
\end{align*}$$

where the last step references the construction of the circulant matrix [above](#circulant-matrix). Therefore $Q \Lambda Q^* = C$.

- We have proved $Q$ is unitary, and by construction $\Lambda$ is diagonal. It follows that the $\lambda_k$ are the eigenvalues of C.

- By previous assertion, the $C$ is positive definite, so has positive, real eigenvalues. Therefore $\Lambda$ so has a real square root matrix constructed efficiently by raising each diagonal entry to the power of $\frac{1}{2}$. Further, all the relevant matrices are symmetric. So, 
$$C = (Q \Lambda^{1/2} Q^*) \cdot (Q \Lambda^{1/2} Q^*) = (Q \Lambda^{1/2} Q^*) \cdot (Q \Lambda^{1/2} Q^*)^*$$

- Let $S := Q \Lambda^{1/2} Q^*$. Therefore, $S$ is real. (TODO)

- Therefore $(\xi_1, ~ \cdots, ~ \xi_{N-1}, ~ \xi_{N-1}, ~ \cdots, ~ \xi_1) = S\zeta$ where $\zeta$ is standard Gaussian. (TODO)

- This can be computed efficiently using DFT: 
$$ \begin{align*}
(\xi_1, ~ \cdots, ~ \xi_{N-1}, ~ \xi_{N-1}, ~ \cdots, ~ \xi_1)

 &= S\zeta & \\

 &= Q \Lambda^{1/2} Q^* \zeta & \\

 &= \left( \frac{1}{\sqrt{M}} Q \right) \Lambda^{1/2} \left( \sqrt{M} Q^* \right) \zeta & \\

 &= \left( \frac{1}{\sqrt{M}} Q \right) \Lambda^{1/2} \text{DFT}^{-1} (\zeta) & \\

 &= \text{DFT} \left( \Lambda^{1/2} \text{DFT}^{-1} (\zeta) \right) & \\
\end{align*}$$

## Simulation results

![Paths of fBm simulated for 9 values of the Hurst parameter](/wood-chan/fBm_9_small.png)