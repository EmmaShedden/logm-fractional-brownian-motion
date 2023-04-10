# Fast Fourier Transformation

See my [GitHub](https://github.com/EmmaShedden/log-m/tree/master/wood-chan#applying-dft) page for a version of this README that renders correctly. Differences between the way GitHub and GitLab render Markdown equations mean that a large fraction of this page is not reader-friendly, but the GitHub version works fine.

## Mathematical Foundation
Source: [_Fast Fourier Transform_, Michel Goemans](https://math.mit.edu/~goemans/18310S15/fft-notes.pdf)

### Discrete Fourier Transformation
- The discrete Fourier transform takes as input the coefficient representation of a polynomial and outputs a point-value representation of the polynomial, where the points lie on the unit circle. In particular:
    - Let $f(x) = \displaystyle{ \sum_{k=0}^{n-1} a_k x^k }$ be a polynomial in $\mathbb{C}$ of degree $n-1$.
    - Then the input is $(a_0, ~ \cdots, ~ a_{n-1}) \in \mathbb{C}^n$.
    - DFT returns the point-value representation where the points are the $n$-th roots of unity, i.e. $(X_0, ~ \cdots, ~ X_{n-1})$ where $X_k = f(x_k)$ for $\displaystyle{ x_k = \exp(-2 \pi i \frac{k}{n}) }$, $k \in \lbrace 0, ~ \cdots, ~ n-1 \rbrace$.


- In other words, for input vector $a = (a_0, ~ \cdots, ~ a_{n-1}) \in \mathbb{C}^n$, DFT outputs $\text{DFT}(a) = X = (X_0, ~ \cdots, ~ X_{n-1}) \in \mathbb{C}^n$ such that:

$$X_j = \sum_{k=0}^{n-1} a_k \exp \left (-2 \pi i \frac j n \right )^k = \sum_{k=0}^{n-1} a_k \exp \left (-2 \pi i \frac{jk}{n} \right ) ~~~~ \forall j \in \lbrace 0, ~ \cdots, ~ n-1 \rbrace$$

- Therefore $\displaystyle{ Qa = \frac{1}{\sqrt{M}} \text{DFT}(a) }$ for any vector $a \in \mathbb{C}^{M}$.

- The inverse DFT simply takes the roots of unity in the opposite direction around the unit circle, i.e. the conjugates of the sequence of points given above; it also scales by a factor of $\frac{1}{n}$.

- Consider $a$, $X$, and $n$ from above.

```math
\begin{aligned}
(\text{DFT}^{-1} \circ \text{DFT})(a) &= \text{DFT}^{-1}(X) & \\
 &= \left ( \frac{1}{n} \sum_{j=0}^{n-1} X_j \exp \left (2 \pi i \frac{r}{n} \right )^j \right )_{r \in \lbrace 0, ~ \cdots, ~ n-1 \rbrace} & \\
 &= \frac{1}{n} \left ( \sum_{j=0}^{n-1} X_j \exp \left (2 \pi i \frac{jr}{n} \right ) \right )_{r \in \lbrace 0, ~ \cdots, ~ n-1 \rbrace} & \\
 &= \frac{1}{n} \left ( \sum_{j=0}^{n-1} \sum_{k=0}^{n-1} a_k \exp \left (-2 \pi i \frac{jk}{n} \right ) \exp \left (2 \pi i \frac{jr}{n} \right ) \right )_{r \in \lbrace 0, ~ \cdots, ~ n-1 \rbrace} & \\
 &= \frac{1}{n} \left ( \sum_{k=0}^{n-1} a_k \sum_{j=0}^{n-1} \exp \left (2 \pi i \frac{j(r-k)}{n} \right ) \right )_{r \in \lbrace 0, ~ \cdots, ~ n-1 \rbrace} & \\
 &= \frac{1}{n} \left ( \sum_{k=0}^{n-1} a_k n \delta_{kr} \right )_{r \in \lbrace 0, ~ \cdots, ~ n-1 \rbrace} & \text{by } (\dagger) \text{ below with } m = r - k \\
 &= \left ( \sum_{k=0}^{n-1} a_k \delta_{kr} \right )_{r \in \lbrace 0, ~ \cdots, ~ n-1 \rbrace} & \\
 &=  ( a_r )_{r \in \lbrace 0, ~ \cdots, ~ n-1 \rbrace} & \\
 &= a & \\
\end{aligned}
```

- Proof of $(\dagger)$:

```math
\begin{aligned}
\sum_{j=0}^{n-1} \exp \left (2 \pi i \frac{jm}{n} \right )
 &= \begin{cases} 
        \displaystyle{\sum_{j=0}^{n-1} 1^{jm / n}}, & n ~ | ~ m, \text{ i.e. } \frac m n \in \mathbb{Z} \\
         &  \\
        \displaystyle{ \frac{1 - \exp \lbrace2 \pi i m\rbrace }{1-\exp \lbrace 2 \pi i \frac{m}{n} \rbrace } }, & n ~ \nmid ~ m \\
    \end{cases} & \text{by the geometric series formula} \\
 &= \begin{cases} 
        n, & n ~ | ~ m, \\
         &  \\
        \displaystyle{ \frac{1 - 1^{m} }{1-\exp \lbrace 2 \pi i \frac{m}{n} \rbrace } }, & n ~ \nmid ~ m \\
    \end{cases} & \text{note the denominator is nonzero by assumption} \\
 &= n \delta_{n ~ | ~ m} & \text{for any } n \neq 0, m \in \mathbb{Z} \hspace2ex (\dagger)
\end{aligned}
```

- So $\displaystyle{ Q^* a = \sqrt{M} \cdot \text{DFT}^{-1}(a) }$ for any vector $a \in \mathbb{C}^{M}$.

- So $\text{DFT} \circ \text{DFT}^{-1} = \text{DFT}^{-1} \circ \text{DFT} = I_n$. This can also be seen using the relationship with $Q$ and $Q^*$ given above.

### Remark
Fast Fourier Transformation is just an improvement to Discrete Fourier Transformation, reducing time complexity from $O(N^2)$ to $O(N\log N)$.
