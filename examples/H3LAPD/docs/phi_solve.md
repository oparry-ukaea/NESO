
potential-vorticity equation for deuterium ions:
$$
\begin{align}
\nabla\cdot(\frac{\over{m_{i} n}}{B^2}\nabla_{\perp}\phi) =& \omega \\
\nabla\cdot(n_{d+}\nabla_{\perp}\phi) =& \frac{B^2 w}{m_{d+}}\\
\end{align}
$$

---

Nektar formulation of the Helmholtz equation with $\lambda=0$ and 'variable coefficients':
$$
\nabla\cdot(\mathcal{C}\nabla\phi) =~\mathrm{rhs} 
$$
where
$$
\mathcal{C}=\left[\begin{array}{ccc}
  d_{00} & d_{01} & d_{02} \\
  d_{01} & d_{11} & d_{12} \\
  d_{02} & d_{12} & d_{22}
\end{array}\right]
$$

---
For $\mathbf{b}=(0,0,1)$, we have
$$
\begin{align}
\nabla_{\perp} =&~\nabla - \mathbf{b}(\mathbf{b}\cdot\nabla)\\
=&~(\frac{\partial}{\partial x},\frac{\partial}{\partial y},0)
\end{align}
$$

To solve (2), set
$$
\begin{align}
d_{00}=d_{00} =&~n_{d+} \\
d_{01}=d_{02}=d_{12} = d_{22}=&~0 \\
\mathrm{rhs} =& \frac{B^2 w}{m_{d+}}
\end{align}
$$
which gives
$$
\begin{align}
\nabla\cdot\left(\left[\begin{array}{ccc}
  n_{d+} & 0 & 0 \\
  0 & n_{d+} & 0 \\
  0 & 0 & 0
\end{array}\right]\left(\begin{array}{c}
  \frac{\partial\phi}{\partial x} \\
  \frac{\partial\phi}{\partial y} \\
  \frac{\partial\phi}{\partial z} 
\end{array}\right)\right) =&\frac{B^2 w}{m_{d+}} \\
\nabla\cdot\left[\left(\begin{array}{c}
  n_{d+}\frac{\partial\phi}{\partial x} \\
  n_{d+}\frac{\partial\phi}{\partial y} \\
  0
\end{array}\right)\right]=&\frac{B^2 w}{m_{d+}} \\
\nabla\cdot\left(n_{d+}\nabla_{\perp}\phi\right)=&\frac{B^2 w}{m_{d+}}
\end{align}
$$

---

[Home](../readme)