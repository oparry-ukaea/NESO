Let $c=2\pi/5$

$\phi={\rm cos}(cx){\rm cos}(cy){\rm sin}(cz/2)$



$$
\begin{align}
\nabla\cdot\nabla_\perp \phi &= \frac{\partial^2\phi}{\partial x^2} + \frac{\partial^2\phi}{\partial y^2} \\
&= -2c^2\phi \\
&= \frac{-8\pi^2}{25}\left[{\rm cos}(2\pi x/5){\rm cos}(2\pi y/5){\rm sin}(2\pi z/10)\right]
\end{align}
$$


---
Neumann BCs:
$$
\begin{align}
&= \nabla\phi\cdot\vec{n} \\
&= \left[c/2*cos(cx)*cos(cy)*cos(cz/2)\right]\cdot\vec{n}
\end{align}
$$
with $\vec{n}=(0,0,-1)$ and $\vec{n}=(0,0,1)$ at low and high-z ends respectively.

