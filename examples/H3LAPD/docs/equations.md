### Equations

Density and parallel momentum ($G_e = m_{e}n_{e}v_{||e}$) equations for electrons
$$
   \begin{align}
   \frac{\partial n_e}{\partial t} =& - \nabla\cdot\left[n_e(\mathbf{v}_{E\times B} + \mathbf{b}v_{||e}\right)]\\
      \frac{\partial G_e}{\partial t} =& - \nabla\cdot\left[G_e(\mathbf{v}_{E\times B} + \mathbf{b}v_{||e}\right)] - \partial_{||}p_e -{e}n_{e}E_{||} + m_{e}n_e\nu_{ei}(v_{||d+} - v_{||e})\\
   \end{align}
$$


Similar equation for ion momentum ($G_{d+} =~m_{d+}n_{d+}v_{||d+}$), with symmetric collision term:

$$
   \begin{align}
      \frac{\partial G_{d+}}{\partial t} =& - \nabla\cdot\left[G_{d+}(\mathbf{v}_{E\times B} + \mathbf{b}v_{||d+}\right)] - \partial_{||}p_{d+} +{e}n_{d+}E_{||} - m_{e}n_{e}\nu_{ei}(v_{||d+} - v_{||e})\\
   \end{align}
$$

Vorticity:

$$
\begin{align}
\nabla\cdot(\frac{\bar{m_{i}}\bar{n}}{B^2}\nabla_{\perp}\phi) =&~\omega \\
\frac{m_{d+}}{B^2}\nabla\cdot(n_{d+}\nabla_{\perp}\phi) =&~\omega \\
    \frac{\partial{\omega}}{\partial{t}} =& -\nabla\cdot(\omega\mathbf{v}_{E\times B}) + \nabla\cdot(n_{d+}v_{\parallel d+} - n_e v_{\parallel e} ) \\
    =& -\nabla\cdot(\omega\mathbf{v}_{E\times B}) - \nabla\cdot\left[n_{e}(v_{\parallel e} - v_{\parallel d+}) \right]
\end{align}
$$

where

$~~\mathbf{v}_{E\times B} = \frac{\mathbf{b}\times\nabla\phi}{B^2}$ is the E-cross-B drift velocity

$~~\nu_{ei}$ is the electron-ion collision rate

$
~~\nabla_{\perp} = \nabla - \mathbf{b}(\mathbf{b}\cdot\nabla) = \mathbf{b}\times\mathbf{b}\times\nabla\\
$

and we assume

$~~n_e = n_{d+}$ (quasi charge neutrality)

---

[Home](../readme)
