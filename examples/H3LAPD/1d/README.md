## 1D SOL for testing in LAPD solver

(Ed Threlfall, 3 November 2023; adapted for NESO by Owen Parry)

### Abstract

A 1D SOL problem for the LAPD solver which includes a source and sonic outflow and uses a heavily simplified equation system.

### 1D SOL solution

One species, no electric field, no transverse dynamics. Isothermal (constant *T*), like LAPD. The equations are the steady-state of equations 10 and 11 from [this report](https://github.com/ExCALIBUR-NEPTUNE/Documents/blob/main/reports/ukaea_reports/CD-EXCALIBUR-FMS0047-M2.2.2.pdf). The third equation (energy balance) isn't required because the system is isothermal.

There are then two ODEs on the domain *x* ∈ \[−1,1\]:

$$
\begin{aligned}
(n u)' &=& n^*,\\
(nu^2)' &=& - T n'.
\end{aligned}
$$

Where prime denotes *x*-derivative.
$n^*$ is a spatially-constant density source (this is what Ben Dudson says he had in Hermes-3 LAPD).
Note RHS of second equation is pressure gradient force.

The solution is obtained easily (high-school maths) as

$$
\begin{aligned}
n &= \frac{c + R}{2T},~~~~~~~(1)\\
u &= \frac{c - R}{2n^* x},~~~~~~~(2)\\
\end{aligned}
$$

where the radical $R=+\sqrt{c^2-4Tn^{*2} x^2}$ and the constant $c=2 n^*\sqrt{T}$ is chosen to give $u=\mp \sqrt{T}$ at the respective boundaries *x* =  − 1 and *x* = 1 (sonic outflow). This condition makes the radical vanish at the boundaries and the density there is $\frac{n^*}{\sqrt{T}}$. Note that it looks a bit like the solution is singular in *u* near *x* = 0 but it’s actually zero there. (Taylor expand radical to see this). Beware of this if trying to initialize the solution in a simulation.

Note the other constant of integration that entered when solving the equations was set to zero, which just corresponds to a choice of origin (the equations are translation-invariant). Note the equations could be expressed in $p ≡ nu$ instead of $u$; note the simple result that $p=n^*x$.

The directly-useable form of the equations is

$$
\begin{aligned}
n = \frac{n^*}{\sqrt{T}} \left ( 1+\sqrt{1-x^2} \right ),\\
u = \sqrt{T} \left ( \frac{1-\sqrt{1-x^2}}{x} \right ).
\end{aligned}
$$

### More general boundary conditions: subsonic outflow

Let the BCs be amended such that $n = \frac{n^*}{\delta \sqrt{T} }$ and $u = \pm \delta \sqrt{T}$ on the boundaries (so $\delta$ is really the Mach number at the boundary), where $\delta$ ∈ \[0<sup>+</sup>,1\]. Note that the momentum $nu$ on the boundaries is just $\pm1$. The solutions described by equations (1) and (2) still apply.

This corresponds to the choice of the constant above as
$c=n^* \sqrt{T} \left ( \delta + \frac{1}{\delta} \right )$.
Note that the radical is simply
$R=n^* \sqrt{T} \left ( \frac{1}{\delta} - \delta \right )$.

The central density maximum is now
$\frac{n^*}{\sqrt{T}} \left ( \delta + \frac{1}{\delta} \right )$.

Note that the difference between the density maxima (centre) and minima
(boundaries) is $\frac{n^* \delta}{\sqrt{T}}$ and so the degree of
compression of the material is less for smaller $\delta$ (i.e. away from the
sonic case).

This is the explanation of the parameter `delta` in the [session file](outflow1d.xml).
Note that $x$ in the equations above corresponds to the $z$ coordinate in our 3D domain,
so expressions for the boundary and initial conditions in the session file are modified
accordingly.

### Running the example

1. Build the H3LAPD solver using spack install (see [top-level readme](../../../README.md)).

1. Generate the mesh with

    `./scripts/geo_to_xml.sh examples/H3LAPD/1d/cuboid.geo`

2. Finally, run the example with

   `./scripts/run_eg.sh H3LAPD 1d` 