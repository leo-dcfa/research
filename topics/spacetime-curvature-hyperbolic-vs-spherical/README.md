# Spacetime curvature — hyperbolic vs spherical

The two ways a space can be curved, side by side: **spherical** ($K > 0$, geodesics converge, angles add up to more than $\pi$) vs **hyperbolic** ($K < 0$, geodesics run apart, angles add up to less).

**The one-sentence version:** the sign of the curvature is the sign of the answer to *"do initially parallel straight lines get closer or further apart?"* — everything else (angle sums, circle circumference, volume growth, cosmological fate) follows from that.

---

## Part 1 — The three model geometries

Constant curvature $K$, in $n$ dimensions. Write $\ell = 1/\sqrt{\lvert K \rvert}$ for the **curvature radius**.

| | Spherical $S^n$ | Euclidean $\mathbb{E}^n$ | Hyperbolic $H^n$ |
|---|---|---|---|
| Curvature | $K = +1/\ell^2$ | $K = 0$ | $K = -1/\ell^2$ |
| Model | round sphere radius $\ell$ | $\mathbb{R}^n$ | upper sheet of $-t^2+\vec x^2 = -\ell^2$ |
| Triangle angles | $> \pi$ | $= \pi$ | $< \pi$ |
| Parallel postulate | *no* parallels | exactly one | infinitely many |
| Geodesics | great circles, all meet twice | straight lines | never meet twice |
| Circumference of radius-$r$ circle | $2\pi\ell\sin(r/\ell)$ | $2\pi r$ | $2\pi\ell\sinh(r/\ell)$ |
| Volume growth | finite, closes up | polynomial $r^n$ | exponential $e^{(n-1)r/\ell}$ |
| Global shape | compact, finite volume | infinite | infinite, "more room than you can fill" |
| Isometry group | $SO(n{+}1)$ | $ISO(n)$ | $SO(1,n)$ |
| Geodesic flow | integrable, refocusing | integrable | **chaotic**, Anosov, mixing |

### One metric, all three cases

$$ds^2 = d\chi^2 + s_k(\chi)^2\, d\Omega_{n-1}^2, \qquad s_k(\chi) = \begin{cases} \ell\sin(\chi/\ell) & k=+1 \\ \chi & k=0 \\ \ell\sinh(\chi/\ell) & k=-1\end{cases}$$

Equivalently, with the areal radius $r = s_k(\chi)$:

$$ds^2 = \frac{dr^2}{1 - k r^2/\ell^2} + r^2 d\Omega_{n-1}^2$$

- Everything below is one identity — $\sin \to \sinh$ — applied over and over. Analytic continuation $\ell \to i\ell$ maps spherical $\leftrightarrow$ hyperbolic.

### Curvature tensors (maximally symmetric spaces)

$$R_{abcd} = K\,(g_{ac}g_{bd} - g_{ad}g_{bc}), \qquad R_{ab} = (n-1)K\,g_{ab}, \qquad R = n(n-1)K$$

- These are the *only* spaces with an isotropy at every point — the curvature tensor is built from the metric alone, with one number $K$.
- In 2D the whole Riemann tensor is one function: $R = 2K$ (Gauss curvature). Curvature only becomes "shape-ful" (Weyl, tidal distortion, gravitational waves) in $n \ge 4$.

---

## Part 2 — What the sign actually does

### Angle excess / defect

$$\alpha + \beta + \gamma - \pi = K \cdot A_{\triangle}$$

- Sphere: a triangle with three right angles has excess $\pi/2$ and area $\pi\ell^2/2$ — an octant.
- Hyperbolic: the defect is *bounded*, so triangle area is bounded by $\pi\ell^2$ no matter how far the vertices are. There is a largest triangle. Similar-but-different-sized triangles don't exist — in $H^n$, shape determines size.
- Global version (2D, Gauss–Bonnet): $\int_M K\, dA + \oint_{\partial M} k_g\, ds = 2\pi\chi(M)$. Curvature integrates to topology.

### Geodesic deviation — the definition of curvature

Separation $\xi$ between neighbouring geodesics obeys the **Jacobi equation**:

$$\frac{D^2\xi}{ds^2} = -K\,\xi \qquad \Rightarrow \qquad \xi(s) \propto \begin{cases}\sin(s/\ell) & K>0 \\ s & K = 0 \\ \sinh(s/\ell) & K<0\end{cases}$$

- **$K>0$ focuses.** $\xi$ returns to zero at $s = \pi\ell$ — a *conjugate point*. Geodesics stop being length-minimising past it. Bonnet–Myers: $K \ge 1/\ell^2$ everywhere forces the space to be compact with diameter $\le \pi\ell$. Positive curvature closes things up.
- **$K<0$ defocuses, exponentially.** No conjugate points ever (Cartan–Hadamard: $H^n$ is diffeomorphic to $\mathbb{R}^n$, every pair of points joined by a unique geodesic). Nearby geodesics separate as $e^{s/\ell}$ — a Lyapunov exponent of $1/\ell$, which is why negatively curved billiards/surfaces are the textbook chaotic systems.

### Holonomy

Parallel transport a vector around a closed loop of area $A$: it comes back rotated by $\Delta\theta = K A$ — counterclockwise on the sphere, clockwise in $H^2$. Curvature *is* that rotation per unit area.

---

## Part 3 — From curved space to curved spacetime

Two things get conflated constantly. Keep them apart:

| | **Spatial** curvature | **Spacetime** curvature |
|---|---|---|
| Object | $^{(3)}R_{ijkl}$ of a constant-time slice | $R_{\mu\nu\rho\sigma}$ of the 4D manifold |
| Slicing-dependent? | **yes** — a different time slicing gives different spatial curvature | no, it's a tensor on spacetime |
| Zero means | slices are flat $\mathbb{E}^3$ | genuinely no gravity (Minkowski) |

- A flat-space FLRW universe ($k=0$) has **flat spatial slices and curved spacetime**. Expansion itself is spacetime curvature.
- Conversely, curvature can hide in the time direction: Newtonian gravity is almost entirely $g_{00} \simeq -(1+2\Phi/c^2)$.

### The Lorentzian sign flip (the big gotcha)

For a maximally symmetric **Lorentzian** spacetime, redo the deviation calculation along a *timelike* geodesic ($u\cdot u = -1$ instead of $+1$):

$$\frac{D^2\xi^\mu}{d\tau^2} = -R^\mu{}_{\nu\rho\sigma}u^\nu \xi^\rho u^\sigma = +K\,\xi^\mu$$

The sign is opposite to the Riemannian case. So the intuition inverts:

| | de Sitter (dS) | Minkowski | anti-de Sitter (AdS) |
|---|---|---|---|
| $\Lambda$ | $> 0$ | $0$ | $< 0$ |
| $K = \Lambda/3$, $R = 4\Lambda = 12K$ | $+1/\ell^2$ | $0$ | $-1/\ell^2$ |
| "Spherical or hyperbolic?" | spherical | flat | hyperbolic |
| Nearby timelike geodesics | $\xi \sim \cosh(\tau/\ell)$ — **fly apart** | $\xi \sim \tau$ | $\xi \sim \cos(\tau/\ell)$ — **refocus at $\pi\ell$** |
| Behaves like | repulsion, accelerated expansion | — | a confining box |
| Boundary | cosmological horizon at $r=\ell$ | none | timelike conformal boundary (null rays reach it in finite $t$) |
| Isometries | $SO(1,4)$ | $ISO(1,3)$ | $SO(2,3)$ |

- Static-patch metrics: $ds^2 = -(1 \mp r^2/\ell^2)dt^2 + (1 \mp r^2/\ell^2)^{-1}dr^2 + r^2 d\Omega^2$, upper sign dS, lower AdS.
- So *positive* curvature (dS) pushes worldlines apart and *negative* curvature (AdS) pulls them back — exactly backwards from the sphere/hyperbolic-plane picture. The culprit is the $-1$ in $u\cdot u$.
- **Raychaudhuri** is the general statement: $\dfrac{d\theta}{d\tau} = -\dfrac{\theta^2}{3} - \sigma_{\mu\nu}\sigma^{\mu\nu} + \omega_{\mu\nu}\omega^{\mu\nu} - R_{\mu\nu}u^\mu u^\nu$. Ordinary matter makes $R_{\mu\nu}u^\mu u^\nu \ge 0$ (strong energy condition) $\Rightarrow$ focusing $\Rightarrow$ singularity theorems. $\Lambda > 0$ violates it, and that's the defocusing.

---

## Part 4 — Which one is our universe?

Friedmann equation, with $k = +1$ spherical / $0$ flat / $-1$ hyperbolic spatial slices:

$$H^2 = \frac{8\pi G}{3}\rho - \frac{kc^2}{a^2} + \frac{\Lambda c^2}{3}, \qquad \Omega_k \equiv -\frac{kc^2}{a^2H^2}, \qquad \sum_i \Omega_i + \Omega_k = 1$$

- **Watch the sign flip:** $\Omega_k > 0$ means $k = -1$ means **open/hyperbolic**; $\Omega_k < 0$ means $k=+1$ means **closed/spherical**. The two conventions are opposite, and this trips everyone.
- Measurement (Planck 2018 TT,TE,EE+lowE+lensing+BAO): $\Omega_k = 0.0007 \pm 0.0019$. Flat to ~0.2%, curvature radius $\gtrsim 100$ Gpc — much bigger than the observable universe.
- **How it's measured:** curvature acts as a lens on the whole sky. The CMB sound horizon $r_s \approx 144$ Mpc is a known ruler at $z_\ast \approx 1090$; its observed angular size $\theta_s = r_s/D_M$ shifts because $D_M = \frac{c}{H_0\sqrt{\Omega_k}}\sinh\!\big(\sqrt{\Omega_k}H_0 D_C/c\big)$ (with $\sinh \to \sin$ for $\Omega_k<0$). Spherical universes **magnify** the CMB spots, hyperbolic ones shrink them.
- **The geometric degeneracy:** $\theta_s$ is a single number, so curvature trades against $H_0$ at fixed physical densities — a closed model with $\Omega_k = -0.05$ fits the same acoustic scale with $H_0 \approx 52$. That's why the tight bound is always "CMB **+ BAO**": you need a second distance probe at low $z$.
- **Why flat:** inflation blows $\ell$ up by $e^{60}$, driving $\Omega_k \to 0$ exponentially. Flatness is a prediction, not a coincidence — and $\lvert\Omega_k\rvert \sim 10^{-4}$ would be the smoking gun if ever detected.
- Curvature $\ne$ topology: a flat universe can still be finite (a 3-torus). Observations bound the curvature, not the global shape.

---

## Part 5 — Where each one shows up

- **Spherical:** closed FLRW; de Sitter (our accelerating late universe, and inflation); the 2-spheres inside Schwarzschild; the interior of a static star.
- **Hyperbolic:** open FLRW; AdS (the backbone of AdS/CFT — the boundary is where the dual CFT lives); constant-negative-curvature slices inside a black hole; the spatial geometry of a bubble universe formed by tunnelling; Regge/hyperbolic lattices in condensed matter.
- **Neither, exactly:** real spacetime, which is only *locally* approximated by constant curvature. Schwarzschild has $R_{\mu\nu}=0$ but $R_{\mu\nu\rho\sigma}\ne 0$ — pure Weyl, tidal stretching along the radius and squeezing transverse to it, with no analogue in constant-curvature models.

---

## Code

| File | What it does |
|---|---|
| `constant_curvature.py` | sympy: metric $\to \Gamma \to$ Riemann $\to$ Ricci $\to R$; verifies $R_{abcd}=K(g_{ac}g_{bd}-g_{ad}g_{bc})$ for $S^2$, $H^2$, and the $k=\pm1,0$ FLRW slices; gets $R = \pm 12/\ell^2$ for dS/AdS |
| `geodesic_deviation.py` | integrates the Jacobi equation both ways ($\sin$/$\sinh$ Riemannian, $\cos$/$\cosh$ Lorentzian), plus circle circumference and triangle angle sums vs radius |
| `cosmology_curvature.py` | $\Omega_k$ $\to$ angular size of the CMB sound horizon; how much curvature the acoustic peak position can hide |

```bash
uv run topics/spacetime-curvature-hyperbolic-vs-spherical/constant_curvature.py
```

## References

- do Carmo, *Riemannian Geometry* — Jacobi fields, comparison theorems (ch. 5, 7).
- Wald, *General Relativity* §5.1–5.2 (FLRW, curvature of slices), §6 (dS/AdS).
- Hawking & Ellis §4.4 — Raychaudhuri, conjugate points, focusing.
- Planck 2018 results VI, *Cosmological parameters* — $\Omega_k$ constraint.
- Thurston, *Three-Dimensional Geometry and Topology* — why hyperbolic dominates 3-manifolds.
