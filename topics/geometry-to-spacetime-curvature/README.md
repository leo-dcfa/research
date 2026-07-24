# From zero geometry to Riemannian geometry and spacetime curvature

A ground-up build: flat-space vectors → curves → surfaces → intrinsic curvature → manifolds → tensors → connections → the Riemann tensor → Lorentzian spacetime → gravity as curvature.

**The one-sentence arc:** *distance* is the only primitive you need; everything else — angles, straight lines, curvature, gravity — is squeezed out of the metric by differentiation.

## Roadmap

| Part | Idea | Key object |
|---|---|---|
| 0 | What "geometry" means | — |
| 1 | Flat space, vectors, index notation | $\delta_{ij}$, $v^i$ |
| 2 | Coordinates ≠ geometry | $ds^2$, Jacobian |
| 3 | Curves | $\kappa$, $\tau$ (Frenet) |
| 4 | Surfaces in $\mathbb{R}^3$ | $\mathrm{I}$, $\mathrm{II}$, $K$, $H$ |
| 5 | **Theorema Egregium** | $K$ is intrinsic |
| 6 | Intrinsic 2D geometry | geodesics, angle excess, Gauss–Bonnet |
| 7 | Manifolds & tensors | $T_pM$, $\partial_\mu$, $dx^\mu$ |
| 8 | Riemannian metric | $g_{\mu\nu}$ |
| 9 | Connection, transport, geodesics | $\Gamma^\lambda_{\mu\nu}$, $\nabla_\mu$ |
| 10 | **Curvature** | $R^\rho{}_{\sigma\mu\nu}$, $R_{\mu\nu}$, $R$, $C_{\mu\nu\rho\sigma}$ |
| 11 | What curvature *does* | geodesic deviation |
| 12 | Lorentzian geometry | $\eta_{\mu\nu}$, light cones, $d\tau$ |
| 13 | Equivalence principle | free fall = geodesic |
| 14 | Einstein field equations | $G_{\mu\nu} = \tfrac{8\pi G}{c^4}T_{\mu\nu}$ |
| 15 | Worked spacetimes | Schwarzschild, FLRW, waves |

---

## Part 0 — What geometry is, in four moves

- **Euclid (~300 BC)** — geometry = axioms about points, lines, circles in *the* plane. Distance is undefined-but-obvious. The 5th (parallel) postulate is the weak joint.
- **Descartes (1637)** — put coordinates on the plane. Geometry becomes algebra: a point is $(x,y)$, distance is $\sqrt{\Delta x^2 + \Delta y^2}$.
- **Gauss (1827)** — a surface has curvature you can measure *from inside it*, with no reference to any surrounding space. Geometry is no longer about the ambient world.
- **Riemann (1854)** — generalise to $n$ dimensions and to a distance rule that varies from point to point. Geometry = a manifold + a metric. (His habilitation lecture explicitly floats the idea that physical space might be curved by matter.)
- **Einstein (1915)** — the manifold is 4D spacetime, the metric is the gravitational field, and matter is what curves it.

Two rival slogans, both useful:

- **Klein's Erlangen program**: a geometry is defined by its *symmetry group* (Euclidean = rotations + translations; affine, projective, Lorentzian…). Good for flat/homogeneous spaces.
- **Riemann's view**: a geometry is defined by an infinitesimal *distance rule* $ds^2$. Good for lumpy, curved, generic spaces — this is the one that carries into GR.

---

## Part 1 — Flat space: vectors, dot products, index notation

### Vectors and the dot product

- A vector $\vec v \in \mathbb{R}^n$ has components $v^i$ ($i = 1,\dots,n$) in some basis $\{\vec e_i\}$: $\;\vec v = \sum_i v^i \vec e_i$.
- **Everything metric** (length, angle, projection, "perpendicular") comes from one bilinear operation, the inner product:

$$\vec u \cdot \vec v = \sum_{i,j} \delta_{ij}\, u^i v^j, \qquad \lVert \vec v \rVert = \sqrt{\vec v \cdot \vec v}, \qquad \cos\theta = \frac{\vec u \cdot \vec v}{\lVert \vec u\rVert\, \lVert\vec v\rVert}$$

- In Cartesian coordinates the matrix of the inner product is the identity $\delta_{ij}$. **That is the special thing about Cartesian coordinates on flat space — nothing else.** Replace $\delta_{ij}$ with a general position-dependent $g_{ij}(x)$ and you have all of Riemannian geometry.

### Index notation (the grammar of everything below)

- **Upper index = vector component** ($v^\mu$, "contravariant"). **Lower index = covector/1-form component** ($\omega_\mu$, "covariant"). They transform oppositely, which is what makes contractions meaningful.
- **Einstein summation**: a repeated index, one up one down, is summed. $\;\omega_\mu v^\mu \equiv \sum_\mu \omega_\mu v^\mu$.
- **Free vs dummy indices**: an index appearing once on each side is *free* (the equation holds for each value); a repeated up/down pair is *dummy* (summed, name irrelevant). Both sides of an equation must have identical free indices — this is the single best typo-detector in the subject.
- Conventions used here: Greek $\mu,\nu,\dots$ run over spacetime $0\ldots3$; Latin $i,j,\dots$ over space $1\ldots3$; signature $(-,+,+,+)$; often $G = c = 1$ ("geometrized units", so mass, length and time all have units of length: $M_\odot \simeq 1.48$ km).

---

## Part 2 — Coordinates are not geometry

The single most common beginner confusion: *curvy coordinates $\neq$ curved space.*

- Flat 2D plane in Cartesian coordinates: $\;ds^2 = dx^2 + dy^2$, so $g_{ij} = \begin{pmatrix}1&0\\0&1\end{pmatrix}$.
- The **same** plane in polar coordinates $x = r\cos\theta,\ y = r\sin\theta$:

$$ds^2 = dr^2 + r^2 d\theta^2, \qquad g_{ij} = \begin{pmatrix}1&0\\0&r^2\end{pmatrix}$$

- The metric components now depend on position, and the Christoffel symbols (Part 9) are nonzero — yet the space is **exactly flat**. Curvature is what survives *after* you quotient out coordinate choices; that is precisely what the Riemann tensor computes.

### The line element

- $ds^2 = g_{\mu\nu}(x)\, dx^\mu dx^\nu$ is the **line element**: the squared distance between neighbouring points. It's the whole geometry, written locally.
- Length of a curve $x^\mu(\lambda)$:

$$L = \int \sqrt{g_{\mu\nu}\frac{dx^\mu}{d\lambda}\frac{dx^\nu}{d\lambda}}\; d\lambda$$

### How components transform

- Under $x^\mu \to x'^{\mu}(x)$, with Jacobian $J^{\mu}{}_{\nu} = \partial x^\mu/\partial x'^{\nu}$:

$$v'^{\mu} = \frac{\partial x'^{\mu}}{\partial x^{\nu}} v^{\nu}, \qquad \omega'_{\mu} = \frac{\partial x^{\nu}}{\partial x'^{\mu}} \omega_{\nu}, \qquad g'_{\mu\nu} = \frac{\partial x^{\alpha}}{\partial x'^{\mu}}\frac{\partial x^{\beta}}{\partial x'^{\nu}} g_{\alpha\beta}$$

- **Why two kinds of index:** components of a *vector* scale like the coordinate increments $dx^\mu$; components of a *gradient* $\partial_\mu f$ scale inversely. Contracting one of each gives a number everyone agrees on — e.g. $\omega_\mu v^\mu$, or the directional derivative $v^\mu \partial_\mu f$.

---

## Part 3 — Curves: the first curvature

Start with the simplest curved thing: a 1D curve in flat $\mathbb{R}^3$.

- Parametrize by **arclength** $s$ (so the tangent $\vec T = d\vec r/ds$ is a unit vector). Arclength is the natural parameter because it's coordinate-independent.
- **Curvature** $\kappa = \lVert d\vec T/ds \rVert$ — how fast the direction turns per unit length. Radius of the best-fitting circle is $1/\kappa$.
- **Torsion** $\tau$ — how fast the curve twists out of its instantaneous plane.

**Frenet–Serret frame** $(\vec T, \vec N, \vec B)$, an orthonormal triad carried along the curve:

$$\frac{d}{ds}\begin{pmatrix}\vec T\\ \vec N\\ \vec B\end{pmatrix} = \begin{pmatrix} 0 & \kappa & 0\\ -\kappa & 0 & \tau \\ 0 & -\tau & 0\end{pmatrix}\begin{pmatrix}\vec T\\ \vec N\\ \vec B\end{pmatrix}$$

For a general parametrization $\vec r(t)$:

$$\kappa = \frac{\lVert \vec r\,' \times \vec r\,''\rVert}{\lVert \vec r\,'\rVert^3}, \qquad \tau = \frac{(\vec r\,' \times \vec r\,'')\cdot \vec r\,'''}{\lVert \vec r\,' \times \vec r\,''\rVert^2}$$

- **Fundamental theorem of curves**: $\kappa(s)$ and $\tau(s)$ determine the curve uniquely up to rigid motion.
- **Crucially, this curvature is extrinsic.** An ant confined to the curve can't feel it: every curve is intrinsically a copy of the real line (arclength is all it has). Bend a wire — no intrinsic change. This is the contrast that makes Gauss's surface result surprising.

---

## Part 4 — Surfaces in $\mathbb{R}^3$: two fundamental forms

Parametrize a surface as $\vec r(u,v)$. Tangent vectors $\vec r_u, \vec r_v$ span the tangent plane; $\vec n = \dfrac{\vec r_u \times \vec r_v}{\lVert \vec r_u \times \vec r_v\rVert}$ is the unit normal.

### First fundamental form $\mathrm{I}$ — the intrinsic part (the metric)

$$\mathrm{I} = E\,du^2 + 2F\,du\,dv + G\,dv^2, \qquad E = \vec r_u\!\cdot\!\vec r_u,\; F = \vec r_u\!\cdot\!\vec r_v,\; G = \vec r_v\!\cdot\!\vec r_v$$

- This *is* $g_{ij}$ in the coordinates $(u,v)$. It gives lengths, angles and areas ($dA = \sqrt{EG - F^2}\, du\, dv$) **measurable by an ant living in the surface**.

### Second fundamental form $\mathrm{II}$ — the extrinsic part (how it bends in $\mathbb{R}^3$)

$$\mathrm{II} = L\,du^2 + 2M\,du\,dv + N\,dv^2, \qquad L = \vec r_{uu}\!\cdot\!\vec n,\; M = \vec r_{uv}\!\cdot\!\vec n,\; N = \vec r_{vv}\!\cdot\!\vec n$$

- $\mathrm{II}$ measures how fast the normal tips over — equivalently the **shape operator** $S = \mathrm{I}^{-1}\mathrm{II}$ (the derivative of the Gauss map $p \mapsto \vec n(p)$).
- **Principal curvatures** $\kappa_1, \kappa_2$ = eigenvalues of $S$ = max/min normal curvature; their eigendirections are orthogonal (principal directions).

### The two curvature scalars

$$K = \kappa_1 \kappa_2 = \frac{LN - M^2}{EG - F^2} \qquad \textbf{(Gaussian)}$$

$$H = \tfrac{1}{2}(\kappa_1 + \kappa_2) = \frac{EN - 2FM + GL}{2(EG - F^2)} \qquad \textbf{(mean)}$$

| Surface | $K$ | Character |
|---|---|---|
| Plane | $0$ | flat |
| Cylinder radius $a$ | $0$ | **flat!** ($\kappa_1 = 1/a,\ \kappa_2 = 0$) |
| Cone (off apex) | $0$ | flat |
| Sphere radius $a$ | $+1/a^2$ | positive, elliptic |
| Torus $(R,a)$ | $\dfrac{\cos v}{a(R + a\cos v)}$ | $>0$ outside, $0$ on top/bottom, $<0$ inside |
| Saddle $z = x^2 - y^2$ (at origin) | $<0$ | negative, hyperbolic |
| Catenoid / helicoid | $<0$ | minimal surfaces: $H = 0$ |
| Pseudosphere | $-1/a^2$ | constant negative |

- **Sign intuition:** $K>0$ ⇒ surface curves the *same* way in all directions (dome); $K<0$ ⇒ opposite ways (saddle); $K=0$ ⇒ flat in at least one direction (developable — can be unrolled onto a plane without stretching).
- $H$ is genuinely extrinsic (it flips sign if you flip $\vec n$; it's what soap films minimise). $K$ is not — see next.

---

## Part 5 — Theorema Egregium: curvature goes intrinsic

> **Gauss (1827):** $K$ depends only on $E, F, G$ and their derivatives — not on $L, M, N$.

$K$ is defined by an *extrinsic* formula involving the embedding, yet it is computable from the metric alone. Consequences:

- **Bending without stretching preserves $K$.** Any isometry (local, distance-preserving map) preserves $K$.
- **A cylinder is intrinsically flat.** Roll a sheet of paper: $\kappa_1$ changes from $0$ to $1/a$, but $K = \kappa_1\kappa_2 = 0$ throughout. An ant on a cylinder measures Euclidean geometry.
- **No perfect flat map of the Earth exists.** Sphere has $K = 1/a^2 \neq 0$, paper has $K = 0$; no isometry can exist. Every map projection must distort area, angle, or distance (Mercator keeps angles, ruins area).
- **The pizza trick.** Fold a slice: you force $\kappa_1 \neq 0$ transversally, so to keep $K = 0$ the longitudinal $\kappa_2$ must stay $0$ — the tip can't flop.
- **Curvature is detectable from inside.** This is the conceptual seed of GR: we can't step outside spacetime, but we can still measure its curvature.

Explicitly, in orthogonal coordinates ($F = 0$):

$$K = -\frac{1}{2\sqrt{EG}}\left[\frac{\partial}{\partial u}\!\left(\frac{G_u}{\sqrt{EG}}\right) + \frac{\partial}{\partial v}\!\left(\frac{E_v}{\sqrt{EG}}\right)\right]$$

The general-metric version is the Brioschi formula; the cleanest route is via the Riemann tensor (Part 10): **in 2D, $R = 2K$**.

---

## Part 6 — Intrinsic 2D geometry: how an ant measures $K$

Once curvature is intrinsic, ask: what experiments does a flatlander run?

### Geodesics — "straight lines"

- Two equivalent definitions: **locally shortest** path, or **straightest** path (tangent parallel-transported along itself, zero geodesic curvature $\kappa_g$).
- Sphere → great circles. Plane → straight lines. Poincaré half-plane → vertical lines and semicircles meeting the boundary at right angles.
- Geodesics are locally, not globally, shortest: the long way round a great circle is still a geodesic.

### Four intrinsic curvature detectors

1. **Angle excess of a geodesic triangle:**

$$\alpha + \beta + \gamma - \pi = \int_\triangle K \, dA$$

   On a unit sphere, a triangle with three right angles has excess $\pi/2$ = its area. On a saddle, angles sum to *less* than $\pi$.

2. **Circumference/area deficit (Bertrand–Puiseux):** draw a circle of geodesic radius $r$,

$$C(r) = 2\pi r\left(1 - \frac{K r^2}{6} + O(r^4)\right), \qquad A(r) = \pi r^2\left(1 - \frac{K r^2}{12} + O(r^4)\right)$$

   $K>0$ ⇒ less circumference than Euclid predicts; $K<0$ ⇒ more. Equivalently $\;K = \lim_{r\to0}\frac{3}{\pi r^3}\left(2\pi r - C(r)\right)$.

3. **Parallel transport around a loop (holonomy):** carry a vector around a closed loop keeping it "as parallel as possible". It comes back **rotated** by

$$\Delta\theta = \oint \text{(connection)} = \iint_{\text{enclosed}} K\, dA$$

   On a sphere, walking a loop at colatitude $\theta_0$ rotates the vector by $2\pi(1 - \cos\theta_0)$ = the enclosed solid angle. **This is curvature's most direct definition and it generalises verbatim to the Riemann tensor.**

4. **Geodesic deviation:** two initially parallel geodesics converge if $K>0$ (great circles meet at the poles), diverge exponentially if $K<0$. This becomes *tidal gravity* in Part 11.

### Gauss–Bonnet — local curvature knows global topology

$$\int_M K\, dA + \oint_{\partial M} \kappa_g\, ds = 2\pi \chi(M)$$

- $\chi$ = Euler characteristic ($2$ for a sphere, $0$ for a torus, $2-2g$ for genus $g$). No matter how you deform a sphere, total curvature is exactly $4\pi$. A torus must have as much negative as positive curvature — exactly what the $\cos v$ in the table says.

### The three constant-curvature model geometries

| | $K$ | Parallels through a point off a line | Triangle angles | Model metric |
|---|---|---|---|---|
| Spherical | $+1$ | none | $>\pi$ | $d\theta^2 + \sin^2\!\theta\, d\phi^2$ |
| Euclidean | $0$ | exactly one | $=\pi$ | $dx^2 + dy^2$ |
| Hyperbolic | $-1$ | infinitely many | $<\pi$ | $\dfrac{dx^2+dy^2}{y^2}$ (half-plane) |

- These are the resolution of the 2000-year parallel-postulate problem: it's *independent*, and denying it gives consistent geometries (Bolyai, Lobachevsky, ~1830).

---

## Part 7 — Manifolds and tensors: geometry without an ambient space

To do this in $n$ dimensions with no $\mathbb{R}^3$ to sit inside, you need the manifold language.

### Manifold

- A **manifold** $M$ of dimension $n$: a space that looks locally like $\mathbb{R}^n$. Formally, an atlas of **charts** $\varphi_\alpha: U_\alpha \to \mathbb{R}^n$ with smooth transition maps $\varphi_\beta \circ \varphi_\alpha^{-1}$ on overlaps.
- Charts = coordinate systems. There is generally **no single global chart** (a sphere needs at least two — hence the coordinate singularity at $\theta = 0$ in spherical coordinates, and the one at $r = 2M$ in Schwarzschild).
- Nothing so far involves distance. A bare manifold has topology and calculus but no geometry.

### Tangent space

- $T_pM$: the vector space of directions at $p$, dimension $n$. Two equivalent definitions:
  - velocities $\dot\gamma(0)$ of curves through $p$;
  - **derivations**: linear maps on functions obeying Leibniz. This is the ambient-free definition, and it's why the natural basis is $\{\partial_\mu\} \equiv \{\partial/\partial x^\mu\}$ — a vector *is* a directional derivative: $v = v^\mu \partial_\mu$.
- **Key point:** tangent spaces at different points are *different vector spaces*. You cannot subtract a vector at $p$ from a vector at $q$ — hence you cannot naively differentiate a vector field. Fixing that requires a connection (Part 9), and the failure of that fix to be path-independent *is* curvature (Part 10).

### Cotangent space and tensors

- $T_p^*M$: linear maps $T_pM \to \mathbb{R}$ (**covectors**, **1-forms**), basis $\{dx^\mu\}$ dual to $\{\partial_\mu\}$: $\;dx^\mu(\partial_\nu) = \delta^\mu_\nu$. The gradient $df = \partial_\mu f\, dx^\mu$ is the canonical 1-form.
- A **$(k,l)$ tensor** is a multilinear map eating $k$ covectors and $l$ vectors:

$$T = T^{\mu_1\dots\mu_k}{}_{\nu_1\dots\nu_l}\; \partial_{\mu_1}\!\otimes\cdots\otimes dx^{\nu_l}$$

  with components transforming by one Jacobian factor per index. **Tensor equations are coordinate-independent statements** — if a tensor vanishes in one frame it vanishes in all. That's the whole reason physics is written in tensors.
- $\Gamma^\lambda_{\mu\nu}$ is famously **not** a tensor (its transformation law has an inhomogeneous second-derivative term) — which is exactly why gravity can be transformed away locally but curvature cannot.

---

## Part 8 — The Riemannian metric

- A **metric** $g$ is a smooth, symmetric, non-degenerate $(0,2)$ tensor field: an inner product on each $T_pM$.

$$ds^2 = g_{\mu\nu}\,dx^\mu dx^\nu, \qquad g(u,v) = g_{\mu\nu}u^\mu v^\nu$$

- **Riemannian**: $g$ positive-definite (signature $+++\dots$) — ordinary lengths and angles.
- **Pseudo-/Lorentzian**: signature $(-,+,+,+)$ — one "time" direction. Vectors come in three flavours (timelike/null/spacelike); this is spacetime.
- The metric supplies:
  - **lengths & angles**: $\lVert v\rVert^2 = g_{\mu\nu}v^\mu v^\nu$;
  - **index gymnastics**: $v_\mu = g_{\mu\nu}v^\nu$, $\;v^\mu = g^{\mu\nu}v_\nu$ with $g^{\mu\nu}$ the matrix inverse. Vectors and covectors become interchangeable *only* once you have a metric;
  - **volume**: $dV = \sqrt{\lvert \det g\rvert}\; d^nx$;
  - **the connection** (Part 9) and hence everything else.
- $g_{\mu\nu}$ has $n(n+1)/2$ independent components (10 in 4D). Of those, $n$ can be killed by coordinate choice, leaving $n(n-1)/2$ physical functions — in 4D, 6, of which 2 are the propagating graviton polarizations.

---

## Part 9 — Connection: differentiating in a curved space

### Why $\partial_\mu v^\nu$ fails

- In polar coordinates the basis vectors $\partial_r, \partial_\theta$ themselves rotate and rescale from point to point. So $\partial_\mu v^\nu$ mixes "the vector changed" with "the basis changed" — and it is **not a tensor**.
- Fix: add a correction encoding how the basis turns.

$$\nabla_\mu v^\nu = \partial_\mu v^\nu + \Gamma^{\nu}_{\mu\lambda} v^{\lambda}, \qquad \nabla_\mu \omega_\nu = \partial_\mu \omega_\nu - \Gamma^{\lambda}_{\mu\nu}\omega_\lambda$$

(one $+\Gamma$ per upper index, one $-\Gamma$ per lower index; $\nabla_\mu f = \partial_\mu f$ on scalars).

### The Levi-Civita connection

Impose two natural conditions:

1. **Metric compatibility** $\nabla_\lambda g_{\mu\nu} = 0$ — transport preserves lengths and angles;
2. **Torsion-free** $\Gamma^\lambda_{\mu\nu} = \Gamma^\lambda_{\nu\mu}$ — infinitesimal parallelograms close.

These determine $\Gamma$ **uniquely** (the fundamental theorem of Riemannian geometry):

$$\boxed{\;\Gamma^{\lambda}_{\mu\nu} = \tfrac{1}{2}\, g^{\lambda\sigma}\left(\partial_\mu g_{\sigma\nu} + \partial_\nu g_{\sigma\mu} - \partial_\sigma g_{\mu\nu}\right)\;}$$

- Interpretation: $\Gamma$ is the "fictitious force" / gravitational-field term. It's built from **first** derivatives of $g$.
- Because it's not a tensor, at any point $p$ you can choose **Riemann normal coordinates** with $g_{\mu\nu}(p) = \eta_{\mu\nu}$ and $\Gamma^\lambda_{\mu\nu}(p) = 0$. Locally, space looks flat and gravity vanishes — the mathematical form of the equivalence principle. What you *cannot* remove is the second derivative:

$$g_{\mu\nu}(x) = \eta_{\mu\nu} - \tfrac{1}{3}R_{\mu\alpha\nu\beta}\,x^\alpha x^\beta + O(x^3)$$

  **Curvature is the leading obstruction to being flat.**

### Parallel transport

- A vector is parallel-transported along a curve $x^\mu(\lambda)$ if its covariant derivative along the curve vanishes:

$$\frac{DV^\mu}{d\lambda} \equiv \frac{dV^\mu}{d\lambda} + \Gamma^{\mu}_{\nu\rho}\frac{dx^\nu}{d\lambda}V^\rho = 0$$

- This is a linear ODE — transport is a linear map between tangent spaces. It's **path-dependent** in curved space; the mismatch around a closed loop (**holonomy**) is curvature.

### Geodesics

- "Transport your own tangent vector along yourself":

$$\boxed{\;\frac{d^2x^\mu}{d\lambda^2} + \Gamma^{\mu}_{\nu\rho}\frac{dx^\nu}{d\lambda}\frac{dx^\rho}{d\lambda} = 0\;}$$

- Equivalently the Euler–Lagrange equations of $S = \int g_{\mu\nu}\dot x^\mu \dot x^\nu\, d\lambda$ — usually the fastest way to get $\Gamma$ by hand.
- $\lambda$ must be an **affine parameter** (proper time/arclength up to $a\lambda + b$); otherwise a $\propto \dot x^\mu$ term appears on the right.
- **Killing vectors** ⇒ conserved quantities: if $\partial_\sigma g_{\mu\nu}=0$ for some coordinate $x^\sigma$, then $p_\sigma = g_{\sigma\nu}\dot x^\nu$ is conserved along geodesics (energy from $t$-independence, angular momentum from $\phi$-independence). This is how you actually solve orbits.

---

## Part 10 — The Riemann curvature tensor

### Definition: curvature = failure of derivatives to commute

$$[\nabla_\mu, \nabla_\nu]\,V^\rho = R^{\rho}{}_{\sigma\mu\nu}\, V^{\sigma}$$

Equivalently, transport $V$ around an infinitesimal parallelogram spanned by $\delta a^\mu, \delta b^\nu$ and it comes back changed by $\delta V^\rho = -R^\rho{}_{\sigma\mu\nu}V^\sigma \delta a^\mu \delta b^\nu$. Same holonomy idea as the sphere in Part 6, now to first order in area.

$$R^{\rho}{}_{\sigma\mu\nu} = \partial_\mu \Gamma^{\rho}_{\nu\sigma} - \partial_\nu \Gamma^{\rho}_{\mu\sigma} + \Gamma^{\rho}_{\mu\lambda}\Gamma^{\lambda}_{\nu\sigma} - \Gamma^{\rho}_{\nu\lambda}\Gamma^{\lambda}_{\mu\sigma}$$

- Schematically $R \sim \partial\Gamma + \Gamma\Gamma \sim \partial^2 g + (\partial g)^2$: **second derivatives of the metric**, which is why it can't be transformed away.
- $R^\rho{}_{\sigma\mu\nu} = 0$ everywhere $\iff$ the space is flat (locally isometric to $\mathbb{R}^n$/Minkowski) $\iff$ parallel transport is path-independent.

### Symmetries and component count

With all indices lowered, $R_{\rho\sigma\mu\nu}$ satisfies

| Symmetry | Statement |
|---|---|
| Antisymmetry (last pair) | $R_{\rho\sigma\mu\nu} = -R_{\rho\sigma\nu\mu}$ |
| Antisymmetry (first pair) | $R_{\rho\sigma\mu\nu} = -R_{\sigma\rho\mu\nu}$ |
| Pair exchange | $R_{\rho\sigma\mu\nu} = R_{\mu\nu\rho\sigma}$ |
| First Bianchi | $R_{\rho[\sigma\mu\nu]} = 0$ |
| Second (differential) Bianchi | $\nabla_{[\lambda}R_{\rho\sigma]\mu\nu} = 0$ |

Independent components: $\;\dfrac{n^2(n^2-1)}{12}$.

| $n$ | components | note |
|---|---|---|
| 1 | 0 | curves are never intrinsically curved |
| 2 | 1 | the single number is $K$: $\;R_{1212} = K(g_{11}g_{22}-g_{12}^2)$, $\;R = 2K$ |
| 3 | 6 | = components of Ricci ⇒ Ricci determines Riemann; **no vacuum gravity in 3D** |
| 4 | 20 | 10 Ricci + 10 Weyl ⇒ vacuum curvature exists: black holes, gravitational waves |

### The contractions

$$R_{\mu\nu} = R^{\lambda}{}_{\mu\lambda\nu} \quad \textbf{(Ricci)}, \qquad R = g^{\mu\nu}R_{\mu\nu} \quad \textbf{(scalar)}$$

$$C_{\rho\sigma\mu\nu} = R_{\rho\sigma\mu\nu} - \big(g_{\rho[\mu}R_{\nu]\sigma} - g_{\sigma[\mu}R_{\nu]\rho}\big) + \tfrac{1}{3}R\,g_{\rho[\mu}g_{\nu]\sigma} \quad \textbf{(Weyl, in 4D)}$$

Geometric meanings — worth memorising:

| Object | Measures | Physics |
|---|---|---|
| **Sectional curvature** $\;\dfrac{R_{\mu\nu\rho\sigma}u^\mu v^\nu u^\rho v^\sigma}{\lVert u\rVert^2\lVert v\rVert^2 - (u\cdot v)^2}$ | Gauss curvature of the 2D surface swept by geodesics in the $u\wedge v$ plane | tidal effect in one plane |
| **Riemann** $R^\rho{}_{\sigma\mu\nu}$ | full tidal information; holonomy per unit area | complete gravitational field |
| **Ricci** $R_{\mu\nu}$ | **volume** change of a small ball of freely-falling dust | sourced directly by local matter |
| **Scalar** $R$ | single-number volume deficit: $\dfrac{V_g(\epsilon)}{V_{\text{flat}}(\epsilon)} = 1 - \dfrac{R\,\epsilon^2}{6(n+2)} + \dots$ | curvature "trace" |
| **Weyl** $C_{\rho\sigma\mu\nu}$ | **shape** distortion at fixed volume (trace-free part) | tidal stretching, gravitational radiation, the free field |

- **Vacuum ($R_{\mu\nu}=0$) $\neq$ flat.** Outside a star, Ricci vanishes but Weyl doesn't: a falling ball of dust keeps its volume while being stretched into a cigar. That's Schwarzschild.

---

## Part 11 — What curvature *does*: geodesic deviation

Take a family of neighbouring geodesics with 4-velocity $u^\mu$ and separation vector $S^\mu$. Then

$$\boxed{\;\frac{D^2 S^{\mu}}{d\tau^2} = -R^{\mu}{}_{\nu\rho\sigma}\, u^{\nu} S^{\rho} u^{\sigma}\;}$$

- **This is the physical definition of curvature.** A single freely-falling observer feels nothing (equivalence principle); *two* of them measure relative acceleration, and that relative acceleration is Riemann.
- Newtonian counterpart: $\;\ddot\xi^i = -\partial_i\partial_j\Phi\;\xi^j$. Matching gives $\;R^{i}{}_{0j0} \approx \dfrac{1}{c^2}\partial_i\partial_j\Phi$: **the Riemann tensor is the relativistic tidal tensor.**
- For a point mass, the tidal tensor in the free-fall frame is $\;\mathrm{diag}\!\left(-\tfrac{2GM}{r^3}, \tfrac{GM}{r^3}, \tfrac{GM}{r^3}\right)$ — stretch along the radius, squeeze transversally: **spaghettification**, and the twice-a-day ocean tides.
- Its trace is $\nabla^2 \Phi = 4\pi G\rho$: **zero in vacuum**. Volume preserved, shape distorted — the Ricci/Weyl split made physical. (In matter the trace is positive and the ball's volume genuinely shrinks — Raychaudhuri's focusing theorem.)
- A passing gravitational wave is the same equation with an oscillating Weyl tensor: a ring of free particles alternately squashed along $+$ and $\times$ polarizations. LIGO measures $\delta L/L \sim 10^{-21}$.

---

## Part 12 — Lorentzian geometry: spacetime

Everything above is unchanged except the metric signature.

$$ds^2 = -c^2dt^2 + dx^2 + dy^2 + dz^2 = \eta_{\mu\nu}dx^\mu dx^\nu, \qquad \eta = \mathrm{diag}(-1,1,1,1)$$

- The minus sign is the whole of special relativity. The invariant interval, not distance or duration separately, is what all observers agree on; the symmetry group is the Lorentz group instead of rotations.
- **Causal classification** of a displacement:

| $ds^2$ | Name | Meaning |
|---|---|---|
| $<0$ | timelike | massive particle can travel it; $d\tau^2 = -ds^2/c^2$ is elapsed proper time |
| $=0$ | null / lightlike | light rays; the **light cone** |
| $>0$ | spacelike | causally disconnected; "simultaneous" for some observer |

- **Worldlines** are curves in spacetime; a particle's **4-velocity** $u^\mu = dx^\mu/d\tau$ always satisfies $g_{\mu\nu}u^\mu u^\nu = -c^2$.
- **Proper time is the arclength of a timelike worldline** — and free fall *maximises* it (the timelike geodesic is the longest path, opposite to Riemannian intuition, thanks to the signature). The twin paradox is just the reverse triangle inequality.
- Light cones are the local causal structure. In curved spacetime they tip over from point to point — that tipping is the picture of a black-hole horizon (the cone tips until all futures point inward).

---

## Part 13 — Equivalence principle: gravity *is* curvature

- **Weak EP**: inertial mass = gravitational mass; all bodies fall identically (Eötvös/MICROSCOPE tests to $\sim10^{-15}$). A universal force is not a force — it's a property of the arena.
- **Einstein EP**: in a small enough freely-falling lab, *all* of physics is that of special relativity. Compare with Part 9: at any point you can choose normal coordinates where $g = \eta$ and $\Gamma = 0$. **The equivalence principle is the statement that spacetime is a Lorentzian manifold.**
- "Small enough" matters: over a bigger lab, tidal effects (Riemann) show up. Curvature is precisely what you cannot transform away.
- So: **a freely-falling particle moves on a timelike geodesic; light moves on a null geodesic.** Gravity is not a force in the equation of motion — it's in $\Gamma$, i.e. in the geometry.

### Recovering Newton

Weak field $g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu}$, $\lVert h \rVert \ll 1$, slow motion, static:

- $\;g_{00} = -\left(1 + \dfrac{2\Phi}{c^2}\right)$, $\;\Gamma^i_{00} \approx \dfrac{1}{c^2}\partial_i \Phi$
- Geodesic equation → $\;\ddot x^i = -\partial_i\Phi$. Newton's law recovered; **gravitational potential = time dilation.**
- **Gravitational time dilation / redshift**: $\;\dfrac{d\tau}{dt} = \sqrt{1 - \dfrac{2GM}{rc^2}} \approx 1 + \dfrac{\Phi}{c^2}$. Clocks tick slower deeper in the well (Pound–Rebka 1959; GPS needs $\sim38\ \mu\text{s/day}$ of correction, or navigation drifts ~10 km/day).
- **Newtonian gravity is almost entirely the curvature of *time*.** For everyday speeds $v \ll c$, the $g_{00}$ term dominates; space curvature contributes the "extra" light bending and Mercury's precession.

---

## Part 14 — The Einstein field equations

Now build the field equation. Requirements: second-order in $g$, tensorial, reduces to $\nabla^2\Phi = 4\pi G\rho$, and consistent with $\nabla^\mu T_{\mu\nu}=0$.

- Newtonian $\nabla^2 \Phi = 4\pi G \rho$ is the trace of the tidal tensor equation ⇒ the relativistic version constrains the trace of Riemann, i.e. **Ricci**.
- Naive guess $R_{\mu\nu} = \kappa T_{\mu\nu}$ fails: $\nabla^\mu R_{\mu\nu} = \tfrac12 \nabla_\nu R \neq 0$, so it would force $\nabla^\mu T_{\mu\nu} \neq 0$ (Einstein's own 1913–15 wrong turn).
- The contracted **Bianchi identity** $\nabla^\mu\!\left(R_{\mu\nu} - \tfrac12 R g_{\mu\nu}\right) = 0$ singles out the Einstein tensor:

$$\boxed{\;G_{\mu\nu} + \Lambda g_{\mu\nu} = \frac{8\pi G}{c^4} T_{\mu\nu}, \qquad G_{\mu\nu} \equiv R_{\mu\nu} - \tfrac12 R\, g_{\mu\nu}\;}$$

- **Lovelock's theorem** makes this a near-uniqueness result: in 4D, $G_{\mu\nu} + \Lambda g_{\mu\nu}$ is the *only* divergence-free symmetric 2-tensor built from $g$ and up to its second derivatives.
- Variational route: extremise the **Einstein–Hilbert action**

$$S = \frac{c^4}{16\pi G}\int (R - 2\Lambda)\sqrt{-g}\; d^4x + S_{\text{matter}}$$

  with $T_{\mu\nu} = -\dfrac{2}{\sqrt{-g}}\dfrac{\delta S_\text{matter}}{\delta g^{\mu\nu}}$. Simplest possible scalar built from curvature ⇒ GR.
- **Structure:** 10 coupled nonlinear PDEs; 4 are constraints, the Bianchi identities remove 4 more, leaving 2 dynamical d.o.f. (the two graviton polarizations). Nonlinear because gravitational energy itself gravitates.
- Full details, solutions and cosmology: see [`topics/einstein-field-equations`](../einstein-field-equations/README.md).

---

## Part 15 — Worked spacetimes

### Schwarzschild (spherical vacuum)

$$ds^2 = -\left(1 - \frac{2GM}{rc^2}\right)c^2 dt^2 + \left(1 - \frac{2GM}{rc^2}\right)^{-1} dr^2 + r^2 d\Omega^2$$

| Feature | Value ($G=c=1$) | Note |
|---|---|---|
| Event horizon | $r_s = 2M$ | coordinate singularity only — $K$ finite; cross it in Eddington–Finkelstein coordinates |
| True singularity | $r = 0$ | Kretschmann $R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma} = 48M^2/r^6 \to \infty$ |
| Photon sphere | $r = 3M$ | unstable light orbits → black hole "shadow" (EHT, 2019) |
| ISCO | $r = 6M$ | inner edge of accretion disks |
| Light deflection | $\Delta\phi = \dfrac{4GM}{c^2 b}$ | $1.75''$ at the solar limb — twice the naive Newtonian value (Eddington 1919) |
| Perihelion precession | $\dfrac{6\pi GM}{c^2 a(1-e^2)}$ per orbit | Mercury: $43''$/century |
| Shapiro delay | $\sim \dfrac{4GM}{c^3}\ln(\dots)$ | radar echo delay, Cassini test to $10^{-5}$ |

- **Birkhoff's theorem**: any spherically symmetric vacuum solution is Schwarzschild. So a pulsating spherical star radiates no gravitational waves, and the exterior field is static.

### FLRW (homogeneous, isotropic — the universe)

$$ds^2 = -c^2dt^2 + a(t)^2\left[\frac{dr^2}{1-kr^2} + r^2 d\Omega^2\right]$$

Plugging into the EFE collapses 10 equations to the two Friedmann equations for $a(t)$ ⇒ expansion, Big Bang, dark energy. Note the two independent curvatures: **spatial** $k$ (flat, per Planck data) and **spacetime** curvature (nonzero — the universe expands).

### Linearized gravity (waves)

$g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu}$ with $\lVert h\rVert \ll 1$; in transverse-traceless gauge the EFE become a wave equation $\;\Box \bar h_{\mu\nu} = -\tfrac{16\pi G}{c^4}T_{\mu\nu}$. Ripples of Weyl curvature at speed $c$, quadrupolar (no monopole/dipole radiation) — GW150914.

---

## Part 16 — Mental models and common traps

| Trap | Reality |
|---|---|
| "Curvy coordinates = curved space" | Polar coordinates on a plane have $\Gamma \neq 0$ but $R = 0$. Only $R^\rho{}_{\sigma\mu\nu}$ decides. |
| "Curved space needs a higher dimension to bend into" | Intrinsic curvature is defined without embedding (Theorema Egregium). Spacetime bends into nothing. |
| "$R_{\mu\nu}=0$ means flat/empty means nothing happens" | Ricci-flat $\neq$ Riemann-flat. Black holes and gravitational waves live in the Weyl part. |
| "The rubber-sheet picture explains gravity" | It smuggles in downward gravity and shows only space curvature. Newtonian gravity is mostly curvature of **time**. |
| "Nothing special happens at the horizon" (as a paradox) | Correct, and it's a coordinate artifact — an infalling observer notices nothing locally; tidal forces $\sim M/r^3$ are *small* at a supermassive horizon. |
| "Space is expanding, so galaxies fly through space" | The metric's $a(t)$ grows; recession is geometry, not motion through space. |
| "Gravity is a force" | It's the $\Gamma$ term — removable at a point. The irremovable part (Riemann) is tides. |
| "$\Gamma$ is a tensor because it has indices" | It isn't; that's exactly why local flatness is possible. |

### Notation cheat sheet

| Symbol | Reads as |
|---|---|
| $g_{\mu\nu}$, $g^{\mu\nu}$ | metric, inverse metric |
| $\partial_\mu$ | $\partial/\partial x^\mu$; also the coordinate basis vector |
| $\Gamma^\lambda_{\mu\nu}$ | Christoffel symbols (connection), $\sim\partial g$ |
| $\nabla_\mu$ | covariant derivative |
| $D/d\lambda$ | covariant derivative along a curve |
| $R^\rho{}_{\sigma\mu\nu}$, $R_{\mu\nu}$, $R$ | Riemann, Ricci, Ricci scalar, $\sim \partial^2 g$ |
| $C_{\rho\sigma\mu\nu}$ | Weyl (trace-free Riemann) |
| $G_{\mu\nu}$ | Einstein tensor |
| $T_{\mu\nu}$ | stress–energy |
| $u^\mu$, $\tau$ | 4-velocity, proper time |
| $[\,\cdot\,]$, $(\,\cdot\,)$ around indices | antisymmetrization, symmetrization |

### The whole subject in one chain

$$g_{\mu\nu} \;\xrightarrow{\ \partial\ }\; \Gamma^{\lambda}_{\mu\nu} \;\xrightarrow{\ \partial\ }\; R^{\rho}{}_{\sigma\mu\nu} \;\xrightarrow{\ \text{contract}\ }\; R_{\mu\nu},\, R \;\longrightarrow\; G_{\mu\nu} \;=\; \frac{8\pi G}{c^4}T_{\mu\nu}$$

with $\Gamma$ giving **motion** (geodesics), Riemann giving **tides**, and the last equality giving the **field equation**.

---

## Code samples

```bash
uv run topics/geometry-to-spacetime-curvature/curves_and_surfaces.py   # Parts 3–5
uv run topics/geometry-to-spacetime-curvature/parallel_transport.py    # Parts 6, 9–10
uv run topics/geometry-to-spacetime-curvature/geodesics.py             # Parts 9, 13, 15
uv run topics/geometry-to-spacetime-curvature/tidal_curvature.py       # Parts 10–11
```

| Script | Shows |
|---|---|
| `curves_and_surfaces.py` | Frenet $\kappa,\tau$ for circle/helix; first and second fundamental forms; $K$ and $H$ for plane/cylinder/sphere/torus/catenoid; **verifies Theorema Egregium** on each by recomputing $K$ from the first fundamental form alone via $R = 2K$ |
| `parallel_transport.py` | RK4 on the transport ODE. Flat plane in polar coordinates: $\Gamma \neq 0$ but holonomy $=0$. Sphere: rotation around a latitude circle $= 2\pi(1-\cos\theta_0)$. Sphere and Poincaré half-plane: holonomy $= \iint K\, dA$ exactly, both signs of $K$. Transport preserves $\lVert V \rVert$ to $10^{-15}$ ($\nabla g = 0$) |
| `geodesics.py` | Generic geodesic integrator from any symbolic metric. Straight lines in polar coordinates; great circles on the sphere; a Schwarzschild orbit conserving $E$, $L$ and $u\cdot u=-1$ to $10^{-13}$, with its precession cross-checked against an exact quadrature; light deflection $4GM/c^2b = 1.751''$; Mercury's $42.97''$/century |
| `tidal_curvature.py` | Riemann → tidal tensor in Schwarzschild: $\mathrm{diag}(-2M/r^3, M/r^3, M/r^3)$, identical for static and infalling observers. Traceless in vacuum → a dust ball's volume is stationary while its shape distorts; inside matter the volume shrinks at $-2\pi\rho$. Kretschmann $=48M^2/r^6$. Real tidal numbers: ocean tides, spaghettification at stellar vs supermassive horizons |

## Further reading

- **Gentle**: Schutz, *A First Course in General Relativity* • Hartle, *Gravity* (physics-first) • Needham, *Visual Differential Geometry and Forms*
- **Standard**: Carroll, *Spacetime and Geometry* (best single text for this arc; his lecture notes are free) • Wald, *General Relativity* (rigorous)
- **Encyclopaedic**: Misner–Thorne–Wheeler, *Gravitation*
- **Pure geometry**: do Carmo, *Differential Geometry of Curves and Surfaces* → *Riemannian Geometry* • Lee, *Introduction to Smooth Manifolds*
