# Einstein field equations

The 10 coupled PDEs at the heart of **general relativity**. They say: **matter/energy tells spacetime how to curve; curvature tells matter how to move.** Gravity isn't a force — it's the geometry of curved 4D spacetime.

$$G_{\mu\nu} + \Lambda\, g_{\mu\nu} = \frac{8\pi G}{c^4}\, T_{\mu\nu}$$

- **Left = geometry** (how spacetime is curved), **right = matter** (what's in it). One equation, but the indices $\mu,\nu$ run 0–3, so it's a $4\times4$ symmetric system → **10 independent equations**.
- Wheeler's summary: *"Spacetime tells matter how to move; matter tells spacetime how to curve."*

## What each symbol is

| Symbol | Name | What it is |
|---|---|---|
| $g_{\mu\nu}$ | **Metric tensor** | The unknown you solve for. Encodes distances/time in curved spacetime. The "shape". |
| $G_{\mu\nu}$ | **Einstein tensor** | $R_{\mu\nu} - \tfrac{1}{2} R\, g_{\mu\nu}$. A specific combination of curvature. Built so $\nabla^\mu G_{\mu\nu} = 0$ automatically. |
| $R_{\mu\nu}$ | **Ricci tensor** | Curvature contracted once — how volumes shrink/grow along geodesics. |
| $R$ | **Ricci scalar** | $g^{\mu\nu} R_{\mu\nu}$. Curvature boiled down to one number per point. |
| $T_{\mu\nu}$ | **Stress–energy tensor** | The source: energy density, momentum, pressure, stress. |
| $\Lambda$ | **Cosmological constant** | Energy of empty space → **dark energy**. Drives accelerating expansion. |
| $\dfrac{8\pi G}{c^4}$ | coupling | Tiny ($\sim 2\times10^{-43}$). Gravity is weak — you need a planet's worth of $T$ to bend spacetime noticeably. |

## The chain of definitions (metric → curvature)

Everything is built from the metric by differentiation. This is the whole computation:

**1. Metric** (given/assumed): $g_{\mu\nu}$

**2. Christoffel symbols** (the connection):

$$\Gamma^{\lambda}_{\mu\nu} = \tfrac{1}{2}\, g^{\lambda\sigma}\left(\partial_\mu g_{\sigma\nu} + \partial_\nu g_{\sigma\mu} - \partial_\sigma g_{\mu\nu}\right)$$

**3. Riemann tensor** ($\partial\Gamma + \Gamma\Gamma$):

$$R^{\rho}{}_{\sigma\mu\nu} = \partial_\mu \Gamma^{\rho}_{\nu\sigma} - \partial_\nu \Gamma^{\rho}_{\mu\sigma} + \Gamma^{\rho}_{\mu\lambda}\Gamma^{\lambda}_{\nu\sigma} - \Gamma^{\rho}_{\nu\lambda}\Gamma^{\lambda}_{\mu\sigma}$$

**4. Ricci tensor** (contract $\rho = \mu$): $\;R_{\sigma\nu} = R^{\mu}{}_{\sigma\mu\nu}$

**5. Ricci scalar** (contract with $g$): $\;R = g^{\sigma\nu} R_{\sigma\nu}$

**6. Einstein tensor** (assemble): $\;G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2} R\, g_{\mu\nu}$

- **$\Gamma$ (Christoffel symbols)** aren't a tensor — they're the "connection" that says how to differentiate vectors in curved space (they encode the fictitious/gravitational acceleration). They vanish in freely-falling (locally inertial) frames.
- **Riemann $R^{\rho}{}_{\sigma\mu\nu}$** is the *full* curvature (20 independent components in 4D). It measures how a vector rotates when parallel-transported around a loop → **tidal forces**. Flat spacetime $\iff$ Riemann $= 0$ everywhere.
- The EFE only constrains the **Ricci** part (10 components). The rest of Riemann — the **Weyl tensor** — is curvature that propagates in vacuum: **tidal distortion and gravitational waves**.

## Why *this* combination — conservation forces it

- Conservation of energy–momentum is $\nabla^\mu T_{\mu\nu} = 0$ (local, automatic). So the geometry side must also be divergence-free.
- $\nabla^\mu R_{\mu\nu} \neq 0$ in general — you **can't** just set $R_{\mu\nu} \propto T_{\mu\nu}$ (that was Einstein's earlier wrong guess).
- The **contracted Bianchi identity** gives $\nabla^\mu \left(R_{\mu\nu} - \tfrac{1}{2} R\, g_{\mu\nu}\right) = 0$ identically. So $G_{\mu\nu}$ is the *unique* (up to $\Lambda$) curvature tensor that's automatically conserved. That's why it's the left-hand side.

## Key solutions (you almost never solve the general case)

The EFE are nonlinear and brutal. Progress comes from imposing **symmetry**:

| Solution | Assumption | Describes |
|---|---|---|
| **Minkowski** | $T = 0$, flat | Special relativity / empty spacetime |
| **Schwarzschild** (1916) | vacuum, spherical, static | Non-spinning black hole / outside any star. Event horizon at $r_s = 2GM/c^2$ |
| **Kerr** (1963) | vacuum, spherical, **spinning** | Rotating black hole (all real ones). Has ergosphere |
| **Reissner–Nordström / Kerr–Newman** | + electric charge | Charged (spinning) black holes |
| **FLRW** | homogeneous, isotropic | The whole universe → **Friedmann equations**, Big Bang, expansion |
| **Gravitational waves** | linearized $g = \eta + h$, $\lVert h\rVert \ll 1$ | Ripples in spacetime (LIGO, 2015) |

- **Vacuum EFE** ($T_{\mu\nu} = 0$, $\Lambda = 0$): contracting shows it reduces to $R_{\mu\nu} = 0$ (**Ricci-flat**) — *not* flat spacetime! Schwarzschild is Ricci-flat but Riemann $\neq 0$ (that's the tidal gravity outside a black hole).
- **Schwarzschild → Newton**: in the weak-field, slow-motion limit the $00$ equation becomes $\nabla^2 \Phi = 4\pi G \rho$ (Poisson) — you recover Newtonian gravity, with $g_{00} \approx -\left(1 + 2\Phi/c^2\right)$.

## FLRW / cosmology (where the EFE meets data)

Plug a homogeneous isotropic metric into the EFE and the 10 equations collapse to **two** (the Friedmann equations) for the scale factor $a(t)$:

$$\left(\frac{\dot a}{a}\right)^2 = \frac{8\pi G}{3}\,\rho \;-\; \frac{k c^2}{a^2} \;+\; \frac{\Lambda c^2}{3} \qquad\text{(expansion rate, } H = \dot a / a\text{)}$$

$$\frac{\ddot a}{a} = -\frac{4\pi G}{3}\left(\rho + \frac{3p}{c^2}\right) \;+\; \frac{\Lambda c^2}{3} \qquad\text{(acceleration)}$$

- $\Lambda > 0$ (or any $p < -\rho c^2/3$) makes $\ddot a > 0$ → **accelerating expansion** (dark energy, discovered 1998).
- Set $k$ (spatial curvature) and the energy contents (matter, radiation, $\Lambda$) → you get the entire expansion history. This is the **ΛCDM** standard model.

## Physical predictions (all confirmed)

- **Mercury's perihelion precession** — the 43″/century Newton couldn't explain (Einstein's first win, 1915).
- **Light bending** — starlight deflected by the Sun (Eddington, 1919). Basis of **gravitational lensing**.
- **Gravitational time dilation** — clocks run slower in stronger gravity. **GPS** must correct for it ($\sim 38\ \mu\text{s/day}$) or navigation drifts ~10 km/day.
- **Black holes & event horizons** — imaged directly (EHT, M87\*, 2019).
- **Gravitational waves** — merging black holes (LIGO, 2015; Nobel 2017).
- **Expanding universe / Big Bang** — from the FLRW solutions.

## Mental models & gotchas

- **It's geometry, not force.** A freely-falling object feels *no* force — it follows a **geodesic** (straightest possible path) through curved spacetime. The "force" of gravity is an artifact of using non-inertial coordinates.
- **The unknown is the metric $g_{\mu\nu}$**, and it appears *nonlinearly* (in $\Gamma$, $R$, and even to invert $g^{\mu\nu}$). Gravity gravitates — that nonlinearity is why the equations are so hard and why GW's don't simply superpose.
- **Coordinates lie.** The Schwarzschild $r = r_s$ "singularity" is a **coordinate** artifact (fixable by changing coordinates); $r = 0$ is a **real** (curvature) singularity. Compute a scalar like the **Kretschmann scalar** $R_{\mu\nu\rho\sigma} R^{\mu\nu\rho\sigma}$ to tell them apart — it blows up only at true singularities.
- **Ricci-flat ≠ flat.** Vacuum means $R_{\mu\nu} = 0$, but Weyl curvature can be huge (tidal forces, waves).
- **10 equations, but not 10 free functions.** 4 are constraints and coordinate freedom (diffeomorphism/gauge) removes 4 more → really ~2 physical degrees of freedom (the two GW polarizations).
- **Index gymnastics.** Raise/lower with $g$, sum repeated indices (Einstein summation). By hand this is error-prone — use a CAS (see samples) or a package like `EinsteinPy` / `sympy.diffgeom`.
- **Units.** Physicists set $c = G = 1$ (geometrized units) to declutter, so you'll see $G_{\mu\nu} = 8\pi T_{\mu\nu}$. Restore factors of $c$, $G$ at the end.

## Samples in this folder

- [`curvature_sympy.py`](curvature_sympy.py) — the full metric→curvature pipeline in sympy. Verifies **Schwarzschild is a vacuum solution** ($G_{\mu\nu} = 0$) and computes the Kretschmann scalar showing $r = 0$ is the real singularity.
- [`friedmann_flrw.py`](friedmann_flrw.py) — feed the FLRW metric to the same machinery and watch the 10 EFE collapse into the **two Friedmann equations**.
