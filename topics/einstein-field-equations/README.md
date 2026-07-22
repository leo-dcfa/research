# Einstein field equations

The 10 coupled PDEs at the heart of **general relativity**. They say: **matter/energy tells spacetime how to curve; curvature tells matter how to move.** Gravity isn't a force — it's the geometry of curved 4D spacetime.

$$G_{\mu\nu} + \Lambda g_{\mu\nu} = \frac{8\pi G}{c^4}\, T_{\mu\nu}$$

- **Left = geometry** (how spacetime is curved), **right = matter** (what's in it). One equation, but the indices `μ,ν` run 0–3, so it's a 4×4 symmetric system → **10 independent equations**.
- Wheeler's summary: *"Spacetime tells matter how to move; matter tells spacetime how to curve."*

## What each symbol is

| Symbol | Name | What it is |
|---|---|---|
| `g_μν` | **Metric tensor** | The unknown you solve for. Encodes distances/time in curved spacetime. The "shape". |
| `G_μν` | **Einstein tensor** | `R_μν − ½ R g_μν`. A specific combination of curvature. Built so `∇^μ G_μν = 0` automatically. |
| `R_μν` | **Ricci tensor** | Curvature contracted once — how volumes shrink/grow along geodesics. |
| `R` | **Ricci scalar** | `g^μν R_μν`. Curvature boiled down to one number per point. |
| `T_μν` | **Stress–energy tensor** | The source: energy density, momentum, pressure, stress. |
| `Λ` | **Cosmological constant** | Energy of empty space → **dark energy**. Drives accelerating expansion. |
| `8πG/c⁴` | coupling | Tiny (~2×10⁻⁴³). Gravity is weak — you need a planet's worth of `T` to bend spacetime noticeably. |

## The chain of definitions (metric → curvature)

Everything is built from the metric by differentiation. This is the whole computation:

```
g_μν                                     # 1. the metric (given/assumed)
  ↓  ∂g, combine
Γ^λ_μν = ½ g^λσ (∂_μ g_σν + ∂_ν g_σμ − ∂_σ g_μν)   # 2. Christoffel symbols (connection)
  ↓  ∂Γ + ΓΓ
R^ρ_σμν = ∂_μ Γ^ρ_νσ − ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ − Γ^ρ_νλ Γ^λ_μσ   # 3. Riemann tensor
  ↓  contract ρ=μ
R_σν = R^μ_σμν                           # 4. Ricci tensor
  ↓  contract with g
R = g^σν R_σν                            # 5. Ricci scalar
  ↓  assemble
G_μν = R_μν − ½ R g_μν                    # 6. Einstein tensor
```

- **Γ (Christoffel symbols)** aren't a tensor — they're the "connection" that says how to differentiate vectors in curved space (they encode the fictitious/gravitational acceleration). They vanish in freely-falling (locally inertial) frames.
- **Riemann `R^ρ_σμν`** is the *full* curvature (20 independent components in 4D). It measures how a vector rotates when parallel-transported around a loop → **tidal forces**. Flat spacetime ⟺ Riemann = 0 everywhere.
- The EFE only constrains the **Ricci** part (10 components). The rest of Riemann — the **Weyl tensor** — is curvature that propagates in vacuum: **tidal distortion and gravitational waves**.

## Why *this* combination (the `−½ R g` and ∇G = 0)

- Conservation of energy–momentum is `∇^μ T_μν = 0` (local, automatic). So the geometry side must also be divergence-free.
- `∇^μ R_μν ≠ 0` in general — you **can't** just set `R_μν ∝ T_μν` (that was Einstein's earlier wrong guess).
- The **contracted Bianchi identity** gives `∇^μ (R_μν − ½ R g_μν) = 0` identically. So `G_μν` is the *unique* (up to `Λ`) curvature tensor that's automatically conserved. That's why it's the left-hand side.

## Key solutions (you almost never solve the general case)

The EFE are nonlinear and brutal. Progress comes from imposing **symmetry**:

| Solution | Assumption | Describes |
|---|---|---|
| **Minkowski** | `T=0`, flat | Special relativity / empty spacetime |
| **Schwarzschild** (1916) | vacuum, spherical, static | Non-spinning black hole / outside any star. Event horizon at `r_s = 2GM/c²` |
| **Kerr** (1963) | vacuum, spherical, **spinning** | Rotating black hole (all real ones). Has ergosphere |
| **Reissner–Nordström / Kerr–Newman** | + electric charge | Charged (spinning) black holes |
| **FLRW** | homogeneous, isotropic | The whole universe → **Friedmann equations**, Big Bang, expansion |
| **Gravitational waves** | linearized `g = η + h`, `\|h\|≪1` | Ripples in spacetime (LIGO, 2015) |

- **Vacuum EFE** (`T_μν = 0`, `Λ=0`): contracting shows it reduces to `R_μν = 0` (**Ricci-flat**) — *not* flat spacetime! Schwarzschild is Ricci-flat but Riemann ≠ 0 (that's the tidal gravity outside a black hole).
- **Schwarzschild → Newton**: in the weak-field, slow-motion limit the `00` equation becomes `∇²Φ = 4πG ρ` (Poisson) — you recover Newtonian gravity, with `g_00 ≈ −(1 + 2Φ/c²)`.

## FLRW / cosmology (where the EFE meets data)

Plug a homogeneous isotropic metric into the EFE and the 10 equations collapse to **two** (the Friedmann equations) for the scale factor `a(t)`:

```
(ȧ/a)² = (8πG/3) ρ  −  kc²/a²  +  Λc²/3        # expansion rate (Hubble H = ȧ/a)
ä/a    = −(4πG/3)(ρ + 3p/c²)   +  Λc²/3        # acceleration
```

- `Λ > 0` (or any `p < −ρc²/3`) makes `ä > 0` → **accelerating expansion** (dark energy, discovered 1998).
- Set `k` (spatial curvature) and the energy contents (matter, radiation, Λ) → you get the entire expansion history. This is the **ΛCDM** standard model.

## Physical predictions (all confirmed)

- **Mercury's perihelion precession** — the 43″/century Newton couldn't explain (Einstein's first win, 1915).
- **Light bending** — starlight deflected by the Sun (Eddington, 1919). Basis of **gravitational lensing**.
- **Gravitational time dilation** — clocks run slower in stronger gravity. **GPS** must correct for it (~38 μs/day) or navigation drifts ~10 km/day.
- **Black holes & event horizons** — imaged directly (EHT, M87\*, 2019).
- **Gravitational waves** — merging black holes (LIGO, 2015; Nobel 2017).
- **Expanding universe / Big Bang** — from the FLRW solutions.

## Mental models & gotchas

- **It's geometry, not force.** A freely-falling object feels *no* force — it follows a **geodesic** (straightest possible path) through curved spacetime. The "force" of gravity is an artifact of using non-inertial coordinates.
- **The unknown is the metric `g_μν`**, and it appears *nonlinearly* (in Γ, R, and even to invert `g^μν`). Gravity gravitates — that nonlinearity is why the equations are so hard and why GW's don't simply superpose.
- **Coordinates lie.** The Schwarzschild `r = r_s` "singularity" is a **coordinate** artifact (fixable by changing coordinates); `r = 0` is a **real** (curvature) singularity. Compute a scalar like the **Kretschmann scalar** `R_μνρσ R^μνρσ` to tell them apart — it blows up only at true singularities.
- **Ricci-flat ≠ flat.** Vacuum means `R_μν = 0`, but Weyl curvature can be huge (tidal forces, waves).
- **10 equations, but not 10 free functions.** 4 are constraints and coordinate freedom (diffeomorphism/gauge) removes 4 more → really ~2 physical degrees of freedom (the two GW polarizations).
- **Index gymnastics.** Raise/lower with `g`, sum repeated indices (Einstein summation). By hand this is error-prone — use a CAS (see samples) or a package like `EinsteinPy` / `sympy.diffgeom`.
- **Units.** Physicists set `c = G = 1` (geometrized units) to declutter, so you'll see `G_μν = 8π T_μν`. Restore factors of `c`, `G` at the end.

## Samples in this folder

- [`curvature_sympy.py`](curvature_sympy.py) — the full metric→curvature pipeline in sympy. Verifies **Schwarzschild is a vacuum solution** (`G_μν = 0`) and computes the Kretschmann scalar showing `r=0` is the real singularity.
- [`friedmann_flrw.py`](friedmann_flrw.py) — feed the FLRW metric to the same machinery and watch the 10 EFE collapse into the **two Friedmann equations**.
