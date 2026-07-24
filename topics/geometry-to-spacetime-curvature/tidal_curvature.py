"""Tidal forces: what the Riemann tensor physically DOES.

Covers Parts 10-11 of the README.

A single freely-falling observer feels nothing at all (equivalence principle).
Two of them, nearby, measure a relative acceleration -- and that relative
acceleration IS the curvature:

    D^2 xi^m / dtau^2  =  - R^m_{n r s} u^n xi^r u^s

Contracting with the observer's 4-velocity defines the TIDAL TENSOR
E_ij = R_(i 0 j 0), the relativistic version of the Newtonian d_i d_j Phi.

What this script does:

  1. Builds the full Riemann tensor of Schwarzschild symbolically, projects it
     onto an observer's orthonormal frame, and gets

         E = diag(-2M/r^3, +M/r^3, +M/r^3)

     -- stretch radially, squeeze transversally: spaghettification.
  2. Shows this is the SAME for a static observer and for one falling in at
     any speed: tidal forces are frame-invariant in the radial plane.
  3. Shows tr(E) = 0 exactly, i.e. Ricci = 0 -- vacuum. A falling ball of dust
     changes SHAPE but not VOLUME. That is the Ricci/Weyl split made physical.
  4. Contrasts with a ball falling THROUGH matter, where tr(E) = 4 pi rho > 0
     and the volume really does shrink.
  5. Real numbers: ocean tides, spaghettification at stellar and supermassive
     black-hole horizons.

Run:  uv run topics/geometry-to-spacetime-curvature/tidal_curvature.py
"""

import numpy as np
import sympy as sp

# ---------------------------------------------------------------------------
# metric -> Christoffel -> Riemann (all indices down)
# ---------------------------------------------------------------------------


def christoffel(g, x):
    n = len(x)
    g_inv = g.inv()
    Gamma = [[[sp.S(0)] * n for _ in range(n)] for _ in range(n)]
    for lam in range(n):
        for m in range(n):
            for nu in range(n):
                Gamma[lam][m][nu] = sp.simplify(
                    sum(
                        g_inv[lam, sig]
                        * (
                            sp.diff(g[sig, nu], x[m])
                            + sp.diff(g[sig, m], x[nu])
                            - sp.diff(g[m, nu], x[sig])
                        )
                        for sig in range(n)
                    )
                    / 2
                )
    return Gamma


def riemann_down(g, x):
    """R_{a b c d}, fully lowered."""
    n = len(x)
    Gamma = christoffel(g, x)

    def up(r, s, m, nu):
        term = sp.diff(Gamma[r][nu][s], x[m]) - sp.diff(Gamma[r][m][s], x[nu])
        term += sum(
            Gamma[r][m][lam] * Gamma[lam][nu][s] - Gamma[r][nu][lam] * Gamma[lam][m][s]
            for lam in range(n)
        )
        return sp.simplify(term)

    R_up = [
        [[[up(a, b, c, d) for d in range(n)] for c in range(n)] for b in range(n)]
        for a in range(n)
    ]
    return [
        [
            [
                [
                    sp.simplify(sum(g[a, e] * R_up[e][b][c][d] for e in range(n)))
                    for d in range(n)
                ]
                for c in range(n)
            ]
            for b in range(n)
        ]
        for a in range(n)
    ], R_up


def ricci_from(R_up, n):
    return sp.Matrix(
        n, n, lambda s, nu: sp.simplify(sum(R_up[m][s][m][nu] for m in range(n)))
    )


# ---------------------------------------------------------------------------
# 1-3. Schwarzschild tidal tensor
# ---------------------------------------------------------------------------


def schwarzschild_tidal():
    print("=" * 72)
    print("1. THE TIDAL TENSOR OF SCHWARZSCHILD  (G = c = 1)")
    print("=" * 72)

    t, r, th, ph, M, v = sp.symbols("t r theta phi M v", positive=True)
    x = [t, r, th, ph]
    f = 1 - 2 * M / r
    g = sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(th) ** 2)

    print("   ds^2 = -(1-2M/r) dt^2 + dr^2/(1-2M/r) + r^2 dOmega^2")
    R, R_up = riemann_down(g, x)
    Rc = ricci_from(R_up, 4)
    print(f"   Ricci tensor R_mn = 0 ?  {Rc == sp.zeros(4, 4)}   (vacuum solution)")

    # Orthonormal frame of a STATIC observer at (r, theta = pi/2):
    #   e_0 = f^(-1/2) d_t,  e_1 = f^(1/2) d_r,  e_2 = (1/r) d_th,  e_3 = (1/r) d_ph
    e0 = sp.Matrix([1 / sp.sqrt(f), 0, 0, 0])
    e1 = sp.Matrix([0, sp.sqrt(f), 0, 0])
    e2 = sp.Matrix([0, 0, 1 / r, 0])
    e3 = sp.Matrix([0, 0, 0, 1 / r])  # at theta = pi/2, sin(theta) = 1

    def nz(vec):
        """Nonzero components of a frame vector -- the frames here are sparse."""
        return [(i, vec[i]) for i in range(4) if vec[i] != 0]

    def tidal(u, spatial):
        """E_ij = R_{a b c d} e_i^a u^b e_j^c u^d, the tidal tensor in a frame."""
        out = sp.zeros(3, 3)
        for i, ei in enumerate(spatial):
            for j, ej in enumerate(spatial):
                out[i, j] = sp.simplify(
                    sum(
                        R[a][b][c][d] * va * vb * vc * vd
                        for a, va in nz(ei)
                        for b, vb in nz(u)
                        for c, vc in nz(ej)
                        for d, vd in nz(u)
                    )
                )
        return sp.simplify(out.subs(th, sp.pi / 2))

    E_static = tidal(e0, [e1, e2, e3])
    print("\n   static observer, tidal tensor E_ij = R_(i 0 j 0):")
    sp.pprint(E_static)
    print(f"\n   trace tr(E) = {sp.simplify(E_static.trace())}   <- ZERO")
    print("   -> E = diag(-2M/r^3, +M/r^3, +M/r^3): stretched along r,")
    print("      squeezed in both transverse directions. Spaghettification.")
    print("   -> tr(E) = R_00 = 0 is exactly the vacuum Einstein equation:")
    print("      the SHAPE of a falling ball of dust distorts, its VOLUME does not.")

    print()
    print("=" * 72)
    print("2. SAME TENSOR FOR AN INFALLING OBSERVER -- boost invariance")
    print("=" * 72)
    gamma = 1 / sp.sqrt(1 - v**2)
    u_boost = sp.simplify(gamma * (e0 + v * e1))  # radial infall at speed v
    e1_boost = sp.simplify(gamma * (e1 + v * e0))
    E_boost = tidal(u_boost, [e1_boost, e2, e3])
    print("   observer falling radially at local speed v, frame boosted by gamma(v):")
    sp.pprint(E_boost)
    print(
        f"\n   identical to the static frame ?  {sp.simplify(E_boost - E_static) == sp.zeros(3, 3)}"
    )
    print("   -> the tidal field does not depend on how fast you cross it, so the")
    print("      Newtonian-looking answer 2GM/r^3 is EXACT here, not an approximation.")

    print()
    print("=" * 72)
    print("3. THE NEWTONIAN LIMIT: E_ij = d_i d_j Phi")
    print("=" * 72)
    rr = sp.Symbol("r", positive=True)
    Phi = -M / rr  # Newtonian potential of a point mass
    d_rr = sp.diff(Phi, rr, 2)
    d_tt = sp.simplify(sp.diff(Phi, rr) / rr)  # transverse piece in spherical coords
    print(f"   Phi = -M/r  ->  d_r d_r Phi = {d_rr},  transverse = {d_tt}")
    print(
        f"   Laplacian  = {sp.simplify(d_rr + 2 * d_tt)}   (= 4 pi rho = 0 in vacuum)"
    )
    print("   -> exactly the eigenvalues of E above. Riemann IS the tidal tensor;")
    print("      Newtonian tides were general relativity all along.")

    # K = R_abcd R^abcd. The metric is diagonal, so raising an index is just a
    # factor of g^{aa} -- no sums needed.
    gi = [g.inv()[a, a] for a in range(4)]
    kretschmann = sp.simplify(
        sum(
            gi[a] * gi[b] * gi[c] * gi[d] * R[a][b][c][d] ** 2
            for a in range(4)
            for b in range(4)
            for c in range(4)
            for d in range(4)
        )
    )
    print(f"\n   Kretschmann scalar R_abcd R^abcd = {kretschmann}")
    print("   -> Ricci = 0 but Riemann != 0. The curvature outside a black hole is")
    print("      pure WEYL: tidal distortion, the part that propagates as waves.")


# ---------------------------------------------------------------------------
# 4. a ball of dust: shape vs volume
# ---------------------------------------------------------------------------


def ball_of_dust():
    print()
    print("=" * 72)
    print("4. A FALLING BALL OF DUST: shape distorts, volume does not")
    print("=" * 72)
    print("   solve  d^2 xi / dtau^2 = -E xi  for an initially spherical cloud.")
    print("   With E constant and diagonal each axis evolves independently:")
    print("       lambda > 0 -> cos(sqrt(lambda) tau)    (squeezed)")
    print("       lambda < 0 -> cosh(sqrt(-lambda) tau)  (stretched)")
    print("   The volume of the cloud is the product of the three axes.\n")

    def axes_of(E, tau):
        vals = np.linalg.eigvalsh(E)
        return np.where(
            vals >= 0,
            np.cos(np.sqrt(np.abs(vals)) * tau),
            np.cosh(np.sqrt(np.abs(vals)) * tau),
        )

    def table(E, taus, label):
        print(
            f"       {'tau':>8} {'stretch':>10} {'squeeze':>10} {'volume':>12}"
            f" {'(V-1)/tau^2':>14}"
        )
        for tau in taus:
            ax = np.sort(axes_of(E, tau))[::-1]
            vol = float(np.prod(ax))
            ratio = "-" if tau == 0 else f"{(vol - 1) / tau**2:>14.3e}"
            print(
                f"       {tau:>8.1f} {ax[0]:>10.6f} {ax[-1]:>10.6f} {vol:>12.8f}"
                f" {ratio:>14}"
            )
        print(
            f"       expected (V-1)/tau^2 -> -tr(E)/2 = {-np.trace(E) / 2:.3e}"
            f"   [{label}]"
        )

    # (a) vacuum: outside a black hole, at r = 20 M, in units G = c = M = 1
    r = 20.0
    E_vac = np.diag([-2 / r**3, 1 / r**3, 1 / r**3])
    print(
        f"   (a) VACUUM at r = {r:.0f}M:  E = diag(-2/r^3, 1/r^3, 1/r^3),  "
        f"tr(E) = {np.trace(E_vac):.1e}"
    )
    table(E_vac, [0, 5, 10, 20, 40], "vacuum")
    print("       -> the cloud is drawn out into a cigar, yet (V-1)/tau^2 -> 0")
    print("          (it falls by 4x every time tau halves, so V - 1 = O(tau^4)):")
    print("          the volume is stationary. Volume change is sourced by RICCI,")
    print("          and Ricci = 0 out here -- the distortion is pure WEYL.")

    # (b) inside matter: tr(E) = R_00 = 4 pi rho (isotropic for uniform dust)
    rho = 3e-6
    E_mat = np.diag([4 * np.pi * rho / 3] * 3)
    print(
        f"\n   (b) INSIDE UNIFORM DUST, rho = {rho:.0e}:  tr(E) = 4 pi rho = "
        f"{np.trace(E_mat):.3e}"
    )
    table(E_mat, [0, 5, 10, 20, 40], "matter")
    print("       -> no shape change at all, but the volume genuinely shrinks, at")
    print("          exactly the rate -tr(E)/2 = -2 pi rho. This is Raychaudhuri's")
    print("          focusing theorem, and it is why matter enters the field")
    print("          equations through the Ricci (trace) part of the curvature.")


# ---------------------------------------------------------------------------
# 5. real tidal numbers
# ---------------------------------------------------------------------------


def real_numbers():
    print()
    print("=" * 72)
    print("5. TIDAL ACCELERATIONS IN THE REAL WORLD")
    print("=" * 72)
    print("   differential acceleration across a separation d:  a = 2 G M d / r^3\n")

    G, c, M_SUN, g_earth = 6.674e-11, 2.998e8, 1.989e30, 9.81

    cases = [
        ("Moon on Earth's oceans", 7.35e22, 3.844e8, 1.274e7),
        ("Sun on Earth's oceans", M_SUN, 1.496e11, 1.274e7),
        ("Earth on a 2 m person (surface)", 5.972e24, 6.371e6, 2.0),
        ("10 Msun black hole, at its horizon", 10 * M_SUN, None, 2.0),
        ("1000 Msun black hole, at its horizon", 1e3 * M_SUN, None, 2.0),
        ("Sgr A* (4.3e6 Msun), at its horizon", 4.3e6 * M_SUN, None, 2.0),
        ("M87* (6.5e9 Msun), at its horizon", 6.5e9 * M_SUN, None, 2.0),
    ]

    print(f"   {'system':<38} {'r [m]':>11} {'tidal a [m/s^2]':>17} {'in g':>12}")
    print(f"   {'-' * 38} {'-' * 11} {'-' * 17} {'-' * 12}")
    accel = {}
    for name, M, r, d in cases:
        if r is None:
            r = 2 * G * M / c**2  # Schwarzschild horizon
        a = 2 * G * M * d / r**3
        accel[name] = a
        print(f"   {name:<38} {r:>11.3e} {a:>17.3e} {a / g_earth:>12.2e}")

    moon = accel["Moon on Earth's oceans"]
    sun = accel["Sun on Earth's oceans"]
    print(f"\n   -> Moon / Sun tidal ratio = {moon / sun:.2f}. The Sun is far more")
    print("      massive but much further away, and tides fall off as 1/r^3, so the")
    print("      Moon wins. Spring tides add the two, neap tides subtract:")
    print(f"      spring / neap = {(moon + sun) / (moon - sun):.1f}.")
    print("   -> horizon tides scale as M/r_s^3 ~ 1/M^2: a stellar-mass black hole")
    print("      shreds you long before the horizon, while at a supermassive one")
    print("      you cross the horizon feeling essentially nothing. Nothing local")
    print("      marks the horizon -- exactly as the equivalence principle demands.")


if __name__ == "__main__":
    schwarzschild_tidal()
    ball_of_dust()
    real_numbers()
