"""Spherical vs hyperbolic, checked symbolically with sympy.

Runs the metric -> Christoffel -> Riemann -> Ricci -> R pipeline on the
constant-curvature spaces and verifies the maximally-symmetric identity

    R_{abcd} = K (g_ac g_bd - g_ad g_bc),  R_ab = (n-1) K g_ab,  R = n(n-1) K

for, in order:

  * S^2 of radius a         -> K = +1/a^2   (spherical)
  * H^2 of radius a         -> K = -1/a^2   (hyperbolic), in two different
                               coordinate charts, to show K is chart-independent
  * the FLRW spatial slice with k = +1/0/-1 kept symbolic  -> K = k
  * 4D de Sitter / anti-de Sitter static patches -> K = +-1/L^2, R = +-12/L^2

Run:  uv run topics/spacetime-curvature-hyperbolic-vs-spherical/constant_curvature.py
"""

import sympy as sp

# ---------------------------------------------------------------------------
# Generic machinery (same conventions as topics/einstein-field-equations).
# ---------------------------------------------------------------------------


def christoffel(g, x):
    """Gamma^l_{m n} = 1/2 g^{l s} (d_m g_{s n} + d_n g_{s m} - d_s g_{m n})."""
    n = len(x)
    g_inv = g.inv()
    Gamma = [[[0] * n for _ in range(n)] for _ in range(n)]
    for lam in range(n):
        for m in range(n):
            for nu in range(n):
                s = 0
                for sig in range(n):
                    s += g_inv[lam, sig] * (
                        sp.diff(g[sig, nu], x[m])
                        + sp.diff(g[sig, m], x[nu])
                        - sp.diff(g[m, nu], x[sig])
                    )
                Gamma[lam][m][nu] = sp.simplify(s / 2)
    return Gamma


def riemann(Gamma, x):
    """R^r_{s m n} = d_m G^r_{n s} - d_n G^r_{m s} + G^r_{m l} G^l_{n s} - G^r_{n l} G^l_{m s}."""
    n = len(x)
    R = [[[[0] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for s in range(n):
            for m in range(n):
                for nu in range(n):
                    term = sp.diff(Gamma[r][nu][s], x[m]) - sp.diff(
                        Gamma[r][m][s], x[nu]
                    )
                    for lam in range(n):
                        term += (
                            Gamma[r][m][lam] * Gamma[lam][nu][s]
                            - Gamma[r][nu][lam] * Gamma[lam][m][s]
                        )
                    R[r][s][m][nu] = sp.simplify(term)
    return R


def ricci(R, x):
    """R_{s n} = R^m_{s m n} (contract first and third index)."""
    n = len(x)
    Rc = sp.zeros(n, n)
    for s in range(n):
        for nu in range(n):
            Rc[s, nu] = sp.simplify(sum(R[m][s][m][nu] for m in range(n)))
    return Rc


def lower_first(Riem, g):
    """R_{a b c d} = g_{a e} R^e_{b c d}."""
    n = g.shape[0]
    return [
        [
            [
                [
                    sp.simplify(sum(g[a, e] * Riem[e][b][c][d] for e in range(n)))
                    for d in range(n)
                ]
                for c in range(n)
            ]
            for b in range(n)
        ]
        for a in range(n)
    ]


def is_zero(expr):
    """Robust 'is this identically zero?'.

    sympy's simplify alone leaves things like -sin(2t)tan(t) - cos(2t) + 1
    (which *is* zero) untouched, so expand the trig first and fall back to
    .equals() before believing a nonzero answer.
    """
    e = sp.simplify(expr)
    if e == 0:
        return True
    e = sp.simplify(sp.expand_trig(sp.expand(e)))
    return e == 0 or bool(e.equals(0))


def curvature_report(name, g, x, K_expected):
    """Full pipeline + the maximally-symmetric checks. Returns the Ricci scalar."""
    n = len(x)
    Gamma = christoffel(g, x)
    Riem = riemann(Gamma, x)
    Rc = ricci(Riem, x)
    g_inv = g.inv()
    R_scalar = sp.simplify(
        sum(g_inv[i, j] * Rc[i, j] for i in range(n) for j in range(n))
    )

    R_low = lower_first(Riem, g)
    # R_{abcd} - K (g_ac g_bd - g_ad g_bc) should vanish identically.
    max_sym = all(
        is_zero(
            R_low[a][b][c][d] - K_expected * (g[a, c] * g[b, d] - g[a, d] * g[b, c])
        )
        for a in range(n)
        for b in range(n)
        for c in range(n)
        for d in range(n)
    )
    ricci_ok = all(
        is_zero(Rc[i, j] - (n - 1) * K_expected * g[i, j])
        for i in range(n)
        for j in range(n)
    )
    scalar_ok = is_zero(R_scalar - n * (n - 1) * K_expected)

    print(f"\n--- {name}  (n = {n}) ---")
    print(f"  ds^2 diag       : {[sp.simplify(g[i, i]) for i in range(n)]}")
    print(f"  sectional K     : {sp.simplify(K_expected)}")
    print(f"  Ricci scalar R  : {R_scalar}")
    print(f"  R_abcd = K(g_ac g_bd - g_ad g_bc) ? {max_sym}")
    print(f"  R_ab   = (n-1) K g_ab             ? {ricci_ok}")
    print(f"  R      = n(n-1) K                 ? {scalar_ok}")
    return R_scalar


def main():
    a, L = sp.symbols("a L", positive=True)

    # -- 2D: the sphere and the hyperbolic plane, same metric with sin -> sinh --
    th, ph = sp.symbols("theta phi", positive=True)
    curvature_report(
        "S^2 radius a:  ds^2 = a^2 (dth^2 + sin^2 th dph^2)",
        sp.diag(a**2, a**2 * sp.sin(th) ** 2),
        [th, ph],
        1 / a**2,
    )

    chi = sp.Symbol("chi", positive=True)
    curvature_report(
        "H^2 radius a:  ds^2 = a^2 (dchi^2 + sinh^2 chi dph^2)",
        sp.diag(a**2, a**2 * sp.sinh(chi) ** 2),
        [chi, ph],
        -1 / a**2,
    )

    # Same H^2, completely different chart (Poincare upper half-plane):
    # curvature is a property of the geometry, not of the coordinates.
    ux, uy = sp.symbols("u_x u_y", positive=True)
    curvature_report(
        "H^2 again, Poincare half-plane:  ds^2 = a^2 (dx^2 + dy^2)/y^2",
        sp.diag(a**2 / uy**2, a**2 / uy**2),
        [ux, uy],
        -1 / a**2,
    )

    # -- 3D FLRW spatial slice, k symbolic: one calculation covers all three --
    k, r = sp.symbols("k r")
    curvature_report(
        "FLRW slice:  ds^2 = dr^2/(1-k r^2) + r^2 dOmega^2  (k = +1/0/-1)",
        sp.diag(1 / (1 - k * r**2), r**2, r**2 * sp.sin(th) ** 2),
        [r, th, ph],
        k,
    )
    print("  -> k = +1 spherical (closed), k = 0 flat, k = -1 hyperbolic (open).")

    # -- 4D Lorentzian: de Sitter and anti-de Sitter static patches --
    t = sp.Symbol("t")
    rr = sp.Symbol("r", positive=True)
    for name, sign in (
        ("de Sitter (Lambda > 0)", +1),
        ("anti-de Sitter (Lambda < 0)", -1),
    ):
        f = 1 - sign * rr**2 / L**2
        R_scalar = curvature_report(
            f"{name}:  f(r) = 1 - ({sign:+d}) r^2/L^2",
            sp.diag(-f, 1 / f, rr**2, rr**2 * sp.sin(th) ** 2),
            [t, rr, th, ph],
            sign / L**2,
        )
        Lam = sp.simplify(R_scalar / 4)  # in 4D vacuum-with-Lambda, R = 4 Lambda
        print(f"  Lambda = R/4    : {Lam}   (= 3K)")

    print(
        "\nTakeaway: one sign flip (sin <-> sinh, +1/L^2 <-> -1/L^2) is the whole\n"
        "difference between spherical and hyperbolic. What the sign *does* to\n"
        "geodesics -- and why it flips again for timelike ones -- is in\n"
        "geodesic_deviation.py."
    )


if __name__ == "__main__":
    main()
