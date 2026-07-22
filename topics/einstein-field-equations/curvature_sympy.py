"""Metric -> curvature pipeline in sympy, applied to Schwarzschild.

Runs the whole chain from the README:
    g -> Christoffel Gamma -> Riemann -> Ricci -> scalar R -> Einstein G

and shows two headline facts about the Schwarzschild black hole:

  1. It is a VACUUM solution: G_mu_nu = 0 everywhere (outside r=0).
     (So Ricci = 0 too -- "Ricci-flat" -- yet spacetime is NOT flat.)
  2. The Kretschmann scalar R_abcd R^abcd = 48 G^2 M^2 / (c^4 r^6)
     blows up only at r=0 -> that is the *real* singularity, while the
     r = r_s = 2GM/c^2 "horizon" is just a coordinate artifact (finite there).

Run:  uv run topics/einstein-field-equations/curvature_sympy.py
"""

import sympy as sp

# ---------------------------------------------------------------------------
# Generic tensor machinery: build everything from a metric g and coords x.
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


def einstein(g, x):
    """Full pipeline -> (Einstein G_mu_nu, Ricci, scalar R, Christoffel, Riemann)."""
    Gamma = christoffel(g, x)
    Riem = riemann(Gamma, x)
    Rc = ricci(Riem, x)
    R_scalar = sp.simplify(
        sum((g.inv())[i, j] * Rc[i, j] for i in range(len(x)) for j in range(len(x)))
    )
    G = sp.simplify(Rc - sp.Rational(1, 2) * R_scalar * g)
    return G, Rc, R_scalar, Gamma, Riem


def kretschmann(Riem, g, x):
    """K = R_{abcd} R^{abcd}. Lower the top index of Riemann, then contract twice."""
    n = len(x)
    # R_{a b c d} = g_{a e} R^e_{b c d}
    R_low = [
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
    g_inv = g.inv()
    K = 0
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    # raise all four indices on the second copy
                    up = 0
                    for e in range(n):
                        for f in range(n):
                            for gg in range(n):
                                for h in range(n):
                                    up += (
                                        g_inv[a, e]
                                        * g_inv[b, f]
                                        * g_inv[c, gg]
                                        * g_inv[d, h]
                                        * R_low[e][f][gg][h]
                                    )
                    K += R_low[a][b][c][d] * up
    return sp.simplify(K)


# ---------------------------------------------------------------------------
# The Schwarzschild metric (geometrized units G = c = 1, so r_s = 2M).
#   ds^2 = -(1 - 2M/r) dt^2 + (1 - 2M/r)^-1 dr^2 + r^2 dtheta^2 + r^2 sin^2(theta) dphi^2
# ---------------------------------------------------------------------------


def main():
    t, r, th, ph, M = sp.symbols("t r theta phi M", positive=True)
    x = [t, r, th, ph]
    f = 1 - 2 * M / r

    g = sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(th) ** 2)

    print("Schwarzschild metric g_mu_nu (diagonal):")
    sp.pprint(g)

    G, Rc, R_scalar, _, Riem = einstein(g, x)

    print("\nEinstein tensor G_mu_nu:")
    sp.pprint(G)
    is_vacuum = G == sp.zeros(4, 4)
    print(f"\n=> G_mu_nu == 0 ?  {is_vacuum}   (vacuum solution: matter-free)")
    print(f"   Ricci scalar R = {R_scalar}")
    print(
        "   Ricci = 0 means 'Ricci-flat', but spacetime is NOT flat -> see Kretschmann below."
    )

    K = kretschmann(Riem, g, x)
    print(f"\nKretschmann scalar  R_abcd R^abcd = {K}")
    print("   -> ~ 1/r^6: finite at the horizon r = 2M, diverges only at r = 0.")
    print(
        "      So r=0 is the true (curvature) singularity; the horizon is a coordinate artifact."
    )


if __name__ == "__main__":
    main()
