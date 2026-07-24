"""Curves and surfaces: Frenet-Serret, the two fundamental forms, and Theorema Egregium.

Covers Parts 3-5 of the README.

  1. Curves in R^3: curvature kappa and torsion tau of a helix.
  2. Surfaces in R^3: first fundamental form (E,F,G)  -> intrinsic,
     second fundamental form (L,M,N)                 -> extrinsic,
     Gaussian K = (LN-M^2)/(EG-F^2), mean H.
  3. THEOREMA EGREGIUM: recompute K from the first fundamental form *alone*
     (metric -> Christoffel -> Riemann -> Ricci scalar, then K = R/2) and
     check it equals the extrinsic answer. The cylinder is the punchline:
     it is bent in R^3 (kappa_1 = 1/a) yet intrinsically flat (K = 0).

Run:  uv run topics/geometry-to-spacetime-curvature/curves_and_surfaces.py
"""

import sympy as sp

# ---------------------------------------------------------------------------
# Part 3 -- curves
# ---------------------------------------------------------------------------


def frenet(r, t):
    """Curvature and torsion of a space curve r(t), for an arbitrary parameter t.

    kappa = |r' x r''| / |r'|^3
    tau   = (r' x r'') . r''' / |r' x r''|^2
    """
    r1, r2, r3 = (sp.diff(r, t, n) for n in (1, 2, 3))
    cross = r1.cross(r2)
    kappa = sp.simplify(cross.norm() / r1.norm() ** 3)
    tau = sp.simplify(cross.dot(r3) / cross.dot(cross))
    return kappa, tau


def curves_demo():
    print("=" * 72)
    print("PART 3 -- curves: Frenet-Serret")
    print("=" * 72)

    t, a, b = sp.symbols("t a b", positive=True)

    circle = sp.Matrix([a * sp.cos(t), a * sp.sin(t), 0])
    kappa, tau = frenet(circle, t)
    print(f"circle radius a : kappa = {kappa}   tau = {tau}")
    print("                  -> constant curvature 1/a, no torsion (planar)")

    helix = sp.Matrix([a * sp.cos(t), a * sp.sin(t), b * t])
    kappa, tau = frenet(helix, t)
    print(f"helix (a, pitch b): kappa = {kappa}   tau = {tau}")
    print("                  -> both constant; b -> 0 recovers the circle")
    print(
        "\nNOTE: this curvature is EXTRINSIC. Intrinsically every curve is just\n"
        "      the real line (arclength is the only measurement available).\n"
        "      Bending a wire changes nothing an ant living on it could detect."
    )


# ---------------------------------------------------------------------------
# Part 4 -- surfaces: the two fundamental forms
# ---------------------------------------------------------------------------


def fundamental_forms(r, u, v, assume=True):
    """(E,F,G), (L,M,N) for a parametrized surface r(u, v) in R^3.

    `assume` is a sympy assumption (e.g. Q.positive(sin(u))) used to resolve the
    Abs() that appears when normalising the normal on a chart of finite range.
    """

    def tidy(e):
        return sp.simplify(sp.refine(sp.simplify(e), assume))

    ru, rv = sp.diff(r, u), sp.diff(r, v)
    E, F, G = tidy(ru.dot(ru)), tidy(ru.dot(rv)), tidy(rv.dot(rv))

    normal = ru.cross(rv)
    n = sp.simplify(normal / normal.norm())

    L = tidy(sp.diff(r, u, 2).dot(n))
    M = tidy(sp.diff(r, u, v).dot(n))
    N = tidy(sp.diff(r, v, 2).dot(n))
    return (E, F, G), (L, M, N)


def gaussian_mean(first, second):
    """K = (LN - M^2)/(EG - F^2)   and   H = (EN - 2FM + GL) / (2(EG - F^2))."""
    E, F, G = first
    L, M, N = second
    det_I = E * G - F**2
    K = sp.simplify((L * N - M**2) / det_I)
    H = sp.simplify((E * N - 2 * F * M + G * L) / (2 * det_I))
    return K, H


# ---------------------------------------------------------------------------
# Part 5 -- the intrinsic route: metric -> Christoffel -> Riemann -> R = 2K
# ---------------------------------------------------------------------------


def christoffel(g, x):
    """Gamma^l_{m n} = 1/2 g^{l s} (d_m g_{s n} + d_n g_{s m} - d_s g_{m n})."""
    n = len(x)
    g_inv = g.inv()
    Gamma = [[[0] * n for _ in range(n)] for _ in range(n)]
    for lam in range(n):
        for m in range(n):
            for nu in range(n):
                s = sum(
                    g_inv[lam, sig]
                    * (
                        sp.diff(g[sig, nu], x[m])
                        + sp.diff(g[sig, m], x[nu])
                        - sp.diff(g[m, nu], x[sig])
                    )
                    for sig in range(n)
                )
                Gamma[lam][m][nu] = sp.simplify(s / 2)
    return Gamma


def ricci_scalar(g, x):
    """Riemann -> Ricci -> R, built from the metric only (no embedding used)."""
    n = len(x)
    Gamma = christoffel(g, x)

    def riem(r, s, m, nu):
        term = sp.diff(Gamma[r][nu][s], x[m]) - sp.diff(Gamma[r][m][s], x[nu])
        term += sum(
            Gamma[r][m][lam] * Gamma[lam][nu][s] - Gamma[r][nu][lam] * Gamma[lam][m][s]
            for lam in range(n)
        )
        return sp.simplify(term)

    Rc = sp.Matrix(
        n, n, lambda s, nu: sp.simplify(sum(riem(m, s, m, nu) for m in range(n)))
    )
    g_inv = g.inv()
    return sp.simplify(sum(g_inv[i, j] * Rc[i, j] for i in range(n) for j in range(n)))


def surfaces_demo():
    print()
    print("=" * 72)
    print("PARTS 4-5 -- surfaces: fundamental forms, K and H, Theorema Egregium")
    print("=" * 72)

    u, v = sp.symbols("u v", real=True)
    a = sp.Symbol("a", positive=True)
    R = sp.Symbol("R", positive=True)

    # name -> (parametrization, chart assumption used to drop Abs())
    surfaces = {
        "plane": (sp.Matrix([u, v, 0]), True),
        "cylinder (radius a)": (sp.Matrix([a * sp.cos(u), a * sp.sin(u), v]), True),
        "sphere (radius a)": (
            sp.Matrix(
                [
                    a * sp.sin(u) * sp.cos(v),
                    a * sp.sin(u) * sp.sin(v),
                    a * sp.cos(u),
                ]
            ),
            sp.Q.positive(sp.sin(u)),  # colatitude u in (0, pi)
        ),
        "torus (R, a)": (
            sp.Matrix(
                [
                    (R + a * sp.cos(v)) * sp.cos(u),
                    (R + a * sp.cos(v)) * sp.sin(u),
                    a * sp.sin(v),
                ]
            ),
            sp.Q.positive(R + a * sp.cos(v)),  # R > a: the tube never self-touches
        ),
        "catenoid (a)": (
            sp.Matrix(
                [
                    a * sp.cosh(v / a) * sp.cos(u),
                    a * sp.cosh(v / a) * sp.sin(u),
                    v,
                ]
            ),
            True,
        ),
    }

    for name, (r, assume) in surfaces.items():
        first, second = fundamental_forms(r, u, v, assume)
        K_ext, H = gaussian_mean(first, second)

        # Theorema Egregium: same K from the first fundamental form alone.
        E, F, G = first
        metric = sp.Matrix([[E, F], [F, G]])
        K_int = sp.simplify(ricci_scalar(metric, [u, v]) / 2)

        agree = sp.simplify(K_ext - K_int) == 0

        print(f"\n{name}")
        print(f"  I  : E={first[0]}, F={first[1]}, G={first[2]}")
        print(f"  II : L={second[0]}, M={second[1]}, N={second[2]}")
        print(f"  K (extrinsic, needs II) = {K_ext}")
        print(f"  K (intrinsic,  R/2 )    = {K_int}    <- agrees: {agree}")
        print(f"  H (mean, genuinely extrinsic) = {H}")

    print(
        "\nTHEOREMA EGREGIUM holds in every case: K never needed the embedding.\n"
        "  * cylinder: K = 0 although it is visibly bent (kappa_1 = 1/a, kappa_2 = 0)\n"
        "      -> a sheet of paper rolls into a cylinder without stretching;\n"
        "         an ant on it measures ordinary Euclidean geometry.\n"
        "  * sphere:   K = 1/a^2 > 0, and no flat map of the Earth can exist.\n"
        "  * torus:    K = cos(v) / (a (R + a cos v)) -- positive on the outer rim,\n"
        "      zero on the top/bottom circles, negative on the inner rim. It must\n"
        "      average to zero: Gauss-Bonnet gives int K dA = 2 pi chi = 0.\n"
        "  * catenoid: H = 0 (a minimal surface -- what a soap film does) but K < 0."
    )


if __name__ == "__main__":
    curves_demo()
    surfaces_demo()
