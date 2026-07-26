"""What the sign of the curvature *does*: focusing vs defocusing.

Three numerical experiments, all in numpy:

  1. Riemannian Jacobi equation  xi'' = -K xi  ->  sin / linear / sinh.
     Spherical curvature refocuses geodesics at a conjugate point (s = pi L);
     hyperbolic curvature blows the separation up like e^{s/L}.

  2. Actual geodesics, integrated with RK4 on the 2D metric
     ds^2 = dchi^2 + s_k(chi)^2 dphi^2, with the separation measured by the
     exact invariant distance formula -> confirms experiment 1 is not a
     linearisation artifact.

  3. The Lorentzian sign flip. Along a *timelike* geodesic, u.u = -1 turns the
     deviation equation into xi'' = +K xi, so de Sitter (K > 0) spreads
     worldlines apart while anti-de Sitter (K < 0) refocuses them at pi L --
     the exact opposite of the Riemannian intuition.

Run:  uv run topics/spacetime-curvature-hyperbolic-vs-spherical/geodesic_deviation.py
"""

import numpy as np

L = 1.0  # curvature radius; K = +-1/L^2


# ---------------------------------------------------------------------------
# 1. Jacobi equation, integrated (not solved analytically).
# ---------------------------------------------------------------------------


def integrate_jacobi(K, s_max, sign, n=200_000):
    """RK4 on xi'' = sign * K * xi with xi(0) = 0, xi'(0) = 1.

    sign = -1 : Riemannian, separation of spacelike geodesics.
    sign = +1 : Lorentzian, separation of nearby timelike geodesics.
    """
    h = s_max / n
    s = np.linspace(0.0, s_max, n + 1)
    y = np.zeros((n + 1, 2))
    y[0] = [0.0, 1.0]  # xi, dxi/ds

    def f(v):
        return np.array([v[1], sign * K * v[0]])

    for i in range(n):
        k1 = f(y[i])
        k2 = f(y[i] + h / 2 * k1)
        k3 = f(y[i] + h / 2 * k2)
        k4 = f(y[i] + h * k3)
        y[i + 1] = y[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return s, y[:, 0]


def first_conjugate_point(s, xi):
    """First s > 0 where the separation returns to zero (geodesics refocus)."""
    sgn = np.sign(xi[1:])
    flips = np.nonzero(np.diff(sgn) != 0)[0]
    return s[flips[0] + 1] if len(flips) else None


def experiment_1():
    print("=" * 72)
    print("1. Riemannian geodesic deviation:  xi'' = -K xi,  xi(0)=0, xi'(0)=1")
    print("=" * 72)
    cases = [
        ("spherical  K=+1", +1.0),
        ("flat       K= 0", 0.0),
        ("hyperbolic K=-1", -1.0),
    ]
    print(f"\n{'geometry':<16}" + "".join(f"{f's={v}':>12}" for v in (0.5, 1, 2, 4, 8)))
    for name, K in cases:
        s, xi = integrate_jacobi(K, 8.0, sign=-1)
        row = "".join(f"{np.interp(v, s, xi):>12.4f}" for v in (0.5, 1, 2, 4, 8))
        print(f"{name:<16}{row}")

    s, xi = integrate_jacobi(1.0, 8.0, sign=-1)
    sc = first_conjugate_point(s, xi)
    print(
        f"\n  spherical : first conjugate point at s = {sc:.5f}   (exact pi L = {np.pi:.5f})"
    )
    print("              -> geodesics from a pole all meet again at the far pole;")
    print("                 past it they no longer minimise length (Bonnet-Myers).")

    s, xi = integrate_jacobi(-1.0, 8.0, sign=-1)
    # asymptotic growth rate = d(log xi)/ds, not log(xi)/s (sinh s ~ e^s / 2,
    # so the latter only converges to 1/L logarithmically slowly)
    lyap = (np.log(xi[-1]) - np.log(np.interp(7.0, s, xi))) / 1.0
    print(f"  hyperbolic: no zeros ever (Cartan-Hadamard); xi(8) = {xi[-1]:.1f}")
    print(f"              growth rate d log(xi)/ds -> {lyap:.5f}  (exact 1/L = 1.0)")
    print("              -> a Lyapunov exponent: geodesic flow is chaotic.")


# ---------------------------------------------------------------------------
# 2. The real thing: integrate geodesics, measure the invariant distance.
# ---------------------------------------------------------------------------


def s_k(chi, K):
    """sin / identity / sinh, unified."""
    if K > 0:
        return np.sin(chi)
    if K < 0:
        return np.sinh(chi)
    return chi


def ds_k(chi, K):
    if K > 0:
        return np.cos(chi)
    if K < 0:
        return np.cosh(chi)
    return np.ones_like(chi)


def geodesic(K, chi0, phi0, dchi0, dphi0, s_max, n=20_000):
    """RK4 on the geodesic equations of ds^2 = dchi^2 + s_k(chi)^2 dphi^2.

    chi'' = s_k s_k' phi'^2
    phi'' = -2 (s_k'/s_k) chi' phi'
    """

    def f(y):
        chi, phi, dchi, dphi = y
        sk, dsk = s_k(chi, K), ds_k(chi, K)
        return np.array([dchi, dphi, sk * dsk * dphi**2, -2 * dsk / sk * dchi * dphi])

    h = s_max / n
    y = np.array([chi0, phi0, dchi0, dphi0], dtype=float)
    for _ in range(n):
        k1 = f(y)
        k2 = f(y + h / 2 * k1)
        k3 = f(y + h / 2 * k2)
        k4 = f(y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return y


def invariant_distance(K, p, q):
    """Exact geodesic distance between (chi, phi) points, from the hyperboloid /
    sphere law of cosines."""
    (c1, f1), (c2, f2) = p, q
    if K > 0:
        return np.arccos(
            np.clip(
                np.cos(c1) * np.cos(c2) + np.sin(c1) * np.sin(c2) * np.cos(f1 - f2),
                -1,
                1,
            )
        )
    if K < 0:
        return np.arccosh(
            max(
                np.cosh(c1) * np.cosh(c2) - np.sinh(c1) * np.sinh(c2) * np.cos(f1 - f2),
                1.0,
            )
        )
    # flat polar coordinates
    return np.hypot(
        c1 * np.cos(f1) - c2 * np.cos(f2), c1 * np.sin(f1) - c2 * np.sin(f2)
    )


def experiment_2():
    print("\n" + "=" * 72)
    print("2. Two geodesics leaving one point at angle delta = 1e-3, integrated")
    print("=" * 72)
    delta = 1e-3
    print(f"\n{'arclength s':>12}{'spherical':>14}{'flat':>12}{'hyperbolic':>14}")
    print(f"{'':>12}{'(sin s)*d':>14}{'s*d':>12}{'(sinh s)*d':>14}")
    for s_max in (0.5, 1.0, 2.0, 3.0, 4.0):
        row = f"{s_max:>12.1f}"
        for K in (1.0, 0.0, -1.0):
            eps = 1e-6  # start just off the coordinate pole (chi = 0 is singular)
            a = geodesic(K, eps, 0.0, 1.0, 0.0, s_max)
            b = geodesic(K, eps, delta, 1.0, 0.0, s_max)
            d = invariant_distance(K, (a[0], a[1]), (b[0], b[1]))
            row += f"{d / delta:>14.5f}" if K != 0.0 else f"{d / delta:>12.5f}"
        print(row)
    print("\n  columns are separation/delta; compare sin(s), s, sinh(s):")
    print(
        "  s=3 ->  sin = {:.5f}   s = {:.5f}   sinh = {:.5f}".format(
            np.sin(3), 3.0, np.sinh(3)
        )
    )
    print("  spherical separation is already shrinking by s=2 (past the equator).")


# ---------------------------------------------------------------------------
# 3. Circles, triangles, and the Lorentzian flip.
# ---------------------------------------------------------------------------


def experiment_3():
    print("\n" + "=" * 72)
    print("3a. How much room there is: circle of geodesic radius r")
    print("=" * 72)
    print(
        f"\n{'r/L':>6}{'C_sph/2pi r':>14}{'C_hyp/2pi r':>14}{'A_sph/pi r^2':>14}{'A_hyp/pi r^2':>14}"
    )
    for r in (0.1, 0.5, 1.0, 2.0, 3.0):
        c_s, c_h = np.sin(r) / r, np.sinh(r) / r
        # A = 2 pi L^2 (1 - cos(r/L))  /  2 pi L^2 (cosh(r/L) - 1)
        a_s = 2 * (1 - np.cos(r)) / r**2
        a_h = 2 * (np.cosh(r) - 1) / r**2
        print(f"{r:>6.1f}{c_s:>14.4f}{c_h:>14.4f}{a_s:>14.4f}{a_h:>14.4f}")
    print("\n  spherical: less perimeter than Euclid allows (ratio < 1), and it")
    print("  vanishes at r = pi L. hyperbolic: exponentially more, forever.")

    print("\n" + "=" * 72)
    print("3b. Triangle angle sum = pi + K * Area")
    print("=" * 72)
    print(f"\n{'area/L^2':>10}{'sum (spherical)':>18}{'sum (hyperbolic)':>18}")
    for area in (0.01, 0.5, 1.0, 2.0, 3.0):
        print(f"{area:>10.2f}{np.pi + area:>18.4f}{np.pi - area:>18.4f}")
    print(f"\n  flat is always {np.pi:.4f}. Hyperbolic triangles have area < pi L^2")
    print("  (the angle sum cannot go below 0), so there is a biggest triangle.")

    print("\n" + "=" * 72)
    print("3c. The Lorentzian flip:  along timelike geodesics  xi'' = +K xi")
    print("=" * 72)
    print(
        f"\n{'proper time tau':>16}{'de Sitter K=+1':>18}{'Minkowski':>13}{'AdS K=-1':>12}"
    )
    s_ds, xi_ds = integrate_jacobi(+1.0, 8.0, sign=+1)  # xi'' = +K xi -> sinh
    s_f, xi_f = integrate_jacobi(0.0, 8.0, sign=+1)
    s_ads, xi_ads = integrate_jacobi(-1.0, 8.0, sign=+1)  # xi'' = -|K| xi -> sin
    for tau in (0.5, 1.0, 2.0, np.pi, 4.0):
        print(
            f"{tau:>16.4f}{np.interp(tau, s_ds, xi_ds):>18.4f}"
            f"{np.interp(tau, s_f, xi_f):>13.4f}{np.interp(tau, s_ads, xi_ads):>12.4f}"
        )
    tc = first_conjugate_point(s_ads, xi_ads)
    print(f"\n  AdS (K<0) refocuses worldlines at tau = {tc:.5f} = pi L -- negative")
    print("  curvature acting like a box. de Sitter (K>0) spreads them apart like")
    print("  cosh/sinh -- that is accelerated expansion. Both signs are the reverse")
    print("  of the Riemannian case in experiment 1, purely because u.u = -1.")


if __name__ == "__main__":
    experiment_1()
    experiment_2()
    experiment_3()
