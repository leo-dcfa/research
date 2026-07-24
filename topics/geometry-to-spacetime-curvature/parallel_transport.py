"""Parallel transport and holonomy: curvature you can measure from the inside.

Covers Parts 6, 9 and 10 of the README.

Carry a vector around a closed loop, keeping it "as parallel as possible":

    dV^m/dl  +  Gamma^m_{n r} (dx^n/dl) V^r  =  0

In flat space it comes back unchanged. In curved space it comes back ROTATED,
and the rotation angle equals the enclosed curvature:

    holonomy angle  =  int_S K dA        (Gauss-Bonnet / the Riemann tensor)

Four experiments, all with the same generic integrator:

  1. Flat plane in POLAR coordinates -- Gamma != 0 everywhere, holonomy = 0.
     Curvy coordinates are not curvature. (Part 2.)
  2. Unit sphere, loops around a circle of latitude -- rotation = enclosed
     solid angle 2 pi (1 - cos theta_0).
  3. Unit sphere, small coordinate loops -- angle -> area, with the sign of K.
  4. Poincare half-plane (K = -1) -- rotation goes the other way.

Only numpy + sympy: sympy builds the Christoffel symbols from the metric,
numpy integrates the transport ODE with RK4.

Run:  uv run topics/geometry-to-spacetime-curvature/parallel_transport.py
"""

import numpy as np
import sympy as sp

# ---------------------------------------------------------------------------
# metric (symbolic)  ->  Christoffel symbols (fast numeric callables)
# ---------------------------------------------------------------------------


def compile_geometry(g, x):
    """Return (g_of, Gamma_of): numeric callables for g_{mn} and Gamma^l_{mn}."""
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

    g_fn = sp.lambdify(x, g, "numpy")
    Gamma_fn = sp.lambdify(x, sp.Array(Gamma), "numpy")

    def g_of(p):
        return np.asarray(g_fn(*p), dtype=float).reshape(n, n)

    def Gamma_of(p):
        return np.asarray(Gamma_fn(*p), dtype=float).reshape(n, n, n)

    return g_of, Gamma_of


# ---------------------------------------------------------------------------
# transport ODE, RK4 along a piecewise-linear loop in coordinate space
# ---------------------------------------------------------------------------


def transport_leg(Gamma_of, start, end, V, steps=4000):
    """Parallel-transport V along the straight coordinate segment start -> end."""
    start, end = np.asarray(start, float), np.asarray(end, float)
    velocity = end - start  # dx/ds with s in [0, 1]
    h = 1.0 / steps

    def dV(s, V):
        p = start + s * velocity
        # dV^m/ds = -Gamma^m_{n r} (dx^n/ds) V^r
        return -np.einsum("mnr,n,r->m", Gamma_of(p), velocity, V)

    for i in range(steps):
        s = i * h
        k1 = dV(s, V)
        k2 = dV(s + h / 2, V + h / 2 * k1)
        k3 = dV(s + h / 2, V + h / 2 * k2)
        k4 = dV(s + h, V + h * k3)
        V = V + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    return V


def transport_loop(Gamma_of, vertices, V0, steps=4000):
    """Transport V0 around the closed loop through `vertices` (returns to the first)."""
    V = np.asarray(V0, float)
    for a, b in zip(vertices, vertices[1:] + [vertices[0]]):
        V = transport_leg(Gamma_of, a, b, V, steps)
    return V


# ---------------------------------------------------------------------------
# measuring the result with the metric
# ---------------------------------------------------------------------------


def orthonormal_frame(g):
    """Gram-Schmidt the coordinate basis into an orthonormal 2D frame."""
    e1 = np.array([1.0, 0.0])
    e1 = e1 / np.sqrt(e1 @ g @ e1)
    e2 = np.array([0.0, 1.0])
    e2 = e2 - (e2 @ g @ e1) * e1
    e2 = e2 / np.sqrt(e2 @ g @ e2)
    return e1, e2


def signed_angle(V, W, g):
    """Angle from V to W in the tangent plane at a point with metric g, in (-pi, pi]."""
    e1, e2 = orthonormal_frame(g)
    av, bv = V @ g @ e1, V @ g @ e2
    aw, bw = W @ g @ e1, W @ g @ e2
    return (np.arctan2(bw, aw) - np.arctan2(bv, av) + np.pi) % (2 * np.pi) - np.pi


def norm(V, g):
    return float(np.sqrt(abs(V @ g @ V)))


def wrap(angle):
    return (angle + np.pi) % (2 * np.pi) - np.pi


# ---------------------------------------------------------------------------
# the geometries
# ---------------------------------------------------------------------------


def flat_polar():
    r, th = sp.symbols("r theta", positive=True)
    g = sp.Matrix([[1, 0], [0, r**2]])  # ds^2 = dr^2 + r^2 dtheta^2
    return g, [r, th]


def unit_sphere():
    th, ph = sp.symbols("theta phi")
    g = sp.Matrix([[1, 0], [0, sp.sin(th) ** 2]])  # ds^2 = dth^2 + sin^2(th) dph^2
    return g, [th, ph]


def poincare_half_plane():
    x, y = sp.symbols("x y", positive=True)
    g = sp.Matrix([[1 / y**2, 0], [0, 1 / y**2]])  # ds^2 = (dx^2 + dy^2)/y^2, K = -1
    return g, [x, y]


# ---------------------------------------------------------------------------
# experiments
# ---------------------------------------------------------------------------


def experiment_flat_polar():
    print("=" * 72)
    print("1. FLAT PLANE IN POLAR COORDINATES -- Gamma != 0, curvature = 0")
    print("=" * 72)

    g, x = flat_polar()
    g_of, Gamma_of = compile_geometry(g, x)

    print("   ds^2 = dr^2 + r^2 dtheta^2")
    print(f"   Gamma^r_(theta theta) at r=2 : {Gamma_of([2.0, 0.0])[0, 1, 1]:+.4f}")
    print(f"   Gamma^theta_(r theta) at r=2 : {Gamma_of([2.0, 0.0])[1, 0, 1]:+.4f}")
    print("   -> the connection is nonzero: the basis vectors rotate and rescale.")

    base = [1.0, 0.0]
    loop = [base, [3.0, 0.0], [3.0, 2.0], [1.0, 2.0]]  # an annular sector
    V0 = np.array([1.0, 0.3])
    Vf = transport_loop(Gamma_of, loop, V0)
    g0 = g_of(base)

    print(f"\n   V initial = {V0}")
    print(f"   V final   = {Vf}")
    print(f"   holonomy angle = {np.degrees(signed_angle(V0, Vf, g0)):+.6f} deg")
    print("   -> ZERO. Nonzero Christoffels, but the space is flat: R = 0.")
    print("      Only the Riemann tensor (not Gamma) decides curvature.")


def experiment_sphere_latitude():
    print()
    print("=" * 72)
    print("2. UNIT SPHERE -- transport around a circle of latitude")
    print("=" * 72)
    print("   loop at colatitude th0, phi: 0 -> 2pi")
    print("   predicted rotation = enclosed solid angle = 2 pi (1 - cos th0)\n")

    g, x = unit_sphere()
    g_of, Gamma_of = compile_geometry(g, x)

    header = f"   {'th0 (deg)':>10} {'cap area':>12} {'mod 2pi':>10} {'measured':>12}"
    print(header + f" {'match':>7} {'|V| drift':>11}")
    print(f"   {'-' * 10} {'-' * 12} {'-' * 10} {'-' * 12} {'-' * 7} {'-' * 11}")
    for th0_deg in (15, 30, 45, 60, 90, 120):
        th0 = np.radians(th0_deg)
        base = [th0, 0.0]
        V0 = np.array([1.0, 0.0])  # points due south along the meridian
        # one leg, not a loop: phi is periodic, so phi = 2pi IS the start point
        Vf = transport_leg(Gamma_of, base, [th0, 2 * np.pi], V0, steps=8000)
        g0 = g_of(base)

        area = 2 * np.pi * (1 - np.cos(th0))  # = int K dA over the polar cap
        measured = signed_angle(V0, Vf, g0)
        drift = abs(norm(Vf, g0) - norm(V0, g0))
        # an angle is only defined mod 2pi, and +pi == -pi
        match = min(abs(wrap(measured - area)), abs(wrap(measured + area))) < 1e-6
        print(
            f"   {th0_deg:>10} {area:>12.6f} {wrap(area):>10.6f} {measured:>12.6f}"
            f" {str(match):>7} {drift:>11.1e}"
        )

    print(
        "\n   measured == area (mod 2 pi; +pi and -pi are the same rotation).\n"
        "   The equator (th0 = 90 deg) encloses 2 pi steradians = 0 mod 2 pi, so a\n"
        "   vector transported around the equator comes back UNCHANGED -- the\n"
        "   equator is a geodesic, and transport along a geodesic keeps the angle\n"
        "   to the tangent fixed. |V| drift ~ 1e-15: transport preserves length,\n"
        "   which is metric compatibility (grad g = 0) showing up numerically."
    )


def experiment_sphere_patch():
    print()
    print("=" * 72)
    print("3. UNIT SPHERE -- small loops: holonomy -> int K dA  (K = +1)")
    print("=" * 72)
    print("   coordinate rectangle th in [th1, th2], phi in [0, dphi]")
    print("   area = dphi (cos th1 - cos th2)\n")

    g, x = unit_sphere()
    g_of, Gamma_of = compile_geometry(g, x)

    print(f"   {'size':>8} {'int K dA':>14} {'holonomy':>14} {'ratio':>10}")
    print(f"   {'-' * 8} {'-' * 14} {'-' * 14} {'-' * 10}")
    for size in (0.8, 0.4, 0.2, 0.1, 0.05):
        th1, th2, dphi = 1.0, 1.0 + size, size
        base = [th1, 0.0]
        loop = [base, [th2, 0.0], [th2, dphi], [th1, dphi]]
        V0 = np.array([1.0, 0.0])
        Vf = transport_loop(Gamma_of, loop, V0, steps=2000)

        area = dphi * (np.cos(th1) - np.cos(th2))
        holo = signed_angle(V0, Vf, g_of(base))
        print(f"   {size:>8} {area:>14.8f} {holo:>14.8f} {holo / area:>10.5f}")

    print(
        "\n   ratio = 1 at every size: holonomy = int K dA is EXACT for any loop\n"
        "   bounding a disk, not just a small one. Shrink the loop and the same\n"
        "   statement becomes the definition of the Riemann tensor,\n"
        "   [grad_m, grad_n] V^r = R^r_(s m n) V^s -- curvature is holonomy per\n"
        "   unit area. In 2D there is one independent component and R = 2K."
    )


def experiment_hyperbolic():
    print()
    print("=" * 72)
    print("4. POINCARE HALF-PLANE -- K = -1, holonomy flips sign")
    print("=" * 72)
    print("   ds^2 = (dx^2 + dy^2) / y^2 ;  area element dx dy / y^2\n")

    g, x = poincare_half_plane()
    g_of, Gamma_of = compile_geometry(g, x)

    print(f"   {'size':>8} {'int K dA':>14} {'holonomy':>14} {'ratio':>10}")
    print(f"   {'-' * 8} {'-' * 14} {'-' * 14} {'-' * 10}")
    for size in (0.8, 0.4, 0.2, 0.1, 0.05):
        x1, y1 = 0.0, 1.0
        x2, y2 = size, 1.0 + size
        base = [x1, y1]
        loop = [base, [x2, y1], [x2, y2], [x1, y2]]
        V0 = np.array([1.0, 0.0])
        Vf = transport_loop(Gamma_of, loop, V0, steps=2000)

        area = (x2 - x1) * (1 / y1 - 1 / y2)  # hyperbolic area of the rectangle
        int_K = -area  # K = -1
        holo = signed_angle(V0, Vf, g_of(base))
        print(f"   {size:>8} {int_K:>14.8f} {holo:>14.8f} {holo / int_K:>10.5f}")

    print(
        "\n   Same law, opposite sign: negative curvature rotates the vector the\n"
        "   other way. Sphere / plane / hyperbolic plane are the three constant-\n"
        "   curvature geometries, and holonomy per unit area is exactly K in each."
    )


if __name__ == "__main__":
    experiment_flat_polar()
    experiment_sphere_latitude()
    experiment_sphere_patch()
    experiment_hyperbolic()
