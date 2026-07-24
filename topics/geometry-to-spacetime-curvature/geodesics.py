"""Geodesics: straight lines in a curved space, and the classic tests of GR.

Covers Parts 9, 13 and 15 of the README.

The geodesic equation is the whole of "how things move" in general relativity:

    d^2 x^m / dl^2  +  Gamma^m_{n r} (dx^n/dl)(dx^r/dl)  =  0

A generic integrator (symbolic metric -> Christoffel symbols -> RK4) is applied to:

  1. Flat plane in POLAR coordinates -> geodesics are straight lines, even
     though Gamma != 0. (Curvy coordinates are not curvature.)
  2. Unit sphere -> geodesics are great circles (checked by the constancy of
     the plane normal in the embedding, which the integrator never sees).
  3. Schwarzschild -> a bound orbit. Checks the two Killing conservation laws
     (E and L) and the normalisation u.u = -1 to machine precision, then
     measures the perihelion precession directly.

Then the two classic solar-system tests, from the exact orbit equation
(the geodesic equation after eliminating t and tau with E and L):

    d^2u/dphi^2 + u = 3 M u^2          (light)
    d^2u/dphi^2 + u = M/L^2 + 3 M u^2  (planets)      with u = 1/r

  4. Light deflection at the solar limb  -> 1.75 arcsec  (Eddington 1919)
  5. Mercury's perihelion precession     -> 43 arcsec/century

Run:  uv run topics/geometry-to-spacetime-curvature/geodesics.py
"""

import numpy as np
import sympy as sp

# ---------------------------------------------------------------------------
# generic machinery: metric -> Christoffel -> geodesic integrator
# ---------------------------------------------------------------------------


def compile_christoffel(g, x):
    """Numeric callable for Gamma^l_{mn}(p) from a symbolic metric."""
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
    fn = sp.lambdify(x, sp.Array(Gamma), "numpy")
    return lambda p: np.asarray(fn(*p), dtype=float).reshape(n, n, n)


def geodesic(Gamma_of, x0, v0, lam_max, steps):
    """RK4-integrate the geodesic equation. Returns (lambdas, positions, velocities)."""

    def rhs(state):
        n = len(state) // 2
        x, v = state[:n], state[n:]
        # dv^m/dl = -Gamma^m_{n r} v^n v^r
        return np.concatenate([v, -np.einsum("mnr,n,r->m", Gamma_of(x), v, v)])

    h = lam_max / steps
    state = np.concatenate([np.asarray(x0, float), np.asarray(v0, float)])
    out = [state]
    for _ in range(steps):
        k1 = rhs(state)
        k2 = rhs(state + h / 2 * k1)
        k3 = rhs(state + h / 2 * k2)
        k4 = rhs(state + h * k3)
        state = state + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        out.append(state)
    out = np.array(out)
    n = len(x0)
    return np.linspace(0, lam_max, steps + 1), out[:, :n], out[:, n:]


# ---------------------------------------------------------------------------
# 1. flat plane in polar coordinates
# ---------------------------------------------------------------------------


def demo_polar_plane():
    print("=" * 72)
    print("1. FLAT PLANE IN POLAR COORDINATES -- geodesics are straight lines")
    print("=" * 72)

    r, th = sp.symbols("r theta", positive=True)
    Gamma_of = compile_christoffel(sp.Matrix([[1, 0], [0, r**2]]), [r, th])

    # start at (x, y) = (1, 0) moving at 60 degrees to the radius, unit speed
    r0, alpha = 1.0, np.radians(60.0)
    x0 = [r0, 0.0]
    v0 = [np.cos(alpha), np.sin(alpha) / r0]  # dr/ds, dtheta/ds

    _, pos, _ = geodesic(Gamma_of, x0, v0, lam_max=4.0, steps=20000)
    xy = np.stack(
        [pos[:, 0] * np.cos(pos[:, 1]), pos[:, 0] * np.sin(pos[:, 1])], axis=1
    )

    start, direction = xy[0], xy[1] - xy[0]
    direction = direction / np.linalg.norm(direction)
    offsets = xy - start
    # distance from the straight line through (start, direction)
    perp = np.abs(offsets[:, 0] * direction[1] - offsets[:, 1] * direction[0])

    print(f"   start (x,y) = ({xy[0, 0]:.4f}, {xy[0, 1]:.4f})")
    print(f"   end   (x,y) = ({xy[-1, 0]:.4f}, {xy[-1, 1]:.4f})")
    print(f"   max deviation from a straight line = {perp.max():.3e}")
    print("   -> exactly straight. Gamma != 0 only because the BASIS rotates.")


# ---------------------------------------------------------------------------
# 2. the sphere
# ---------------------------------------------------------------------------


def demo_sphere():
    print()
    print("=" * 72)
    print("2. UNIT SPHERE -- geodesics are great circles")
    print("=" * 72)

    th, ph = sp.symbols("theta phi")
    Gamma_of = compile_christoffel(sp.Matrix([[1, 0], [0, sp.sin(th) ** 2]]), [th, ph])

    beta = np.radians(30.0)  # inclination of the launch direction to the equator
    x0 = [np.pi / 2, 0.0]
    v0 = [np.sin(beta), np.cos(beta)]  # unit speed at theta = pi/2

    s, pos, vel = geodesic(Gamma_of, x0, v0, lam_max=2 * np.pi, steps=20000)
    t_, p_ = pos[:, 0], pos[:, 1]

    # embed in R^3 -- the integrator never used this
    X = np.stack([np.sin(t_) * np.cos(p_), np.sin(t_) * np.sin(p_), np.cos(t_)], axis=1)
    dX = np.gradient(X, s, axis=0)
    n = np.cross(X, dX)
    n /= np.linalg.norm(n, axis=1, keepdims=True)

    print(f"   launched at colatitude 90 deg, inclined {np.degrees(beta):.0f} deg")
    print(
        f"   colatitude range visited: "
        f"{np.degrees(t_.min()):.4f} .. {np.degrees(t_.max()):.4f} deg"
        f"   (great circle predicts {90 - np.degrees(beta):.0f} .. "
        f"{90 + np.degrees(beta):.0f})"
    )
    print(f"   wobble of the orbital-plane normal = {np.abs(n - n[0]).max():.3e}")
    print("   -> the path lies in a fixed plane through the centre: a GREAT CIRCLE.")
    print(
        f"   closure after arclength 2pi: |X(2pi) - X(0)| = "
        f"{np.linalg.norm(X[-1] - X[0]):.3e}"
    )
    print("   -> every geodesic on the unit sphere closes at length 2pi, and all")
    print("      of them intersect: positive curvature focuses geodesics.")


# ---------------------------------------------------------------------------
# 3. Schwarzschild: a bound orbit, conservation laws, precession
# ---------------------------------------------------------------------------


def demo_schwarzschild_orbit():
    print()
    print("=" * 72)
    print("3. SCHWARZSCHILD -- bound orbit, conserved E and L, precession")
    print("=" * 72)

    # equatorial slice (theta = pi/2), geometrized units G = c = M = 1
    t, r, ph = sp.symbols("t r phi", positive=True)
    f = 1 - 2 / r
    g = sp.diag(-f, 1 / f, r**2)
    Gamma_of = compile_christoffel(g, [t, r, ph])

    # orbit with semi-latus rectum p and eccentricity e (exact Schwarzschild ICs)
    p, e = 20.0, 0.3
    L = np.sqrt(p**2 / (p - 3 - e**2))
    E = np.sqrt((p - 2 - 2 * e) * (p - 2 + 2 * e) / (p * (p - 3 - e**2)))
    r0 = p / (1 + e)  # start at periapsis
    x0 = [0.0, r0, 0.0]
    v0 = [E / (1 - 2 / r0), 0.0, L / r0**2]  # dt/dtau, dr/dtau, dphi/dtau

    print(f"   p = {p}M, e = {e}  ->  r_peri = {r0:.4f}M, r_apo = {p / (1 - e):.4f}M")
    print(f"   E = {E:.10f}   L = {L:.10f} M")

    # ~4 radial periods (Keplerian estimate 2 pi (p/(1-e^2))^{3/2} ~ 650 M)
    tau, pos, vel = geodesic(Gamma_of, x0, v0, lam_max=2600.0, steps=260000)
    rr, pp = pos[:, 1], pos[:, 2]
    dt, dr, dp = vel[:, 0], vel[:, 1], vel[:, 2]

    # the two Killing conservation laws + the mass-shell normalisation
    E_num = (1 - 2 / rr) * dt
    L_num = rr**2 * dp
    norm = -(1 - 2 / rr) * dt**2 + dr**2 / (1 - 2 / rr) + rr**2 * dp**2

    print("\n   drift over the whole integration:")
    print(f"     E   (from d/dt symmetry)  : {np.ptp(E_num):.3e}")
    print(f"     L   (from d/dphi symmetry): {np.ptp(L_num):.3e}")
    print(f"     u.u (should be exactly -1): {np.abs(norm + 1).max():.3e}")
    print("   -> Killing vectors give exact constants of motion; the integrator")
    print("      reproduces them without ever being told about them.")

    # successive periapsis passages: minima of r, refined by parabola fit
    peri = []
    for i in range(1, len(rr) - 1):
        if rr[i] <= rr[i - 1] and rr[i] < rr[i + 1]:
            y0, y1, y2 = rr[i - 1], rr[i], rr[i + 1]
            shift = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)  # sub-step vertex
            phis = pp[i - 1 : i + 2]
            peri.append(np.interp(shift, [-1, 0, 1], phis))

    measured = float(np.mean(np.diff(peri) - 2 * np.pi))
    weak = 6 * np.pi / p  # 6 pi G M / (c^2 a (1 - e^2)), leading order in M/p
    exact = exact_precession(p, e)
    print(f"\n   periapsis passages found: {len(peri)}")
    print(
        f"   advance per orbit, RK4 geodesic : {measured:.9f} rad "
        f"= {np.degrees(measured):.4f} deg"
    )
    print(
        f"   advance per orbit, exact quadrature: {exact:.9f} rad "
        f"(check: {abs(measured - exact):.2e})"
    )
    print(
        f"   advance per orbit, weak field 6piM/p: {weak:.9f} rad "
        f"= {np.degrees(weak):.4f} deg"
    )
    print(
        f"   -> the leading-order formula is {abs(measured - weak) / measured:.0%} "
        f"low at p = 20M: this orbit is deep in the strong field,"
    )
    print("      precessing 70 deg per revolution. In the solar system p/M ~ 1e7")
    print("      and 6 pi M/p is exact for all practical purposes (test 5).")


def exact_precession(p, e):
    """Exact Schwarzschild periapsis advance, by quadrature (M = 1).

    (du/dphi)^2 = 2(u - u1)(u2 - u)(u3 - u) with u = 1/r; the roots are the
    apoapsis, the periapsis and a third root fixed by  u1 + u2 + u3 = 1/2.
    Substituting u = u1 + (u2 - u1) sin^2(psi) removes both endpoint
    singularities and leaves a smooth integrand for Simpson's rule.
    """
    u1, u2 = (1 - e) / p, (1 + e) / p
    u3 = 0.5 - u1 - u2
    n = 200001
    psi = np.linspace(0, np.pi / 2, n)
    f = 2.0 / np.sqrt(2 * (u3 - (u1 + (u2 - u1) * np.sin(psi) ** 2)))
    w = np.ones(n)
    w[1:-1:2], w[2:-1:2] = 4, 2
    return 2 * ((np.pi / 2) / (n - 1)) / 3 * np.sum(w * f) - 2 * np.pi


# ---------------------------------------------------------------------------
# 4-5. the classic tests, from the orbit equation in u = 1/r
# ---------------------------------------------------------------------------


def rk4_step(rhs, s, y, h):
    k1 = rhs(s, y)
    k2 = rhs(s + h / 2, y + h / 2 * k1)
    k3 = rhs(s + h / 2, y + h / 2 * k2)
    k4 = rhs(s + h, y + h * k3)
    return y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def advance(rhs, s0, y0, h, s_target):
    """Plain RK4 from s0 to (at least) s_target. Returns (s, y)."""
    s, y = s0, np.asarray(y0, float)
    while s < s_target:
        y = rk4_step(rhs, s, y, h)
        s += h
    return s, y


def integrate_to_root(rhs, s0, y0, event, h, s_max):
    """Integrate until `event(y)` changes sign; bisect the last step for the root."""
    s, y = s0, np.asarray(y0, float)
    sign0 = np.sign(event(y))
    while s < s_max:
        s_new, y_new = s + h, rk4_step(rhs, s, y, h)
        if np.sign(event(y_new)) != sign0:
            lo, hi = 0.0, h
            for _ in range(80):  # bisect with a single accurate RK4 step
                mid = 0.5 * (lo + hi)
                if np.sign(event(rk4_step(rhs, s, y, mid))) != sign0:
                    hi = mid
                else:
                    lo = mid
            return s + hi, rk4_step(rhs, s, y, hi)
        s, y = s_new, y_new
    raise RuntimeError("no root found")


ARCSEC = 180 * 3600 / np.pi  # radians -> arcseconds
GM_SUN = 1476.6250  # G M_sun / c^2, in metres
R_SUN = 6.957e8  # metres


def demo_light_deflection():
    print()
    print("=" * 72)
    print("4. LIGHT DEFLECTION BY THE SUN")
    print("=" * 72)

    # null orbit equation, non-dimensionalised with w = b/r :
    #     d^2 w / dphi^2 + w = 3 (M/b) w^2
    # start far away (w = 0) moving inward; the ray is straight if M = 0.
    b = R_SUN
    eps = GM_SUN / b

    def rhs(_phi, y):
        w, wp = y
        return np.array([wp, -w + 3 * eps * w**2])

    # w = 0 at the start too, so run past the point of closest approach first
    h = 1e-5
    phi_mid, y_mid = advance(rhs, 0.0, [0.0, 1.0], h, np.pi / 2)
    phi_end, _ = integrate_to_root(
        rhs, phi_mid, y_mid, event=lambda y: y[0], h=h, s_max=10.0
    )
    deflection = phi_end - np.pi

    print(f"   impact parameter b = R_sun = {b:.4e} m,  GM/c^2 = {GM_SUN} m")
    print(f"   M/b = {eps:.4e}   (the whole relativistic correction)")
    print(f"\n   integrated total swing  = {phi_end:.12f} rad  (pi + deflection)")
    print(
        f"   deflection (numeric)    = {deflection:.6e} rad = "
        f"{deflection * ARCSEC:.4f} arcsec"
    )
    print(
        f"   deflection (4GM/c^2 b)  = {4 * eps:.6e} rad = "
        f"{4 * eps * ARCSEC:.4f} arcsec"
    )
    print(
        f"   Newtonian 'photon as a particle' (2GM/c^2 b) = "
        f"{2 * eps * ARCSEC:.4f} arcsec"
    )
    print("\n   -> the factor of 2 is the point: half the bending comes from the")
    print("      curvature of TIME (the Newtonian potential), half from the")
    print("      curvature of SPACE. Eddington's 1919 eclipse measured 1.75\".")


def demo_mercury():
    print()
    print("=" * 72)
    print("5. MERCURY'S PERIHELION PRECESSION")
    print("=" * 72)

    a = 5.790905e10  # semi-major axis, metres
    e = 0.205630
    period_days = 87.9691
    p = a * (1 - e**2)  # semi-latus rectum
    eps = GM_SUN / p

    # timelike orbit equation with w = p/r :  w'' + w = 1 + 3 (M/p) w^2
    def rhs(_phi, y):
        w, wp = y
        return np.array([wp, 1 - w + 3 * eps * w**2])

    # start at perihelion (w = 1 + e, w' = 0), step past aphelion so that the
    # next w' = 0 event is the *next* perihelion
    h = 1e-5
    phi, y = advance(rhs, 0.0, [1 + e, 0.0], h, np.pi)
    phi_peri, _ = integrate_to_root(
        rhs, phi, y, event=lambda y: y[1], h=h, s_max=4 * np.pi
    )

    per_orbit = phi_peri - 2 * np.pi
    predicted = 6 * np.pi * eps
    orbits_per_century = 36525.0 / period_days

    print(f"   a = {a:.6e} m, e = {e},  p = a(1-e^2) = {p:.6e} m")
    print(f"   M/p = {eps:.6e}")
    print(f"\n   perihelion-to-perihelion angle = {phi_peri:.12f} rad")
    print(f"   2 pi                           = {2 * np.pi:.12f} rad")
    print(
        f"   advance per orbit (numeric)    = {per_orbit:.6e} rad = "
        f"{per_orbit * ARCSEC:.6f} arcsec"
    )
    print(
        f"   advance per orbit (6 pi M/p)   = {predicted:.6e} rad = "
        f"{predicted * ARCSEC:.6f} arcsec"
    )
    print(f"\n   orbits per century = {orbits_per_century:.2f}")
    print(f"   => {per_orbit * ARCSEC * orbits_per_century:.2f} arcsec/century")
    print('\n   -> the famous 43"/century. Newtonian perturbations from the other')
    print('      planets account for ~532"/century of the observed ~575"; the')
    print('      leftover 43" was unexplained from 1859 until Einstein, in 1915,')
    print("      computed it from pure geometry with no adjustable parameters.")


if __name__ == "__main__":
    demo_polar_plane()
    demo_sphere()
    demo_schwarzschild_orbit()
    demo_light_deflection()
    demo_mercury()
