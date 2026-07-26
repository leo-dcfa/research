"""Is our universe spherical or hyperbolic? Turning Omega_k into an observable.

Spatial curvature acts as a lens on the whole sky: the CMB sound horizon is a
standard ruler of known length r_s ~ 144 Mpc (comoving) sitting at z* ~ 1090,
so its observed angular size depends on whether light travelled through a
spherical, flat, or hyperbolic geometry to reach us.

    D_C = (c/H0) * int_0^z dz' / E(z')                       comoving distance
    D_M = (c/H0) / sqrt(Om_k) * sinh( sqrt(Om_k) H0 D_C/c )   Om_k > 0, open
        = D_C                                                 Om_k = 0, flat
        = (c/H0) / sqrt(|Om_k|) * sin( ... )                  Om_k < 0, closed
    theta_s = r_s / D_M      (comoving ruler with comoving distance; the two
                              factors of (1+z) in r_s^phys / D_A cancel)

Sign convention warning: Omega_k > 0 <=> k = -1 <=> HYPERBOLIC / open.

Run:  uv run topics/spacetime-curvature-hyperbolic-vs-spherical/cosmology_curvature.py
"""

import numpy as np

C_KMS = 299_792.458  # km/s
H0 = 67.4  # km/s/Mpc, so h = 0.674   (Planck 2018)
OM_M_H2 = 0.1430  # physical matter density omega_m = Omega_m h^2 -- this is what
#                   the CMB damping/peak-ratio structure actually pins down, and
#                   it also fixes r_s. Omega_m = omega_m / h^2 then follows.
R_S = 144.4  # Mpc, comoving sound horizon at the drag epoch
Z_STAR = 1089.9  # redshift of last scattering


def E(z, om_m, om_k, om_lam):
    """H(z)/H0 for matter + curvature + Lambda (radiation ignored: <0.1% here)."""
    return np.sqrt(om_m * (1 + z) ** 3 + om_k * (1 + z) ** 2 + om_lam)


def comoving_distance(z, om_m, om_k, h0, n=200_000):
    """D_C = (c/H0) int_0^z dz'/E(z'), Simpson's rule."""
    om_lam = 1.0 - om_m - om_k
    zs = np.linspace(0.0, z, n + 1)
    integrand = 1.0 / E(zs, om_m, om_k, om_lam)
    step = z / n
    total = (
        (
            integrand[0]
            + integrand[-1]
            + 4 * integrand[1:-1:2].sum()
            + 2 * integrand[2:-1:2].sum()
        )
        * step
        / 3
    )
    return (C_KMS / h0) * total


def transverse_distance(z, om_m, om_k, h0):
    """D_M -- the sin/sinh that turns a comoving distance into an angular one."""
    d_c = comoving_distance(z, om_m, om_k, h0)
    hubble_dist = C_KMS / h0
    if om_k > 1e-12:  # hyperbolic / open
        rt = np.sqrt(om_k)
        return hubble_dist / rt * np.sinh(rt * d_c / hubble_dist)
    if om_k < -1e-12:  # spherical / closed
        rt = np.sqrt(-om_k)
        return hubble_dist / rt * np.sin(rt * d_c / hubble_dist)
    return d_c


def acoustic_scale(om_k, h0=H0):
    """theta_s in degrees, the acoustic multipole l_A = pi/theta_s, and D_M.

    Omega_m is derived from the fixed physical density omega_m = Omega_m h^2,
    so varying h0 here moves along the direction the CMB leaves open.
    """
    om_m = OM_M_H2 / (h0 / 100.0) ** 2
    d_m = transverse_distance(Z_STAR, om_m, om_k, h0)
    theta = R_S / d_m  # r_s is comoving, so it pairs with the comoving D_M
    return np.degrees(theta), np.pi / theta, d_m


def h0_matching_flat(om_k, target_theta, lo=30.0, hi=200.0, tol=1e-8):
    """Bisect for the H0 that reproduces the flat model's theta_s at this Omega_k.

    This traces the CMB's 'geometric degeneracy' direction: the acoustic scale is
    one number, so curvature can always be traded against expansion rate.
    """

    def f(h):
        return acoustic_scale(om_k, h)[0] - target_theta

    if f(lo) * f(hi) > 0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def main():
    om_m_flat = OM_M_H2 / (H0 / 100) ** 2
    print("=" * 78)
    print("Curvature as a lens: the CMB acoustic scale vs Omega_k")
    print(
        f"(H0 = {H0}, omega_m = Omega_m h^2 = {OM_M_H2} -> Omega_m = {om_m_flat:.3f},"
        f" r_s = {R_S} Mpc, z* = {Z_STAR})"
    )
    print("=" * 78)
    print(
        f"\n{'Omega_k':>9}{'k':>4}{'geometry':>13}{'D_M(z*) /Mpc':>15}"
        f"{'theta_s /deg':>14}{'l_A':>8}"
    )
    flat_theta, flat_l, _ = acoustic_scale(0.0)
    for om_k in (-0.10, -0.05, -0.01, 0.0, 0.01, 0.05, 0.10):
        theta, ell, d_m = acoustic_scale(om_k)
        if om_k < 0:
            k, geom = "+1", "spherical"
        elif om_k > 0:
            k, geom = "-1", "hyperbolic"
        else:
            k, geom = " 0", "flat"
        print(f"{om_k:>9.2f}{k:>4}{geom:>13}{d_m:>15.0f}{theta:>14.4f}{ell:>8.1f}")

    print("\n  A closed (spherical) universe MAGNIFIES the CMB spots -- bigger theta,")
    print("  smaller l_A. An open (hyperbolic) one shrinks them.")
    print(f"  l_A = pi/theta_s is the acoustic *scale* ({flat_l:.0f}; Planck: 301.7),")
    print("  not the first peak -- driving effects shift that one down to l ~ 220.")

    # How well can we actually pin it down?
    print("\n" + "=" * 78)
    print("Sensitivity, and the geometric degeneracy")
    print("=" * 78)
    d_theta = (acoustic_scale(0.01)[0] - acoustic_scale(-0.01)[0]) / 0.02
    print(f"\n  flat prediction     : theta_s = {flat_theta:.4f} deg")
    print(f"  d(theta_s)/d(Om_k)  = {d_theta:.4f} deg per unit Omega_k")
    print(
        f"  Planck+BAO Omega_k = 0.0007 +- 0.0019  ->  theta_s shifts by only "
        f"{abs(d_theta) * 0.0019:.5f} deg"
    )

    # The catch: at fixed physical densities, H0 can imitate curvature exactly.
    print("\n  The catch: theta_s is ONE number, and H0 can absorb the curvature.")
    print("  Holding omega_m (and hence r_s) fixed, solve for the H0 that gives")
    print("  the flat model's theta_s exactly:")
    print(f"\n{'Omega_k':>9}{'H0 needed':>12}{'Omega_m':>10}{'theta_s /deg':>14}")
    for om_k in (-0.05, -0.02, 0.0, 0.02, 0.05):
        h = h0_matching_flat(om_k, flat_theta)
        if h is None:
            print(f"{om_k:>9.2f}{'none in 30-200':>16}")
            continue
        om_m = OM_M_H2 / (h / 100) ** 2
        print(f"{om_k:>9.2f}{h:>12.1f}{om_m:>10.3f}{acoustic_scale(om_k, h)[0]:>14.4f}")
    print("\n  Identical acoustic scale, wildly different geometry AND expansion rate:")
    print("  that is the geometric degeneracy, and it runs towards closed + low H0.")
    print("  (Omega_k = +0.05 has no partner at all: at fixed omega_m the biggest")
    print("  theta_s an open model can reach is 0.5735 deg, at H0 ~ 130 -- still")
    print("  short of the flat value. Open universes are the harder ones to hide.)")
    print("\n  It is broken by a second distance probe -- low-z BAO, lensing,")
    print("  supernovae -- which is why the tight Omega_k bound is always quoted as")
    print("  'CMB + BAO'. This is also exactly why the mild closed-universe pull in")
    print("  Planck's CMB-only likelihood drags H0 down to ~54: it is sliding along")
    print("  this line, straight away from every direct H0 measurement.")

    # Curvature radius implied by the bound.
    print("\n" + "=" * 78)
    print("How big is the curvature radius, if there is one?")
    print("=" * 78)
    hubble_dist = C_KMS / H0
    for om_k in (0.0019, 0.02, 0.10):
        r_curv = hubble_dist / np.sqrt(om_k)
        print(
            f"  |Omega_k| = {om_k:<7.4f} -> curvature radius = {r_curv / 1000:>6.1f} Gpc"
            f"   ({r_curv / hubble_dist:>5.1f} Hubble radii)"
        )
    print("\n  The 1-sigma bound already pushes the curvature radius past 100 Gpc --")
    print("  ~20x the radius of the observable universe. Space is flat to |Omega_k|")
    print("  < 0.002; inflation predicts exactly this (it drives Omega_k -> 0 like")
    print("  e^{-2N}, so 60 e-folds leaves nothing measurable behind).")


if __name__ == "__main__":
    main()
