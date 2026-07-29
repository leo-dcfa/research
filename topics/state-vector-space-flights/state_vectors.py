"""State vectors in space flight: conversions, invariants, and propagation.

Run with:  uv run topics/state-vector-space-flights/state_vectors.py
"""

import numpy as np

MU_EARTH = 398_600.4418  # km^3 / s^2
R_EARTH = 6378.137  # km
J2 = 1.08262668e-3
DEG = np.pi / 180.0


# --------------------------------------------------------------------------
# state vector  ->  classical orbital elements
# --------------------------------------------------------------------------
def rv_to_coe(r, v, mu=MU_EARTH):
    """(r, v) in km, km/s  ->  dict of a, e, i, raan, argp, nu (angles in rad)."""
    r = np.asarray(r, dtype=float)
    v = np.asarray(v, dtype=float)
    r_mag = np.linalg.norm(r)
    v_mag = np.linalg.norm(v)

    h = np.cross(r, v)  # specific angular momentum
    h_mag = np.linalg.norm(h)
    n = np.cross([0.0, 0.0, 1.0], h)  # node vector
    n_mag = np.linalg.norm(n)

    e_vec = ((v_mag**2 - mu / r_mag) * r - np.dot(r, v) * v) / mu
    e = np.linalg.norm(e_vec)

    energy = v_mag**2 / 2 - mu / r_mag
    a = -mu / (2 * energy)

    i = np.arccos(h[2] / h_mag)

    # RAAN: undefined for equatorial orbits, take 0 by convention
    if n_mag > 1e-12:
        raan = np.arctan2(n[1], n[0]) % (2 * np.pi)
    else:
        raan = 0.0

    # argument of periapsis: undefined for circular orbits
    if n_mag > 1e-12 and e > 1e-12:
        argp = np.arccos(np.clip(np.dot(n, e_vec) / (n_mag * e), -1, 1))
        if e_vec[2] < 0:
            argp = 2 * np.pi - argp
    else:
        argp = 0.0

    if e > 1e-12:
        nu = np.arccos(np.clip(np.dot(e_vec, r) / (e * r_mag), -1, 1))
        if np.dot(r, v) < 0:  # descending -> second half of the orbit
            nu = 2 * np.pi - nu
    else:  # circular: fall back on the argument of latitude
        nu = np.arctan2(np.dot(r, np.cross(h, n)) / h_mag, np.dot(r, n)) % (2 * np.pi)

    return {"a": a, "e": e, "i": i, "raan": raan, "argp": argp, "nu": nu}


# --------------------------------------------------------------------------
# classical orbital elements  ->  state vector
# --------------------------------------------------------------------------
def coe_to_rv(a, e, i, raan, argp, nu, mu=MU_EARTH):
    """Angles in radians. Returns (r, v) in the same inertial frame."""
    p = a * (1 - e**2)
    r_mag = p / (1 + e * np.cos(nu))

    # perifocal (PQW) frame: x towards periapsis, z along angular momentum
    r_pqw = r_mag * np.array([np.cos(nu), np.sin(nu), 0.0])
    v_pqw = np.sqrt(mu / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])

    # PQW -> ECI is R3(-raan) R1(-i) R3(-argp)
    rot = _r3(-raan) @ _r1(-i) @ _r3(-argp)
    return rot @ r_pqw, rot @ v_pqw


def _r1(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])


def _r3(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])


# --------------------------------------------------------------------------
# invariants — cheap sanity checks on any state
# --------------------------------------------------------------------------
def invariants(r, v, mu=MU_EARTH):
    r = np.asarray(r, dtype=float)
    v = np.asarray(v, dtype=float)
    r_mag = np.linalg.norm(r)
    energy = np.dot(v, v) / 2 - mu / r_mag
    a = -mu / (2 * energy)
    return {
        "energy": energy,  # km^2/s^2, conserved in two-body
        "h": np.linalg.norm(np.cross(r, v)),  # km^2/s, conserved in two-body
        "a": a,
        "period": 2 * np.pi * np.sqrt(a**3 / mu) if a > 0 else np.inf,
    }


# --------------------------------------------------------------------------
# dynamics and propagation
# --------------------------------------------------------------------------
def accel(r, mu=MU_EARTH, with_j2=False):
    """Two-body acceleration, optionally plus the J2 oblateness term."""
    r_mag = np.linalg.norm(r)
    a = -mu * r / r_mag**3
    if with_j2:
        x, y, z = r
        k = 1.5 * J2 * mu * R_EARTH**2 / r_mag**5
        zr2 = 5 * z**2 / r_mag**2
        a = a - k * np.array([x * (1 - zr2), y * (1 - zr2), z * (3 - zr2)])
    return a


def derivative(state, mu=MU_EARTH, with_j2=False):
    """f(x) for the 6-dimensional state x = [r, v]."""
    return np.concatenate([state[3:], accel(state[:3], mu, with_j2)])


def propagate(state, t_end, dt=10.0, mu=MU_EARTH, with_j2=False):
    """Fixed-step RK4. Fine for a demo; real work uses RKF7(8) or Gauss-Jackson."""
    x = np.asarray(state, dtype=float).copy()
    t = 0.0
    while t < t_end:
        h = min(dt, t_end - t)  # exact landing on t_end
        k1 = derivative(x, mu, with_j2)
        k2 = derivative(x + h / 2 * k1, mu, with_j2)
        k3 = derivative(x + h / 2 * k2, mu, with_j2)
        k4 = derivative(x + h * k3, mu, with_j2)
        x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += h
    return x


def j2_secular_rates(a, e, i, mu=MU_EARTH):
    """Closed-form secular drift of RAAN and argument of periapsis (rad/s)."""
    p = a * (1 - e**2)
    n = np.sqrt(mu / a**3)
    factor = J2 * (R_EARTH / p) ** 2 * n
    raan_dot = -1.5 * factor * np.cos(i)
    argp_dot = 0.75 * factor * (5 * np.cos(i) ** 2 - 1)
    return raan_dot, argp_dot


def main():
    # A typical sun-synchronous LEO: 700 km altitude, i = 98.2 deg
    coe0 = {
        "a": R_EARTH + 700.0,
        "e": 0.001,
        "i": 98.2 * DEG,
        "raan": 40.0 * DEG,
        "argp": 30.0 * DEG,
        "nu": 0.0,
    }
    r0, v0 = coe_to_rv(**coe0)
    print("state vector at epoch (ECI, J2000)")
    print(f"  r = {np.array2string(r0, precision=3)} km")
    print(f"  v = {np.array2string(v0, precision=6)} km/s")

    inv = invariants(r0, v0)
    print(f"  |h| = {inv['h']:.3f} km^2/s   energy = {inv['energy']:.6f} km^2/s^2")
    print(f"  a = {inv['a']:.3f} km   period = {inv['period'] / 60:.3f} min")

    # round trip: state -> elements -> state
    back = rv_to_coe(r0, v0)
    print("\nround trip state -> elements -> state")
    for key, value in coe0.items():
        got = back[key]
        unit = "km" if key == "a" else ("" if key == "e" else "deg")
        scale = 1.0 if key in ("a", "e") else 1 / DEG
        print(f"  {key:>4}: in {value * scale:11.6f} {unit:3}  out {got * scale:11.6f}")
    r1, v1 = coe_to_rv(**back)
    print(f"  max position error: {np.max(np.abs(r1 - r0)):.3e} km")

    # two-body: invariants must be conserved over one orbit
    state0 = np.concatenate([r0, v0])
    one_orbit = inv["period"]
    state_kep = propagate(state0, one_orbit, dt=1.0)
    inv_kep = invariants(state_kep[:3], state_kep[3:])
    print("\ntwo-body propagation over one period (RK4, dt = 1 s)")
    print(f"  closure error: {np.linalg.norm(state_kep[:3] - r0):.3e} km")
    print(f"  d(energy)    : {inv_kep['energy'] - inv['energy']:+.3e} km^2/s^2")
    print(f"  d(|h|)       : {inv_kep['h'] - inv['h']:+.3e} km^2/s")

    # with J2: the plane precesses — compare integrated vs closed-form rate
    day = 86_400.0
    state_j2 = propagate(state0, day, dt=5.0, with_j2=True)
    coe_j2 = rv_to_coe(state_j2[:3], state_j2[3:])
    d_raan = (coe_j2["raan"] - coe0["raan"] + np.pi) % (2 * np.pi) - np.pi
    raan_dot, argp_dot = j2_secular_rates(coe0["a"], coe0["e"], coe0["i"])
    print("\nJ2 nodal regression over 1 day")
    print(f"  integrated : {d_raan / DEG:+.4f} deg/day")
    print(f"  closed form: {raan_dot * day / DEG:+.4f} deg/day")
    print(f"  sun-sync target: +{360 / 365.2422:.4f} deg/day")
    print(f"  apsidal drift  : {argp_dot * day / DEG:+.4f} deg/day")

    # an impulsive burn only touches the velocity block
    dv = 0.05 * v0 / np.linalg.norm(v0)  # 50 m/s prograde
    state_burn = state0 + np.concatenate([np.zeros(3), dv])
    coe_burn = rv_to_coe(state_burn[:3], state_burn[3:])
    print("\n50 m/s prograde burn at periapsis")
    print(
        f"  a: {coe0['a']:.3f} -> {coe_burn['a']:.3f} km  (+{coe_burn['a'] - coe0['a']:.3f})"
    )
    print(f"  e: {coe0['e']:.6f} -> {coe_burn['e']:.6f}")
    apo = coe_burn["a"] * (1 + coe_burn["e"]) - R_EARTH
    print(f"  apoapsis altitude: {apo:.3f} km")


if __name__ == "__main__":
    main()
