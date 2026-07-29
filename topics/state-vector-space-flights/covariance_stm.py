"""State transition matrix and covariance propagation for a 6D state vector.

Shows the two things every navigation filter needs:
  1. Phi(t, t0) = d x(t) / d x(t0), obtained by integrating the variational
     equation alongside the trajectory, and validated against finite differences.
  2. P(t) = Phi P0 Phi^T + Q, and why the error ellipsoid grows along-track.

Run with:  uv run topics/state-vector-space-flights/covariance_stm.py
"""

import numpy as np

from state_vectors import MU_EARTH, R_EARTH, coe_to_rv, derivative, propagate

DEG = np.pi / 180.0


def gravity_gradient(r, mu=MU_EARTH):
    """d a / d r for two-body gravity — the tidal tensor."""
    r_mag = np.linalg.norm(r)
    r_hat = r / r_mag
    return -mu / r_mag**3 * (np.eye(3) - 3 * np.outer(r_hat, r_hat))


def jacobian(state, mu=MU_EARTH):
    """A = d f / d x for x = [r, v], two-body dynamics."""
    a = np.zeros((6, 6))
    a[:3, 3:] = np.eye(3)
    a[3:, :3] = gravity_gradient(state[:3], mu)
    return a


def augmented_derivative(y, mu=MU_EARTH):
    """Derivative of the 6 + 36 augmented vector [x, vec(Phi)]."""
    state = y[:6]
    phi = y[6:].reshape(6, 6)
    return np.concatenate([derivative(state, mu), (jacobian(state, mu) @ phi).ravel()])


def propagate_with_stm(state, t_end, dt=1.0, mu=MU_EARTH):
    """RK4 on the trajectory and the variational equation together."""
    y = np.concatenate([np.asarray(state, dtype=float), np.eye(6).ravel()])
    t = 0.0
    while t < t_end:
        h = min(dt, t_end - t)
        k1 = augmented_derivative(y, mu)
        k2 = augmented_derivative(y + h / 2 * k1, mu)
        k3 = augmented_derivative(y + h / 2 * k2, mu)
        k4 = augmented_derivative(y + h * k3, mu)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += h
    return y[:6], y[6:].reshape(6, 6)


def stm_finite_difference(state, t_end, dt=1.0, eps=1e-4):
    """Brute-force Phi: perturb each initial component and re-propagate."""
    phi = np.zeros((6, 6))
    for j in range(6):
        step = np.zeros(6)
        step[j] = eps
        plus = propagate(state + step, t_end, dt)
        minus = propagate(state - step, t_end, dt)
        phi[:, j] = (plus - minus) / (2 * eps)
    return phi


def ric_rotation(state):
    """Rotation from inertial to the RIC (radial / in-track / cross-track) frame."""
    r, v = state[:3], state[3:]
    radial = r / np.linalg.norm(r)
    cross = np.cross(r, v)
    cross = cross / np.linalg.norm(cross)
    in_track = np.cross(cross, radial)
    return np.vstack([radial, in_track, cross])


def main():
    a0 = R_EARTH + 700.0
    r0, v0 = coe_to_rv(
        a=a0, e=0.001, i=98.2 * DEG, raan=40 * DEG, argp=30 * DEG, nu=0.0
    )
    state0 = np.concatenate([r0, v0])
    period = 2 * np.pi * np.sqrt(a0**3 / MU_EARTH)

    # --- the STM, two ways -------------------------------------------------
    quarter = period / 4
    _, phi = propagate_with_stm(state0, quarter, dt=1.0)
    phi_fd = stm_finite_difference(state0, quarter, dt=1.0)
    print(f"STM over a quarter orbit ({quarter / 60:.2f} min)")
    print(
        f"  max |variational - finite difference| = {np.max(np.abs(phi - phi_fd)):.3e}"
    )
    print(f"  det(Phi) = {np.linalg.det(phi):.9f}   (symplectic flow -> exactly 1)")

    # symplecticity: Phi^T J Phi = J for a Hamiltonian system
    j_mat = np.block([[np.zeros((3, 3)), np.eye(3)], [-np.eye(3), np.zeros((3, 3))]])
    print(
        f"  max |Phi^T J Phi - J|  = {np.max(np.abs(phi.T @ j_mat @ phi - j_mat)):.3e}"
    )

    # --- covariance propagation -------------------------------------------
    # 100 m position, 0.1 m/s velocity, uncorrelated at epoch
    p0 = np.diag([0.1**2] * 3 + [1e-4**2] * 3)
    print("\ncovariance growth (1-sigma, RIC frame)")
    print(
        f"  {'orbits':>7}  {'radial [km]':>12}  {'in-track [km]':>14}  {'cross [km]':>11}"
    )
    for n_orbits in (0, 1, 5, 20, 50):
        t = n_orbits * period
        if t == 0:
            state_t, phi_t = state0, np.eye(6)
        else:
            state_t, phi_t = propagate_with_stm(state0, t, dt=5.0)
        p_t = phi_t @ p0 @ phi_t.T  # no process noise: pure dynamics
        rot = ric_rotation(state_t)
        p_ric = rot @ p_t[:3, :3] @ rot.T
        sig = np.sqrt(np.diag(p_ric))
        print(f"  {n_orbits:>7}  {sig[0]:>12.4f}  {sig[1]:>14.4f}  {sig[2]:>11.4f}")

    print("\n  -> in-track dominates: an error in a becomes an error in the period,")
    print("     which integrates into a secular along-track drift.")

    # --- what a filter measurement does to that covariance -----------------
    # one perfectly-modelled range measurement from a ground station
    t = 3 * period
    state_t, phi_t = propagate_with_stm(state0, t, dt=5.0)
    p_minus = phi_t @ p0 @ phi_t.T
    station = np.array([-2516.7, 4653.5, 3551.9])  # ECEF-ish, treated as inertial here
    los = state_t[:3] - station
    h_row = np.concatenate([los / np.linalg.norm(los), np.zeros(3)])  # d(range)/dx
    sigma_rho = 0.005  # 5 m ranging noise
    gain = p_minus @ h_row / (h_row @ p_minus @ h_row + sigma_rho**2)
    p_plus = (np.eye(6) - np.outer(gain, h_row)) @ p_minus
    trace_before = np.sqrt(np.trace(p_minus[:3, :3]))
    trace_after = np.sqrt(np.trace(p_plus[:3, :3]))
    print(f"\none 5 m range update after {t / period:.0f} orbits")
    print(f"  position RSS sigma: {trace_before:.4f} km -> {trace_after:.4f} km")
    print(
        "  (one scalar measurement only squeezes the line-of-sight direction. It bites"
    )
    print(
        "   hard here because the propagated error is nearly rank-1 along-track and the"
    )
    print(
        "   station sees a chunk of it — full 6D observability still needs a tracking arc.)"
    )


if __name__ == "__main__":
    main()
