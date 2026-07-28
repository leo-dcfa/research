"""Amplitudes and qubits: normalization, Born rule, phase, interference,
tensor products, entanglement and (partial) measurement.

Run with:  uv run python topics/amplitudes-and-qubits/amplitudes.py
"""

import numpy as np

KET0 = np.array([1, 0], dtype=complex)
KET1 = np.array([0, 1], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)


def probabilities(psi):
    """Born rule: Pr(x) = |alpha_x|^2."""
    return np.abs(psi) ** 2


def is_normalized(psi):
    return np.isclose(probabilities(psi).sum(), 1.0)


def bloch_angles(psi):
    """(theta, phi) on the Bloch sphere, with the global phase divided out."""
    alpha, beta = psi
    theta = 2 * np.arccos(np.clip(np.abs(alpha), 0.0, 1.0))
    # rotate so alpha is real and non-negative -> only the relative phase is left
    phi = np.angle(beta) - np.angle(alpha)
    return theta, np.mod(phi, 2 * np.pi)


def from_bloch_angles(theta, phi):
    return np.array(
        [np.cos(theta / 2), np.exp(1j * phi) * np.sin(theta / 2)], dtype=complex
    )


def measure_qubit(psi, qubit, n_qubits, rng):
    """Measure one qubit of an n-qubit state (qubit 0 = most significant).

    Returns (outcome, collapsed_state).
    """
    amps = psi.reshape((2,) * n_qubits)
    idx = [slice(None)] * n_qubits
    branches = []
    for b in (0, 1):
        idx[qubit] = b
        branch = np.zeros_like(amps)
        branch[tuple(idx)] = amps[tuple(idx)]  # projector P_b applied
        branches.append(branch.reshape(-1))
    weights = np.array([np.vdot(b, b).real for b in branches])
    outcome = rng.choice(2, p=weights / weights.sum())
    return outcome, branches[outcome] / np.sqrt(weights[outcome])


def is_entangled(psi2):
    """Two-qubit product-state test: a00*a11 - a01*a10 == 0 iff separable."""
    a00, a01, a10, a11 = psi2
    return not np.isclose(a00 * a11 - a01 * a10, 0)


def main():
    rng = np.random.default_rng(0)

    # --- a single qubit ---------------------------------------------------
    plus = (KET0 + KET1) / np.sqrt(2)
    minus = (KET0 - KET1) / np.sqrt(2)
    ket_i = (KET0 + 1j * KET1) / np.sqrt(2)

    print("state      amplitudes                 Pr(0), Pr(1)")
    for name, psi in [("|0>", KET0), ("|+>", plus), ("|->", minus), ("|i>", ket_i)]:
        assert is_normalized(psi)
        print(f"{name:<10} {np.round(psi, 3)}   {np.round(probabilities(psi), 3)}")
    print("-> |+>, |->, |i> are indistinguishable in one measurement\n")

    # --- global vs relative phase ----------------------------------------
    shifted = np.exp(1j * 0.7) * plus
    assert np.allclose(probabilities(shifted), probabilities(plus))
    assert np.allclose(bloch_angles(shifted), bloch_angles(plus))
    print(f"global phase is invisible: |+> and e^(0.7i)|+> -> {bloch_angles(plus)}")
    print(f"relative phase is not:     |-> -> {bloch_angles(minus)}")
    print(f"                           |i> -> {bloch_angles(ket_i)}\n")

    # Bloch round-trip
    theta, phi = bloch_angles(ket_i)
    assert np.allclose(from_bloch_angles(theta, phi), ket_i)

    # --- interference: H twice is the identity ---------------------------
    once = H @ KET0
    twice = H @ once
    print(f"|0> --H--> {np.round(once, 3)}  (Pr = {np.round(probabilities(once), 3)})")
    print(f"    --H--> {np.round(twice, 3)}  <- |1> amplitudes +1/2 and -1/2 cancel")
    assert np.allclose(twice, KET0)

    # the same two "paths" with probabilities instead of amplitudes: no cancelling
    stochastic = np.array([[0.5, 0.5], [0.5, 0.5]])
    print(f"classical coin flipped twice: {stochastic @ (stochastic @ [1.0, 0.0])}\n")

    # --- two qubits -------------------------------------------------------
    phi_plus = (np.kron(KET0, KET0) + np.kron(KET1, KET1)) / np.sqrt(2)
    product = np.kron(plus, KET0)
    print(f"|Phi+>   = {np.round(phi_plus, 3)}  entangled={is_entangled(phi_plus)}")
    print(f"|+>x|0>  = {np.round(product, 3)}  entangled={is_entangled(product)}\n")

    # --- measurement ------------------------------------------------------
    counts = np.zeros(4, dtype=int)
    for _ in range(1000):
        outcomes = rng.choice(4, p=probabilities(phi_plus))
        counts[outcomes] += 1
    labels = ["|00>", "|01>", "|10>", "|11>"]
    print("sampling |Phi+> 1000x:", dict(zip(labels, counts)))

    outcome, collapsed = measure_qubit(phi_plus, qubit=0, n_qubits=2, rng=rng)
    print(f"measured qubit 0 -> {outcome}; state collapses to {np.round(collapsed, 3)}")
    print("-> the second qubit is now forced to match: that's entanglement")


if __name__ == "__main__":
    main()
