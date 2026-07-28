"""Hermitian adjoint (conjugate transpose): worked examples."""

import numpy as np


def dagger(a: np.ndarray) -> np.ndarray:
    """Hermitian adjoint A† = (A*)^T = (A^T)*."""
    return a.conj().T


def main() -> None:
    # --- complex conjugate reminder -------------------------------------
    for z in (1 + 1j, 2 + 3j, 1j, 5 + 0j):
        print(f"({z})* = {np.conj(z)}   |z|^2 = z z* = {(z * np.conj(z)).real}")

    # --- vector adjoint: column (ket) -> row (bra) -----------------------
    v = np.array([[1], [1j], [2 + 3j]])
    print("\nv =\n", v)
    print("v† =", dagger(v))

    # --- matrix adjoint --------------------------------------------------
    a = np.array([[1, 2 - 1j], [3j, 4]])
    print("\nA =\n", a)
    print("A† =\n", dagger(a))
    print("(A†)† == A:", np.allclose(dagger(dagger(a)), a))

    # --- bras of the computational basis ---------------------------------
    ket0 = np.array([[1], [0]])
    ket1 = np.array([[0], [1]])
    print("\n<0| =", dagger(ket0), "  <1| =", dagger(ket1))

    # a ket with a complex amplitude: the i flips sign in the bra
    psi = np.array([[1], [1j]]) / np.sqrt(2)  # |psi> = (|0> + i|1>)/sqrt(2)
    print("|psi> =", psi.ravel(), " ->  <psi| =", dagger(psi).ravel())

    # --- outer products: |v><w| builds operators --------------------------
    print("\n|0><0| =\n", ket0 @ dagger(ket0))  # projector onto |0>
    x = ket0 @ dagger(ket1) + ket1 @ dagger(ket0)  # X = |0><1| + |1><0|
    print("X = |0><1| + |1><0| =\n", x)
    print(
        "completeness sum_k |k><k| == I:",
        np.allclose(ket0 @ dagger(ket0) + ket1 @ dagger(ket1), np.eye(2)),
    )
    print(
        "(|0><1|)† == |1><0| :",
        np.allclose(dagger(ket0 @ dagger(ket1)), ket1 @ dagger(ket0)),
    )

    # --- Born rule and expectation values ---------------------------------
    amp0 = (dagger(ket0) @ psi).item()  # <0|psi>
    print("\n<0|psi> =", amp0, " -> Pr(0) = |<0|psi>|^2 =", abs(amp0) ** 2)
    print(
        "<psi|X|psi> =", (dagger(psi) @ x @ psi).item().real, "(real: X is Hermitian)"
    )

    # --- properties ------------------------------------------------------
    b = np.array([[0, 1 + 2j], [1, -1j]])
    alpha = 2 + 5j
    print(
        "\n(alpha A)† == alpha* A†:",
        np.allclose(dagger(alpha * a), np.conj(alpha) * dagger(a)),
    )
    print("(AB)†     == B† A†     :", np.allclose(dagger(a @ b), dagger(b) @ dagger(a)))

    # Hadamard: Hermitian *and* unitary
    h = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    print("H Hermitian (H† == H) :", np.allclose(dagger(h), h))
    print("H unitary (H† H == I) :", np.allclose(dagger(h) @ h, np.eye(2)))

    # --- inner product <v|w> = v† w --------------------------------------
    w = np.array([[0], [1], [1j]])
    print("\n<v|w> =", (dagger(v) @ w).item())
    print("<v|v> =", (dagger(v) @ v).item().real, "(real, >= 0)")
    print(
        "<w|v> == <v|w>* :",
        np.allclose((dagger(w) @ v).item(), np.conj((dagger(v) @ w).item())),
    )


if __name__ == "__main__":
    main()
