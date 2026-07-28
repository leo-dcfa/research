# Hermitian Adjoint (†)

The complex-vector-space replacement for "transpose". Written $A^\dagger$ ("A dagger"), and it's just **transpose + complex conjugate**. In quantum computing it's how you turn a ket into a bra.

## Reminder: the complex conjugate

For $z = a + bi$ with $a, b \in \mathbb{R}$, the conjugate is

$$z^* = \overline{z} = a - bi$$

i.e. **flip the sign of the imaginary part**. Geometrically: reflect across the real axis.

- $(1 + i)^* = 1 - i$, $(2 + 3i)^* = 2 - 3i$, $\;i^* = -i$, $\;5^* = 5$ (reals are their own conjugate).
- $z z^* = a^2 + b^2 = \lvert z \rvert^2$ — always real and $\geq 0$. This is why it shows up in inner products/probabilities.
- In polar form $z = re^{i\theta} \Rightarrow z^* = re^{-i\theta}$ (negate the phase).
- Rules: $(z + w)^* = z^* + w^*$, $(zw)^* = z^* w^*$, $(z^*)^* = z$.

## Definition

For a matrix or vector $A$, the Hermitian adjoint (a.k.a. **conjugate transpose**, adjoint) is

$$A^\dagger \equiv (A^*)^T = (A^T)^*$$

Conjugate every entry, then transpose — order doesn't matter. Entry-wise: $(A^\dagger)_{ij} = (A_{ji})^*$.

For a column vector $v \in \mathbb{C}^n$:

$$v = \begin{bmatrix} v_1 \\ v_2 \\ \vdots \\ v_n \end{bmatrix}
\quad\Longrightarrow\quad
v^\dagger = \begin{bmatrix} v_1^* & v_2^* & \cdots & v_n^* \end{bmatrix}$$

A column becomes a row, and every entry gets conjugated.

## Worked examples

**Vector** (Exercise 6.3):

$$v = \begin{bmatrix} 1 \\ i \\ 2 + 3i \end{bmatrix}
\quad\Longrightarrow\quad
v^\dagger = \begin{bmatrix} 1 & -i & 2 - 3i \end{bmatrix}$$

**Row → column**:

$$\left( \begin{bmatrix} 1+i \\ 1-i \end{bmatrix}^{\,T} \right)^{*} = \begin{bmatrix} 1+i & 1-i \end{bmatrix}^{*} = \begin{bmatrix} 1-i & 1+i \end{bmatrix}$$

**Matrix**:

$$A = \begin{bmatrix} 1 & 2-i \\ 3i & 4 \end{bmatrix}
\quad\Longrightarrow\quad
A^\dagger = \begin{bmatrix} 1 & -3i \\ 2+i & 4 \end{bmatrix}$$

## Bra–ket connection

A **ket** $\lvert v \rangle$ is a column vector; the **bra** $\langle v \rvert$ is its adjoint:

$$\langle v \rvert \equiv \lvert v \rangle^\dagger
\qquad\qquad
\lvert v \rangle \equiv \langle v \rvert^\dagger$$

(Definition 6.7. The relationship goes both ways since $(A^\dagger)^\dagger = A$.)

So for the computational basis:

$$\langle 0 \rvert = \lvert 0 \rangle^\dagger = \begin{bmatrix} 1 \\ 0 \end{bmatrix}^\dagger = \begin{bmatrix} 1 & 0 \end{bmatrix}
\qquad
\langle 1 \rvert = \lvert 1 \rangle^\dagger = \begin{bmatrix} 0 \\ 1 \end{bmatrix}^\dagger = \begin{bmatrix} 0 & 1 \end{bmatrix}$$

Both are real, so here the dagger is *just* a transpose. But for e.g. $\lvert v\rangle = \tfrac{1}{\sqrt2}(\lvert 0\rangle + i\lvert 1\rangle)$ you get $\langle v \rvert = \tfrac{1}{\sqrt2}\begin{bmatrix} 1 & -i \end{bmatrix}$ — the $i$ flips sign.

This is what makes the inner product $\langle v \rvert w \rangle = v^\dagger w = \sum_k v_k^* w_k$ behave: $\langle v \rvert v \rangle = \sum_k \lvert v_k \rvert^2 \geq 0$ is a real norm, so amplitudes square to probabilities.

## Properties worth memorising

| Property | Statement |
|---|---|
| Involution | $(A^\dagger)^\dagger = A$ |
| Anti-linearity | $(\alpha A)^\dagger = \alpha^* A^\dagger$ — the scalar conjugates! |
| Additivity | $(A + B)^\dagger = A^\dagger + B^\dagger$ |
| Reverses products | $(AB)^\dagger = B^\dagger A^\dagger$ (so $(A\lvert v\rangle)^\dagger = \langle v \rvert A^\dagger$) |
| Inverse | $(A^{-1})^\dagger = (A^\dagger)^{-1}$ |
| Real matrices | $A^\dagger = A^T$ |

Two classes defined by the dagger:

- **Hermitian / self-adjoint**: $A^\dagger = A$ → real eigenvalues. These are the **observables**.
- **Unitary**: $U^\dagger U = U U^\dagger = I$, so $U^{-1} = U^\dagger$ → preserves norms. These are the **quantum gates**.

## In Python

```python
import numpy as np

v = np.array([[1], [1j], [2 + 3j]])          # column vector (ket)
v_dag = v.conj().T                            # bra  ->  [[1.-0.j, 0.-1.j, 2.-3.j]]

A = np.array([[1, 2 - 1j], [3j, 4]])
A_dag = A.conj().T                            # == A.T.conj(); order doesn't matter

# Hadamard is both Hermitian and unitary
H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
np.allclose(H.conj().T, H)                    # True  -> Hermitian
np.allclose(H.conj().T @ H, np.eye(2))        # True  -> unitary

# Inner product <v|w> = v† w
w = np.array([[0], [1], [1j]])
braket = (v.conj().T @ w).item()               # complex scalar
norm_sq = (v.conj().T @ v).item().real         # always real >= 0
```

See [`adjoint.py`](adjoint.py) for a runnable version with the examples above.

## Gotchas

- `.T` alone is **not** the adjoint for complex arrays — a silent, very common bug. Always `.conj().T`.
- Anti-linearity: pulling a scalar out of a bra conjugates it. $\langle \alpha v \rvert = \alpha^* \langle v \rvert$, while $\lvert \alpha v\rangle = \alpha \lvert v \rangle$.
- $(AB)^\dagger$ **reverses** the order. Same as transpose/inverse.
- "Adjoint" also means the classical adjugate (cofactor matrix) in older linear-algebra texts. Different thing — in QC it always means conjugate transpose.
