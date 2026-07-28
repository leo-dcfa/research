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

## Dirac (bra–ket) notation

The name comes from splitting "bra-c-ket" $\langle \cdot \rvert \cdot \rangle$ in two. The dagger is exactly what maps one half to the other.

- **Ket** $\lvert v \rangle$ — a column vector, $n \times 1$. The state.
- **Bra** $\langle v \rvert$ — a row vector, $1 \times n$. Its adjoint (Definition 6.7):

$$\langle v \rvert \equiv \lvert v \rangle^\dagger
\qquad\qquad
\lvert v \rangle \equiv \langle v \rvert^\dagger$$

Both directions hold because $(A^\dagger)^\dagger = A$. A bra is a *functional*: feed it a ket, get a number out.

### The computational basis

$$\lvert 0 \rangle = \begin{bmatrix} 1 \\ 0 \end{bmatrix},\quad
\lvert 1 \rangle = \begin{bmatrix} 0 \\ 1 \end{bmatrix}
\qquad\Longrightarrow\qquad
\langle 0 \rvert = \begin{bmatrix} 1 & 0 \end{bmatrix},\quad
\langle 1 \rvert = \begin{bmatrix} 0 & 1 \end{bmatrix}$$

(Example 6.2.) Both are real, so here the dagger is *just* a transpose. But for $\lvert v\rangle = \tfrac{1}{\sqrt2}\big(\lvert 0\rangle + i\lvert 1\rangle\big)$ you get $\langle v \rvert = \tfrac{1}{\sqrt2}\begin{bmatrix} 1 & -i \end{bmatrix}$ — the $i$ flips sign.

### Products — it's all just matrix shapes

| Expression | Shapes | Result | Name |
|---|---|---|---|
| $\langle v \rvert w \rangle$ | $(1{\times}n)(n{\times}1)$ | scalar | **inner product** / bracket |
| $\lvert v \rangle \langle w \rvert$ | $(n{\times}1)(1{\times}n)$ | $n \times n$ matrix | **outer product** |
| $A \lvert v \rangle$ | $(n{\times}n)(n{\times}1)$ | ket | operator acting on a state |
| $\langle v \rvert A \lvert w \rangle$ | $(1{\times}n)(n{\times}n)(n{\times}1)$ | scalar | **matrix element** / expectation |

**Inner product** — conjugate the left slot:

$$\langle v \rvert w \rangle = v^\dagger w = \sum_k v_k^* w_k
\qquad
\langle v \rvert v \rangle = \sum_k \lvert v_k \rvert^2 \geq 0$$

So the norm is real and non-negative — that's the dagger earning its keep, and why amplitudes square to probabilities. Conjugate symmetry: $\langle w \rvert v \rangle = \langle v \rvert w \rangle^*$. Orthonormal basis: $\langle i \rvert j \rangle = \delta_{ij}$.

**Outer product** — builds operators out of states:

$$\lvert 0 \rangle \langle 0 \rvert = \begin{bmatrix} 1 \\ 0 \end{bmatrix}\begin{bmatrix} 1 & 0 \end{bmatrix} = \begin{bmatrix} 1 & 0 \\ 0 & 0 \end{bmatrix}
\qquad
\lvert 0 \rangle \langle 1 \rvert = \begin{bmatrix} 0 & 1 \\ 0 & 0 \end{bmatrix}$$

$\lvert v\rangle\langle v\rvert$ (unit $v$) is the **projector** onto $v$; $\big(\lvert v\rangle\langle w\rvert\big)^\dagger = \lvert w\rangle\langle v\rvert$. **Completeness**: $\sum_k \lvert k \rangle \langle k \rvert = I$, which is how you insert a resolution of the identity anywhere. Gates get written this way too, e.g. $X = \lvert 0\rangle\langle 1\rvert + \lvert 1\rangle\langle 0\rvert$.

### Reading a state

For $\lvert \psi \rangle = \alpha \lvert 0 \rangle + \beta \lvert 1 \rangle$ with $\lvert\alpha\rvert^2 + \lvert\beta\rvert^2 = 1$:

- Amplitude of outcome $k$: $\alpha = \langle 0 \rvert \psi \rangle$, $\beta = \langle 1 \rvert \psi \rangle$.
- **Born rule**: $\Pr(k) = \lvert \langle k \rvert \psi \rangle \rvert^2 = \langle \psi \rvert k \rangle \langle k \rvert \psi \rangle$.
- Expectation of an observable: $\langle A \rangle = \langle \psi \rvert A \lvert \psi \rangle$ — real precisely when $A^\dagger = A$.
- Gate then measure: $(U\lvert\psi\rangle)^\dagger = \langle\psi\rvert U^\dagger$, so $\lVert U \lvert \psi\rangle \rVert^2 = \langle \psi\rvert U^\dagger U \lvert\psi\rangle = 1$ — unitarity preserves total probability.

### Daggering an expression

Reverse the order, flip every bra $\leftrightarrow$ ket, conjugate the scalars, dagger the operators:

$$\big(\alpha\, \langle v \rvert A \lvert w \rangle\big)^\dagger = \alpha^* \langle w \rvert A^\dagger \lvert v \rangle$$

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

# Bra-ket: kets are columns, bras are their adjoint
ket0, ket1 = np.array([[1], [0]]), np.array([[0], [1]])
bra0 = ket0.conj().T                           # <0| = [[1, 0]]

# Outer product |0><1| -> matrix;  X = |0><1| + |1><0|
X = ket0 @ ket1.conj().T + ket1 @ ket0.conj().T

psi = np.array([[1], [1j]]) / np.sqrt(2)       # |psi> = (|0> + i|1>)/sqrt(2)
amp0 = (bra0 @ psi).item()                     # <0|psi>  -> amplitude
prob0 = abs(amp0) ** 2                         # Born rule -> 0.5
exp_X = (psi.conj().T @ X @ psi).item().real   # <psi|X|psi>, real since X† = X
```

See [`adjoint.py`](adjoint.py) for a runnable version with the examples above.

## Gotchas

- `.T` alone is **not** the adjoint for complex arrays — a silent, very common bug. Always `.conj().T`.
- Anti-linearity: pulling a scalar out of a bra conjugates it. $\langle \alpha v \rvert = \alpha^* \langle v \rvert$, while $\lvert \alpha v\rangle = \alpha \lvert v \rangle$.
- $(AB)^\dagger$ **reverses** the order. Same as transpose/inverse.
- "Adjoint" also means the classical adjugate (cofactor matrix) in older linear-algebra texts. Different thing — in QC it always means conjugate transpose.
