# Amplitudes and qubits

A qubit is a unit vector in $\mathbb{C}^2$. The coefficients are **amplitudes** — complex numbers that are *not* probabilities, but square to them. Everything quantum computing does that classical probability can't comes from amplitudes being complex and therefore able to **cancel**.

Builds on [`hermitian-adjoint`](../hermitian-adjoint/README.md) (bra–ket, inner products).

## The qubit

$$\lvert \psi \rangle = \alpha \lvert 0 \rangle + \beta \lvert 1 \rangle = \begin{bmatrix} \alpha \\ \beta \end{bmatrix},
\qquad \alpha, \beta \in \mathbb{C}$$

- $\alpha, \beta$ are the **amplitudes** of the basis states $\lvert 0 \rangle, \lvert 1 \rangle$.
- **Normalization** (it's a *unit* vector): $\lvert\alpha\rvert^2 + \lvert\beta\rvert^2 = 1$.
- **Born rule** — measuring in the computational basis:

$$\Pr(0) = \lvert \alpha \rvert^2, \qquad \Pr(1) = \lvert \beta \rvert^2$$

- Extract an amplitude by projecting: $\alpha = \langle 0 \rvert \psi \rangle$, $\beta = \langle 1 \rvert \psi \rangle$.
- A **bit** is one of $\{0, 1\}$. A **probabilistic bit** is a point on the line $p_0 + p_1 = 1$. A **qubit** is a point on the unit sphere in $\mathbb{C}^2$ — strictly bigger, and the extra room is where phase lives.

### Named states worth knowing

| State | Vector | $\Pr(0), \Pr(1)$ |
|---|---|---|
| $\lvert 0 \rangle$ | $\begin{bmatrix} 1 & 0\end{bmatrix}^T$ | $1, 0$ |
| $\lvert 1 \rangle$ | $\begin{bmatrix} 0 & 1\end{bmatrix}^T$ | $0, 1$ |
| $\lvert + \rangle = \tfrac{1}{\sqrt2}(\lvert 0\rangle + \lvert 1\rangle)$ | $\tfrac{1}{\sqrt2}\begin{bmatrix} 1 & 1\end{bmatrix}^T$ | $\tfrac12, \tfrac12$ |
| $\lvert - \rangle = \tfrac{1}{\sqrt2}(\lvert 0\rangle - \lvert 1\rangle)$ | $\tfrac{1}{\sqrt2}\begin{bmatrix} 1 & -1\end{bmatrix}^T$ | $\tfrac12, \tfrac12$ |
| $\lvert i \rangle = \tfrac{1}{\sqrt2}(\lvert 0\rangle + i\lvert 1\rangle)$ | $\tfrac{1}{\sqrt2}\begin{bmatrix} 1 & i\end{bmatrix}^T$ | $\tfrac12, \tfrac12$ |

The last three are **indistinguishable by a single measurement** — same probabilities, different states. The difference is phase, and it shows up only after further gates.

## Amplitudes ≠ probabilities

This is the whole point.

| | Probabilities $p_k$ | Amplitudes $\alpha_k$ |
|---|---|---|
| Field | $\mathbb{R}_{\geq 0}$ | $\mathbb{C}$ |
| Constraint | $\sum_k p_k = 1$ ($L^1$) | $\sum_k \lvert \alpha_k \rvert^2 = 1$ ($L^2$) |
| Combining paths | add, always grow | add, **can cancel** |
| Evolution | stochastic matrix | unitary matrix |

Two paths to the same outcome with amplitudes $\alpha_1, \alpha_2$ give

$$\Pr = \lvert \alpha_1 + \alpha_2 \rvert^2 = \lvert\alpha_1\rvert^2 + \lvert\alpha_2\rvert^2 + 2\,\mathrm{Re}(\alpha_1^\ast \alpha_2)$$

That last **interference term** is the quantum surplus. It's positive (constructive), negative (destructive), or zero. Classically it's always zero.

### Interference in one line: $H H = I$

$$\lvert 0 \rangle \;\xrightarrow{\;H\;}\; \tfrac{1}{\sqrt2}\big(\lvert 0\rangle + \lvert 1\rangle\big) \;\xrightarrow{\;H\;}\; \tfrac12\big(\lvert 0\rangle + \lvert 1\rangle\big) + \tfrac12\big(\lvert 0\rangle - \lvert 1\rangle\big) = \lvert 0 \rangle$$

The $\lvert 1 \rangle$ amplitudes are $+\tfrac12$ and $-\tfrac12$ — they **destructively cancel**. A classical coin flipped twice is still random; a qubit "flipped" twice is back to certainty. Algorithms (Deutsch–Jozsa, Grover, Shor) are all engineered interference: arrange phases so wrong answers cancel and right ones add.

## Phase

Write amplitudes in polar form, $\alpha = r_\alpha e^{i\phi_\alpha}$.

- **Global phase** — $e^{i\gamma}\lvert\psi\rangle$ is *physically the same state*. Every probability $\lvert e^{i\gamma}\alpha_k\rvert^2 = \lvert\alpha_k\rvert^2$ is untouched, and it stays untouched under any later gate. Unobservable, always.
- **Relative phase** — the *difference* $\phi_\beta - \phi_\alpha$ between amplitudes. Very observable: it's what separates $\lvert+\rangle$ from $\lvert-\rangle$ from $\lvert i\rangle$.

So a qubit has 4 real parameters, minus 1 for normalization, minus 1 for global phase = **2 real degrees of freedom**:

$$\lvert \psi \rangle = \cos\tfrac{\theta}{2}\,\lvert 0 \rangle + e^{i\varphi}\sin\tfrac{\theta}{2}\,\lvert 1 \rangle$$

which is exactly a point $(\theta, \varphi)$ on the **Bloch sphere**. $\lvert 0\rangle$ at the north pole, $\lvert 1\rangle$ at the south, $\lvert \pm \rangle$ and $\lvert \pm i\rangle$ around the equator. Antipodal $\Leftrightarrow$ orthogonal. Gates are rotations of the sphere.

## $n$ qubits: $2^n$ amplitudes

Compose with the **tensor product** (Kronecker product):

$$\lvert 0 \rangle \otimes \lvert 1 \rangle \equiv \lvert 01 \rangle = \begin{bmatrix} 0 & 1 & 0 & 0 \end{bmatrix}^T$$

An $n$-qubit state is one unit vector in $\mathbb{C}^{2^n}$, indexed by bitstrings:

$$\lvert \psi \rangle = \sum_{x \in \{0,1\}^n} \alpha_x \lvert x \rangle, \qquad \sum_x \lvert \alpha_x \rvert^2 = 1$$

- Amplitudes grow **exponentially**: 50 qubits = $2^{50} \approx 10^{15}$ complex numbers. That's the simulation memory wall (see [`quantum-simulation-on-gpu`](../quantum-simulation-on-gpu/README.md)) and the source of the hoped-for speedup.
- Ordering convention matters: `|q1 q0>` big-endian (textbooks) vs little-endian (Qiskit). A constant source of transposed results.

### Product vs entangled

$\lvert \psi \rangle$ is a **product state** if it factors as $\lvert a \rangle \otimes \lvert b \rangle$; otherwise it's **entangled**.

$$\lvert \Phi^+ \rangle = \tfrac{1}{\sqrt2}\big(\lvert 00 \rangle + \lvert 11 \rangle\big)$$

For 2 qubits the test is a determinant: with amplitudes $\alpha_{00}, \alpha_{01}, \alpha_{10}, \alpha_{11}$, the state is a product **iff**

$$\alpha_{00}\,\alpha_{11} - \alpha_{01}\,\alpha_{10} = 0$$

For $\lvert \Phi^+ \rangle$ that's $\tfrac12 - 0 \neq 0$ — entangled. Note $n$ qubits need $2^n$ amplitudes but only $2n$ if they were all separable; the gap is entanglement, and it's most of the space.

## Measurement

Measuring is not a gate — it's the one non-unitary, irreversible step.

- Outcome $x$ with $\Pr(x) = \lvert \alpha_x \rvert^2$, then the state **collapses** to $\lvert x \rangle$. Re-measuring gives $x$ again.
- **Partial measurement**: measuring qubit 0 of $\lvert \psi\rangle$ and getting $b$ leaves the projected-and-renormalized remainder,

$$\lvert \psi' \rangle = \frac{P_b \lvert \psi \rangle}{\lVert P_b \lvert \psi \rangle \rVert}, \qquad P_b = \lvert b \rangle\langle b \rvert \otimes I$$

  For $\lvert \Phi^+\rangle$, measuring the first qubit as $0$ forces the second to $0$ — that's the correlation people call "spooky".
- **Other bases**: measuring in $\{\lvert+\rangle, \lvert-\rangle\}$ means using those amplitudes instead, $\Pr(+) = \lvert\langle + \rvert \psi\rangle\rvert^2$. Equivalently: apply $H$, then measure in the computational basis. The basis is a choice; the state isn't.
- You never see an amplitude. One shot gives one bitstring. Estimating $\lvert\alpha_x\rvert^2$ to precision $\epsilon$ needs $O(1/\epsilon^2)$ shots, and phases need interference tricks to reveal at all.

## In Python

```python
import numpy as np

ket0, ket1 = np.array([1, 0], dtype=complex), np.array([0, 1], dtype=complex)

psi = (ket0 + 1j * ket1) / np.sqrt(2)        # |i> = (|0> + i|1>)/sqrt(2)
probs = np.abs(psi) ** 2                     # Born rule -> [0.5, 0.5]
assert np.isclose(probs.sum(), 1.0)          # normalization

# global phase is invisible: same probabilities, forever
assert np.allclose(np.abs(np.exp(1j * 0.7) * psi) ** 2, probs)

# interference: H twice returns |0> exactly
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
assert np.allclose(H @ (H @ ket0), ket0)     # |1> amplitudes +1/2 and -1/2 cancel

# two qubits: tensor product, 2^n amplitudes
phi_plus = (np.kron(ket0, ket0) + np.kron(ket1, ket1)) / np.sqrt(2)
a00, a01, a10, a11 = phi_plus
entangled = not np.isclose(a00 * a11 - a01 * a10, 0)   # True

# sampling: you observe bitstrings, never amplitudes
rng = np.random.default_rng(0)
outcomes = rng.choice(4, size=1000, p=np.abs(phi_plus) ** 2)
# -> only 0 (|00>) and 3 (|11>), roughly 50/50
```

See [`amplitudes.py`](amplitudes.py) for a runnable version, plus partial measurement and a Bloch-angle round-trip.

## Gotchas

- **Amplitudes can be negative or complex; probabilities can't.** Forgetting this is forgetting why quantum computing works.
- Normalize with $\sum \lvert \alpha_x \rvert^2 = 1$, not $\sum \alpha_x = 1$.
- $\lvert\alpha\rvert^2 = \alpha^\ast\alpha$, not $\alpha^2$ — the latter is complex and meaningless as a probability.
- Global phase is unobservable; **relative** phase is the entire game. Don't "normalize away" a relative phase.
- Superposition is not "the qubit is secretly 0 or 1 with some probability". That's a classical mixture, has no interference, and is a different object (a density matrix, not a state vector).
- A qubit in superposition holds $2^n$ amplitudes, but measurement gives you $n$ bits. No free lunch — you need interference to concentrate amplitude on the answer first.
- Endianness: $\lvert 01 \rangle$ means different vectors in different frameworks. Check before trusting a state vector dump.
