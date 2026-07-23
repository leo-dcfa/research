# Toy Models of Superposition

Anthropic's [*Toy Models of Superposition*](https://transformer-circuits.pub/2022/toy_model/index.html) (Elhage et al., 2022). A minimal, fully-understood model that reproduces **superposition** — a network representing **more features than it has dimensions** by packing them into nearly-orthogonal directions and tolerating a little interference.

## The core question

Real networks seem to have more "features" (interpretable directions) than neurons. Either features get a dedicated neuron each, or many features share dimensions — **superposition**. The toy model isolates *when and why* superposition happens.

## Setup

- **Synthetic features** $x \in \mathbb{R}^n$. Each feature is:
  - **sparse** — active (nonzero) with prob $1 - S$, where $S$ is the *sparsity*. When active, $x_i \sim U[0,1]$; else $0$.
  - **importance-weighted** — feature $i$ has importance $I_i$ (typically decaying, e.g. $I_i = r^i$).
- **Model** — squeeze $n$ features through an $m$-dim bottleneck ($m < n$) and reconstruct:

$$h = W x \in \mathbb{R}^m, \qquad x' = \operatorname{ReLU}\!\big(W^\top W\, x + b\big)\in\mathbb{R}^n$$

  A single weight matrix $W \in \mathbb{R}^{m\times n}$ ties encoder and decoder (weight tying), plus a bias $b$ and an output ReLU.
- **Loss** — importance-weighted reconstruction error:

$$L = \mathbb{E}_{x}\sum_{i=1}^{n} I_i\,\big(x_i - x'_i\big)^2$$

## $W^\top W$ is the whole story

Column $W_i \in \mathbb{R}^m$ is feature $i$'s **embedding direction** in the bottleneck. The Gram matrix $W^\top W$ captures everything:

| Entry | Meaning |
|---|---|
| $\lVert W_i\rVert^2 = (W^\top W)_{ii}$ | how strongly feature $i$ is represented ($0$ = dropped) |
| $W_i\cdot W_j = (W^\top W)_{ij}$ | **interference** between features $i$ and $j$ |

- **No superposition** $\iff$ columns orthogonal $\iff$ $W^\top W$ diagonal. Then at most $m$ features can be stored.
- **Superposition** $\iff$ nonzero off-diagonals — $>m$ features stored, at the cost of interference. The output ReLU + bias let the model *clean up* small negative interference, and **sparsity** means interfering features rarely fire together.

## What happens: the phase change

Sweep sparsity $S$ from dense to sparse:

- **Dense ($S\to 0$)** — behaves like **PCA**: keep the $m$ most important features orthogonally, throw the rest away. No superposition.
- **Sparse ($S\to 1$)** — **superposition**: represent *all* $n$ features (even with $m \ll n$) as non-orthogonal directions. Rare co-activation makes interference cheap.
- The transition is a **sharp phase change**, not gradual. Per feature, as sparsity rises it jumps: *not represented* → *in superposition* → *dedicated dimension*.

## Geometry: features form uniform polytopes

Superposed features don't scatter randomly — they snap into symmetric arrangements that spread directions apart (minimising interference):

- 2 features in 2D → **antipodal pair** (180°)
- 3 → **triangle** (120°), 4 → **square/tetrahedron**, 5 → **pentagon**, …
- Higher dims → tetrahedra, and **tegum products** (orthogonal stacks) of these.

**Feature dimensionality** — fraction of a bottleneck dimension a feature "owns":

$$D_i = \frac{\lVert W_i\rVert^2}{\sum_{j}\big(\hat W_i \cdot W_j\big)^2}, \qquad \hat W_i = \frac{W_i}{\lVert W_i\rVert}$$

The model gets "stuck" at **sticky fractions** — $D_i = 1$ (own dimension), $\tfrac12$ (antipodal pair), $\tfrac23$ (triangle), $\tfrac34$ (tetrahedron)… — with unstable gaps between them.

## Why it matters (interpretability)

- **Explains polysemanticity** — a neuron fires for several unrelated concepts *because* features are in superposition, not one-per-neuron.
- **Motivates sparse dictionary learning / SAEs** — if features are a sparse overcomplete code crammed into activations, an overcomplete sparse autoencoder can pull them back out.
- Connects to **adversarial fragility** (interference is an attack surface) and **computation in superposition** (a ReLU-hidden variant computes $|x|$-like functions on more features than it has neurons).

## Key quantities at a glance

| Symbol | Meaning |
|---|---|
| $n, m$ | \# features, bottleneck width ($m < n$) |
| $S$ | sparsity — $P(\text{feature} = 0)$ |
| $I_i$ | importance weight of feature $i$ in the loss |
| $W_i$ | embedding direction of feature $i$ (column of $W$) |
| $\lVert W_i\rVert^2$ | representation strength ($0$ ⇒ feature dropped) |
| $W_i\cdot W_j$ | interference between two features |
| $D_i$ | feature dimensionality (fraction of a dim owned) |

## Code

- `toy_model.py` — trains the ReLU-output model in pure NumPy (manual backprop + Adam) and sweeps sparsity for $n=8$ features in an $m=2$ bottleneck. Reproduces the phase change — features kept grows $2 \to 4 \to 5 \to 6$ as $S$ rises while $\sum_i D_i$ stays $\approx m$ — and prints $W^\top W$: near-diagonal when dense (PCA), full of $\approx\!-1$ antipodal interference when sparse (superposition).

```bash
uv run python topics/toy-models-of-superposition/toy_model.py
```
