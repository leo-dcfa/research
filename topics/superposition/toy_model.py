"""Toy model of superposition (Anthropic, 2022) from scratch in NumPy.

A tiny linear-ish autoencoder is forced to squeeze `n` sparse features through
an `m`-dimensional bottleneck (m < n) and reconstruct them:

    h    = W x                     # (m,)  compress into the bottleneck
    x_'  = ReLU(W^T h + b)         # (n,)  decompress (tied weights)
    L    = sum_i  I_i (x_i - x'_i)^2

When features are **dense**, the model can only afford to represent the `m` most
important ones (an orthogonal basis) and ignores the rest. When features are
**sparse**, interference between features is rare, so the model packs *more than
m* features into the m dimensions at non-orthogonal angles — **superposition**.

Manual gradients (the model is tiny, so no autograd needed) + Adam.

Run:  uv run python topics/superposition/toy_model.py
"""

import numpy as np

rng = np.random.default_rng(0)


def make_batch(batch, n_features, sparsity, importance_decay=0.9):
    """Sparse feature vectors. Each feature is 0 w.p. `sparsity`, else U[0,1].

    Returns (X, importance) with X shape (batch, n_features).
    """
    active = rng.random((batch, n_features)) > sparsity
    values = rng.random((batch, n_features))
    x = active * values
    importance = importance_decay ** np.arange(n_features)  # I_0 > I_1 > ...
    return x, importance


def relu(z):
    return np.maximum(z, 0.0)


def forward(x, w, b):
    """x:(B,n)  w:(m,n)  b:(n,) -> (out, hidden, pre)."""
    hidden = x @ w.T  # (B, m)
    pre = hidden @ w + b  # (B, n)
    return relu(pre), hidden, pre


def loss_and_grads(x, imp, w, b):
    batch = x.shape[0]
    out, hidden, pre = forward(x, w, b)
    diff = out - x  # (B, n)
    loss = np.mean(np.sum(imp * diff**2, axis=1))

    d_out = 2.0 * imp * diff / batch  # (B, n)
    g = d_out * (pre > 0)  # ReLU backward, (B, n)
    grad_b = g.sum(axis=0)  # (n,)
    # W appears twice: pre = W^T (W x) + b
    grad_w = hidden.T @ g + (g @ w.T).T @ x  # (m, n)
    return loss, grad_w, grad_b


def train(n_features=20, m_hidden=5, sparsity=0.7, steps=8000, batch=1024, lr=1e-2):
    """Train one model; returns weights W (m, n) and bias b (n,)."""
    w = rng.normal(0, 1.0, size=(m_hidden, n_features)) * 0.1
    b = np.zeros(n_features)

    # Adam state
    mw, vw = np.zeros_like(w), np.zeros_like(w)
    mb, vb = np.zeros_like(b), np.zeros_like(b)
    b1, b2, eps = 0.9, 0.999, 1e-8

    for t in range(1, steps + 1):
        x, imp = make_batch(batch, n_features, sparsity)
        _, gw, gb = loss_and_grads(x, imp, w, b)

        for p, g, mm, vv in ((w, gw, mw, vw), (b, gb, mb, vb)):
            mm[:] = b1 * mm + (1 - b1) * g
            vv[:] = b2 * vv + (1 - b2) * g * g
            mhat = mm / (1 - b1**t)
            vhat = vv / (1 - b2**t)
            p -= lr * mhat / (np.sqrt(vhat) + eps)

    return w, b


if __name__ == "__main__":
    w, b = train(n_features=20, m_hidden=5, sparsity=0.7)
    gram = w.T @ w  # (n, n): diagonal = ||W_i||^2, off-diag = interference
    norms = np.sqrt(np.diag(gram))
    represented = int(np.sum(norms > 0.1))
    print(f"features represented (||W_i|| > 0.1): {represented} / {w.shape[1]}")
    print("feature norms:", np.round(norms, 2))
