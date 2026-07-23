"""Toy model of superposition (Elhage et al., 2022), from scratch in NumPy.

We squeeze n sparse features through an m-dim bottleneck and reconstruct them:

    h  = W x                          # encode  (m-dim, m < n)
    x' = ReLU(W^T W x + b)            # decode
    L  = E_x  sum_i  I_i (x_i - x'_i)^2

The single matrix W ties encoder and decoder. Column W_i is feature i's
direction in the bottleneck; the Gram matrix W^T W is the whole story:
its diagonal is how strongly each feature is represented, its off-diagonal
is interference between features.

The demo sweeps the sparsity of the same model and shows the phase change:
a DENSE model keeps only m features (PCA-like), while sparser models pack
progressively more of the n features into the m-dim plane (superposition).

Run:  uv run python topics/toy-models-of-superposition/toy_model.py
"""

import numpy as np


def make_batch(rng, n, batch, sparsity):
    """n features, each active w.p. (1 - sparsity); value ~ U[0,1] when active."""
    x = rng.uniform(0.0, 1.0, size=(batch, n))
    mask = rng.uniform(size=(batch, n)) > sparsity
    return x * mask


def train(n, m, sparsity, importance, steps=8000, batch=1024, lr=1e-2, seed=0):
    """Train the ReLU-output model with manual backprop + Adam. Returns (W, b)."""
    rng = np.random.default_rng(seed)
    W = rng.normal(scale=0.1, size=(m, n))
    b = np.zeros(n)
    imp = importance  # (n,) importance weights I_i

    # Adam state
    mW = np.zeros_like(W)
    vW = np.zeros_like(W)
    mb = np.zeros_like(b)
    vb = np.zeros_like(b)
    b1, b2, eps = 0.9, 0.999, 1e-8

    for t in range(1, steps + 1):
        x = make_batch(rng, n, batch, sparsity)  # (B, n)

        # forward
        h = x @ W.T  # (B, m)
        z = h @ W + b  # (B, n)
        xp = np.maximum(z, 0.0)  # (B, n)

        # backward (mean over batch, importance-weighted)
        dxp = (2.0 / batch) * imp * (xp - x)  # (B, n)
        dz = dxp * (z > 0.0)  # ReLU grad
        db = dz.sum(axis=0)  # (n,)
        dh = dz @ W.T  # (B, m)
        dW = h.T @ dz + dh.T @ x  # outer + inner W

        # Adam update
        for p, g, mm, vv in ((W, dW, mW, vW), (b, db, mb, vb)):
            mm[...] = b1 * mm + (1 - b1) * g
            vv[...] = b2 * vv + (1 - b2) * g * g
            mhat = mm / (1 - b1**t)
            vhat = vv / (1 - b2**t)
            p -= lr * mhat / (np.sqrt(vhat) + eps)

    return W, b


def feature_dimensionality(W):
    """D_i = ||W_i||^2 / sum_j (Wi_hat . W_j)^2  — fraction of a dim feature i owns."""
    norms = np.linalg.norm(W, axis=0)  # ||W_i||
    what = W / np.clip(norms, 1e-9, None)  # unit columns
    interference = (what.T @ W) ** 2  # (n, n): (Wi_hat . W_j)^2
    return norms**2 / interference.sum(axis=0)


def report(W):
    """Return (# represented features, sum of feature dimensionalities)."""
    norms2 = (W**2).sum(axis=0)
    represented = int((norms2 > 0.1).sum())
    return represented, feature_dimensionality(W).sum()


def main():
    n, m = 8, 2  # 8 features into a 2-D bottleneck
    importance = 0.9 ** np.arange(n)  # feature i matters 0.9^i  (mild decay)

    print(f"n = {n} features, m = {m} bottleneck dims, importance = 0.9^i")
    print("Same model, same importances -- only the sparsity S changes.\n")
    print(f"{'sparsity S':>11} | {'features kept':>13} | {'sum of D_i':>10}")
    print("-" * 41)

    W_at = {}
    for S in (0.0, 0.7, 0.9, 0.99):
        W, _ = train(n, m, sparsity=S, importance=importance)
        W_at[S] = W
        rep, sumD = report(W)
        print(f"{S:>11.2f} | {rep:>10}/{n} | {sumD:>10.2f}")

    print("\nPhase change: as features get sparser, the model represents MORE")
    print(f"than m={m} of them -- superposition. Yet sum(D_i) stays ~= m: the")
    print("bottleneck still holds only m dimensions' worth of stuff, now")
    print("spread thinly across many non-orthogonal, mutually-interfering features.\n")

    # Show the Gram matrix at the two extremes.
    for S in (0.0, 0.99):
        G = W_at[S].T @ W_at[S]
        tag = (
            "dense -> near-diagonal (orthogonal, PCA)"
            if S == 0.0
            else "sparse -> off-diagonals appear (interference / superposition)"
        )
        print(f"W^T W at S={S:.2f}  [{tag}]:")
        print(np.array2string(G, precision=2, suppress_small=True), "\n")


if __name__ == "__main__":
    main()
