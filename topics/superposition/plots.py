"""Generate the figures for the superposition notes.

Run:  uv run python topics/superposition/plots.py
Writes PNGs into this folder.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from toy_model import train

HERE = os.path.dirname(__file__)


def fig_geometry():
    """m=2 bottleneck: as features get sparser, columns of W arrange into
    regular polygons — the visual signature of superposition."""
    sparsities = [0.0, 0.7, 0.9, 0.97]
    n_features = 5
    fig, axes = plt.subplots(1, len(sparsities), figsize=(14, 3.6))

    for ax, s in zip(axes, sparsities):
        # equal importance -> clean, symmetric geometry
        w, _ = train(n_features=n_features, m_hidden=2, sparsity=s, steps=6000)
        colors = plt.cm.viridis(np.linspace(0, 1, n_features))
        for i in range(n_features):
            ax.plot([0, w[0, i]], [0, w[1, i]], "-", color=colors[i], lw=2)
            ax.plot(w[0, i], w[1, i], "o", color=colors[i], ms=7)
        lim = 1.3
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.axhline(0, color="0.85", lw=0.8, zorder=0)
        ax.axvline(0, color="0.85", lw=0.8, zorder=0)
        n_rep = int(np.sum(np.linalg.norm(w, axis=0) > 0.1))
        ax.set_title(f"sparsity={s}\n{n_rep}/{n_features} features in 2D")

    fig.suptitle(
        "Feature geometry in a 2D bottleneck (each arrow = one feature's column of W)",
        y=1.02,
    )
    fig.tight_layout()
    path = os.path.join(HERE, "geometry.png")
    fig.savefig(path, dpi=110, bbox_inches="tight")
    print("wrote", path)


def fig_gram():
    """W^T W heatmaps. Diagonal = how strongly a feature is represented;
    off-diagonal = interference between features sharing directions."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, s in zip(axes, [0.0, 0.8, 0.95]):
        w, _ = train(n_features=20, m_hidden=5, sparsity=s, steps=8000)
        gram = w.T @ w
        im = ax.imshow(gram, cmap="RdBu_r", vmin=-1, vmax=1)
        n_rep = int(np.sum(np.sqrt(np.diag(gram)) > 0.1))
        ax.set_title(f"sparsity={s}\n{n_rep}/20 represented")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(im, ax=axes, shrink=0.8, label="W$^T$W entry")
    fig.suptitle(
        "W$^T$W as sparsity rises: dense -> 5 orthogonal features; "
        "sparse -> >5 features with off-diagonal interference (superposition)",
        y=1.04,
    )
    path = os.path.join(HERE, "gram.png")
    fig.savefig(path, dpi=110, bbox_inches="tight")
    print("wrote", path)


def fig_phase():
    """Number of features represented vs sparsity: the phase transition from
    'store the m most important' to 'pack them all in via superposition'."""
    sparsities = np.linspace(0.0, 0.98, 12)
    m_hidden = 5
    n_features = 20
    counts = []
    for s in sparsities:
        w, _ = train(n_features=n_features, m_hidden=m_hidden, sparsity=s, steps=6000)
        counts.append(int(np.sum(np.linalg.norm(w, axis=0) > 0.1)))

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(sparsities, counts, "o-", color="#c0392b", lw=2)
    ax.axhline(m_hidden, ls="--", color="0.5", label=f"m = {m_hidden} (bottleneck dim)")
    ax.set_xlabel("feature sparsity  (P[feature = 0])")
    ax.set_ylabel("# features represented  (||W$_i$|| > 0.1)")
    ax.set_title("More sparsity -> more features squeezed into m dimensions")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    path = os.path.join(HERE, "phase.png")
    fig.savefig(path, dpi=110, bbox_inches="tight")
    print("wrote", path)


if __name__ == "__main__":
    fig_geometry()
    fig_gram()
    fig_phase()
