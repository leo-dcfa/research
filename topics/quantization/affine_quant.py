"""Affine (uniform) quantization from scratch.

Shows the core math for symmetric vs asymmetric INT8, then demonstrates why
*granularity* (per-tensor vs per-channel) is the single biggest lever when the
data has outliers.

Run:  uv run python topics/quantization/affine_quant.py
"""

import numpy as np


def quantize_affine(x, bits=8, symmetric=False):
    """Quantize array x to `bits` ints. Returns (x_q, scale, zero_point)."""
    qmax = 2**bits - 1
    if symmetric:
        # range centered on 0; zero-point is fixed at the midpoint
        a = np.max(np.abs(x))
        scale = (2 * a) / qmax
        zero = (qmax) // 2  # integer that maps to real 0
    else:
        lo, hi = np.min(x), np.max(x)
        scale = (hi - lo) / qmax
        zero = np.round(-lo / scale)

    x_q = np.clip(np.round(x / scale) + zero, 0, qmax).astype(np.int64)
    return x_q, scale, zero


def dequantize_affine(x_q, scale, zero):
    return scale * (x_q.astype(np.float64) - zero)


def roundtrip_error(x, **kw):
    x_q, s, z = quantize_affine(x, **kw)
    x_hat = dequantize_affine(x_q, s, z)
    return np.sqrt(np.mean((x - x_hat) ** 2))  # RMSE


def main():
    rng = np.random.default_rng(0)

    # --- 1. symmetric vs asymmetric on skewed (all-positive) data --------
    print("=== symmetric vs asymmetric (post-ReLU-like, all >= 0) ===")
    acts = np.abs(rng.normal(size=10_000))  # skewed, one-sided
    print(f"symmetric  RMSE: {roundtrip_error(acts, symmetric=True):.5f}")
    print(f"asymmetric RMSE: {roundtrip_error(acts, symmetric=False):.5f}")
    print("-> asymmetric wins on one-sided data: it doesn't waste half the range.\n")

    # --- 2. per-tensor vs per-channel with an outlier --------------------
    print("=== per-tensor vs per-channel (one row has a big outlier) ===")
    W = rng.normal(scale=0.1, size=(8, 256))  # weight matrix, 8 output channels
    W[3, 0] = 12.0  # single fat outlier in channel 3

    # per-tensor: ONE scale for the whole matrix -> outlier coarsens everything
    per_tensor = roundtrip_error(W, symmetric=True)

    # per-channel: one scale PER ROW -> outlier only hurts its own channel
    per_channel = np.sqrt(
        np.mean([roundtrip_error(W[i], symmetric=True) ** 2 for i in range(W.shape[0])])
    )

    print(f"per-tensor  RMSE: {per_tensor:.5f}")
    print(f"per-channel RMSE: {per_channel:.5f}")
    print(f"-> per-channel is {per_tensor / per_channel:.1f}x more accurate here.")
    print("   This is why 4-bit LLM weights are always quantized per-channel/group.")


if __name__ == "__main__":
    main()
