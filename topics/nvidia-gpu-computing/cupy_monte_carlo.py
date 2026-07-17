"""Estimate pi by Monte Carlo on an NVIDIA GPU using CuPy.

Throw N random points into the unit square; the fraction landing inside the
quarter circle approximates pi/4. This is embarrassingly parallel — a perfect
fit for the GPU, and CuPy makes it read exactly like NumPy.

Run on a CUDA machine or Google Colab (GPU runtime):
    pip install cupy-cuda12x   # match your CUDA major version
    python cupy_monte_carlo.py
"""

import cupy as cp


def estimate_pi(n: int) -> float:
    # Random points generated *on the GPU* — no CPU->GPU transfer.
    x = cp.random.random(n, dtype=cp.float32)
    y = cp.random.random(n, dtype=cp.float32)
    inside = (x * x + y * y) <= 1.0
    # .mean() reduces on-device; only the single scalar comes back to the CPU.
    return float(4.0 * inside.mean())


if __name__ == "__main__":
    n = 100_000_000  # 100M points — big enough for the GPU to shine

    # Warm up: the first call compiles kernels, so it is not representative.
    estimate_pi(1000)

    start = cp.cuda.Event()
    end = cp.cuda.Event()
    start.record()
    pi = estimate_pi(n)
    end.record()
    end.synchronize()  # wait for the async GPU work to actually finish

    ms = cp.cuda.get_elapsed_time(start, end)
    print(f"pi ~= {pi:.6f}  ({n:,} points in {ms:.1f} ms)")
