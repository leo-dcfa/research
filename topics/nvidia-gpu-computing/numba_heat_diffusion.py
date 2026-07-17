"""2D heat diffusion with a hand-written CUDA kernel via Numba.

Shows the level below CuPy: when your simulation needs custom per-cell logic,
you write the kernel yourself. Each GPU thread updates one grid cell per step
using the 5-point stencil (average of its 4 neighbours), simulating heat
spreading out from a hot square in the middle of a cold plate.

Run on a CUDA machine or Google Colab (GPU runtime):
    pip install numba
    python numba_heat_diffusion.py
"""

import numpy as np
from numba import cuda

N = 512  # grid is N x N
STEPS = 1000  # simulation steps
ALPHA = 0.2  # diffusion rate (keep < 0.25 for stability)


@cuda.jit
def diffuse_step(grid, out):
    """One time step: every thread computes one cell from its neighbours."""
    i, j = cuda.grid(2)  # this thread's (row, col) in the grid
    if 1 <= i < grid.shape[0] - 1 and 1 <= j < grid.shape[1] - 1:
        lap = (
            grid[i + 1, j]
            + grid[i - 1, j]
            + grid[i, j + 1]
            + grid[i, j - 1]
            - 4.0 * grid[i, j]
        )
        out[i, j] = grid[i, j] + ALPHA * lap
    elif i < grid.shape[0] and j < grid.shape[1]:
        out[i, j] = grid[i, j]  # hold the fixed boundary


def main() -> None:
    grid = np.zeros((N, N), dtype=np.float32)
    grid[N // 2 - 20 : N // 2 + 20, N // 2 - 20 : N // 2 + 20] = 100.0  # hot square

    # Move both buffers to the GPU once; ping-pong between them on-device.
    d_a = cuda.to_device(grid)
    d_b = cuda.device_array_like(d_a)

    threads = (16, 16)  # 256 threads per block
    blocks = ((N + 15) // 16, (N + 15) // 16)

    for _ in range(STEPS):
        diffuse_step[blocks, threads](d_a, d_b)
        d_a, d_b = d_b, d_a  # swap: output becomes next input

    result = d_a.copy_to_host()  # single transfer back at the very end
    print(f"after {STEPS} steps: centre temp = {result[N // 2, N // 2]:.2f}")
    print(f"grid min/max = {result.min():.2f} / {result.max():.2f}")


if __name__ == "__main__":
    main()
