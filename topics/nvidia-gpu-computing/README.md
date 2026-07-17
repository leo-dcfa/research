# Running GPU simulations on NVIDIA GPUs

How to run numerical/scientific simulations (GPGPU computing) on NVIDIA hardware.

## The stack (bottom → top)

- **NVIDIA driver** — talks to the physical GPU. `nvidia-smi` shows GPUs, memory, driver + CUDA version.
- **CUDA Toolkit** — the compute platform (compiler `nvcc`, runtime, libs like cuBLAS/cuFFT/cuRAND). Everything above is built on it.
- **Python libraries** — where you'll actually work. Pick by how much control you need.

## Which Python library to use

| Library | Use it for | Feel |
|---|---|---|
| **CuPy** | Array math on GPU | Drop-in NumPy (`import cupy as cp`) |
| **Numba** (`@cuda.jit`) | Custom kernels, per-element sim logic | Write CUDA in Python |
| **PyTorch / JAX** | Anything needing autodiff or tensor ops | ML-style, but great for physics too |
| **NVIDIA Warp** | Particle / physics / spatial sims | Purpose-built for simulation kernels |
| **RAPIDS** (cuDF/cuML) | GPU dataframes + classic ML | pandas/scikit-learn on GPU |
| **PyCUDA** | Full low-level control | Raw CUDA C in Python strings |

**Rule of thumb:** start with **CuPy** (vectorised math) → drop to **Numba** when you need a custom kernel → use **Warp** for real particle/physics sims.

## The core mental model

- GPU is fast because it runs **thousands of threads in parallel** — good for *data-parallel* work (same op over a huge array).
- Data lives in **GPU memory**. Moving CPU↔GPU (`cp.asarray` / `.get()`) is the usual bottleneck — do it rarely, keep the sim loop entirely on-device.
- Speedups show up on **big arrays**; tiny problems are slower on GPU than CPU (transfer + launch overhead dominates).

## Where to actually run it

- **Local NVIDIA GPU** — install matching driver + CUDA, then `pip install cupy-cuda12x` (match your CUDA major version).
- **Google Colab** — free NVIDIA GPU, CuPy/PyTorch preinstalled. Fastest way to try the samples. *(Runtime → Change runtime type → GPU.)*
- **Cloud VMs** — AWS `g5`/`p4`, GCP `a2`, Lambda, RunPod, etc. Rent by the hour.
- ⚠️ **Not on this Mac** — Apple Silicon has no NVIDIA GPU/CUDA. The `.py` samples here run on a CUDA machine or Colab.

## Verify the GPU is there

```bash
nvidia-smi                       # driver + GPU list
python -c "import cupy; print(cupy.cuda.runtime.getDeviceProperties(0)['name'])"
```

## Samples in this folder

- [`cupy_monte_carlo.py`](cupy_monte_carlo.py) — estimate π with millions of random points (CuPy, vectorised).
- [`numba_heat_diffusion.py`](numba_heat_diffusion.py) — 2D heat diffusion with a custom CUDA kernel (Numba).

## Gotchas

- **CuPy wheel must match your CUDA version** — `cupy-cuda12x` vs `cupy-cuda11x`. Wrong one = import error.
- **First kernel call is slow** (JIT/compile). Time the *second* run.
- **Async by default** — CuPy/CUDA return before work finishes. Call `cp.cuda.Stream.null.synchronize()` before timing.
- **Out-of-memory** — GPU RAM is small (e.g. 24 GB). Watch it in `nvidia-smi`; use `float32` over `float64`.
