# Simulating quantum circuits on NVIDIA GPUs

Running **quantum computing simulations** (classically simulating qubits/circuits) on NVIDIA
hardware. You don't have a real quantum computer — you simulate one, and GPUs make the big
simulations fast.

## Why a GPU helps

- Simulating `n` qubits means tracking a **state vector of 2ⁿ complex amplitudes**, and every gate
  is a big matrix–vector multiply over that vector.
- That's dense, data-parallel linear algebra — exactly what GPUs are built for. Speedups of
  **10–100×** over CPU are typical for 20+ qubit circuits.

## The two simulation methods

| Method | Cost | Good for | Qubit reach |
|---|---|---|---|
| **State vector** | Exact, memory = 2ⁿ amplitudes | Any circuit | ~30–36 qubits (memory-bound) |
| **Tensor network** | Depends on entanglement | Shallow / low-entanglement circuits | 100s–1000s of qubits |

**The memory wall (state vector):** amplitudes are complex.
- `complex64` (8 bytes): **30 qubits ≈ 8.6 GB**, 33 ≈ 69 GB — each extra qubit **doubles** memory.
- This is why past ~33 qubits you need multi-GPU, or you switch to tensor networks.

## NVIDIA's stack

- **cuQuantum** — the low-level SDK that does the actual GPU work. Two engines:
  - **cuStateVec** → accelerates *state-vector* simulation.
  - **cuTensorNet** → accelerates *tensor-network* simulation.
  - You rarely call it directly — it plugs in under the frameworks below.
- **CUDA-Q** (`cudaq`) — NVIDIA's native quantum programming platform. Write kernels once, run on GPU
  simulators *or* real QPUs; built for hybrid quantum-classical and multi-GPU scale-out.

## Frameworks (what you actually write), all GPU-backed via cuQuantum

| Framework | GPU backend | Flip to GPU with |
|---|---|---|
| **PennyLane** | `lightning.gpu` (cuStateVec) | `qml.device("lightning.gpu", wires=n)` |
| **Qiskit** (Aer) | Aer-GPU (cuQuantum) | `AerSimulator(device="GPU")` |
| **Cirq** | qsim / qsimcirq | `qsimcirq.QSimSimulator(... gpu ...)` |
| **CUDA-Q** | native | `cudaq.set_target("nvidia")` |

**Rule of thumb:** already using PennyLane/Qiskit? Just switch the device string. Starting fresh and
want NVIDIA-native + real-QPU path later? Use **CUDA-Q**.

## Where to run it

- **Local NVIDIA GPU** — install CUDA, then the GPU plugin (`pip install pennylane-lightning-gpu`,
  `qiskit-aer-gpu`, or `cudaq`).
- **Google Colab** — free NVIDIA GPU; quickest way to try the samples *(Runtime → GPU)*.
- **Cloud** — AWS `g5`/`p4`, NVIDIA DGX Cloud, or the **NVIDIA Quantum Appliance** container (cuQuantum preinstalled).
- ⚠️ **Not on this Mac** — no NVIDIA GPU/CUDA. The `.py` samples run on a CUDA machine or Colab
  (the CPU device strings, e.g. `lightning.qubit`, work locally for tiny circuits).

## Samples in this folder

- [`pennylane_gpu_bell.py`](pennylane_gpu_bell.py) — Bell state + a scaling GHZ circuit on `lightning.gpu`.
- [`qiskit_aer_gpu.py`](qiskit_aer_gpu.py) — build a circuit, run 30 qubits on `AerSimulator(device="GPU")`.

CUDA-Q, for reference, looks like this:

```python
import cudaq
cudaq.set_target("nvidia")           # GPU state-vector backend

@cudaq.kernel
def bell():
    q = cudaq.qvector(2)
    h(q[0])
    x.ctrl(q[0], q[1])

print(cudaq.sample(bell))            # {'00': ~500, '11': ~500}
```

## Gotchas

- **GPU only wins on big circuits.** Below ~18–20 qubits the CPU is faster (launch overhead). Don't
  benchmark a 3-qubit Bell state and conclude the GPU is slow.
- **You hit RAM, not time.** The wall is memory (2ⁿ), so simulation dies with an OOM long before it
  gets "too slow". Check the memory table before picking a qubit count.
- **`complex64` vs `complex128`** halves memory — use single precision unless you need the accuracy.
- **Match the GPU plugin to your CUDA version**, same as any cuQuantum install.
