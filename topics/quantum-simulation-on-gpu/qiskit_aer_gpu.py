"""Quantum circuit simulation on an NVIDIA GPU with Qiskit Aer.

Aer's GPU backend uses cuQuantum under the hood. Flipping to the GPU is a single
argument: AerSimulator(device="GPU"). Everything else is ordinary Qiskit.

Run on a CUDA machine or Google Colab (GPU runtime):
    pip install qiskit qiskit-aer-gpu
    python qiskit_aer_gpu.py

No GPU? Use device="CPU" (the default) to run smaller circuits locally.
"""

from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

DEVICE = "GPU"  # -> "CPU" to fall back to CPU simulation


def ghz(n_qubits: int) -> QuantumCircuit:
    """Build a GHZ circuit: H on qubit 0, then a ladder of CNOTs."""
    qc = QuantumCircuit(n_qubits, n_qubits)
    qc.h(0)
    for i in range(n_qubits - 1):
        qc.cx(i, i + 1)
    qc.measure(range(n_qubits), range(n_qubits))
    return qc


def main() -> None:
    n_qubits = 30  # ~8.6 GB of state vector at complex64 — needs a real GPU
    sim = AerSimulator(device=DEVICE, precision="single")  # single = complex64

    qc = ghz(n_qubits)
    result = sim.run(qc, shots=1024).result()
    counts = result.get_counts()

    # A perfect GHZ state only ever collapses to all-zeros or all-ones.
    print(f"GHZ({n_qubits}) on {DEVICE}: {len(counts)} distinct outcomes")
    for bitstring, n in sorted(counts.items(), key=lambda kv: -kv[1])[:2]:
        print(f"  {bitstring[:8]}... : {n} shots")


if __name__ == "__main__":
    main()
