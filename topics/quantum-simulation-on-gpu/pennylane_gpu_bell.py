"""Quantum circuit simulation on an NVIDIA GPU with PennyLane + lightning.gpu.

lightning.gpu runs state-vector simulation on cuQuantum (cuStateVec). The only
change from CPU simulation is the device string: "lightning.gpu" instead of
"lightning.qubit" / "default.qubit".

Run on a CUDA machine or Google Colab (GPU runtime):
    pip install pennylane pennylane-lightning-gpu
    python pennylane_gpu_bell.py

No GPU? Set DEVICE = "lightning.qubit" to run the small cases on CPU.
"""

import pennylane as qml
from pennylane import numpy as np

DEVICE = "lightning.gpu"  # -> "lightning.qubit" to fall back to CPU


def bell_state() -> None:
    """The simplest entangled state: measures 00 and 11 with equal probability."""
    dev = qml.device(DEVICE, wires=2)

    @qml.qnode(dev)
    def circuit():
        qml.Hadamard(wires=0)
        qml.CNOT(wires=[0, 1])
        return qml.probs(wires=[0, 1])

    probs = circuit()
    print(f"Bell state probabilities [00, 01, 10, 11]: {np.round(probs, 3)}")


def ghz_scaling(n_qubits: int) -> None:
    """GHZ state on n qubits — the state vector has 2**n amplitudes.

    This is where the GPU earns its keep: push n_qubits up (28, 30, 32...) and
    the GPU stays fast while a CPU crawls. Watch memory: 2**n complex numbers.
    """
    dev = qml.device(DEVICE, wires=n_qubits)

    @qml.qnode(dev)
    def circuit():
        qml.Hadamard(wires=0)
        for i in range(n_qubits - 1):
            qml.CNOT(wires=[i, i + 1])
        # Expectation of Z on the last qubit (0 for a perfect GHZ state).
        return qml.expval(qml.PauliZ(n_qubits - 1))

    print(f"GHZ({n_qubits}) <Z_last> = {float(circuit()):.3f}")


if __name__ == "__main__":
    bell_state()
    ghz_scaling(24)  # bump toward 30+ on a real GPU to see it shine
