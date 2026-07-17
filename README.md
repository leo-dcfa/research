# research

A personal research repo. I ask Claude to research topics for me and drop the results here as
markdown notes, code samples, and the occasional visualisation.

## Layout

Each topic lives in its own folder under `topics/`:

```
topics/
  <topic-name>/
    README.md      # the notes — bullet points, short read
    *.py           # code samples (Python)
    *.py           # Manim scenes for video/animated examples
```

## Topics

- [Quantization](topics/quantization/README.md) — storing/computing numbers with fewer bits (FP32 → INT8/INT4) to shrink models and speed up inference.
- [Simulating quantum circuits on NVIDIA GPUs](topics/quantum-simulation-on-gpu/README.md) — classically simulating qubits/circuits on NVIDIA hardware with cuQuantum, CUDA-Q, and friends.

## Setup

This is a [uv](https://docs.astral.sh/uv/) project.

```bash
uv sync              # install dependencies
uv run <script.py>   # run a code sample
pre-commit install   # enable ruff on commit
```
