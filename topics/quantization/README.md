# Quantization

Storing/computing numbers with fewer bits (e.g. FP32 → INT8) to shrink models and speed up inference. Mostly used on neural-net **weights** (and sometimes activations).

## Why

- **Memory**: FP16 → INT4 is a 4× smaller model. A 70B model goes from ~140 GB to ~35 GB — the difference between "needs a cluster" and "fits on one GPU".
- **Speed**: integer math + less memory traffic. Inference is usually *memory-bandwidth bound*, so smaller weights = faster.
- **Energy**: fewer bits moved = less power.

## The core math (affine / uniform quantization)

Map a real value `x` in range `[β, α]` to a `b`-bit integer. This is **affine** (a.k.a. uniform — equal step sizes):

```
scale   s = (α − β) / (2^b − 1)          # size of one quantization step
zero    z = round(−β / s)                # integer that maps to real 0

quantize:    x_q = clip(round(x / s) + z,  0, 2^b − 1)
dequantize:  x̂  = s · (x_q − z)
```

- `s` (scale) and `z` (zero-point) are the **quantization parameters** — you store them alongside the ints.
- The error `x − x̂` is the **quantization error**, bounded by `±s/2` from rounding (plus clipping error for out-of-range values).
- **Symmetric** (`z = 0`, range `[−α, α]`): simpler/faster, wastes a code, good for weights (roughly zero-centered).
- **Asymmetric** (`z ≠ 0`): uses the full range, better for skewed data like post-ReLU activations (all ≥ 0).

## Key design axes

| Axis | Options | Notes |
|---|---|---|
| **Uniformity** | Uniform (affine) vs non-uniform | Non-uniform (e.g. NF4, log, k-means codebooks) puts more levels where data is dense |
| **Symmetry** | Symmetric vs asymmetric | Weights → symmetric; activations → asymmetric |
| **Granularity** | Per-tensor → per-channel → per-group | Finer = one `(s,z)` per fewer weights = less error, more overhead |
| **When** | PTQ vs QAT | Post-training (fast, no data/little data) vs quant-aware training (retrain, best quality) |
| **What** | Weight-only vs weight+activation | W-only (W4A16) is easiest & most common for LLMs; W8A8 needs activation calibration |

**Granularity matters most.** A single outlier weight blows up the per-tensor range and coarsens every step. Per-channel (one scale per output channel) or **per-group** (one scale per block of e.g. 32/64/128 weights) is why 4-bit LLMs work at all.

## Calibration — picking the range `[β, α]`

The clip range trades off two errors: too wide → coarse steps (rounding error); too narrow → clipping error. Strategies:

- **Min/max**: `α = max(x)`, `β = min(x)`. Simple, but one outlier ruins it.
- **Percentile**: clip at e.g. 99.9th percentile — deliberately clip outliers to keep steps fine.
- **MSE / KL-divergence**: search for the range minimizing reconstruction error or output distribution shift. What most PTQ toolkits do.

## PTQ vs QAT

- **PTQ (Post-Training Quantization)** — quantize an already-trained model. Needs a small *calibration set* (or none for weight-only). Fast, cheap. Modern LLM methods (below) are PTQ.
- **QAT (Quant-Aware Training)** — simulate quantization during training so the model *learns* to be robust. Uses the **Straight-Through Estimator (STE)**: `round()` has zero gradient everywhere, so on the backward pass you pretend it's the identity (gradient = 1). Best accuracy, but you need the training pipeline + compute.

## The LLM quantization zoo (what people actually run)

| Method | Idea | Typical |
|---|---|---|
| **GPTQ** | Layer-wise; uses inverse-Hessian info to compensate rounding error weight-by-weight | W4, one-shot |
| **AWQ** | Activation-aware: scale up the ~1% "salient" weight channels (spotted via activation magnitude) before quantizing | W4, W-only |
| **SmoothQuant** | Migrate activation outliers into weights (per-channel scaling) so both quantize to INT8 cleanly | W8A8 |
| **bitsandbytes NF4** | 4-bit **NormalFloat** — non-uniform codebook whose levels are quantiles of a normal dist (matches weight stats). Core of **QLoRA** | W4 |
| **GGUF k-quants** | llama.cpp's block formats (`Q4_K`, `Q6_K`, …): per-block scales, mixed precision per layer | CPU/edge |
| **FP8 / MXFP4** | Low-bit *floating* point (E4M3/E5M2, micro-scaled block formats). Hardware-native on H100/B200 | Training + inference |

**FP vs INT at low bit-width**: floating-point formats (FP8, FP4) spend bits on an exponent, giving a wider dynamic range per code — better for tensors with outliers. Integer formats give uniform precision. Newer GPUs accelerate FP8/FP4 directly.

## "Lossless" quantization — two very different meanings

**1. Near-lossless (colloquial).** "No *measurable* quality drop" on benchmarks. This is what the community usually means by a "lossless quant":
- INT8 / W8A8 is near-lossless for most models.
- 6-bit and good 5-bit GGUF quants (`Q6_K`, `Q5_K_M`) are near-lossless in practice for many LLMs.
- ⚠️ It's still **lossy** — outputs differ bit-for-bit, and perplexity/eval gaps show up if you look hard enough. "Lossless" here means "good enough that you can't tell", not zero error.

**2. Truly lossless (information-theoretic).** Bit-**identical** outputs — the decompressed weights reconstruct the originals exactly. This is not quantization in the lossy sense; it's **entropy coding** (compression):
- BF16 weights don't use their 16 bits uniformly — the **exponent** field is heavily concentrated (weights cluster near 0). Entropy of the actual distribution < 16 bits.
- **DFloat11 / lossless float compression**: Huffman-code the exponent bits, store mantissa/sign raw. Gets BF16 down to ~11 bits/param with **exact** reconstruction (GPU kernel decompresses on the fly). ~30% smaller, zero accuracy change.
- General principle: `bits_needed ≥ entropy H(X)` (Shannon). You can losslessly compress *only* down to the data's entropy — below that, information is destroyed. So true-lossless gains are modest (~30%); the big 4× wins are inherently lossy.

**Bottom line:** if someone says "lossless 4-bit", they mean meaning #1 (near-lossless benchmarks). Genuine bit-exact compression (#2) can't reach 4 bits for BF16 weights — that would violate the entropy bound.

## Gotchas

- **Outliers dominate.** A handful of large-magnitude weights/activations set the range and wreck everything else. Every good LLM method (AWQ, SmoothQuant, GPTQ) is fundamentally an outlier story.
- **Activations are harder than weights** — dynamic, input-dependent, outlier-prone. Weight-only (W4A16) sidesteps this and is the default for LLMs.
- **Measure the right thing.** Low weight-MSE ≠ good model. Track end-task metrics (perplexity, eval accuracy), not just reconstruction error.
- **Overhead is real.** Per-group scales add bytes back. `Q4_K` isn't 4.0 bits/weight — it's ~4.5 once you count scales.
- **Dequant cost.** W4A16 stores 4-bit but computes in FP16 — you dequantize on the fly. The win is bandwidth (loading fewer bytes), not the matmul itself.

## Samples in this folder

- [`affine_quant.py`](affine_quant.py) — symmetric & asymmetric INT8 from scratch; per-tensor vs per-channel error on a matrix with an outlier.
- [`lossless_entropy.py`](lossless_entropy.py) — why BF16 weights compress losslessly: measure exponent entropy, show the ~11-bit bound.
