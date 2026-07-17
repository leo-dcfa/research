"""Why 'truly lossless' quantization of BF16 weights lands near ~11 bits.

Lossless compression can't beat the data's Shannon entropy. BF16 weights waste
bits because the *exponent* field is heavily concentrated (weights cluster near
0), so entropy-coding the exponents (a la DFloat11) shrinks the model with
BIT-EXACT reconstruction -- but only by ~30%, never down to 4 bits.

Run:  uv run python topics/quantization/lossless_entropy.py
"""

import numpy as np


def bf16_fields(x):
    """Split float32 values into BF16 (sign, 8-bit exponent, 7-bit mantissa).

    BF16 = top 16 bits of a float32, so we can read the fields off the float32
    bit pattern directly.
    """
    bits = x.astype(np.float32).view(np.uint32)
    bf16 = (bits >> 16).astype(np.uint16)  # truncate to BF16
    sign = (bf16 >> 15) & 0x1
    exponent = (bf16 >> 7) & 0xFF  # 8 bits
    mantissa = bf16 & 0x7F  # 7 bits
    return sign, exponent, mantissa


def entropy_bits(symbols):
    """Shannon entropy H(X) in bits -- the lossless coding lower bound."""
    _, counts = np.unique(symbols, return_counts=True)
    p = counts / counts.sum()
    return float(-np.sum(p * np.log2(p)))


def main():
    rng = np.random.default_rng(0)

    # realistic LLM weights: ~N(0, 0.02), heavy mass near zero
    weights = rng.normal(scale=0.02, size=2_000_000).astype(np.float32)

    sign, exponent, mantissa = bf16_fields(weights)

    h_sign = entropy_bits(sign)
    h_exp = entropy_bits(exponent)
    h_man = entropy_bits(mantissa)

    print("Per-field Shannon entropy of BF16 weights (bits):")
    print(f"  sign     : {h_sign:5.2f}  (of 1 stored bit)")
    print(f"  exponent : {h_exp:5.2f}  (of 8 stored bits)  <- concentrated!")
    print(
        f"  mantissa : {h_man:5.2f}  (of 7 stored bits)  <- ~uniform, near-incompressible"
    )

    # A DFloat11-style scheme: Huffman-code the exponent, store sign+mantissa raw.
    lossless_bits = h_exp + 1 + 7  # entropy-coded exp + raw sign + raw mantissa
    print("\nStored per weight   : 16.00 bits (BF16)")
    print(f"Lossless floor      : {h_sign + h_exp + h_man:5.2f} bits (full entropy)")
    print(f"DFloat11-style      : {lossless_bits:5.2f} bits (code exponent only)")
    print(f"  -> {16 / lossless_bits:.2f}x smaller, BIT-EXACT reconstruction.")

    print(
        "\nKey point: the mantissa is ~uniform, so you CANNOT losslessly reach 4\n"
        "bits -- that would drop below the entropy bound and destroy information.\n"
        "4-bit quantization is fundamentally LOSSY; ~11-bit lossless is the real floor."
    )


if __name__ == "__main__":
    main()
