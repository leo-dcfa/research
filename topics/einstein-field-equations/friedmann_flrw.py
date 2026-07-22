"""Feed the FLRW metric to the same curvature machinery -> Friedmann equations.

The 10 Einstein field equations look terrifying, but if you *assume* the
universe is homogeneous and isotropic (the cosmological principle), the metric
has a single unknown function a(t) (the scale factor), and the 10 equations
collapse to just TWO -- the Friedmann equations that run modern cosmology.

FLRW metric (flat, k=0, geometrized units G = c = 1):
    ds^2 = -dt^2 + a(t)^2 (dx^2 + dy^2 + dz^2)

With a perfect-fluid source T^mu_nu = diag(-rho, p, p, p), the EFE
G_mu_nu = 8 pi T_mu_nu give:

    tt:  3 (a'/a)^2            = 8 pi rho          (Friedmann I)
    ii:  -2 a''/a - (a'/a)^2   = 8 pi p             (Friedmann II)

Rearranged, these are the two equations in the README. Here we let sympy
derive G_mu_nu and read them straight off.

Run:  uv run topics/einstein-field-equations/friedmann_flrw.py
"""

import sympy as sp

from curvature_sympy import einstein


def main():
    t, x, y, z = sp.symbols("t x y z")
    a = sp.Function("a")(t)  # scale factor a(t): the ONE unknown
    coords = [t, x, y, z]

    # flat FLRW metric
    g = sp.diag(-1, a**2, a**2, a**2)

    print("Flat FLRW metric g_mu_nu:")
    sp.pprint(g)

    G, _, _, _, _ = einstein(g, coords)

    print("\nEinstein tensor G_mu_nu (mixed lower indices):")
    sp.pprint(G)

    # Mixed-index G^mu_nu is nicer -- multiply by g^{-1}. Diagonal entries are
    # the two independent equations.
    G_mixed = sp.simplify(g.inv() * G)

    G_tt = sp.simplify(-G_mixed[0, 0])  # = 3 (a'/a)^2
    G_ii = sp.simplify(G_mixed[1, 1])  # = -(2 a''/a + (a'/a)^2)

    rho, p = sp.symbols("rho p")
    print("\n--- Friedmann I  (from the tt component,  G^t_t = -8 pi rho) ---")
    print("   3 (a'/a)^2 = 8 pi rho")
    print("   sympy G gives -G^t_t =", G_tt)
    sp.pprint(sp.Eq(G_tt, 8 * sp.pi * rho))

    print("\n--- Friedmann II (from the spatial component, G^i_i = 8 pi p) ---")
    print("   -2 a''/a - (a'/a)^2 = 8 pi p")
    print("   sympy G gives  G^i_i =", G_ii)
    sp.pprint(sp.Eq(G_ii, 8 * sp.pi * p))

    print(
        "\nTakeaway: 10 coupled PDEs -> 2 ODEs for a(t). Plug in matter/radiation/Lambda\n"
        "and you get the entire expansion history of the universe (the LCDM model)."
    )


if __name__ == "__main__":
    main()
