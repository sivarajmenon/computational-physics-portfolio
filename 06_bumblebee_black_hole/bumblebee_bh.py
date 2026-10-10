"""Thermodynamics of a black hole with cosmological constant in bumblebee gravity.

Horizon condition f(r)=0 with f = -(1 - 2M/r - (1+l) Lambda r^2 / 3) gives M(r). With
S = pi r^2 (area law) everything follows from M(S):

    T = dM/dS,   F = M - T S,   C = T / (dT/dS)

l is the Lorentz-violating parameter, Lambda the effective cosmological constant.
The symbolic results are derived with sympy and checked against the thesis equations.
"""
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

S = sp.Symbol("S", positive=True)
l, Lam = sp.symbols("l Lambda", real=True)
a = (1 + l) * Lam  # only this combination appears

M = sp.sqrt(S) / (2 * sp.sqrt(sp.pi)) - a * S ** sp.Rational(3, 2) / (6 * sp.pi ** sp.Rational(3, 2))
T = sp.simplify(sp.diff(M, S))
F = sp.simplify(M - T * S)
C = sp.simplify(T / sp.diff(T, S))


def check_against_thesis():
    T_thesis = 1 / (4 * sp.sqrt(sp.pi * S)) - a * sp.sqrt(S) / (4 * sp.pi ** sp.Rational(3, 2))
    F_thesis = (3 * sp.pi * sp.sqrt(S) - a * S ** sp.Rational(3, 2)) / (12 * sp.pi ** sp.Rational(3, 2))
    print("T  =", T)
    print("F  =", F)
    print("C  =", sp.factor(C))
    print("T matches thesis Eq. 5.5:", sp.simplify(T - T_thesis) == 0)
    print("F matches thesis Eq. 5.6:", sp.simplify(F - F_thesis) == 0,
          "(difference:", sp.simplify(F - F_thesis), ")")


def curves(lam, lorentz, s_max, s_min=0.01, n=800):
    s = np.linspace(s_min, s_max, n)
    sub = {Lam: lam, l: lorentz}
    f = {name: sp.lambdify(S, expr.subs(sub), "numpy") for name, expr in
         dict(M=M, T=T, F=F, C=C).items()}
    return s, {k: np.broadcast_to(fn(s), s.shape).astype(float) for k, fn in f.items()}


def main():
    check_against_thesis()
    lorentz = 1.0
    cases = {"Anti-de Sitter (Lambda=-1, l=1)": (-1.0, 12.0), "de Sitter (Lambda=+1, l=1)": (1.0, 6.0)}
    fig, axes = plt.subplots(2, 4, figsize=(16, 7))
    for row, (title, (lam, s_max)) in zip(axes, cases.items()):
        s, y = curves(lam, lorentz, s_max)
        a_val = (1 + lorentz) * lam
        s_pole = np.pi / abs(a_val)  # T=0 (dS, Nariai limit) or C pole (AdS)
        for ax, (name, label) in zip(row, [("M", "mass"), ("T", "temperature"),
                                           ("F", "free energy"), ("C", "heat capacity")]):
            values = y[name].copy()
            if name == "C":
                values[np.abs(values) > 40] = np.nan  # hide the pole
                ax.set_ylim(-40, 40)
            ax.plot(s, values)
            ax.axvline(s_pole, color="r", ls="--", lw=1)
            ax.axhline(0, color="gray", lw=0.5)
            ax.set(xlabel="entropy S", title=label)
        row[0].set_ylabel(title)
        if lam > 0:
            print(f"dS: M>0 and T>0 only for S < pi/((1+l)Lambda) = {s_pole:.3f} (Nariai limit)")
        else:
            print(f"AdS: C diverges at S = pi/|(1+l)Lambda| = {s_pole:.3f}; C<0 below (unstable), C>0 above (stable)")
    fig.tight_layout()
    fig.savefig("bumblebee_bh.png", dpi=140)


if __name__ == "__main__":
    main()
