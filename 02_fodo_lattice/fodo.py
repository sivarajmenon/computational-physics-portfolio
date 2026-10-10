"""Stability of a FODO lattice in a ring, using 6x6 transfer matrices.

Phase-space vector: (x, x', y, y', z, delta). One cell = drift -> QD -> dipole -> QF
(matrices multiply right-to-left). A ring of N cells is stable if |Tr(M_cell)/2| < 1
in both planes. The "most stable" gradient k maximises the stability margin, i.e.
minimises max(|Tr_x|, |Tr_y|)/2.
"""
import matplotlib.pyplot as plt
import numpy as np

L_QUAD, L_DIPOLE, L_DRIFT = 0.5, 1.0, 1.0  # m
N_CELLS = 6
R_BEND = N_CELLS * L_DIPOLE / (2 * np.pi)  # bending radius so the dipoles close the ring


def drift(length):
    m = np.eye(6)
    m[0, 1] = m[2, 3] = length
    return m


def dipole(length, radius):
    a = length / radius
    m = np.eye(6)
    m[0, :2] = [np.cos(a), radius * np.sin(a)]
    m[0, 4] = radius * (1 - np.cos(a))
    m[1, :2] = [-np.sin(a) / radius, np.cos(a)]
    m[1, 4] = np.sin(a)
    m[2, 3] = length
    return m


def quad(length, k):
    """Quadrupole of normalised gradient k [1/m^2]; k>0 focuses in x, k<0 in y."""
    ks = np.sqrt(abs(k))
    foc = np.array([[np.cos(ks * length), np.sin(ks * length) / ks],
                    [-ks * np.sin(ks * length), np.cos(ks * length)]])
    dfo = np.array([[np.cosh(ks * length), np.sinh(ks * length) / ks],
                    [ks * np.sinh(ks * length), np.cosh(ks * length)]])
    m = np.eye(6)
    m[0:2, 0:2], m[2:4, 2:4] = (foc, dfo) if k > 0 else (dfo, foc)
    return m


def cell_matrix(k):
    return quad(L_QUAD, k) @ dipole(L_DIPOLE, R_BEND) @ quad(L_QUAD, -k) @ drift(L_DRIFT)


def stability(k):
    m = cell_matrix(k)
    return abs(np.trace(m[0:2, 0:2]) / 2), abs(np.trace(m[2:4, 2:4]) / 2)


def track(k, x0, turns):
    ring = np.linalg.matrix_power(cell_matrix(k), N_CELLS)
    out = np.empty((turns + 1, 6))
    out[0] = x0
    for i in range(turns):
        out[i + 1] = ring @ out[i]
    return out


def main():
    ks = np.arange(0.05, 3.0, 0.01)
    margin = np.array([max(stability(k)) for k in ks])
    stable = margin < 1
    k_best = ks[np.argmin(margin)]
    print(f"bending radius R = {R_BEND:.3f} m")
    print(f"stable for k in [{ks[stable].min():.2f}, {ks[stable].max():.2f}] 1/m^2")
    print(f"most stable gradient k = {k_best:.2f} 1/m^2  (max|Tr/2| = {margin.min():.3f})")

    x0 = np.array([0.1, 0.0, 0.5, 0.0, 0.0, 0.0])
    turns = 200
    good, bad = track(k_best, x0, turns), track(2.5, x0, turns)
    print(f"max |x|, |y| over {turns} turns at k={k_best:.2f}: {np.abs(good[:, [0, 2]]).max():.3f} m")
    print(f"max |x|, |y| over {turns} turns at k=2.50 (unstable): {np.abs(bad[:, [0, 2]]).max():.3g} m")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.2))
    ax1.plot(ks, [stability(k)[0] for k in ks], label="x plane")
    ax1.plot(ks, [stability(k)[1] for k in ks], label="y plane")
    ax1.axhline(1, color="r", ls="--", label="stability limit")
    ax1.axvline(k_best, color="k", ls=":", label=f"k = {k_best:.2f}")
    ax1.set(xlabel="quadrupole gradient k [1/m$^2$]", ylabel="|Tr(M)/2|", ylim=(0, 1.5),
            title="FODO stability diagram")
    ax1.legend()
    n = np.arange(turns + 1)
    ax2.plot(n, good[:, 0], label="x")
    ax2.plot(n, good[:, 2], label="y")
    ax2.set(xlabel="turn", ylabel="displacement [m]", title=f"Tracking at k = {k_best:.2f} (bounded)")
    ax2.legend()
    fig.tight_layout()
    fig.savefig("fodo.png", dpi=150)


if __name__ == "__main__":
    main()
