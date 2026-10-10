"""Lane-Emden equation solved with RK4 and checked against Chandrasekhar's table.

    (1/xi^2) d/dxi ( xi^2 dtheta/dxi ) = -theta^n

Written as a first-order system y = (theta, phi = dtheta/dxi):
    theta' = phi,   phi' = -theta^n - 2 phi / xi
The singular point xi = 0 is avoided by starting at small xi0 from the series
    theta = 1 - xi^2/6 + n xi^4/120 - ...
The first zero xi_1 (stellar surface) is located by linear interpolation.
"""
import matplotlib.pyplot as plt
import numpy as np

# xi_1 for the first zero of theta_n, Chandrasekhar (1939), Table 4
CHANDRASEKHAR_XI1 = {0: 2.4494, 1: 3.14159, 1.5: 3.65375, 2: 4.35287, 3: 6.89685, 4: 14.97155}


def rhs(xi, y, n):
    theta, phi = y
    return np.array([phi, -max(theta, 0.0) ** n - 2.0 * phi / xi])


def rk4_step(xi, y, h, n):
    k1 = rhs(xi, y, n)
    k2 = rhs(xi + h / 2, y + h * k1 / 2, n)
    k3 = rhs(xi + h / 2, y + h * k2 / 2, n)
    k4 = rhs(xi + h, y + h * k3, n)
    return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6


def solve(n, h=1e-3, xi0=1e-4, xi_max=30.0):
    """Return (xi, theta, xi_1) with xi_1 the first zero of theta."""
    y = np.array([1 - xi0**2 / 6 + n * xi0**4 / 120, -xi0 / 3 + n * xi0**3 / 30])
    xi, xs, ths = xi0, [xi0], [y[0]]
    while xi < xi_max:
        y_new = rk4_step(xi, y, h, n)
        if y_new[0] <= 0:  # linear interpolation to the zero crossing
            xi1 = xi + h * y[0] / (y[0] - y_new[0])
            xs.append(xi1)
            ths.append(0.0)
            return np.array(xs), np.array(ths), xi1
        xi, y = xi + h, y_new
        xs.append(xi)
        ths.append(y[0])
    return np.array(xs), np.array(ths), np.nan  # no surface (n >= 5)


def main():
    print(f"{'n':>4} {'xi_1 (RK4)':>12} {'xi_1 (table)':>13} {'rel. error':>11}")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.2))
    for n in [0, 1, 1.5, 2, 3, 4, 5]:
        xi, theta, xi1 = solve(n)
        ax1.plot(xi, theta, label=f"n = {n}")
        if n in CHANDRASEKHAR_XI1:
            ref = CHANDRASEKHAR_XI1[n]
            print(f"{n:>4} {xi1:12.5f} {ref:13.5f} {abs(xi1 - ref) / ref:11.2e}")
    ax1.set(xlabel=r"$\xi$", ylabel=r"$\theta_n(\xi)$", xlim=(0, 10), ylim=(0, 1.05),
            title="Lane-Emden solutions (RK4)")
    ax1.legend()

    # Exact solutions exist for n = 0, 1, 5: check the whole curve, not just the zero
    xi, th1, _ = solve(1)
    ax2.semilogy(xi, np.abs(th1 - np.sin(xi) / xi) + 1e-18, label="n = 1 vs sin(xi)/xi")
    xi, th0, _ = solve(0)
    ax2.semilogy(xi, np.abs(th0 - (1 - xi**2 / 6)) + 1e-18, label="n = 0 vs 1 - xi^2/6")
    xi, th5, _ = solve(5, xi_max=10)
    ax2.semilogy(xi, np.abs(th5 - (1 + xi**2 / 3) ** -0.5) + 1e-18, label="n = 5 vs (1+xi^2/3)^(-1/2)")
    ax2.set(xlabel=r"$\xi$", ylabel="absolute error", title="Error against exact solutions")
    ax2.legend()
    fig.tight_layout()
    fig.savefig("lane_emden.png", dpi=150)


if __name__ == "__main__":
    main()
