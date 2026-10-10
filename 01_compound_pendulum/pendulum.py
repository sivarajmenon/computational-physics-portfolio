"""Compound (physical) pendulum: nonlinear vs small-angle solution.

Equation of motion (from the Lagrangian L = 1/2 I thetadot^2 + m g d cos(theta)):
    I * theta'' = -m g d sin(theta)

Usage:
    python pendulum.py                       # uniform rod, pivoted at one end
    python pendulum.py --tracker data.csv    # overlay Tracker data (columns: t [s], theta [deg])
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipk

G = 9.81  # m/s^2


def rod_about_end(mass, length):
    """Return (d, I) for a uniform rod pivoted at one end."""
    return length / 2.0, mass * length**2 / 3.0


def solve(theta0, mass, d, inertia, t_end, linear=False, n_points=4000):
    """Integrate the pendulum from rest at angle theta0 (radians)."""
    restoring = (lambda th: th) if linear else np.sin

    def rhs(_, y):
        return [y[1], -mass * G * d / inertia * restoring(y[0])]

    t = np.linspace(0.0, t_end, n_points)
    sol = solve_ivp(rhs, (0, t_end), [theta0, 0.0], t_eval=t, rtol=1e-10, atol=1e-12)
    return sol.t, sol.y[0]


def small_angle_period(mass, d, inertia):
    return 2 * np.pi * np.sqrt(inertia / (mass * G * d))


def exact_period(theta0, mass, d, inertia):
    """Exact period for amplitude theta0, via the complete elliptic integral K."""
    return small_angle_period(mass, d, inertia) * 2 / np.pi * ellipk(np.sin(theta0 / 2) ** 2)


def measured_period(t, theta):
    """Period from successive upward zero crossings (linear interpolation)."""
    idx = np.where((theta[:-1] < 0) & (theta[1:] >= 0))[0]
    crossings = t[idx] - theta[idx] * (t[idx + 1] - t[idx]) / (theta[idx + 1] - theta[idx])
    return np.mean(np.diff(crossings))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mass", type=float, default=1.175, help="kg")
    ap.add_argument("--length", type=float, default=1.0, help="rod length, m")
    ap.add_argument("--amplitude", type=float, default=5.0, help="initial angle, deg")
    ap.add_argument("--t-end", type=float, default=20.0)
    ap.add_argument("--tracker", help="CSV with columns t [s], theta [deg]")
    ap.add_argument("--out", default="pendulum.png")
    a = ap.parse_args()

    d, inertia = rod_about_end(a.mass, a.length)
    th0 = np.radians(a.amplitude)
    t, th_nl = solve(th0, a.mass, d, inertia, a.t_end)
    _, th_lin = solve(th0, a.mass, d, inertia, a.t_end, linear=True)

    t_small = small_angle_period(a.mass, d, inertia)
    t_exact = exact_period(th0, a.mass, d, inertia)
    t_num = measured_period(t, th_nl)
    print(f"small-angle period : {t_small:.5f} s")
    print(f"exact (elliptic)   : {t_exact:.5f} s")
    print(f"numerical (nonlin) : {t_num:.5f} s   (error vs exact: {abs(t_num - t_exact) / t_exact:.2e})")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.2))
    ax1.plot(t, np.degrees(th_nl), label="nonlinear")
    ax1.plot(t, np.degrees(th_lin), "--", label="small-angle")
    if a.tracker:
        data = np.loadtxt(a.tracker, delimiter=",", skiprows=1)
        ax1.plot(data[:, 0], data[:, 1], "k.", ms=3, label="Tracker")
        model = np.degrees(np.interp(data[:, 0], t, th_nl))
        print(f"RMS error vs Tracker: {np.sqrt(np.mean((model - data[:, 1]) ** 2)):.3f} deg")
    ax1.set(xlabel="time [s]", ylabel="theta [deg]", title=f"Amplitude {a.amplitude:g} deg")
    ax1.legend()

    amps = np.radians(np.linspace(1, 120, 60))
    ax2.plot(np.degrees(amps), [exact_period(x, a.mass, d, inertia) for x in amps], label="exact (elliptic K)")
    ax2.axhline(t_small, color="gray", ls="--", label="small-angle")
    num = [measured_period(*solve(x, a.mass, d, inertia, 40.0)) for x in amps[::6]]
    ax2.plot(np.degrees(amps[::6]), num, "o", label="numerical")
    ax2.set(xlabel="amplitude [deg]", ylabel="period [s]", title="Period vs amplitude")
    ax2.legend()
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)


if __name__ == "__main__":
    main()
