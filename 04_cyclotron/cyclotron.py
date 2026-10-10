"""Relativistic proton in an ideal cyclotron, integrated with RK4.

State: position r and momentum p. Equation of motion:
    dp/dt = q (E + v x B),   v = p / (gamma m)
B is uniform along z. E (along x) acts only in the gap between the dees and its
sign follows the particle's x-velocity (ideal, perfectly phase-locked square wave).
Each gap crossing then adds exactly q*E0*gap of kinetic energy, which gives an
analytic check on the integrator.
"""
import matplotlib.pyplot as plt
import numpy as np

C = 299_792_458.0
Q, M = 1.602176634e-19, 1.67262192e-27
B0 = 1.5          # T
E0 = 5.0e6        # V/m
N_CROSSINGS = 46
V0 = 0.05 * C     # injection speed


def gamma_of(p):
    return np.sqrt(1.0 + (np.linalg.norm(p) / (M * C)) ** 2)


def derivs(state, gap):
    r, p = state[:3], state[3:]
    v = p / (gamma_of(p) * M)
    e = np.zeros(3)
    if -gap / 2 < r[0] < gap / 2:
        e[0] = E0 * np.sign(v[0])
    force = Q * (e + np.cross(v, [0.0, 0.0, B0]))
    return np.concatenate([v, force])


def rk4(state, dt, gap):
    k1 = derivs(state, gap)
    k2 = derivs(state + dt * k1 / 2, gap)
    k3 = derivs(state + dt * k2 / 2, gap)
    k4 = derivs(state + dt * k3, gap)
    return state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6


def kinetic_energy(p):
    return (gamma_of(p) - 1.0) * M * C**2


def main():
    r0 = M * V0 / (Q * B0)          # initial gyroradius (non-relativistic)
    gap = 0.5 * r0
    period = 2 * np.pi * M / (Q * B0)
    dt = period / 2000
    state = np.array([-gap / 2, 0, 0, M * V0, 0, 0], dtype=float)

    traj, times, speeds, radii = [state[:3].copy()], [0.0], [V0], [r0]
    crossings, inside, t = 0, True, 0.0
    while crossings < N_CROSSINGS:
        state = rk4(state, dt, gap)
        t += dt
        now_inside = abs(state[0]) < gap / 2
        if inside and not now_inside:
            crossings += 1
            p = state[3:]
            speeds.append(np.linalg.norm(p) / (gamma_of(p) * M))
            radii.append(np.linalg.norm(p) / (Q * B0))
            times.append(t)
        inside = now_inside
        traj.append(state[:3].copy())

    traj = np.array(traj)
    ke0 = 0.5 * M * V0**2
    ke_final = kinetic_energy(state[3:])
    ke_expected = ke0 + N_CROSSINGS * Q * E0 * gap
    v_final = speeds[-1]
    print(f"initial radius {r0:.4f} m, gap {gap:.4f} m, period {period * 1e9:.2f} ns")
    print(f"final speed  {v_final:.4e} m/s = {v_final / C * 100:.2f}% of c, gamma = {gamma_of(state[3:]):.4f}")
    print(f"final radius {radii[-1]:.3f} m, elapsed time {t * 1e6:.3f} us")
    print(f"final KE {ke_final / Q / 1e6:.3f} MeV vs expected {ke_expected / Q / 1e6:.3f} MeV "
          f"(relative difference {abs(ke_final - ke_expected) / ke_expected:.1e})")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    axes[0].plot(traj[:, 0], traj[:, 1], lw=0.6)
    axes[0].set(aspect="equal", xlabel="x [m]", ylabel="y [m]", title="Spiral trajectory")
    axes[1].step(np.array(times) * 1e6, np.array(speeds) / C, where="post")
    axes[1].set(xlabel="time [us]", ylabel="v / c", title="Speed after each gap crossing")
    axes[2].plot(radii, np.array(speeds) / C, "o-", ms=3)
    axes[2].set(xlabel="radius [m]", ylabel="v / c", title="Speed vs radius")
    fig.tight_layout()
    fig.savefig("cyclotron.png", dpi=150)


if __name__ == "__main__":
    main()
