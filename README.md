# Computational physics portfolio

Numerical and modelling projects in Python (M.Sc. Physics, St. Xavier's College, Mumbai). Each folder has a script, a README with the method and a key figure, and runs with `python <script>.py`.

| Project | What it shows |
|---|---|
| [01 Compound pendulum](01_compound_pendulum) | Lagrangian model, nonlinear ODE vs small-angle, exact elliptic-integral check, hook for Tracker video data |
| [02 FODO lattice](02_fodo_lattice) | Transfer-matrix model of a ring, stability analysis, turn-by-turn tracking |
| [03 Lane-Emden](03_lane_emden) | RK4 for a singular ODE, validated against exact solutions and published tables |
| [04 Cyclotron](04_cyclotron) | Relativistic RK4 integration of a driven charged particle, energy-balance check |
| [05 Spatial filtering](05_spatial_filtering) | Noise reduction (mean / median / Gaussian) with a quantitative SNR comparison |
| [06 Bumblebee black hole](06_bumblebee_black_hole) | Symbolic derivation (sympy) of black-hole thermodynamics, with checks of published formulas |

Setup: `pip install -r requirements.txt`
