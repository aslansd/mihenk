"""mihenk 0.3 checks on synthetic systems with known answers.

    python experiments/physics_synthetic/scripts/sweep.py 0 10     # calibration seeds 0-9
    python experiments/physics_synthetic/scripts/sweep.py 10 20    # test seeds 10-19 (see ../PROTOCOL.md)

Each row: one case at one noise level over the seeds. The count is for the
case's target check (error / warning / clean), plus a number that check
reports. Noise is a fraction of each variable's standard deviation.
"""

from __future__ import annotations

import sys
import warnings

import numpy as np

from mihenk import Model, scan, synth
from mihenk.adapters import pysindy as A

H_PEND = "0.5*omega**2 - 64*cos(theta)"           # energy of the synthetic pendulum (theta = 0 down)
H_LV = "0.5*x - log(x) + 0.5*y - log(y)"           # conserved by x' = x - 0.5xy, y' = -y + 0.5xy


def hand(text, b):
    return Model.from_text(text, variables=b.variables)


def sindy(b, degree, threshold, fourier=0):
    return A.fit_bundle(b, degree=degree, threshold=threshold, fourier=fourier)[1]


# name -> (system, sensor, model builder, scan kwargs, target check, metric)
CASES = {
    # PHYSICS_VIOLATION, dissipates
    "P1 true pendulum law": ("pendulum", None, lambda b: hand("theta' = omega\nomega' = -64 sin(theta) - 0.07 omega", b),
                             {"properties": {"dissipates": H_PEND}}, "PHYSICS_VIOLATION", "ratio"),
    "P2 SINDy, clean pendulum": ("pendulum", None, lambda b: sindy(b, 1, 0.02, 1),
                                 {"properties": {"dissipates": H_PEND}}, "PHYSICS_VIOLATION", "ratio"),
    "P3 energy-injecting law": ("pendulum", None, lambda b: hand("theta' = omega\nomega' = -64 sin(theta) + 0.03 omega", b),
                                {"properties": {"dissipates": H_PEND}}, "PHYSICS_VIOLATION", "ratio"),
    "P4 SINDy, omega 0.7 ms late": ("pendulum", {"omega": {"delay": 7e-4}}, lambda b: sindy(b, 1, 0.02, 1),
                                    {"properties": {"dissipates": H_PEND}, "kinematics": {}}, "PHYSICS_VIOLATION", "ratio"),
    "P5 over-damped law (x3)": ("pendulum", None, lambda b: hand("theta' = omega\nomega' = -64 sin(theta) - 0.21 omega", b),
                                {"properties": {"dissipates": H_PEND}}, "PHYSICS_VIOLATION", "ratio"),
    # PHYSICS_VIOLATION, conserves
    "C1 true Lotka-Volterra law": ("lotka_volterra", None, lambda b: hand("x' = x - 0.5 x y\ny' = -y + 0.5 x y", b),
                                   {"properties": {"conserves": H_LV}}, "PHYSICS_VIOLATION", "drift"),
    "C2 SINDy, clean Lotka-Volterra": ("lotka_volterra", None, lambda b: sindy(b, 2, 0.05),
                                       {"properties": {"conserves": H_LV}}, "PHYSICS_VIOLATION", "drift"),
    "C3 law with 5% damping on x": ("lotka_volterra", None, lambda b: hand("x' = 0.95 x - 0.5 x y\ny' = -y + 0.5 x y", b),
                                    {"properties": {"conserves": H_LV}}, "PHYSICS_VIOLATION", "drift"),
    # ROLLOUT_TRAP, long rollouts
    "L1 undamped pendulum law": ("pendulum", None, lambda b: hand("theta' = omega\nomega' = -64 sin(theta)", b),
                                 {}, "ROLLOUT_TRAP", "long_r2"),
    "L2 SINDy, clean pendulum": ("pendulum", None, lambda b: sindy(b, 1, 0.02, 1), {}, "ROLLOUT_TRAP", "long_r2"),
    "L3 SINDy, clean Van der Pol": ("van_der_pol", None, lambda b: sindy(b, 3, 0.05), {}, "ROLLOUT_TRAP", "long_r2"),
    "L4 SINDy, clean Duffing": ("duffing", None, lambda b: sindy(b, 3, 0.05), {}, "ROLLOUT_TRAP", "long_r2"),
    "L5 SINDy, Lorenz (chaotic)": ("lorenz", None, lambda b: sindy(b, 2, 0.5), {}, "ROLLOUT_TRAP", "chaotic"),
    # UNDERFIT two degrees out
    "U1 cubic library, large swings": ("pendulum", None, lambda b: sindy(b, 3, 0.02), {"kinematics": {}},
                                       "UNDERFIT", "best_term"),
    # centred relations
    "R1 small swings about pi": ("pendulum_small", {"theta": {"offset": float(np.pi)}}, lambda b: sindy(b, 1, 0.02, 1),
                                 {"kinematics": {}}, "RELATION", "relation"),
}


def run_case(name, seed, noise):
    system, sensor, build, kw, target, metric = CASES[name]
    b = synth.simulate(system, seed=seed, noise=noise, sensor=sensor)
    b.meta.pop("true_equations", None)
    r = scan(b, build(b), **kw)
    if target == "RELATION":
        rows = r.details.get("Relations among candidate terms on your data", {}).get("rows", [])
        hit = any("sin(theta)" in row[0] and "theta" in row[0].replace("sin(theta)", "") for row in rows)
        return ("found" if hit else "missing"), (rows[0][0] if rows else "-")
    f = r.by_id(target)
    sev = f[0].severity if f else "clean"
    if metric == "ratio":
        p = r.details["Declared physical properties"]["rows"][0]
        val = p[3] / p[2] if p[2] else float("nan")           # equation rate / data rate (1 = matched)
    elif metric == "drift":
        val = r.details["Declared physical properties"]["rows"][0][6]
    elif metric in ("long_r2", "chaotic"):
        full = [row for row in r.details["Rollouts (integrating the equation)"]["rows"] if row[1] == "full"]
        val = full[0][4] if metric == "long_r2" else ("chaotic" in full[0][5])
    else:
        rows = r.details["Best left-out term per equation"]["rows"]
        val = max(rows, key=lambda row: row[3])[1] if rows else "-"
    if metric == "chaotic":
        full = [row for row in r.details["Rollouts (integrating the equation)"]["rows"] if row[1] == "full"]
        val = "chaotic" if "chaotic" in full[0][5] and "not" not in full[0][5] else "not chaotic"
        if any("drifts away" in x.what for x in f):
            val += " + drift warning"
    return sev, val


def main(argv):
    warnings.simplefilter("ignore")
    lo, hi = int(argv[0]), int(argv[1])
    for noise in (0.0, 0.005, 0.01):
        for name in CASES:
            res = [run_case(name, seed, noise) for seed in range(lo, hi)]
            sev = [s for s, _ in res]
            vals = [v for _, v in res]
            if isinstance(vals[0], float):
                summ = f"median {np.nanmedian(vals):.3g} (range {np.nanmin(vals):.3g} to {np.nanmax(vals):.3g})"
            else:
                from collections import Counter
                summ = ", ".join(f"{k}: {n}" for k, n in Counter(map(str, vals)).most_common(2))
            counts = ("found {:>2}  missing {:>2}".format(sev.count("found"), sev.count("missing"))
                      if CASES[name][4] == "RELATION" else
                      "error {:>2}  warning {:>2}  clean {:>2}".format(sev.count("error"), sev.count("warning"),
                                                                     sev.count("clean")))
            print(f"noise {noise:<6} {name:<32} {counts} | {summ}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["0", "10"]))
