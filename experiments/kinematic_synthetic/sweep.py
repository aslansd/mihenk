"""KINEMATIC_MISMATCH on synthetic pendulums with known sensor imperfections.

    python experiments/kinematic_synthetic/sweep.py 0 10     # calibration seeds 0-9
    python experiments/kinematic_synthetic/sweep.py 10 20    # test seeds 10-19 (see PROTOCOL.md)

Each row: one imperfection at one noise level over the seeds. The pendulum is
mihenk's `pendulum` system (theta' = omega, omega' = -64 sin(theta) - 0.07 omega).
It is sampled at 2 ms and fitted with PySINDy (polynomial degree 1 + Fourier 1,
STLSQ threshold 0.02). The imperfection is applied to the recorded omega only.
"""
import warnings, sys; warnings.filterwarnings("ignore")
import numpy as np
from mihenk import synth, scan
from mihenk.adapters import pysindy as A
seeds = range(int(sys.argv[1]), int(sys.argv[2]))
conds = [("none", None), ("lag0.7ms", {"omega": {"delay": 7e-4}}), ("lag0.3ms", {"omega": {"delay": 3e-4}}),
         ("gain1.02", {"omega": {"gain": 1.02}}), ("offset0.05", {"omega": {"offset": 0.05}})]
for noise in (0.0, 0.002, 0.005, 0.01):
    for name, sensor in conds:
        res = []
        for seed in seeds:
            b = synth.simulate("pendulum", seed=seed, sensor=sensor, noise=noise)
            s, m = A.fit_bundle(b, degree=1, threshold=0.02, fourier=1)
            r = scan(b, m)
            f = [x for x in r.findings if x.id == "KINEMATIC_MISMATCH"]
            row = r.details["Kinematic pairs (derivative-free)"]["rows"][0]
            res.append((f[0].severity if f else "clean", row[1], row[5], row[6]))
        sev = [x[0] for x in res]
        print(f"noise {noise:<6} {name:<11} error {sev.count('error'):>2}  warning {sev.count('warning'):>2}  clean {sev.count('clean'):>2} | "
              f"lag est {np.median([x[1] for x in res]):+.3f} samp, share {np.median([x[2] for x in res]):.4f}, F median {np.median([x[3] for x in res]):.3g}")
