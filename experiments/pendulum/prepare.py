"""Build mihenk bundles from the MultiArm-Pendulum repository (single pendulum).

    git clone --depth 1 https://github.com/dynamicslab/MultiArm-Pendulum
    python experiments/pendulum/prepare.py path/to/MultiArm-Pendulum

Writes pendulum_large.npz and pendulum_small.npz next to this script. See
PROTOCOL.md for what each one is. Needs numpy and scipy only.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.io as sio
from scipy.signal import savgol_filter

from mihenk import Bundle

HERE = Path(__file__).resolve().parent
SOURCE = ("Kaheman et al. (2022), The Experimental Multi-Arm Pendulum on a Cart; "
          "https://github.com/dynamicslab/MultiArm-Pendulum (MIT license)")


def reference_law(repo: Path) -> tuple[dict, dict]:
    e = sio.loadmat(repo / "ParameterEstimation/SinglePendulum/EstimatedParameters_SinglePendulum.mat",
                    squeeze_me=True)["EstimatedParameters"]
    a1, m1, I1, k1, g = map(float, e)
    J = m1 * a1 ** 2 + I1
    c_sin, c_damp = a1 * g * m1 / J, -k1 / J
    params = {"a1": a1, "m1": m1, "I1": I1, "k1": k1, "g": g, "c_sin": c_sin, "c_damp": c_damp}
    law = {"theta": "omega", "omega": f"{c_sin!r}*sin(theta) + {c_damp!r}*omega"}
    return law, params


def large(repo: Path, law, params) -> Bundle:
    p = sio.loadmat(repo / "ParameterEstimation/SinglePendulum/SinglePendulumDataForParameterEstimation.mat",
                    squeeze_me=True)
    xs, ts, val = [], [], []
    for key, held_out in (("Y_id", False), ("Y_vad", True)):
        for y, t in zip(p[key], p["Time_" + key[2:]]):
            y = np.asarray(y, float)
            xs.append(y.T.copy())                 # (T, 2): theta, omega
            ts.append(np.asarray(t, float).ravel())
            val.append(held_out)
    meta = {"name": "pendulum_large", "source": SOURCE, "reference_law": law, "reference_params": params,
            "note": "theta and omega measured by encoder; omega denoised by the lab (median + wavelet)"}
    return Bundle(xs, ts, ["theta", "omega"], val=np.array(val), meta=meta)


def small(repo: Path, law, params, t0=60.0, t1=83.0, window=101, step=10) -> Bundle:
    d = sio.loadmat(repo / "Datas/SinglePendulum/SingleDataFreeSwing_1_Dt_0_0001.mat",
                    squeeze_me=True, struct_as_record=False)
    t = np.asarray(d["Theta1"].time, float)
    theta = np.asarray(d["Theta1"].signals.values, float)
    dt = float(np.median(np.diff(t)))
    omega = savgol_filter(theta, window, 3, deriv=1, delta=dt)
    keep = (t >= t0) & (t <= t1)
    idx = np.where(keep)[0][::step]
    x = np.column_stack([theta[idx], omega[idx]])
    meta = {"name": "pendulum_small", "source": SOURCE, "reference_law": law, "reference_params": params,
            "note": f"theta measured at 10 kHz; omega = Savitzky-Golay derivative ({window} samples, cubic); "
                    f"then every {step}th sample, {t0}-{t1} s"}
    return Bundle([x], [t[idx]], ["theta", "omega"], meta=meta)


def main(argv=None) -> int:
    argv = argv or sys.argv[1:]
    if not argv:
        print(__doc__)
        return 2
    repo = Path(argv[0])
    law, params = reference_law(repo)
    for b in (large(repo, law, params), small(repo, law, params)):
        out = b.save(HERE / f"{b.name}.npz")
        amp = [float(np.ptp(x[:, 0]) / 2) for x in b.xs]
        print(f"wrote {out.name}: {b.n_traj} trajectories, {sum(len(t) for t in b.ts)} samples, "
              f"theta half-ranges {', '.join(f'{a:.2f}' for a in amp)} rad")
    print(f"reference law: omega' = {params['c_sin']:.4f} sin(theta) {params['c_damp']:+.4f} omega")
    return 0


if __name__ == "__main__":
    sys.exit(main())
