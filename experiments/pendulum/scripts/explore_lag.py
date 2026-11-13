"""Exploratory, not pre-registered: is the lab's filtered velocity late?

Analysis A found theta' = 0.134 - 0.0427 theta + omega instead of theta' = omega.
A velocity signal that lags the angle by D seconds predicts exactly that:
theta' = omega_meas - 64 D (theta - pi), so D ~ 0.0427 / 64 ~ 0.67 ms.

This script (1) measures the lag directly, and (2) refits analysis A with omega
taken from the angle itself (Savitzky-Golay, 11 samples at 1 kHz, cubic).

    python experiments/pendulum/scripts/explore_lag.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import savgol_filter

import mihenk
from mihenk import Bundle

HERE = Path(__file__).resolve().parent          # experiments/pendulum/scripts
ROOT = HERE.parent                              # experiments/pendulum
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FORMATS = ("html", "json")                      # html to read, json for the numbers
sys.path.insert(0, str(HERE))
import run  # noqa: E402


def measured_lag(x, t) -> float:
    """Shift of the measured omega that best matches d(theta)/dt, in seconds."""
    wd = savgol_filter(x[:, 0], 11, 3, deriv=1, delta=np.median(np.diff(t)))
    cs = CubicSpline(t, x[:, 1])
    shifts = np.linspace(-3e-3, 3e-3, 601)
    err = [np.mean((cs(t[100:-100] + s) - wd[100:-100]) ** 2) for s in shifts]
    return float(shifts[int(np.argmin(err))])


def main() -> int:
    warnings.simplefilter("ignore")
    b = Bundle.load(DATA / "pendulum_large.npz")
    lags = [measured_lag(x, t) for x, t in zip(b.xs, b.ts)]
    print("lag of measured omega behind d(theta)/dt, per segment (ms): "
          + ", ".join(f"{1e3 * s:.2f}" for s in lags))
    xs = [np.column_stack([x[:, 0], savgol_filter(x[:, 0], 11, 3, deriv=1, delta=np.median(np.diff(t)))])
          for x, t in zip(b.xs, b.ts)]
    b1 = Bundle(xs, b.ts, b.variables, val=b.val, meta={"name": "pendulum_large_omega_from_theta"})
    sindy = run.fit(b1, run.informed_library())
    report = mihenk.scan_pysindy(sindy, b1.xs, b1.ts, val=b1.val, seed=0,
                                 name="pendulum E1: omega from the angle (exploratory)")
    report.save(RESULTS / "0.1-analyses" / "E1", formats=FORMATS)
    print(report.render_text())
    return 0


if __name__ == "__main__":
    sys.exit(main())
