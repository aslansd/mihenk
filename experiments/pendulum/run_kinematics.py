"""Development case for KINEMATIC_MISMATCH (mihenk 0.2): analyses A, B, C and E1 with theta' = omega declared.

    python experiments/pendulum/run_kinematics.py

Not a test: the check was designed after these data showed the velocity lag.
The test is experiments/kinematic_synthetic. Reports go to results/kinematics/<analysis>/.
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.signal import savgol_filter

import mihenk
from mihenk import Bundle

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run  # noqa: E402

KIN = {"theta": "omega"}


def main() -> int:
    warnings.simplefilter("ignore")
    out = {}
    large = Bundle.load(HERE / "pendulum_large.npz")
    small = Bundle.load(HERE / "pendulum_small.npz")
    e1_xs = [np.column_stack([x[:, 0], savgol_filter(x[:, 0], 11, 3, deriv=1, delta=np.median(np.diff(t)))])
             for x, t in zip(large.xs, large.ts)]
    e1 = Bundle(e1_xs, large.ts, large.variables, val=large.val, meta={"name": "pendulum_E1"})
    cases = {
        "A": (large, run.informed_library), "B": (large, run.ANALYSES["B"][1]),
        "C": (small, run.informed_library), "E1": (e1, run.informed_library),
    }
    for key, (b, lib) in cases.items():
        sindy = run.fit(b, lib())
        r = mihenk.scan_pysindy(sindy, b.xs, b.ts, val=b.val, seed=0, kinematics=KIN,
                                name=f"pendulum {key} with theta' = omega declared")
        r.save(HERE / "results" / "kinematics" / key)
        f = [x for x in r.findings if x.id == "KINEMATIC_MISMATCH"]
        out[key] = {"kinematic": f[0].severity if f else "clean", "what": f[0].what if f else "",
                    "evidence": f[0].evidence if f else {},
                    "per_trajectory": r.details.get("Kinematic fit per trajectory", {}).get("rows"),
                    "impact": r.details.get("Coefficients after correcting the kinematic mismatch", {}).get("rows")}
        print(f"===== {key}: {out[key]['kinematic']}\n{out[key]['what']}\n")
    (HERE / "results" / "kinematics" / "summary.json").write_text(json.dumps(out, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
