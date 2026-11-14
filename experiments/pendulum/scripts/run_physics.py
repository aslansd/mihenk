"""Development case for mihenk 0.3: analyses A, B, C and E1 with physics declared.

    python experiments/pendulum/scripts/run_physics.py

Declares theta' = omega and dissipates H = 0.5 omega^2 + 64.22 cos(theta)
(theta = pi hanging down, the lab's estimate of a1 g m1 / J). Not a test: the
0.3 checks were designed after these data showed what mihenk missed. The test
is experiments/physics_synthetic. Reports go to results/0.3-physics/<analysis>/.
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

HERE = Path(__file__).resolve().parent          # experiments/pendulum/scripts
ROOT = HERE.parent                              # experiments/pendulum
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FORMATS = ("html", "json")
sys.path.insert(0, str(HERE))
import run  # noqa: E402

KIN = {"theta": "omega"}
PROPS = {"dissipates": "0.5*omega**2 + 64.22*cos(theta)"}


def main() -> int:
    warnings.simplefilter("ignore")
    large = Bundle.load(DATA / "pendulum_large.npz")
    small = Bundle.load(DATA / "pendulum_small.npz")
    e1_xs = [np.column_stack([x[:, 0], savgol_filter(x[:, 0], 11, 3, deriv=1, delta=np.median(np.diff(t)))])
             for x, t in zip(large.xs, large.ts)]
    e1 = Bundle(e1_xs, large.ts, large.variables, val=large.val, meta={"name": "pendulum_E1"})
    cases = {"A": (large, run.informed_library), "B": (large, run.ANALYSES["B"][1]),
             "C": (small, run.informed_library), "E1": (e1, run.informed_library)}
    out = {}
    for key, (b, lib) in cases.items():
        sindy = run.fit(b, lib())
        r = mihenk.scan_pysindy(sindy, b.xs, b.ts, val=b.val, seed=0, kinematics=KIN, properties=PROPS,
                                name=f"pendulum {key}, with theta' = omega and dissipation declared")
        r.save(RESULTS / "0.3-physics" / key, formats=FORMATS)
        phys = r.details.get("Declared physical properties", {}).get("rows", [[None] * 7])[0]
        rolls = [row for row in r.details["Rollouts (integrating the equation)"]["rows"] if row[1] == "full"]
        under = r.details.get("Best left-out term per equation", {}).get("rows", [])
        rel = r.details.get("Relations among candidate terms on your data", {}).get("rows", [])
        out[key] = {
            "equation": r.equation,
            "flags": {f.id: f.severity for f in r.findings if f.tier < 2},
            "physics": {"data dH/dt": phys[2], "equation dH/dt": phys[3], "rollout gain": phys[5]},
            "physics_finding": next((f.what for f in r.findings if f.id == "PHYSICS_VIOLATION"), ""),
            "long_rollouts": [[row[0], row[4], row[5], row[6]] for row in rolls],
            "rollout_finding": next((f.what for f in r.findings if f.id == "ROLLOUT_TRAP"), ""),
            "underfit": [[row[0], row[1], row[3], row[4]] for row in under],
            "underfit_finding": next((f.what for f in r.findings if f.id == "UNDERFIT"), ""),
            "relations": [row[0] for row in rel],
        }
        print(f"===== {key}\n{json.dumps(out[key], indent=1, default=str)}\n")
    (RESULTS / "0.3-physics" / "summary.json").write_text(json.dumps(out, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
