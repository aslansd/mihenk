"""Run analyses A, B and C of PROTOCOL.md and save the mihenk reports.

    python experiments/pendulum/run.py

Needs pysindy. Reports go to experiments/pendulum/results/<analysis>/.
"""

from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path

import pysindy as ps

import mihenk
from mihenk import Bundle

HERE = Path(__file__).resolve().parent
THRESHOLD = 0.02  # fixed in PROTOCOL.md


def informed_library():
    return ps.PolynomialLibrary(degree=1) + ps.FourierLibrary(n_frequencies=1)


ANALYSES = {
    "A": ("pendulum_large.npz", informed_library, "informed library, large swings"),
    "B": ("pendulum_large.npz", lambda: ps.PolynomialLibrary(degree=3), "polynomial library, large swings"),
    "C": ("pendulum_small.npz", informed_library, "informed library, small swings"),
}


def fit(bundle: Bundle, library):
    pairs = list(zip(bundle.xs, bundle.ts))
    if bundle.val is not None:  # fit on the fitting segments only, as the lab did
        pairs = [p for p, v in zip(pairs, bundle.val) if not v]
    sindy = ps.SINDy(optimizer=ps.STLSQ(threshold=THRESHOLD), feature_library=library)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sindy.fit([x for x, _ in pairs], t=[t for _, t in pairs], feature_names=bundle.variables)
    return sindy


def main() -> int:
    summary = {}
    for key, (data, lib, label) in ANALYSES.items():
        b = Bundle.load(HERE / data)
        t0 = time.time()
        sindy = fit(b, lib())
        report = mihenk.scan_pysindy(sindy, b.xs, b.ts, name=f"pendulum {key}: {label}", val=b.val,
                                     meta={"source": b.meta["source"]}, seed=0)
        out = report.save(HERE / "results" / key)
        summary[key] = {
            "label": label,
            "pysindy": sindy.equations(precision=4),
            "flags": sorted(report.flagged()),
            "errors": sorted({f.id for f in report.findings if f.severity == "error" and f.tier < 2}),
            "metrics": report.metrics,
            "seconds": round(time.time() - t0, 1),
        }
        print(f"===== {key}: {label}\n{report.render_text()}")
    (HERE / "results" / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
