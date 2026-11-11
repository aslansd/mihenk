"""mihenk command line.

    mihenk scan BUNDLE --model MODEL [--out DIR] [--json] [--set name=value ...] [--daftar]
    mihenk fit BUNDLE -o MODEL.json [--degree D] [--threshold T] [--diff fd|smoothed]
    mihenk synth SYSTEM -o BUNDLE.npz [--noise F] [--n-traj K] [--seed S]
    mihenk validate [--quick]
    mihenk info

``scan`` exits with status 1 when it finds errors, so it can gate a pipeline.
"""

from __future__ import annotations

import argparse
import json
import sys
from contextlib import nullcontext
from pathlib import Path

from . import __version__, provenance
from .bundle import Bundle
from .model import Model


def _parse_sets(items):
    out = {}
    for item in items or []:
        if "=" not in item:
            raise SystemExit(f"--set expects name=value, got {item!r}")
        k, v = item.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def cmd_scan(args) -> int:
    from .scan import Thresholds, scan

    bundle = Bundle.load(args.bundle)
    model = Model.load(args.model)
    th = Thresholds().update(_parse_sets(args.set))
    refit = variants = None
    cfg = model.meta.get("pysindy_config")
    if cfg and not args.no_refit:
        try:
            from .adapters import pysindy as A

            sindy, _ = A.fit_bundle(bundle, **cfg)
            refit = A.make_refit(sindy, bundle.variables)
            variants = A.make_variants(sindy, bundle.variables)
        except ImportError:
            print("note: PySINDy is not installed, so tier-1 checks are skipped "
                  '(pip install "mihenk[pysindy]")', file=sys.stderr)
    if args.daftar and not provenance.available():
        print('daftar is not installed; run `pip install "mihenk[provenance]"`', file=sys.stderr)
        return 2
    params = {"thresholds": _parse_sets(args.set), "seed": args.seed, "bundle": str(args.bundle),
              "model": str(args.model)}
    ctx = provenance.tracked("mihenk scan", params, inputs=[args.bundle, args.model], seed=args.seed) \
        if args.daftar else nullcontext(None)
    with ctx as run:
        report = scan(bundle, model, refit=refit, variants=variants, thresholds=th, seed=args.seed,
                      kinematics=args.kinematic or None)
        out = None
        if args.out:
            out = report.save(args.out)
        provenance.log_report(run, report, out)
    if args.json:
        print(json.dumps(report.to_dict(), indent=2, default=str))
    else:
        print(report.render_text(), end="")
        if out:
            print(f"\nreport written to {out}/report.html (also .json, .md, .txt)")
    return 1 if report.has_errors else 0


def cmd_fit(args) -> int:
    from .adapters import pysindy as A

    bundle = Bundle.load(args.bundle)
    _, model = A.fit_bundle(bundle, degree=args.degree, threshold=args.threshold, diff=args.diff,
                            fourier=args.fourier)
    model.save(args.output)
    print(model.pretty())
    print(f"\nmodel written to {args.output}; next: mihenk scan {args.bundle} --model {args.output}")
    return 0


def cmd_synth(args) -> int:
    from . import synth

    if args.list or not args.system:
        for name, s in synth.SYSTEMS.items():
            print(f"{name:<24} {', '.join(s.variables):<10} {s.note}")
        return 0
    b = synth.simulate(args.system, noise=args.noise, n_traj=args.n_traj, dt=args.dt,
                       duration=args.duration, seed=args.seed)
    out = args.output or f"{b.name}.npz"
    b.save(out)
    print(f"wrote {out}: {b.n_traj} trajectory(ies), {len(b.ts[0])} samples each, variables "
          f"{', '.join(b.variables)}")
    print("true law: " + "; ".join(f"{v}' = {e}" for v, e in b.true_equations.items()))
    return 0


def cmd_validate(args) -> int:
    from . import validate

    if args.seeds > 1:
        text, ok = validate.sweep(range(args.seed, args.seed + args.seeds), quick=args.quick)
        print(text)
        return 0 if ok else 1
    results = validate.run_all(quick=args.quick, seed=args.seed)
    print(validate.render(results))
    return 0 if all(r.passed for r in results) else 1


def cmd_info(_args) -> int:
    print(f"mihenk {__version__}")
    for mod, why in (("pysindy", 'extra [pysindy]: fit, tier-1 refits'),
                     ("daftar", "extra [provenance]: scan --daftar")):
        try:
            m = __import__(mod)
            status = getattr(m, "__version__", "installed")
        except ImportError:
            status = "not installed"
        print(f"  {mod:<10} {status:<16} {why}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="mihenk", description="A touchstone for discovered equations.")
    p.add_argument("--version", action="version", version=f"mihenk {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("scan", help="check a discovered equation against its data")
    s.add_argument("bundle", help="data: .npz bundle or .csv with a 't' column")
    s.add_argument("--model", required=True, help="equation: .json from `mihenk fit`/Model.save, or a text "
                                                   "file with lines like  x' = -10 x + 10 y")
    s.add_argument("--out", help="directory for report.html/.json/.md/.txt")
    s.add_argument("--json", action="store_true", help="print the report as JSON")
    s.add_argument("--set", action="append", metavar="NAME=VALUE", help="override a threshold")
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--kinematic", action="append", metavar="P=V",
                   help="declare P' = V (e.g. theta=omega); repeatable; overrides the bundle's own pairs")
    s.add_argument("--no-refit", action="store_true", help="skip tier-1 refits even if possible")
    s.add_argument("--daftar", action="store_true", help="record the scan with daftar")
    s.set_defaults(fn=cmd_scan)

    f = sub.add_parser("fit", help="fit PySINDy (polynomial library, STLSQ) and save the equation")
    f.add_argument("bundle")
    f.add_argument("-o", "--output", default="model.json")
    f.add_argument("--degree", type=int, default=2)
    f.add_argument("--threshold", type=float, default=0.1)
    f.add_argument("--diff", choices=["fd", "smoothed"], default="fd")
    f.add_argument("--fourier", type=int, default=0, help="add FourierLibrary(n_frequencies=N)")
    f.set_defaults(fn=cmd_fit)

    y = sub.add_parser("synth", help="simulate a system with a known law")
    y.add_argument("system", nargs="?")
    y.add_argument("-o", "--output")
    y.add_argument("--noise", type=float, default=0.0, help="noise SD as a fraction of each variable's SD")
    y.add_argument("--n-traj", type=int, default=1)
    y.add_argument("--dt", type=float)
    y.add_argument("--duration", type=float)
    y.add_argument("--seed", type=int, default=0)
    y.add_argument("--list", action="store_true")
    y.set_defaults(fn=cmd_synth)

    v = sub.add_parser("validate", help="run mihenk's self-validation cases")
    v.add_argument("--quick", action="store_true")
    v.add_argument("--seed", type=int, default=0)
    v.add_argument("--seeds", type=int, default=1, help="run over this many seeds and summarize rates")
    v.set_defaults(fn=cmd_validate)

    i = sub.add_parser("info", help="version and optional dependencies")
    i.set_defaults(fn=cmd_info)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.fn(args)
    except (ValueError, KeyError, FileNotFoundError, NotImplementedError) as exc:
        print(f"mihenk: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
